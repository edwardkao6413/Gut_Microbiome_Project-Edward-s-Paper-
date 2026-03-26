# title: Supporting functions
# Authors: Joyce Hu, Edward Kao
# Date: 16 Sep 2025

### paste functions below here




## install packages
install_if_missing <- function(packages) {
  new_pkgs <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(new_pkgs)) install.packages(new_pkgs)
  invisible(lapply(packages, library, character.only = TRUE))
}


## filter metadata to derive wanted diseases
filter_meta <- function(meta, diseases, cols_to_keep = NULL) {
  
  # 1. Check input
  if (!"disease" %in% colnames(meta)) {
    stop("`meta` must contain a column named 'disease'.")
  }
  if (!is.null(cols_to_keep)) {
    missing_cols <- setdiff(cols_to_keep, colnames(meta))
    if (length(missing_cols) > 0) {
      stop("The following columns are missing from metadata: ",
           paste(missing_cols, collapse = ", "))
    }
  }
  
  # 2. Normalize disease names (case-insensitive)
  meta <- meta %>%
    mutate(disease = map_chr(disease, function(x) {
      matched <- NA
      for (pattern in diseases) {
        if (str_detect(x, regex(pattern, ignore_case = TRUE))) {
          matched <- pattern
          break
        }
      }
      if (is.na(matched)) x else matched
    }))
  
  # 3. Filter metadata
  filtered_meta <- meta %>% filter(disease %in% diseases)
  
  # 4. Select columns
  if (!is.null(cols_to_keep)) {
    filtered_meta <- filtered_meta %>% dplyr::select(disease, all_of(cols_to_keep))
  } else {
    filtered_meta <- filtered_meta %>% dplyr::select(disease)
  }
  
  # 5. Return (row names preserved automatically by R)
  return(filtered_meta)
}





## Split Count Filtering by Disease
split_count_filtering <- function(labelled_abund, diseases, percent, filtered_meta) {
  
  # Case: less than 2 diseases
  if (length(diseases) < 2) {
    message("Not enough diseases to split; running low count removal on full dataset")
    res <- low_count_removal(labelled_abund, percent)
    abundance_filtered <- res$data.filter
  } else {
    all_keep_otu <- c()
    
    for (disease_name in diseases) {
      
      # Subset samples matching this disease, preserving rownames
      disease_sids <- rownames(filtered_meta)[filtered_meta$disease == disease_name]
      disease_df <- labelled_abund[disease_sids, , drop = FALSE]
      
      
      # Apply low count removal
      res <- low_count_removal(disease_df, percent)
      
      all_keep_otu <- union(all_keep_otu, res$keep.otu)
    }
    
    # Filter original abundance table by OTUs that survived across all groups
    
    abundance_filtered <- labelled_abund[, all_keep_otu, drop = FALSE]
    
    # Ensure numeric matrix, no lists, preserve rownames
    abundance_filtered <- as.data.frame(lapply(abundance_filtered, function(col) {
      if (is.list(col)) col <- unlist(col)
      as.numeric(col)
    }), row.names = rownames(abundance_filtered), check.names = FALSE)
  }
  #View(abundance_filtered)
  return(abundance_filtered)
}



## low counts removal
low_count_removal <- function(
    data,    # OTU count dataframe of size n (samples) x p (OTUs)
    percent = 0.01  # cutoff percentage threshold
) {
  
  # Check input type
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a dataframe or matrix of counts.", call. = FALSE)
  }
  
  if (!all(sapply(data, is.numeric))) stop("All columns must be numeric.")
  
  
  # Identify OTUs above the relative abundance cutoff
  keep.otu <- which(
    colSums(data) * 100 / sum(colSums(data)) > percent
  )
  
  # Filter the dataframe to keep only those OTUs
  data.filter <- data[, keep.otu, drop = FALSE]
  
  # Return both filtered data and indices
  return(list(
    data.filter = data.filter,
    keep.otu = keep.otu
  ))
}


## CLR Transformation
clr.transformation <- function(abundance.filtered, offset_value) {
  # Apply CLR transformation using MixOmics
  # logratio.transfo requires a matrix input
  data.clr <- logratio.transfo(
    as.matrix(abundance.filtered),
    logratio = "CLR",
    offset = offset_value
  )
  
  # Convert result to a regular dataframe
  clr_mat <- unclass(data.clr)
  clr_df <- as.data.frame(clr_mat, check.names = FALSE)
  
  return(clr_df)
}



## processing single dataset
preprocess_single <- function(meta, abund, diseases, keep_cols, percent, offset_value, CLR=TRUE) {
  # 1. Filter metadata
  filtered_meta <- filter_meta(meta, diseases, keep_cols)
  #View(filtered_meta)
  # 2. Select abundance rows that match filtered_meta
  selected_abund <- as.data.frame(abund[rownames(abund) %in% rownames(filtered_meta), ], check.names = FALSE)
  #View(selected_abund)
  # 3. Apply split count filtering
  abundance_filtered <- split_count_filtering(selected_abund, diseases, percent, filtered_meta)
  
  if (!CLR) return(abundance_filtered)
  # 5. Convert to matrix for CLR
  abundance_matrix <- as.matrix(abundance_filtered)
  
  # 6. Apply CLR transformation
  clr_transformed <- clr.transformation(abundance_matrix, offset_value)
  #View(clr_transformed)
  # 7. Return results
  clr_transformeddf <- as.data.frame(clr_transformed, check.names = FALSE)
  return(list(
    data = clr_transformeddf,
    meta = filtered_meta
  ))
  
}


## Process multiple datasets
preprocess_multiple <- function(abundance_list, meta_list, diseases, keep_cols, percent, offset_value, threshold = NULL) {
  results <- list()
  # iterate through abund and meta lists to proprocess individually
  for (nm in names(abundance_list)) {
    abund <- abundance_list[[nm]]
    meta  <- meta_list[[nm]]
    result <- preprocess_single(meta, abund, diseases, keep_cols, percent, offset_value)
    
    # store it inside results under the dataset named results accessed by name
    results[[nm]] <- list(
      data = result$data,
      meta = result$meta
    )
    cat("\n--- Dataset:", nm, "---\n")
    print(dim(result$data))
    # View(result$data)
    # View(result$meta)
  }
  
  combined <- integrate_filter(results, threshold)
  return(combined)
}



## Integrate Filtered Abundance Across Studies
integrate_filter <- function(results, threshold) {
  
  # Step 1: Determine taxa present above threshold across studies
  taxa_present_list <- purrr::map(results, function(study) {
    taxa_cols <- colnames(study$data)
    present_taxa <- taxa_cols[
      purrr::map_lgl(taxa_cols, ~ sum(study$data[[.x]], na.rm = TRUE) > 0)
    ]
    return(present_taxa)
  })
  
  taxa_counts <- table(unlist(taxa_present_list))
  surviving_taxa <- names(taxa_counts)[taxa_counts >= threshold]
  
  # Step 2: Filter each study
  comb_data_list <- list()
  meta_list <- list()
  
  
  
  
  final <- purrr::imap(results, function(study, study_name) {
    
    # Filter abundance to surviving taxa
    abund_filtered <- study$data %>%
      dplyr::select(any_of(surviving_taxa))
    
    # Preserve rownames
    rownames(abund_filtered) <- rownames(study$data)
    
    # Filter metadata to match abundance rows
    meta_filtered <- study$meta[rownames(abund_filtered), , drop = FALSE]
    
    # Return a nested list, automatically named by imap using study_name
    list(
      abund = abund_filtered,
      meta = meta_filtered
    )
  })
  
  
  return(final)
}



## Run Preprocessing on Single or Multiple Datasets
run_preprocessing <- function(abundance_list, meta_list, diseases,
                              keep_cols, percent, offset_value, threshold = NULL) {
  # Check how many datasets are in the abundance list
  if (length(abundance_list) == 1) {
    print("single")
    # Only one dataset → run preprocess_single
    #how to access this type of database
    # abund <-  as.data.frame((abundance_list)[[1]])
    # meta  <- as.data.frame((meta_list)[[1]])
    abund <-  as.data.frame(abundance_list[[1]], check.names = FALSE)
    meta  <- as.data.frame(meta_list[[1]], check.names = FALSE)
    result <- preprocess_single(meta, abund, diseases, keep_cols, percent, offset_value)
    
  } else {
    print("Multiple")
    # Multiple datasets → run preprocess_combine (preprocess_multiple)
    result <- preprocess_multiple(abundance_list, meta_list,
                                  diseases, keep_cols, percent, offset_value, threshold)
  }
  return(result)
}


## Individual filtering
individual_filtering <- function(data, err_thres = 0.55, ncomp = 2,
                                 validation = "Mfold", folds = 5, nrepeat = 10, dist = "centroids.dist") {
  
  # Helper to check a single study
  ## parameter study here means a single dataset, along with its abundance matrix and metadata
  check_single_dataset <- function(study) {
    abund <- study$abund   # renamed from bund
    meta  <- study$meta
    Y <- factor(meta$disease)
    
    plsda.model <- plsda(abund, Y, ncomp = ncomp)
    perf.model <- perf(
      object = plsda.model,
      validation = validation,
      folds = folds,
      progressBar = FALSE,
      dist = dist, 
      nrepeat = nrepeat
    )
    
    #class_errors <- perf.model$error.rate.class[[dist]]
    #exceed_threshold <- apply(class_errors, 1, function(x) any(x > err_thres))
    class_errors <- perf.model$error.rate.class[[dist]][, ncomp]
    exceed_threshold <- any(class_errors > err_thres)
    
    keep_dataset <- !any(exceed_threshold)
    return(keep_dataset)
  }
  
  # Case 1: single study (has both abund and meta)
  if (is.list(data) && all(c("abund","meta") %in% names(data))) {
    return(check_single_dataset(data))
  }
  
  # Case 2: nested list of studies
  if (is.list(data) && !all(c("abund","meta") %in% names(data))) {
    results <- purrr::imap(data, function(study, study_name) {
      if (!all(c("abund","meta") %in% names(study))) {
        stop(paste("Study", study_name, "is missing abund or meta"))
      }
      keep <- check_single_dataset(study)
      if (keep) return(study) else return(NULL)
    })
    results <- results[!purrr::map_lgl(results, is.null)]
    return(results)
  }
  
  stop("Input must be a single study (abund + meta) or a nested list of studies")
}



## Calculate the BER for single study, with 5 folds for CV and 10 times repeats
# output the BER value
calculate_BER = function(single_study, ncomp=2, validation="Mfold", 
                         folds=5, nrepeat=10, dist = "centroids.dist") {
  
  abund = single_study$abund
  meta = single_study$meta
  Y = factor(meta$disease)
  
  plsda.model = plsda(X=abund, Y=Y, ncomp=ncomp)
  perf.model = perf(
    object = plsda.model,
    validation = validation,
    folds = folds,
    progressBar = FALSE,
    dist = dist, 
    nrepeat = nrepeat
  )
  ## return the error rate of all the components, two classes
  return(perf.model$error.rate$BER[ncomp, ])
}


## Eliminate the taxa that only exist in single study
# study_col: let user specify the column name for study name; n_taxa. eliminate the taxa that only exist in n studies across all studies
eliminate_single_taxa = function(lst_df, study_col="study", n_taxa=1) {
  
  dff = data.frame()
  df = lst_df$abund
  df$study = lst_df$meta[[study_col]]
  # df$disease = lst_df$meta$disease
  taxa_lst = setdiff(colnames(df), c(study_col))
  for (species in taxa_lst) {
    # Initialize a temporary data frame for the current species
    df_taxa <- data.frame(species = species, stringsAsFactors = FALSE)
    
    # Check for presence in each study
    for (study in unique(df$study)) {
      # Filter data for the current study and species
      study_data <- df[df$study == study, species]
      
      # Determine if the species is detected or not
      if (length(unique(study_data)) == n_taxa) {
        df_taxa[[study]] <- "x"
      } else {
        df_taxa[[study]] <- "Exist"
      }
    }
    
    # Append to the main data frame
    dff <- rbind(dff, df_taxa)
  }
  
  dff1 = dff
  dff1[, unique(df$study)] <- lapply(dff1[, unique(df$study)], function(x) ifelse(x == "Exist", 1, 0))
  dff$Exist_count <- rowSums(dff1[, unique(df$study)])
  one_taxa <- dff[dff$Exist_count == 1, "species"]
  return(df[, !(colnames(df) %in% one_taxa)])
}


## Calculate BER for integrative dataset including multiple studies, eliminate_single_taxa function is embedded.
calculate_multipleDB_BER = function(lst, ncomp=2, study_col="study", validation="Mfold", 
                                    folds=5, nrepeat=10, dist="centroids.dist", n_taxa=1) {
  
  ## add study column in each metadata matrix
  for (name in names(lst)) {
    lst[[name]]$meta$study = strsplit(name, '\\.')[[1]][1]
    # lst[[name]]$meta$study = name
  }
  # binding together
  lst_df = bind_rows(lst)
  
  # eliminate taxa only appear in one study
  comb_df = eliminate_single_taxa(lst_df, study_col = study_col, n_taxa=n_taxa)
  abund_X = comb_df[ ,setdiff(names(comb_df), study_col), drop = FALSE]
  Y = lst_df$meta$disease
  corresponding_study = comb_df$study
  model = mint.splsda(X=abund_X, Y=Y, study=corresponding_study, ncomp=ncomp)
  perf.model = perf(model, validation=validation, folds=folds, dist=dist)
  
  # return the global overall BER
  return(perf.model$global.error$BER[ncomp, ])
}


## Forward Selection
# data: study pool
forward_selection = function(data, ncomp=2, study_col="study", validation='Mfold', folds=5, nrepeat=20, 
                             dist="centroids.dist", n_taxa=1) {
  remain_studies = list()   # studies included in combination
  i = 0  # calculate round
  
  ## tell which study should be put in the combination as first study
  cat("########## ROUND ", i+1, " ##########", "\n")
  cat("Total study: ", length(data), "\n")
  cat("Processing... Testing studies: ", names(data), "\n")
  first_study = ""
  ERR = 1
  for (study in names(data)) {
    study_ERR = calculate_BER(data[[study]])
    
    if ( study_ERR < ERR ) {
      first_study = study
      ERR = study_ERR
    }
  }
  # input the first study
  cat("Introduce ", first_study, " into combination.", "\n")
  remain_studies[[first_study]] = data[[first_study]]   # plugin the first study into combination permanently
  cat("########## Current combination:", names(remain_studies), " round ", i+1, "\n") 
  data[[first_study]] = NULL     # delete first study from the study pool permanently
  for (j in 1:2) {
    cat("                                     ", "\n")
  }
  
  
  # start formal forward selection
  dynamic_BER = 1  # lowest BER of combination
  current_BER = 1  # BER for current combination
  while ( (length(data) > 0) | (current_BER > dynamic_BER) ) {
    ## searching for the datasets in all available studies
    i = i + 1
    cat("########## ROUND ", i+1, " ##########", "\n")
    
    combin_err = 1
    selected_study = ""
    for (study in names(data)) {
      ## add study into combination iteratively
      remain_studies[[study]] = data[[study]]; # cat("adding ", study, " in combination", "\n")
      err = calculate_multipleDB_BER(remain_studies, ncomp=ncomp, study_col=study_col, 
                                     validation=validation, dist=dist, n_taxa=n_taxa)
      cat( "Processing... Test new combination: ", names(remain_studies), "\n" )
      
      if (err < combin_err) {
        combin_err = err
        selected_study = study
      }
      # delete the inserted studies, and go calculating BER for adding next study
      remain_studies[[study]] = NULL; #  cat("removing ", study, " out of combination", "\n")
    }
    
    # After calculating all the combinations, calculate the combination of BER that generates the lowest BER
    #cat("grab studies ", "(", selected_study, ")", " into round ", i , "\n")
    #cat("Its BER is ", combin_err, "\n")
    remain_studies[[selected_study]] = data[[selected_study]]  # formally add a new study into current combination
    data[[selected_study]] = NULL    # delete the study which has been added into current combination from the study pool
    current_BER = combin_err  # assign the lowest BER produced in current number of combination
    #current_BER = calculate_multipleDB_BER(remain_studies, ncomp=ncomp, study_col=study_col, 
    #                                 validation=validation, dist=dist)
    
    if (current_BER < dynamic_BER) {
      #cat("current_BER: ", current_BER, "\n")
      #cat("dynamic_BER: ", dynamic_BER, "\n")
      #cat("current BER lower than dynamic global BER, assign it as new dynamic global BER", "\n")
      dynamic_BER = current_BER
      cat("Introduce ", selected_study, " into combination.", "\n")
      cat("########## Current combination:", names(remain_studies), " round ", i+1, "\n")
      # remain_studies[[selected_study]] = data[[selected_study]]
    } else {
      cat("No improvement in BER after introducing new dataset into combination", "\n")
      cat("########## Current combination:", names(remain_studies), " round ", i+1, "\n")
      #cat("current_BER: ", current_BER, "\n")
      #cat("dynamic_BER: ", dynamic_BER, "\n")
      # Bcz adding new study cannot bring us lower BER, we kick it out
      remain_studies[[selected_study]] = NULL
      break
    }
    
    for (j in 1:2) {
      cat("                                     ", "\n")
    }
  }
  
  if (length(data) == 0) {
    print("All studies have been selected since BER is improved round to round.")
  }
  cat("                                     ", "\n")
  cat("Number of study in the final combination: ", length(remain_studies), "\n")
  cat("Studies include: ", names(remain_studies) )
  return(remain_studies)
}



# draw testing set from HGMA
draw_TestSet = function(threshold_for_LowCountsRm=0.01, disease_lst=c("T2D", "Healthy", "NGT", "healthy"), dataset_study) {
  # ingest rawdata
  print("Reading HGMA rawdata...")
  sample_md = read_excel("sample_metadata.xlsx")
  abund_data = read.csv("vect_atlas.csv\\vect_atlas.csv")
  corr_taxa = read.csv("corresponding_taxa.csv")
  
  ## Preprocessing 1st round
  print("Processing metadata and rawdata...")
  corr_taxa <- corr_taxa %>%
    filter(!str_detect(name, "unclassified"))
  
  multi_taxa <- corr_taxa %>%
    filter(str_detect(name, "\\s\\d$"))
  
  multi_taxa <- multi_taxa %>%
    mutate(name = str_sub(name, 1, -3))
  
  multi_taxa_lst <- unique(multi_taxa$name)
  
  # Process metadata and rawdata
  metadata = sample_md[sample_md['BioProject'] == dataset_study, ]
  metadata = metadata %>% 
    filter(Disease %in% disease_lst)
  sp_id = metadata$sample.ID
  rawdata = abund_data[, c(sp_id, "X")]
  rawdata <- merge(
    rawdata,
    corr_taxa,
    by.x = "X",
    by.y = "id",
    all.x = TRUE  # all.x=TRUE: Keep all rows from x even if no matching key in y (left join)
  )
  rawdata <- rawdata[, !(names(rawdata) %in% c("X", "id"))]  # drop 'X' and 'id' column
  rawdata <- rawdata[!is.na(rawdata$name), ]  # drop the NaN value in 'name' column
  
  
  taxa_lst = rawdata$name
  rawdata_v1 = data.frame(t(rawdata[, !(names(rawdata) %in% c("name"))]))
  colnames(rawdata_v1) = taxa_lst
  names(rawdata_v1)[names(rawdata_v1) == "index"] = "sample_id"
  rawdata_v1$sample_id = rownames(rawdata_v1)
  
  rawdata_v1 = merge(rawdata_v1, 
                     metadata[, c("sample.ID", "Disease")], 
                     by.x="sample_id", 
                     by.y="sample.ID",
                     all.x=TRUE)
  names(rawdata_v1)[names(rawdata_v1) == "Blautia coccoides == Blautia producta"] = "Blautia producta"
  
  check_main_species <- function(lst) {
    integers <- as.character(1:9)
    for (sub in lst) {
      last_char <- substr(sub, nchar(sub), nchar(sub))
      if (!(last_char %in% integers)) {
        return(TRUE)
      }
    }
    return(FALSE)
  }
  
  for (taxa in multi_taxa_lst) {
    if ( !grepl("subtype", tolower(taxa), fixed = T) ) {
      cols = colnames(rawdata_v1)
      rawdata_v1$taxa = 0
      alltaxa_columns <- cols[grepl(taxa, cols, fixed = T)]
      
      if ( length(alltaxa_columns) == 1 ) {
        rawdata_v1[[taxa]] <- rawdata_v1[[alltaxa_columns[1]]]
      } else if ( check_main_species(alltaxa_columns) == T ) {
        alltaxa_columns <- setdiff(alltaxa_columns, taxa)
        rawdata_v1 <- rawdata_v1[, !(names(rawdata_v1) %in% alltaxa_columns), drop = FALSE]
        
      } else {
        rawdata_v1[[taxa]] <- rowSums(rawdata_v1[, alltaxa_columns, drop = FALSE], na.rm = TRUE)
        rawdata_v1 <- rawdata_v1[, !(names(rawdata_v1) %in% alltaxa_columns), drop = FALSE]
      }
    }
  }
  
  # rename the columns
  rename_map <- c(
    "Lachnospiraceae bacterium hRUG904 / Clostridiales bacterium RUG495 / [Clostridium] sp. CAG:127" = "[Clostridium] sp. CAG:127",
    "Firmicutes bacterium CAG:65 / [Clostridium] sp. 2789STDY5608883" = "Firmicutes bacterium CAG:65",
    "Lachnospiraceae bacterium CIM:MAG 844 / Firmicutes bacterium CAG:95" = "Firmicutes bacterium CAG:95",
    "[Ruminococcus] sp. CAG:17 / Blautia sp. 2789STDY5608880" = "[Ruminococcus] sp. CAG:17",
    "Roseburia sp. CAG:100 / Lachnospiraceae bacterium CIM:MAG 54" = "Roseburia sp. CAG:100",
    "Clostridiales bacterium UBA7273 / Firmicutes bacterium CAG:240" = "Firmicutes bacterium CAG:240",
    "Firmicutes bacterium CAG:24 / [Clostridium] sp. 2789STDY5834869 & 2789STDY5834922" = "Firmicutes bacterium CAG:24",
    "Coprococcus sp. 2789STDY5608819 / [Clostridium] sp. CAG:264" = "[Clostridium] sp. CAG:264",
    "[Ruminococcus] sp. 2789STDY5608794 & sp. 2789STDY5834890 / Firmicutes bacterium CAG:56" = "Firmicutes bacterium CAG:56",
    "Lachnospiraceae bacterium TF01-11 / [Clostridium] sp. CAG:122" = "[Clostridium] sp. CAG:122",
    "[Acinetobacter] sp. N54.MGS-139 / Proteobacteria bacterium CAG:139" = "Proteobacteria bacterium CAG:139",
    "[Ruminococcus] sp. CAG:60 / Blautia sp. 2789STDY5608836" = "[Ruminococcus] sp. CAG:60",
    "[Eubacterium] sp. CAG:251 / [Clostridium] sp. A254.MGS-251" = "[Eubacterium] sp. CAG:251",
    "[Clostridium] sp. CAG:217 / Clostridiales bacterium CIM:MAG 317_1" = "[Clostridium] sp. CAG:217",
    "Lachnoclostridium sp. SNUG30386 / [Clostridium] sp. CAG:43" = "[Clostridium] sp. CAG:43",
    "[Ruminococcaceae] bacterium D16 / Clostridiales bacterium Choco116" = "[Ruminococcaceae] bacterium D16",
    "Firmicutes bacterium CAG:41 / [Clostridium] sp. 2789STDY5834935 & sp. 2789STDY5608853" = "Firmicutes bacterium CAG:41",
    "[Clostridiaceae] bacterium CIM:MAG 755 / [Clostridium] sp. CAG:230" = "[Clostridium] sp. CAG:230",
    "Firmicutes bacterium CAG:212 / [Clostridium] sp. 2789STDY5834871" = "Firmicutes bacterium CAG:212",
    "Bacteroidales bacterium H2 / [Prevotella] sp. CAG:251" = "[Prevotella] sp. CAG:251",
    "Ruminococcaceae bacterium UBA1821 / [Clostridium] sp. CAG:678" = "[Clostridium] sp. CAG:678",
    "Oscillibacter sp. ER4 / Firmicutes bacterium CAG:129_59_24" = "Oscillibacter sp. ER4",
    "Clostridiales bacterium RUG303 & RUG439 & RUG453 / Firmicutes bacterium CAG:129" = "Firmicutes bacterium CAG:129",
    "Candidatus Methanomethylophilus alvus / Methanoculleus sp. CAG:1088" = "Methanoculleus sp. CAG:1088",
    "[Clostridiaceae] bacterium CIM:MAG 987 / [Clostridium] sp. CAG:533" = "[Clostridium] sp. CAG:533",
    "[Firmicutes] bacterium CIM:MAG 721 / [Bacillus] sp. CAG:988" = "[Bacillus] sp. CAG:988", 
    "Disease" = "disease"
  )
  
  cols <- colnames(rawdata_v1)
  hits <- cols %in% names(rename_map)
  cols[hits] <- unname(rename_map[cols[hits]])
  colnames(rawdata_v1) <- cols
  
  print("Organizing the metadata and the rawdata into proper output form")
  hgma_lst = list(abund=subset(rawdata_v1, select=-c(sample_id, disease)), meta=subset(rawdata_v1, select=c(sample_id, disease)) )
  # calculate the second smallest value
  x <- unlist(hgma_lst$abund)
  x <- x[x > 0]
  second_smallest <- sort(unique(x), na.last = NA)[2]
  processed_lst = preprocess_single(meta=hgma_lst$meta, abund=hgma_lst$abund, diseases=disease_lst, 
                                    keep_cols=c("sample_id"), offset_value=second_smallest, percent=0.01, CLR=T )
  names(processed_lst)[names(processed_lst) == "data"] = "abund"
  processed_lst$meta$study = dataset_study
  
  return(processed_lst)
}


## Feature alignment: CMD (training) vs HGMA (testing) column names
# Promoted from Test_Framework.Rmd — Edward Kao, 2025
# Strips 'species:' prefix from both matrices, then resolves HGMA '/' alias columns
# against CMD training column names. Returns aligned training and testing matrices
# containing only their common taxa, plus a change_log of all renames performed.
align_columns <- function(p_test, p_selected) {

  # Remove 'species:' prefix from column names
  colnames(p_selected) <- sub("^species:\\s*", "", colnames(p_selected), ignore.case = TRUE)
  colnames(p_test)     <- sub("^species:\\s*", "", colnames(p_test),     ignore.case = TRUE)

  change_log  <- list()
  message_log <- character()

  for (i in seq_along(colnames(p_test))) {

    col_test <- colnames(p_test)[i]

    # Only process test columns that contain '/'
    if (grepl("/", col_test, fixed = TRUE)) {

      taxa_test         <- trimws(unlist(strsplit(col_test, "/", fixed = TRUE)))
      matched_train_cols <- c()

      # Find matching training columns
      for (col_train in colnames(p_selected)) {
        taxa_train <- trimws(unlist(strsplit(col_train, "/", fixed = TRUE)))
        if (all(taxa_train %in% taxa_test)) {
          matched_train_cols <- c(matched_train_cols, col_train)
        }
      }

      # Single match → rename p_test column to the training name
      if (length(matched_train_cols) == 1) {
        old_name            <- col_test
        colnames(p_test)[i] <- matched_train_cols[1]
        change_log[[old_name]] <- paste0("Single match: renamed test column to '", matched_train_cols[1], "'")
        message_log <- c(message_log,
          paste0("Single match found: test column '", old_name, "' renamed to '", matched_train_cols[1], "'"))

      # Multiple matches → sum matched training columns into one combined column
      } else if (length(matched_train_cols) > 1) {
        new_train_col <- paste(matched_train_cols, collapse = "/")

        if (!(new_train_col %in% colnames(p_selected))) {
          p_selected[[new_train_col]] <- rowSums(p_selected[, matched_train_cols, drop = FALSE])
          change_log[[paste0("p_selected: ", new_train_col)]] <- paste0(
            "Created new combined column from: ", paste(matched_train_cols, collapse = ", "))
        }

        old_name            <- col_test
        colnames(p_test)[i] <- new_train_col
        change_log[[old_name]] <- paste0("Multiple matches: renamed test column to '", new_train_col, "'")
      }
      # No match → leave column name unchanged
    }
  }

  # Subset both matrices to their common taxa
  common_cols <- intersect(colnames(p_selected), colnames(p_test))
  training    <- p_selected[, common_cols, drop = FALSE]
  testing     <- p_test[,     common_cols, drop = FALSE]

  return(list(
    training   = training,
    testing    = testing,
    change_log = change_log,
    messages   = paste(message_log, collapse = "\n")
  ))
}


## MINT sPLS-DA modelling: sequential hyperparameter tuning + model fitting + evaluation
# Promoted from Joyce/modelling.Rmd — Joyce Hu, 2025
# Tunes comp 1 first, then comp 2 using already.tested.X (mixOmics best-practice).
# study must exist in the calling environment as a factor vector (one entry per sample).
# Returns: model, optimal.keepX, performance, selected features per component, all unique features.
run_mint <- function(X, Y, ncomp = 2, dist = "centroids.dist", nrepeat_perf=20, fold_perf=5) {

  # STEP 1: Tune component 1
  tune.comp1 <- tune(method = "mint.splsda",
                     X = X, Y = Y, study = study,
                     ncomp = 1,
                     test.keepX = seq(10, ncol(X), 1),
                     dist = dist)

  optimal.comp1 <- tune.comp1$choice.keepX
  cat("Optimal keepX for Component 1:", optimal.comp1, "\n\n")
  plot(tune.comp1, main = "Component 1 Tuning")

  # STEP 2: Tune component 2, conditioning on optimal comp 1
  cat("Step 2: Tuning Component 2...\n")
  tune.comp2 <- tune(method = "mint.splsda",
                     X = X, Y = Y, study = study,
                     ncomp = 2,
                     test.keepX = seq(10, ncol(X), 1),
                     already.tested.X = optimal.comp1,
                     dist = dist)

  optimal.comp2 <- tune.comp2$choice.keepX[2]
  cat("Optimal keepX for Component 2:", optimal.comp2, "\n\n")
  plot(tune.comp2, main = "Component 2 Tuning")

  optimal.keepX <- c(optimal.comp1, optimal.comp2)

  # STEP 3: Build final model
  cat("Step 3: Building final model...\n")
  cat("Using keepX =", paste(optimal.keepX, collapse = ", "), "\n\n")

  model <- mint.splsda(X = X, Y = Y, study = study,
                       ncomp = ncomp,
                       keepX = optimal.keepX)

  # STEP 4: Evaluate model (leave-one-study-out CV — default for MINT)
  cat("Step 4: Evaluating model...\n")
  perf <- perf(model, validation="Mfolds", folds=fold_perf, nrepeat=nrepeat_perf, dist=dist)

  cat("\nBalanced Error Rate (BER):\n")
  print(perf$global.error$BER)

  cat("\nError Rate per Class:\n")
  print(perf$global.error$error.rate.class)

  # Extract selected taxa for each component
  selected.comp1 <- selectVar(model, comp = 1)$name
  selected.comp2 <- selectVar(model, comp = 2)$name
  ttl_taxa <- tolower(unique(c(selected.comp1, selected.comp2)))

  results <- list(
    model                  = model,
    optimal.keepX          = optimal.keepX,
    performance            = perf,
    selected.features.comp1 = selected.comp1,
    selected.features.comp2 = selected.comp2,
    features.all           = ttl_taxa
  )

  return(results)
}


## MINT sPLS-DA visualisation: score plot + loadings for both components
# Promoted from Joyce/modelling.Rmd — Joyce Hu, 2025
# results: output list from run_mint(); Y: disease factor; study: study factor
plot_mint_results <- function(results, Y, study) {

  model <- results$model

  # Score plot (samples coloured by disease, shaped by study)
  plotIndiv(model,
            group = Y,
            pch = as.numeric(factor(study)) + 13,
            pch.levels = study,
            title = 'MINT sPLS-DA',
            legend = TRUE,
            legend.title = 'Disease',
            legend.title.pch = 'Study',
            ellipse = TRUE)

  # Loadings — Component 1
  plotLoadings(model,
               comp = 1,
               legend = TRUE,
               contrib = 'max',
               title = 'Feature Loadings - Component 1',
               method = 'median',
               ndisplay = 20)
  legend("topright", legend = levels(Y), fill = c("orange", "blue"))

  # Loadings — Component 2
  plotLoadings(model,
               comp = 2,
               legend = TRUE,
               contrib = 'max',
               title = 'Feature Loadings - Component 2',
               method = 'median',
               ndisplay = 20)
  legend("topright", legend = levels(Y), fill = c("orange", "blue"))
}


