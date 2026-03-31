# title: Supporting functions
# Authors: Joyce Hu, Edward Kao
# Date: 16 Sep 2025
#
# OVERVIEW
# --------
# This file defines all reusable helper functions for the gut microbiome MINT pipeline.
# Source this file at the top of every .Rmd notebook:  source("Function_Helpers.R")
#
# PIPELINE LAYERS (in execution order)
# 1. install_if_missing        – package management
# 2. filter_meta               – metadata filtering by disease label
# 3. low_count_removal         – per-disease OTU prevalence filtering
# 4. split_count_filtering     – wrapper: per-disease low-count removal + union
# 5. clr.transformation        – centered log-ratio (CLR) normalisation
# 6. preprocess_single         – full preprocessing for one dataset
# 7. preprocess_multiple       – iterates preprocess_single across multiple datasets
# 8. integrate_filter          – cross-study taxa intersection
# 9. run_preprocessing         – entry point: routes to single or multiple
# 10. individual_filtering     – Layer 2: per-study PLS-DA BER check
# 11. calculate_BER            – helper: BER for one study
# 12. eliminate_single_taxa    – helper: drop taxa present in only one study
# 13. calculate_multipleDB_BER – helper: BER for a multi-study combination
# 14. forward_selection        – Layer 3: greedy study selection by BER
# 15. draw_TestSet             – load and preprocess HGMA external test data
# 16. align_columns            – align training (CMD) vs testing (HGMA) column names
# 17. run_mint                 – Layer 4/5: MINT sPLS-DA tuning + fitting + evaluation
# 18. plot_mint_results        – visualise MINT score plot and feature loadings

### paste functions below here




## ============================================================
## 1. PACKAGE MANAGEMENT
## ============================================================

# install_if_missing(packages)
# ----------------------------
# Installs any packages from `packages` that are not yet installed,
# then loads all of them into the current R session.
# Parameters:
#   packages – character vector of CRAN/Bioconductor package names
install_if_missing <- function(packages) {
  # Identify packages not yet installed
  new_pkgs <- packages[!(packages %in% installed.packages()[,"Package"])]
  # Install any missing packages
  if (length(new_pkgs)) install.packages(new_pkgs)
  # Load every package (returns invisibly to suppress the list output)
  invisible(lapply(packages, library, character.only = TRUE))
}




## ============================================================
## 2. METADATA FILTERING
## ============================================================

# filter_meta(meta, diseases, cols_to_keep)
# ------------------------------------------
# Filters a metadata data frame to retain only rows belonging to the
# specified disease groups, then selects the `disease` column plus any
# additional columns requested via cols_to_keep.
#
# Parameters:
#   meta         – data frame; rows = samples, must contain a column named "disease"
#   diseases     – character vector of disease labels to keep (e.g. c("Healthy","T2D"))
#   cols_to_keep – optional character vector of extra columns to retain (e.g. c("gender","body_site"))
#
# Returns: filtered metadata data frame with `disease` + requested columns;
#          row names are sample IDs (preserved from the input).
filter_meta <- function(meta, diseases, cols_to_keep = NULL) {

  # 1. Validate that the required 'disease' column is present
  if (!"disease" %in% colnames(meta)) {
    stop("`meta` must contain a column named 'disease'.")
  }
  # 1b. Validate that all requested extra columns exist
  if (!is.null(cols_to_keep)) {
    missing_cols <- setdiff(cols_to_keep, colnames(meta))
    if (length(missing_cols) > 0) {
      stop("The following columns are missing from metadata: ",
           paste(missing_cols, collapse = ", "))
    }
  }

  # 2. Normalize disease labels: map raw CMD disease strings to canonical labels
  #    (case-insensitive regex match against the user-supplied disease list)
  meta <- meta %>%
    mutate(disease = map_chr(disease, function(x) {
      matched <- NA
      for (pattern in diseases) {
        # Try matching each canonical label against the raw value
        if (str_detect(x, regex(pattern, ignore_case = TRUE))) {
          matched <- pattern
          break
        }
      }
      # If no pattern matched, keep the original value (will be filtered out next)
      if (is.na(matched)) x else matched
    }))

  # 3. Keep only samples whose disease label is in the requested set
  filtered_meta <- meta %>% filter(disease %in% diseases)

  # 4. Select disease + optional extra columns
  if (!is.null(cols_to_keep)) {
    filtered_meta <- filtered_meta %>% dplyr::select(disease, all_of(cols_to_keep))
  } else {
    filtered_meta <- filtered_meta %>% dplyr::select(disease)
  }

  # 5. Return the filtered metadata; row names (sample IDs) are preserved by R automatically
  return(filtered_meta)
}




## ============================================================
## 3 & 4. PREVALENCE FILTERING
## ============================================================

# split_count_filtering(labelled_abund, diseases, percent, filtered_meta)
# -----------------------------------------------------------------------
# Applies low-count (prevalence) filtering separately within each disease
# group, then takes the UNION of surviving OTUs across groups. This ensures
# taxa present in at least one disease class are not accidentally dropped.
#
# Parameters:
#   labelled_abund – data frame; rows = samples, columns = taxa (OTU counts)
#   diseases       – character vector of disease labels used to subset samples
#   percent        – minimum relative-abundance threshold (passed to low_count_removal)
#   filtered_meta  – filtered metadata; row names must match row names in labelled_abund
#
# Returns: filtered abundance data frame (samples × surviving taxa)
split_count_filtering <- function(labelled_abund, diseases, percent, filtered_meta) {

  # Edge case: fewer than 2 disease groups → skip per-group splitting
  if (length(diseases) < 2) {
    message("Not enough diseases to split; running low count removal on full dataset")
    res <- low_count_removal(labelled_abund, percent)
    abundance_filtered <- res$data.filter
  } else {
    all_keep_otu <- c()

    for (disease_name in diseases) {

      # Extract sample IDs belonging to this disease group (using metadata row names)
      disease_sids <- rownames(filtered_meta)[filtered_meta$disease == disease_name]
      # Subset the abundance table to those samples only
      disease_df <- labelled_abund[disease_sids, , drop = FALSE]

      # Run low-count removal on the disease-specific subset
      res <- low_count_removal(disease_df, percent)

      # Accumulate OTUs that survive in ANY disease group (union)
      all_keep_otu <- union(all_keep_otu, res$keep.otu)
    }

    # Retain the union of OTUs across all groups from the full abundance table
    abundance_filtered <- labelled_abund[, all_keep_otu, drop = FALSE]

    # Coerce each column to numeric (guard against list columns from bind_rows edge cases)
    # check.names=FALSE preserves taxa names with spaces/special characters
    abundance_filtered <- as.data.frame(lapply(abundance_filtered, function(col) {
      if (is.list(col)) col <- unlist(col)
      as.numeric(col)
    }), row.names = rownames(abundance_filtered), check.names = FALSE)
  }
  #View(abundance_filtered)
  return(abundance_filtered)
}



# low_count_removal(data, percent)
# ---------------------------------
# Removes taxa (columns) whose summed relative abundance across all samples
# falls below `percent`% of the total. This removes very rare taxa that are
# unlikely to be biologically meaningful.
#
# Parameters:
#   data    – data frame or matrix; rows = samples, columns = OTUs (must be numeric)
#   percent – minimum % of total abundance an OTU must have to be kept (default 0.01)
#
# Returns: list with
#   $data.filter – filtered data frame
#   $keep.otu    – named integer vector of surviving column indices
low_count_removal <- function(
    data,    # OTU count dataframe of size n (samples) x p (OTUs)
    percent = 0.01  # cutoff percentage threshold
) {

  # Validate that input is a data frame or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a dataframe or matrix of counts.", call. = FALSE)
  }

  if (!all(sapply(data, is.numeric))) stop("All columns must be numeric.")


  # Identify OTUs whose relative contribution exceeds the percent threshold
  # colSums(data) * 100 / sum(colSums(data)) = per-OTU % of total abundance
  keep.otu <- which(
    colSums(data) * 100 / sum(colSums(data)) > percent
  )

  # Subset the abundance table to surviving OTUs
  data.filter <- data[, keep.otu, drop = FALSE]

  # Return both the filtered data and the indices of kept OTUs
  return(list(
    data.filter = data.filter,
    keep.otu = keep.otu
  ))
}




## ============================================================
## 5. CLR TRANSFORMATION
## ============================================================

# clr.transformation(abundance.filtered, offset_value)
# -----------------------------------------------------
# Applies the Centered Log-Ratio (CLR) transformation to an abundance matrix.
# CLR converts compositional data to real-valued data suitable for multivariate
# methods. An offset is added before log to handle zero values.
#
# Parameters:
#   abundance.filtered – numeric matrix or data frame; rows = samples, columns = taxa
#   offset_value       – small positive constant added to all values before CLR
#                        (e.g. 1 for count data; second-smallest for relative abundances)
#
# Returns: data frame of CLR-transformed values (same dimensions as input)
clr.transformation <- function(abundance.filtered, offset_value) {
  # logratio.transfo (mixOmics) requires a matrix input; offset handles zeros
  data.clr <- logratio.transfo(
    as.matrix(abundance.filtered),
    logratio = "CLR",
    offset = offset_value
  )

  # unclass() strips the S3 class wrapper returned by logratio.transfo
  clr_mat <- unclass(data.clr)
  # check.names=FALSE preserves taxa names with spaces/special characters
  clr_df <- as.data.frame(clr_mat, check.names = FALSE)

  return(clr_df)
}




## ============================================================
## 6. SINGLE DATASET PREPROCESSING
## ============================================================

# preprocess_single(meta, abund, diseases, keep_cols, percent, offset_value, CLR)
# -------------------------------------------------------------------------------
# Full preprocessing pipeline for ONE study dataset:
#   metadata filtering → abundance alignment → prevalence filtering → CLR transformation
#
# Parameters:
#   meta         – metadata data frame for this study (rows = samples)
#   abund        – abundance data frame for this study (rows = samples, cols = taxa)
#   diseases     – disease labels to retain (passed to filter_meta)
#   keep_cols    – extra metadata columns to retain (passed to filter_meta)
#   percent      – low-count removal threshold (passed to split_count_filtering)
#   offset_value – CLR offset for zeros (passed to clr.transformation)
#   CLR          – logical; if FALSE returns filtered abundance without CLR (default TRUE)
#
# Returns (when CLR=TRUE): list with
#   $data – CLR-transformed abundance data frame (samples × taxa)
#   $meta – filtered metadata data frame (samples × requested columns)
preprocess_single <- function(meta, abund, diseases, keep_cols, percent, offset_value, CLR=TRUE) {
  # 1. Filter metadata: keep only samples in target disease groups + requested columns
  filtered_meta <- filter_meta(meta, diseases, keep_cols)
  #View(filtered_meta)

  # 2. Align abundance rows to the filtered metadata sample list
  #    (removes samples that were excluded during metadata filtering)
  selected_abund <- as.data.frame(abund[rownames(abund) %in% rownames(filtered_meta), ], check.names = FALSE)
  #View(selected_abund)

  # 3. Apply per-disease prevalence filtering (union of surviving taxa across disease groups)
  abundance_filtered <- split_count_filtering(selected_abund, diseases, percent, filtered_meta)

  # Early return: skip CLR if caller only wants filtered counts (e.g. HGMA preprocessing check)
  if (!CLR) return(abundance_filtered)

  # 4. (Step 4 intentionally numbered as 5/6 in original; CLR=TRUE path continues below)
  # Convert to matrix for logratio.transfo (requires matrix, not data.frame)
  abundance_matrix <- as.matrix(abundance_filtered)

  # 5. Apply CLR transformation (adds offset to handle zeros, then computes log-ratios)
  clr_transformed <- clr.transformation(abundance_matrix, offset_value)
  #View(clr_transformed)

  # 6. Return as data frame with taxa column names preserved (check.names=FALSE)
  clr_transformeddf <- as.data.frame(clr_transformed, check.names = FALSE)
  return(list(
    data = clr_transformeddf,  # CLR-transformed abundance; rows = samples, cols = taxa
    meta = filtered_meta       # filtered metadata aligned to the same samples
  ))

}




## ============================================================
## 7. MULTIPLE DATASET PREPROCESSING
## ============================================================

# preprocess_multiple(abundance_list, meta_list, diseases, keep_cols, percent, offset_value, threshold)
# -----------------------------------------------------------------------------------------------------
# Runs preprocess_single on each dataset in abundance_list / meta_list, then calls
# integrate_filter to find taxa present across at least `threshold` studies.
#
# Parameters:
#   abundance_list – named list of abundance data frames; one per study
#   meta_list      – named list of metadata data frames; must have the same names as abundance_list
#   diseases       – disease labels to retain
#   keep_cols      – extra metadata columns to retain
#   percent        – low-count removal threshold
#   offset_value   – CLR offset
#   threshold      – minimum number of studies a taxon must appear in to be kept (NULL = keep all)
#
# Returns: named nested list (same structure as integrate_filter output):
#   result[[study_name]]$abund – CLR-transformed abundance (taxa filtered to cross-study survivors)
#   result[[study_name]]$meta  – filtered metadata
preprocess_multiple <- function(abundance_list, meta_list, diseases, keep_cols, percent, offset_value, threshold = NULL) {
  results <- list()

  # Preprocess each study independently
  for (nm in names(abundance_list)) {
    abund <- abundance_list[[nm]]
    meta  <- meta_list[[nm]]
    # Run the single-study pipeline: filter → prevalence filter → CLR
    result <- preprocess_single(meta, abund, diseases, keep_cols, percent, offset_value)

    # Store as a named entry in results; accessed by study name downstream
    results[[nm]] <- list(
      data = result$data,  # CLR-transformed abundance for this study
      meta = result$meta   # filtered metadata for this study
    )
    cat("\n--- Dataset:", nm, "---\n")
    print(dim(result$data))
    # View(result$data)
    # View(result$meta)
  }

  # Intersect taxa across all studies (keeps only taxa present in >= threshold studies)
  combined <- integrate_filter(results, threshold)
  return(combined)
}




## ============================================================
## 8. CROSS-STUDY TAXA INTEGRATION
## ============================================================

# integrate_filter(results, threshold)
# -------------------------------------
# After per-study preprocessing, this function enforces cross-study taxon consistency.
# It keeps only taxa that appear (have non-zero total abundance) in at least `threshold`
# studies. This is critical so that all studies share a common feature space before
# integration in MINT.
#
# Parameters:
#   results   – named list from preprocess_multiple; each element has $data and $meta
#   threshold – minimum number of studies a taxon must appear in (integer, e.g. 1 = keep all;
#               setting equal to length(results) = keep only universal taxa)
#
# Returns: named nested list:
#   final[[study_name]]$abund – abundance filtered to surviving taxa
#   final[[study_name]]$meta  – metadata aligned to the same samples
integrate_filter <- function(results, threshold) {

  # Step 1: For each study, identify taxa with any non-zero abundance
  taxa_present_list <- purrr::map(results, function(study) {
    taxa_cols <- colnames(study$data)
    # A taxon "is present" if its column sum > 0 across all samples in this study
    present_taxa <- taxa_cols[
      purrr::map_lgl(taxa_cols, ~ sum(study$data[[.x]], na.rm = TRUE) > 0)
    ]
    return(present_taxa)
  })

  # Count how many studies each taxon appears in
  taxa_counts <- table(unlist(taxa_present_list))
  # Keep only taxa that meet the threshold (present in at least `threshold` studies)
  surviving_taxa <- names(taxa_counts)[taxa_counts >= threshold]

  # Step 2: Apply the surviving-taxa filter to each study independently
  comb_data_list <- list()
  meta_list <- list()


  # imap provides both the study object and its name; result is automatically named
  final <- purrr::imap(results, function(study, study_name) {

    # Subset abundance columns to surviving taxa (any_of handles taxa missing in this study)
    abund_filtered <- study$data %>%
      dplyr::select(any_of(surviving_taxa))

    # Re-attach row names (sample IDs) — dplyr select drops them
    rownames(abund_filtered) <- rownames(study$data)

    # Align metadata rows to the same samples as the filtered abundance
    meta_filtered <- study$meta[rownames(abund_filtered), , drop = FALSE]

    # Return as abund/meta pair (renamed from data/meta for consistency downstream)
    list(
      abund = abund_filtered,  # filtered CLR abundance; rows = samples, cols = surviving taxa
      meta  = meta_filtered    # filtered metadata; rows = same samples
    )
  })


  return(final)
}




## ============================================================
## 9. PREPROCESSING ENTRY POINT
## ============================================================

# run_preprocessing(abundance_list, meta_list, diseases, keep_cols, percent, offset_value, threshold)
# ----------------------------------------------------------------------------------------------------
# Top-level dispatcher: routes to preprocess_single (one dataset) or
# preprocess_multiple (many datasets) based on the size of abundance_list.
# This is the function called in the .Rmd notebooks.
#
# Parameters: same as preprocess_multiple (threshold only used for multiple datasets)
#
# Returns:
#   - If one dataset:  list(data = CLR_df, meta = filtered_meta)
#   - If >1 datasets:  named nested list from integrate_filter
run_preprocessing <- function(abundance_list, meta_list, diseases,
                              keep_cols, percent, offset_value, threshold = NULL) {
  # Check how many datasets are in the abundance list
  if (length(abundance_list) == 1) {
    print("single")
    # Only one dataset → run preprocess_single directly
    abund <-  as.data.frame(abundance_list[[1]], check.names = FALSE)
    meta  <- as.data.frame(meta_list[[1]], check.names = FALSE)
    result <- preprocess_single(meta, abund, diseases, keep_cols, percent, offset_value)

  } else {
    print("Multiple")
    # Multiple datasets → run preprocess_multiple (includes integrate_filter internally)
    result <- preprocess_multiple(abundance_list, meta_list,
                                  diseases, keep_cols, percent, offset_value, threshold)
  }
  return(result)
}




## ============================================================
## 10. LAYER 2 — INDIVIDUAL FILTERING
## ============================================================

# individual_filtering(data, err_thres, ncomp, validation, folds, nrepeat, dist)
# -------------------------------------------------------------------------------
# Layer 2 of the pipeline. For each study, fits a single-study PLS-DA model
# and checks whether any per-class error rate at `ncomp` exceeds `err_thres`.
# Studies that FAIL the check (i.e. any class BER > threshold) are dropped.
#
# Rationale: removes studies that are too noisy for integrative analysis,
# preventing them from harming the MINT model's discriminative ability.
#
# Parameters:
#   data      – output from run_preprocessing: named list of study objects (abund + meta),
#               OR a single study object with $abund and $meta
#   err_thres – maximum tolerated per-class error rate (default 0.55; thesis specifies 0.7)
#               NOTE: currently called with 0.5 in sandbox — reconcile before full_framework.Rmd
#   ncomp     – number of PLS-DA components to fit and evaluate (default 2)
#   validation– CV strategy (default "Mfold")
#   folds     – number of CV folds (default 5)
#   nrepeat   – number of CV repeats (default 10)
#   dist      – classification rule (default "centroids.dist")
#
# Returns:
#   - Single study input: logical TRUE/FALSE (keep or drop)
#   - Nested list input:  named list of study objects that passed the filter
individual_filtering <- function(data, err_thres = 0.55, ncomp = 2,
                                 validation = "Mfold", folds = 5, nrepeat = 10, dist = "centroids.dist") {

  # Inner helper: evaluates one study object and returns TRUE if it passes
  ## parameter study here means a single dataset, along with its abundance matrix and metadata
  check_single_dataset <- function(study) {
    abund <- study$abund   # CLR-transformed abundance matrix for this study
    meta  <- study$meta    # metadata with disease labels
    Y <- factor(meta$disease)  # disease labels as factor for PLS-DA

    # Fit a PLS-DA model on this single study
    plsda.model <- plsda(abund, Y, ncomp = ncomp)
    # Evaluate with repeated k-fold CV
    perf.model <- perf(
      object = plsda.model,
      validation = validation,
      folds = folds,
      progressBar = FALSE,
      dist = dist,
      nrepeat = nrepeat
    )

    # Extract per-class error rates at the final component
    # error.rate.class[[dist]] is a matrix: rows = components, cols = classes
    #class_errors <- perf.model$error.rate.class[[dist]]
    #exceed_threshold <- apply(class_errors, 1, function(x) any(x > err_thres))
    class_errors <- perf.model$error.rate.class[[dist]][, ncomp]
    # Study fails if ANY class error rate exceeds the threshold at ncomp
    exceed_threshold <- any(class_errors > err_thres)

    # keep_dataset = TRUE means the study passed (no class error exceeded threshold)
    keep_dataset <- !any(exceed_threshold)
    return(keep_dataset)
  }

  # Case 1: caller passed a single study directly (has both $abund and $meta)
  if (is.list(data) && all(c("abund","meta") %in% names(data))) {
    return(check_single_dataset(data))
  }

  # Case 2: caller passed a named list of multiple studies
  if (is.list(data) && !all(c("abund","meta") %in% names(data))) {
    results <- purrr::imap(data, function(study, study_name) {
      if (!all(c("abund","meta") %in% names(study))) {
        stop(paste("Study", study_name, "is missing abund or meta"))
      }
      keep <- check_single_dataset(study)
      # Return study if it passed; NULL if it should be dropped
      if (keep) return(study) else return(NULL)
    })
    # Remove NULL entries (dropped studies)
    results <- results[!purrr::map_lgl(results, is.null)]
    return(results)
  }

  stop("Input must be a single study (abund + meta) or a nested list of studies")
}




## ============================================================
## 11. BER HELPER
## ============================================================

## Calculate the BER for single study, with 5 folds for CV and 10 times repeats
# output the BER value

# calculate_BER(single_study, ncomp, validation, folds, nrepeat, dist)
# ---------------------------------------------------------------------
# Fits a PLS-DA model on a single study and returns the Balanced Error Rate (BER)
# at component `ncomp`. Used internally by forward_selection to rank studies.
#
# BER = (1/K) * Σ(n_k_misclassified / n_k)  where K = number of classes.
# BER is robust to class imbalance unlike standard error rate.
#
# Parameters:
#   single_study – list with $abund (abundance matrix) and $meta (metadata with $disease)
#   ncomp        – component at which to read BER (default 2)
#   (other params passed directly to perf())
#
# Returns: numeric vector of BER values at ncomp (one value per CV repeat, usually length 1 after averaging)
calculate_BER = function(single_study, ncomp=2, validation="Mfold",
                         folds=5, nrepeat=10, dist = "centroids.dist") {

  abund = single_study$abund  # CLR-transformed abundance for this study
  meta = single_study$meta    # metadata with disease labels
  Y = factor(meta$disease)    # disease vector as factor

  # Fit PLS-DA and evaluate with repeated k-fold CV
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
  # BER matrix rows = components, cols = CV repeats; return the row for ncomp
  return(perf.model$error.rate$BER[ncomp, ])
}




## ============================================================
## 12. SINGLE-STUDY TAXA REMOVAL HELPER
## ============================================================

## Eliminate the taxa that only exist in single study
# study_col: let user specify the column name for study name; n_taxa. eliminate the taxa that only exist in n studies across all studies

# eliminate_single_taxa(lst_df, study_col, n_taxa)
# ------------------------------------------------
# Removes taxa (columns) that appear in only `n_taxa` studies (default: 1).
# This is important for MINT because taxa unique to one study cannot be
# integrated reliably across the cohort.
#
# Parameters:
#   lst_df    – data frame produced by bind_rows(result); must have an $abund sub-df
#               and $meta with a study column
#   study_col – column name in metadata that holds the study label (default "study")
#   n_taxa    – a taxon is "single-study" if its number of unique values within a study
#               equals n_taxa; these are removed (default 1 = constant within one study)
#
# Returns: filtered abundance data frame including the study column
eliminate_single_taxa = function(lst_df, study_col="study", n_taxa=1) {

  dff = data.frame()
  df = lst_df$abund
  # Attach study label from metadata to the abundance data frame
  df$study = lst_df$meta[[study_col]]
  # Taxa columns = all columns except the study label column
  taxa_lst = setdiff(colnames(df), c(study_col))

  for (species in taxa_lst) {
    # Build a row for this species: record "Exist" or "x" for each study
    df_taxa <- data.frame(species = species, stringsAsFactors = FALSE)

    # Check presence/absence of this species in each study
    for (study in unique(df$study)) {
      # Get all values of this species within this study
      study_data <- df[df$study == study, species]

      # "x" = effectively absent (only n_taxa unique value across study — usually all zeros)
      # "Exist" = has more than n_taxa unique values (i.e. some variation / presence)
      if (length(unique(study_data)) == n_taxa) {
        df_taxa[[study]] <- "x"
      } else {
        df_taxa[[study]] <- "Exist"
      }
    }

    # Accumulate species rows into the full summary data frame
    dff <- rbind(dff, df_taxa)
  }

  # Convert "Exist"/"x" to 1/0 for counting
  dff1 = dff
  dff1[, unique(df$study)] <- lapply(dff1[, unique(df$study)], function(x) ifelse(x == "Exist", 1, 0))
  # Exist_count = number of studies this taxon appears in
  dff$Exist_count <- rowSums(dff1[, unique(df$study)])
  # Identify taxa appearing in only 1 study
  one_taxa <- dff[dff$Exist_count == 1, "species"]
  # Return the abundance data frame without those single-study taxa
  return(df[, !(colnames(df) %in% one_taxa)])
}




## ============================================================
## 13. MULTI-STUDY BER HELPER
## ============================================================

## Calculate BER for integrative dataset including multiple studies, eliminate_single_taxa function is embedded.

# calculate_multipleDB_BER(lst, ncomp, study_col, ...)
# -----------------------------------------------------
# Computes the Balanced Error Rate for a combination of multiple studies
# using MINT sPLS-DA. Called iteratively by forward_selection to evaluate
# candidate study combinations.
#
# Parameters:
#   lst       – named list of study objects (abund + meta); the candidate combination
#   ncomp     – number of MINT components (default 2)
#   study_col – metadata column holding study label (default "study")
#   (other params forwarded to perf)
#   n_taxa    – passed to eliminate_single_taxa (default 1)
#
# Returns: numeric BER value at ncomp (scalar) for this study combination
calculate_multipleDB_BER = function(lst, ncomp=2, study_col="study", validation="Mfold",
                                    folds=5, nrepeat=10, dist="centroids.dist", n_taxa=1) {

  ## Add study column to each study's metadata so bind_rows can label samples
  for (name in names(lst)) {
    # Extract short study name (strip trailing ".relative_abundance" suffix)
    lst[[name]]$meta$study = strsplit(name, '\\.')[[1]][1]
    # lst[[name]]$meta$study = name
  }
  # Combine all studies into a single long data frame (abund rows stacked)
  lst_df = bind_rows(lst)

  # Remove taxa that appear in only one study (uninformative for integration)
  comb_df = eliminate_single_taxa(lst_df, study_col = study_col, n_taxa=n_taxa)
  # Separate the abundance matrix from the study label column
  abund_X = comb_df[ ,setdiff(names(comb_df), study_col), drop = FALSE]
  Y = lst_df$meta$disease             # disease outcome vector (all studies combined)
  corresponding_study = comb_df$study  # study membership vector for MINT

  # Fit MINT sPLS-DA on the combined data
  model = mint.splsda(X=abund_X, Y=Y, study=corresponding_study, ncomp=ncomp)
  perf.model = perf(model, validation=validation, folds=folds, dist=dist)

  # return the global overall BER (averaged across studies via leave-one-study-out)
  return(perf.model$global.error$BER[ncomp, ])
}




## ============================================================
## 14. LAYER 3 — FORWARD SELECTION
## ============================================================

## Forward Selection
# data: study pool

# forward_selection(data, ncomp, study_col, validation, folds, nrepeat, dist, n_taxa)
# ------------------------------------------------------------------------------------
# Layer 3 of the pipeline. Greedy algorithm that selects an optimal SUBSET of
# studies to integrate into the MINT model, maximising discriminative performance
# (minimising BER).
#
# Algorithm:
#   Round 0  – Seed: find the single study with the lowest individual BER.
#   Round i+ – Greedy expansion: test each remaining study added to the current
#              combination; keep the one that reduces multi-study BER the most.
#              Stop if adding ANY study fails to improve BER.
#
# Parameters:
#   data      – output from run_preprocessing (or individual_filtering); named list of study objects
#   ncomp     – number of MINT components (default 2)
#   study_col – metadata column for study label (default "study")
#   nrepeat   – CV repeats for calculate_BER / calculate_multipleDB_BER (default 20)
#   dist      – classification rule (default "centroids.dist")
#   n_taxa    – taxa-prevalence threshold passed to eliminate_single_taxa (default 1)
#
# Returns: named list of selected study objects (abund + meta) — the optimal combination
forward_selection = function(data, ncomp=2, study_col="study", validation='Mfold', folds=5, nrepeat=20,
                             dist="centroids.dist", n_taxa=1) {
  remain_studies = list()   # studies included in the current combination
  i = 0  # round counter

  ## ---- ROUND 0: seed selection ----
  ## Find the single study with the lowest individual BER to start the combination
  cat("########## ROUND ", i+1, " ##########", "\n")
  cat("Total study: ", length(data), "\n")
  cat("Processing... Testing studies: ", names(data), "\n")
  first_study = ""
  ERR = 1  # initialise to worst possible BER (100%)
  for (study in names(data)) {
    study_ERR = calculate_BER(data[[study]])
    # Update seed study if this study has a better (lower) BER
    if ( study_ERR < ERR ) {
      first_study = study
      ERR = study_ERR
    }
  }
  # Add the best single study to the combination permanently
  cat("Introduce ", first_study, " into combination.", "\n")
  remain_studies[[first_study]] = data[[first_study]]   # plugin the first study into combination permanently
  cat("########## Current combination:", names(remain_studies), " round ", i+1, "\n")
  data[[first_study]] = NULL     # remove seed study from the available pool
  for (j in 1:2) {
    cat("                                     ", "\n")
  }


  ## ---- ROUNDS 1+: greedy expansion ----
  dynamic_BER = 1  # best (lowest) BER achieved so far across all rounds
  current_BER = 1  # BER of the current candidate combination
  # Continue while there are studies left to try and the last addition improved BER
  while ( (length(data) > 0) | (current_BER > dynamic_BER) ) {
    i = i + 1
    cat("########## ROUND ", i+1, " ##########", "\n")

    combin_err = 1     # best BER found in this round
    selected_study = ""  # candidate study to add in this round

    for (study in names(data)) {
      ## Temporarily add this study to the combination and evaluate BER
      remain_studies[[study]] = data[[study]]; # cat("adding ", study, " in combination", "\n")
      err = calculate_multipleDB_BER(remain_studies, ncomp=ncomp, study_col=study_col,
                                     validation=validation, dist=dist, n_taxa=n_taxa)
      cat( "Processing... Test new combination: ", names(remain_studies), "\n" )

      # Track the study that gives the lowest BER when added
      if (err < combin_err) {
        combin_err = err
        selected_study = study
      }
      # Remove temporary study to test the next one
      remain_studies[[study]] = NULL; #  cat("removing ", study, " out of combination", "\n")
    }

    # After scanning all remaining studies, add the best candidate permanently
    #cat("grab studies ", "(", selected_study, ")", " into round ", i , "\n")
    #cat("Its BER is ", combin_err, "\n")
    remain_studies[[selected_study]] = data[[selected_study]]  # formally add best study
    data[[selected_study]] = NULL    # remove it from the pool
    current_BER = combin_err  # BER of the new combination

    if (current_BER < dynamic_BER) {
      # Adding selected_study improved BER → accept and continue
      #cat("current_BER: ", current_BER, "\n")
      #cat("dynamic_BER: ", dynamic_BER, "\n")
      #cat("current BER lower than dynamic global BER, assign it as new dynamic global BER", "\n")
      dynamic_BER = current_BER  # update the running best BER
      cat("Introduce ", selected_study, " into combination.", "\n")
      cat("########## Current combination:", names(remain_studies), " round ", i+1, "\n")
    } else {
      # Adding the best available study FAILED to improve BER → stop selection
      cat("No improvement in BER after introducing new dataset into combination", "\n")
      cat("########## Current combination:", names(remain_studies), " round ", i+1, "\n")
      # Remove the last-added study (it didn't help)
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




## ============================================================
## 15. HGMA TEST SET BUILDER
## ============================================================

# draw_TestSet(threshold_for_LowCountsRm, disease_lst, dataset_study)
# -------------------------------------------------------------------
# Loads and preprocesses ONE study from the Human Gut Microbiome Atlas (HGMA)
# to serve as an external test set. HGMA uses numeric taxon IDs that must be
# mapped to species names via corresponding_taxa.csv, and ambiguous multi-alias
# taxon names (containing "/") are resolved by aggregation or renaming.
#
# Required files (must be present in the working directory):
#   sample_metadata.xlsx    – HGMA sample metadata; columns include BioProject, sample.ID, Disease
#   vect_atlas.csv/vect_atlas.csv – HGMA abundance matrix; rows = taxa (indexed by numeric ID)
#   corresponding_taxa.csv  – lookup table mapping HGMA numeric IDs to species names
#
# Parameters:
#   threshold_for_LowCountsRm – low-count removal threshold (default 0.01; unused directly,
#                                used implicitly via preprocess_single)
#   disease_lst               – disease labels to retain in HGMA (e.g. c("T2D","Healthy","NGT","healthy"))
#   dataset_study             – BioProject ID to select from HGMA (e.g. "PRJEB1786" or "PRJNA422434")
#
# Returns: list with
#   $abund – CLR-transformed abundance data frame (rows = samples, cols = taxa)
#   $meta  – metadata data frame with $disease and $study columns
draw_TestSet = function(threshold_for_LowCountsRm=0.01, disease_lst=c("T2D", "Healthy", "NGT", "healthy"), dataset_study) {
  # --- Load raw HGMA files ---
  print("Reading HGMA rawdata...")
  sample_md = read_excel("sample_metadata.xlsx")        # sample-level metadata
  abund_data = read.csv("vect_atlas.csv\\vect_atlas.csv") # abundance matrix (rows = taxa IDs)
  corr_taxa = read.csv("corresponding_taxa.csv")          # ID-to-name mapping table

  ## --- Pre-process the taxa name lookup table ---
  print("Processing metadata and rawdata...")
  # Remove rows where the species name contains "unclassified" (ambiguous taxa)
  corr_taxa <- corr_taxa %>%
    filter(!str_detect(name, "unclassified"))

  # Multi-taxa entries: rows where name ends in " 1", " 2", etc. (subtypes of same species)
  multi_taxa <- corr_taxa %>%
    filter(str_detect(name, "\\s\\d$"))

  # Strip trailing digit + space from multi-taxa names to get the base species name
  multi_taxa <- multi_taxa %>%
    mutate(name = str_sub(name, 1, -3))

  multi_taxa_lst <- unique(multi_taxa$name)  # unique base names with multiple subtypes

  ## --- Filter metadata and abundance to the requested study + diseases ---
  # Select samples from the target BioProject
  metadata = sample_md[sample_md['BioProject'] == dataset_study, ]
  # Keep only samples with the requested disease labels
  metadata = metadata %>%
    filter(Disease %in% disease_lst)
  sp_id = metadata$sample.ID  # sample IDs to extract from abundance matrix
  # Subset abundance matrix to matching samples + keep the row ID column ("X")
  rawdata = abund_data[, c(sp_id, "X")]

  # Map numeric taxon IDs ("X" column) to species names via corresponding_taxa.csv
  rawdata <- merge(
    rawdata,
    corr_taxa,
    by.x = "X",
    by.y = "id",
    all.x = TRUE  # all.x=TRUE: Keep all rows from x even if no matching key in y (left join)
  )
  rawdata <- rawdata[, !(names(rawdata) %in% c("X", "id"))]  # drop numeric ID columns
  rawdata <- rawdata[!is.na(rawdata$name), ]  # drop rows with no name match

  ## --- Transpose: make rows = samples, cols = taxa ---
  taxa_lst = rawdata$name  # taxon names (will become column headers)
  rawdata_v1 = data.frame(t(rawdata[, !(names(rawdata) %in% c("name"))]))
  colnames(rawdata_v1) = taxa_lst  # assign taxon names as column headers
  names(rawdata_v1)[names(rawdata_v1) == "index"] = "sample_id"
  rawdata_v1$sample_id = rownames(rawdata_v1)  # sample IDs as a column

  # Attach disease label from metadata
  rawdata_v1 = merge(rawdata_v1,
                     metadata[, c("sample.ID", "Disease")],
                     by.x="sample_id",
                     by.y="sample.ID",
                     all.x=TRUE)
  # Fix one known alias: two equivalent species names merged to the preferred name
  names(rawdata_v1)[names(rawdata_v1) == "Blautia coccoides == Blautia producta"] = "Blautia producta"

  ## --- Resolve multi-subtype taxa (aggregation logic) ---
  # Helper: check whether any alias in a "/" list lacks a trailing digit
  # (i.e. a "main" species entry exists, not just subtypes)
  check_main_species <- function(lst) {
    integers <- as.character(1:9)
    for (sub in lst) {
      last_char <- substr(sub, nchar(sub), nchar(sub))
      if (!(last_char %in% integers)) {
        return(TRUE)  # found a column without a trailing digit → main species present
      }
    }
    return(FALSE)
  }

  # For each multi-subtype base species, consolidate the subtype columns
  for (taxa in multi_taxa_lst) {
    if ( !grepl("subtype", tolower(taxa), fixed = T) ) {
      cols = colnames(rawdata_v1)
      rawdata_v1$taxa = 0
      # Find all columns that match this base species (includes subtypes)
      alltaxa_columns <- cols[grepl(taxa, cols, fixed = T)]

      if ( length(alltaxa_columns) == 1 ) {
        # Only one column found → directly assign it as the canonical species column
        rawdata_v1[[taxa]] <- rawdata_v1[[alltaxa_columns[1]]]
      } else if ( check_main_species(alltaxa_columns) == T ) {
        # A main-species column exists → remove the subtype columns, keep the main
        alltaxa_columns <- setdiff(alltaxa_columns, taxa)
        rawdata_v1 <- rawdata_v1[, !(names(rawdata_v1) %in% alltaxa_columns), drop = FALSE]

      } else {
        # All columns are subtypes → sum them into one combined column, drop subtypes
        rawdata_v1[[taxa]] <- rowSums(rawdata_v1[, alltaxa_columns, drop = FALSE], na.rm = TRUE)
        rawdata_v1 <- rawdata_v1[, !(names(rawdata_v1) %in% alltaxa_columns), drop = FALSE]
      }
    }
  }

  ## --- Rename ambiguous "/" alias columns to the canonical CMD-compatible name ---
  # These mappings resolve HGMA taxa names with "/" (multiple possible identities)
  # to the single name used in the CMD training set
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

  # Apply rename_map to matching column names
  cols <- colnames(rawdata_v1)
  hits <- cols %in% names(rename_map)
  cols[hits] <- unname(rename_map[cols[hits]])
  colnames(rawdata_v1) <- cols

  ## --- Final packaging: separate abundance from metadata, then CLR-transform ---
  print("Organizing the metadata and the rawdata into proper output form")
  # Split into abundance (all taxa cols) and metadata (sample_id + disease) sub-lists
  hgma_lst = list(abund=subset(rawdata_v1, select=-c(sample_id, disease)), meta=subset(rawdata_v1, select=c(sample_id, disease)) )

  # Compute the second-smallest non-zero value as the CLR offset
  # (offset = second-smallest gives a conservative near-zero value specific to this dataset)
  x <- unlist(hgma_lst$abund)
  x <- x[x > 0]
  second_smallest <- sort(unique(x), na.last = NA)[2]

  # Run full preprocessing (filter_meta → split_count_filtering → CLR)
  processed_lst = preprocess_single(meta=hgma_lst$meta, abund=hgma_lst$abund, diseases=disease_lst,
                                    keep_cols=c("sample_id"), offset_value=second_smallest, percent=0.01, CLR=T )
  # Rename $data → $abund for consistency with the rest of the pipeline
  names(processed_lst)[names(processed_lst) == "data"] = "abund"
  # Add study label to metadata (used by align_columns and MINT)
  processed_lst$meta$study = dataset_study

  return(processed_lst)
}




## ============================================================
## 16. COLUMN ALIGNMENT: TRAINING vs TESTING
## ============================================================

## Feature alignment: CMD (training) vs HGMA (testing) column names
# Promoted from Test_Framework.Rmd — Edward Kao, 2025
# Strips 'species:' prefix from both matrices, then resolves HGMA '/' alias columns
# against CMD training column names. Returns aligned training and testing matrices
# containing only their common taxa, plus a change_log of all renames performed.

# align_columns(p_test, p_selected)
# ----------------------------------
# Reconciles the taxon column names between the CMD training set and the
# HGMA testing set so they can be compared directly.
#
# Problem: CMD taxa names have a "species: " prefix (from curatedMetagenomicData),
# while HGMA taxa may use "/" to list multiple equivalent species names.
#
# Steps:
#   1. Strip "species: " prefix from both matrices.
#   2. For each HGMA "/" alias column, find matching training column(s):
#      - Single match  → rename the test column to the training name.
#      - Multiple match → sum those training columns into a new combined column,
#                         then rename the test column to match.
#      - No match       → leave unchanged (will be excluded at intersection step).
#   3. Subset both matrices to their intersection of column names.
#
# Parameters:
#   p_test     – HGMA test abundance matrix (rows = test samples, cols = taxa)
#   p_selected – CMD training abundance matrix (rows = training samples, cols = taxa)
#
# Returns: list with
#   $training  – training matrix subset to common taxa
#   $testing   – testing matrix subset to common taxa
#   $change_log – named list recording every rename/merge performed
#   $messages   – single string summarising all changes (for logging/printing)
align_columns <- function(p_test, p_selected) {

  # Step 1: Remove 'species:' prefix from column names in both matrices
  # This prefix is added by curatedMetagenomicData but is absent in HGMA
  colnames(p_selected) <- sub("^species:\\s*", "", colnames(p_selected), ignore.case = TRUE)
  colnames(p_test)     <- sub("^species:\\s*", "", colnames(p_test),     ignore.case = TRUE)

  change_log  <- list()    # records each rename/merge operation
  message_log <- character()  # human-readable summary lines

  # Step 2: Iterate test columns and resolve "/" alias names
  for (i in seq_along(colnames(p_test))) {

    col_test <- colnames(p_test)[i]

    # Only process test columns that contain '/' (i.e. HGMA alias entries)
    if (grepl("/", col_test, fixed = TRUE)) {

      # Split the "/" alias into individual candidate species names
      taxa_test         <- trimws(unlist(strsplit(col_test, "/", fixed = TRUE)))
      matched_train_cols <- c()

      # Search training columns for any that are a subset of the test aliases
      for (col_train in colnames(p_selected)) {
        taxa_train <- trimws(unlist(strsplit(col_train, "/", fixed = TRUE)))
        # A training column matches if ALL its taxa are present in the test alias list
        if (all(taxa_train %in% taxa_test)) {
          matched_train_cols <- c(matched_train_cols, col_train)
        }
      }

      # Case A: exactly one training column matched → simple rename
      if (length(matched_train_cols) == 1) {
        old_name            <- col_test
        colnames(p_test)[i] <- matched_train_cols[1]
        change_log[[old_name]] <- paste0("Single match: renamed test column to '", matched_train_cols[1], "'")
        message_log <- c(message_log,
          paste0("Single match found: test column '", old_name, "' renamed to '", matched_train_cols[1], "'"))

      # Case B: multiple training columns matched → sum them into a new combined column
      } else if (length(matched_train_cols) > 1) {
        new_train_col <- paste(matched_train_cols, collapse = "/")

        # Create the combined column in training if it doesn't already exist
        if (!(new_train_col %in% colnames(p_selected))) {
          p_selected[[new_train_col]] <- rowSums(p_selected[, matched_train_cols, drop = FALSE])
          change_log[[paste0("p_selected: ", new_train_col)]] <- paste0(
            "Created new combined column from: ", paste(matched_train_cols, collapse = ", "))
        }

        old_name            <- col_test
        colnames(p_test)[i] <- new_train_col
        change_log[[old_name]] <- paste0("Multiple matches: renamed test column to '", new_train_col, "'")
      }
      # Case C: no match → column name left unchanged; excluded at intersection step
    }
  }

  # Step 3: Subset both matrices to their common taxa only
  # taxa not shared between training and testing are dropped from both
  common_cols <- intersect(colnames(p_selected), colnames(p_test))
  training    <- p_selected[, common_cols, drop = FALSE]
  testing     <- p_test[,     common_cols, drop = FALSE]

  return(list(
    training   = training,    # training abundance restricted to common taxa
    testing    = testing,     # testing abundance restricted to common taxa
    change_log = change_log,  # list of all rename/merge operations
    messages   = paste(message_log, collapse = "\n")  # printable summary
  ))
}




## ============================================================
## 17. LAYER 4/5 — MINT sPLS-DA MODELLING
## ============================================================

## MINT sPLS-DA modelling: sequential hyperparameter tuning + model fitting + evaluation
# Promoted from Joyce/modelling.Rmd — Joyce Hu, 2025
# Tunes comp 1 first, then comp 2 using already.tested.X (mixOmics best-practice).
# study must exist in the calling environment as a factor vector (one entry per sample).
# Returns: model, optimal.keepX, performance, selected features per component, all unique features.

# run_mint(X, Y, ncomp, dist, nrepeat_perf, fold_perf, label)
# -------------------------------------------------------------
# Tunes, fits, and evaluates a MINT sPLS-DA model on the forward-selected
# training studies. Uses sequential tuning (comp 1 first, then comp 2 conditioned
# on comp 1) which is the mixOmics-recommended approach.
#
# IMPORTANT: `study` must exist as a factor in the CALLING ENVIRONMENT.
#            It is a factor vector of study membership (one entry per sample row in X).
#
# Parameters:
#   X            – numeric matrix; rows = samples (all forward-selected studies combined),
#                  cols = taxa (CLR-transformed)
#   Y            – factor vector; disease label per sample (same row order as X)
#   ncomp        – number of MINT components to tune and fit (default 2)
#   dist         – classification rule for tuning and evaluation (default "centroids.dist")
#   nrepeat_perf – number of repeats for final model evaluation via perf() (default 20)
#   fold_perf    – number of folds for final perf() evaluation (default 5)
#   label        – short string identifying the data integration strategy used to produce X;
#                  used as a filename prefix for the output Excel file (default "mint").
#                  Use "forward" for forward_selection() results, "individual" for
#                  individual_filtering() results.
#
# Returns: named list with
#   $model                   – fitted mint.splsda object
#   $optimal.keepX           – integer vector of selected keepX for each component
#   $performance             – perf() output object (contains BER, error.rate.class, etc.)
#   $selected.features.comp1 – character vector of taxa selected on component 1
#   $selected.features.comp2 – character vector of taxa selected on component 2
#   $features.all            – character vector of all unique selected taxa (lowercased)
run_mint <- function(X, Y, ncomp = 2, dist = "centroids.dist", nrepeat_perf=20, fold_perf=5, label = "mint") {

  # STEP 1: Tune component 1
  # test.keepX tries keeping 10 to all taxa; leave-one-study-out CV is automatic for MINT
  tune.comp1 <- tune(method = "mint.splsda",
                     X = X, Y = Y, study = study,
                     ncomp = 1,
                     test.keepX = seq(10, ncol(X), 1),
                     dist = dist)

  # Extract the optimal number of taxa to keep for component 1
  optimal.comp1 <- tune.comp1$choice.keepX
  cat("Optimal keepX for Component 1:", optimal.comp1, "\n\n")
  plot(tune.comp1, main = "Component 1 Tuning")

  # STEP 2: Tune component 2, conditioning on the locked optimal comp 1 keepX
  # already.tested.X = optimal.comp1 fixes comp 1 while searching for comp 2 optimum
  cat("Step 2: Tuning Component 2...\n")
  tune.comp2 <- tune(method = "mint.splsda",
                     X = X, Y = Y, study = study,
                     ncomp = 2,
                     test.keepX = seq(10, ncol(X), 1),
                     already.tested.X = optimal.comp1,
                     dist = dist)

  # Extract optimal keepX for component 2 (index [2] = comp 2 entry of choice.keepX)
  optimal.comp2 <- tune.comp2$choice.keepX[2]
  cat("Optimal keepX for Component 2:", optimal.comp2, "\n\n")
  plot(tune.comp2, main = "Component 2 Tuning")

  # Combine both optimal keepX values into the vector required by mint.splsda
  optimal.keepX <- c(optimal.comp1, optimal.comp2)

  # STEP 3: Build final MINT sPLS-DA model with tuned hyperparameters
  cat("Step 3: Building final model...\n")
  cat("Using keepX =", paste(optimal.keepX, collapse = ", "), "\n\n")

  model <- mint.splsda(X = X, Y = Y, study = study,
                       ncomp = ncomp,
                       keepX = optimal.keepX)

  # STEP 4: Evaluate model using repeated k-fold CV
  # validation="Mfolds" with leave-one-study-out is the MINT standard
  cat("Step 4: Evaluating model...\n")
  perf <- perf(model, validation="Mfolds", folds=fold_perf, nrepeat=nrepeat_perf, dist=dist)

  cat("\nBalanced Error Rate (BER):\n")
  print(perf$global.error$BER)

  cat("\nError Rate per Class:\n")
  print(perf$global.error$error.rate.class)

  # Extract the names of taxa selected on each component
  selected.comp1 <- selectVar(model, comp = 1)$name
  selected.comp2 <- selectVar(model, comp = 2)$name
  # Combine into a single de-duplicated list (lowercased for case-insensitive comparisons)
  ttl_taxa <- tolower(unique(c(selected.comp1, selected.comp2)))

  # STEP 5: Print a human-readable summary of selected taxa counts and names
  cat("\n========================================\n")
  cat("SELECTED TAXA SUMMARY\n")
  cat("========================================\n")
  cat("Component 1:", length(selected.comp1), "taxa selected (keepX =", optimal.comp1, ")\n")
  for (i in seq_along(selected.comp1)) {
    cat("  ", i, ".", selected.comp1[i], "\n")
  }
  cat("\nComponent 2:", length(selected.comp2), "taxa selected (keepX =", optimal.comp2, ")\n")
  for (i in seq_along(selected.comp2)) {
    cat("  ", i, ".", selected.comp2[i], "\n")
  }
  cat("\nTotal unique taxa across both components:", length(ttl_taxa), "\n")
  cat("========================================\n\n")

  # STEP 6: Save selected taxa + loading values and per-class error rates to Excel
  # -------------------------------------------------------------------------------
  # selectVar()$value returns a data.frame: row names = taxa names, col 1 = loading value
  # We extract col 1 explicitly to be robust to version differences in column naming.
  loadings.comp1 <- data.frame(
    taxa    = rownames(selectVar(model, comp = 1)$value),
    loading = selectVar(model, comp = 1)$value[, 1],
    row.names = NULL
  )
  loadings.comp2 <- data.frame(
    taxa    = rownames(selectVar(model, comp = 2)$value),
    loading = selectVar(model, comp = 2)$value[, 1],
    row.names = NULL
  )

  # perf$global.error$error.rate.class[[dist]] is a matrix:
  # rows = disease classes (e.g. Healthy, T2D), cols = components (1, 2, ...)
  err_class_mat <- perf$global.error$error.rate.class[[dist]]
  err_comp1 <- data.frame(
    class      = rownames(err_class_mat),
    error_rate = err_class_mat[, 1],  # column 1 = component 1 error rates
    row.names  = NULL
  )
  err_comp2 <- data.frame(
    class      = rownames(err_class_mat),
    error_rate = err_class_mat[, 2],  # column 2 = component 2 error rates
    row.names  = NULL
  )

  # Create output directory if it doesn't exist
  if (!dir.exists("output_files")) dir.create("output_files")

  # Install openxlsx if needed, then write 4-sheet Excel workbook
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
  library(openxlsx)

  wb <- createWorkbook()
  addWorksheet(wb, "Comp1_Taxa_Loadings")          # Sheet 1: comp 1 taxa + loading values
  addWorksheet(wb, "Comp2_Taxa_Loadings")          # Sheet 2: comp 2 taxa + loading values
  addWorksheet(wb, "Comp1_Error_Rate_per_Class")   # Sheet 3: comp 1 per-class error rate
  addWorksheet(wb, "Comp2_Error_Rate_per_Class")   # Sheet 4: comp 2 per-class error rate
  writeData(wb, "Comp1_Taxa_Loadings",        loadings.comp1)
  writeData(wb, "Comp2_Taxa_Loadings",        loadings.comp2)
  writeData(wb, "Comp1_Error_Rate_per_Class", err_comp1)
  writeData(wb, "Comp2_Error_Rate_per_Class", err_comp2)
  # Prefix the filename with `label` so outputs from different integration strategies
  # (e.g. "forward", "individual") are saved as separate files and never overwrite each other.
  xlsx_path <- paste0("output_files/", label, "_selected_taxa.xlsx")
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat("Excel saved →", xlsx_path, "\n\n")

  results <- list(
    model                  = model,          # fitted MINT sPLS-DA object
    optimal.keepX          = optimal.keepX,  # c(comp1_keepX, comp2_keepX)
    performance            = perf,           # perf() output; access $global.error$BER etc.
    selected.features.comp1 = selected.comp1, # taxa names selected on component 1
    selected.features.comp2 = selected.comp2, # taxa names selected on component 2
    features.all           = ttl_taxa        # all unique selected taxa (lowercased)
  )

  return(results)
}




## ============================================================
## 18. MINT RESULT VISUALISATION
## ============================================================

## MINT sPLS-DA visualisation: score plot + loadings for both components
# Promoted from Joyce/modelling.Rmd — Joyce Hu, 2025
# Updated: removed duplicate upper legends, fixed title truncation, added dpi=300 export — Edward Kao, 2025
# results: output list from run_mint(); Y: disease factor; study: study factor
# Figures are saved as PNG (dpi=300) in the working directory (same as the calling script).

# plot_mint_results(results, Y, study, label)
# --------------------------------------------
# Produces three plots from a fitted MINT model:
#   1. Score plot (samples coloured by disease, shaped by study, with ellipses)
#   2. Feature loadings – Component 1 (top 20 taxa by median loading)
#   3. Feature loadings – Component 2 (top 20 taxa by median loading)
#
# All plots are saved as PNG (dpi=300) under output_files/, prefixed with `label`.
#
# Parameters:
#   results – output list from run_mint() (must contain $model)
#   Y       – factor vector of disease labels (same order as samples in the model)
#   study   – factor vector of study labels (same order as samples)
#   label   – short string identifying the data integration strategy; used as a filename
#             prefix so outputs from different strategies don't overwrite each other
#             (default "mint"). Use "forward" or "individual" to match run_mint() calls.
plot_mint_results <- function(results, Y, study, label = "mint") {

  model <- results$model  # extract the fitted mint.splsda object

  # Ensure output directory exists before saving any figures
  if (!dir.exists("output_files")) dir.create("output_files")

  # --- Plot 1: Score plot ---
  # plotIndiv() shows each sample as a point in the 2D MINT latent space.
  # Colour = disease group; shape = study (pch offset +13 gives distinct symbols).
  # ellipse = TRUE draws 95% confidence ellipses per disease group.
  # plotIndiv() returns a list; $graph is the ggplot2 object for ggsave.
  indiv_plot <- plotIndiv(model,
                          group = Y,
                          pch = as.numeric(factor(study)) + 13,
                          pch.levels = study,
                          title = 'MINT sPLS-DA',
                          legend = TRUE,
                          legend.title = 'Disease',
                          legend.title.pch = 'Study',
                          ellipse = TRUE)
  ggsave(paste0("output_files/", label, "_scoreplot.png"), plot = indiv_plot$graph, dpi = 300, width = 8, height = 6)

  # --- Plot 2: Feature loadings – Component 1 ---
  # plotLoadings() shows the top ndisplay=20 taxa by their median loading magnitude.
  # contrib='max' colours taxa by the disease class they are most associated with.
  #
  # [ORIGINAL] Assumed plotLoadings() returns a ggplot object — this fails in newer mixOmics
  #            because plotLoadings() returns a data.frame (loading table), not a ggplot.
  #            Calling print() on a data.frame triggers the grid.draw error; ggsave() also fails.
  # p_load1 <- plotLoadings(model, comp=1, legend=TRUE, contrib='max',
  #                         title='Feature Loadings - Component 1', method='median', ndisplay=20)
  # print(p_load1)
  # ggsave("MINT_loadings_comp1.png", plot=p_load1, dpi=300, width=10, height=7)
  #
  # [UPDATED] Call plotLoadings() without capturing return value — it renders the plot inline
  #           as a side effect. Save separately using png()/dev.off() (works for base R graphics).
  plotLoadings(model,
               comp = 1,
               legend = TRUE,
               contrib = 'max',
               title = 'Feature Loadings - Component 1',
               method = 'median',
               ndisplay = 20)
  png(paste0("output_files/", label, "_loadings_comp1.png"), width = 10 * 300, height = 7 * 300, res = 300)
  plotLoadings(model,
               comp = 1,
               legend = TRUE,
               contrib = 'max',
               title = 'Feature Loadings - Component 1',
               method = 'median',
               ndisplay = 20)
  dev.off()

  # --- Plot 3: Feature loadings – Component 2 ---
  # [ORIGINAL] Same issue as comp 1 above.
  # p_load2 <- plotLoadings(model, comp=2, legend=TRUE, contrib='max',
  #                         title='Feature Loadings - Component 2', method='median', ndisplay=20)
  # print(p_load2)
  # ggsave("MINT_loadings_comp2.png", plot=p_load2, dpi=300, width=10, height=7)
  #
  # [UPDATED] Same fix as comp 1: render inline, then save via png()/dev.off().
  plotLoadings(model,
               comp = 2,
               legend = TRUE,
               contrib = 'max',
               title = 'Feature Loadings - Component 2',
               method = 'median',
               ndisplay = 20)
  png(paste0("output_files/", label, "_loadings_comp2.png"), width = 10 * 300, height = 7 * 300, res = 300)
  plotLoadings(model,
               comp = 2,
               legend = TRUE,
               contrib = 'max',
               title = 'Feature Loadings - Component 2',
               method = 'median',
               ndisplay = 20)
  dev.off()
}
