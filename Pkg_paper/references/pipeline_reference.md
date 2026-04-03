# Pipeline Reference — Gut Microbiome T2D MINT Workflow

**Project**: Gut Microbiome T2D — Integrative MINT Workflow
**Researcher**: Edward Kao (Pei-Hsiu Kao), MSc Bioinformatics, University of Melbourne
**Supervisors**: Prof Kim-Anh Lê Cao, Dr Saritha Kodikara, Dr Yiwen Wang
**Collaborator**: Joyce Hu (co-developer, helper functions and preprocessing)
**Last updated**: 2026-04-04

---

## How to use this document

This file is the primary reference for the pipeline. Read it at the start of every session to recall what every variable, function, and script does.

**Two copies are maintained in sync — always update both when changes are made:**
- `ClaudeHelp/Gut_Microbiome_Project-Edward-s-Paper-/Pkg_paper/references/pipeline_reference.md` (project-level, alongside the scripts)
- `ClaudeHelp/references/pipeline_reference.md` (top-level, for quick access from any session)

It covers:

1. [Pipeline overview and architecture](#1-pipeline-overview)
2. [Data sources](#2-data-sources)
3. [Script map](#3-script-map)
4. [Variable glossary — Test_Framework(BioInterp).Rmd](#4-variable-glossary)
5. [Function reference — Function_Helpers.R](#5-function-reference)
6. [Key data structures](#6-key-data-structures)
7. [Known issues and pending decisions](#7-known-issues-and-pending-decisions)

---

## 1. Pipeline Overview

**Goal**: Identify gut microbial species that robustly discriminate Type 2 Diabetes (T2D) from Normal Glucose Tolerance (NGT/Healthy) individuals across multiple independent cohort studies.

**Key challenge**: Each study has batch effects (genetic, geographic, methodological differences). Taxa names may differ, and not all taxa are shared across all studies.

**Solution**: The MINT (Multivariate INTegrative) framework accounts for batch effects by modelling study as a blocking factor inside the multivariate model, without requiring identical taxa across studies.

### Pipeline flow (from Research_Structure.jpg)

```
CMD (5 studies: Feng2015, Karl2013, LiJ2014, YuJ2015, Sank2015)
       ↓  [HMP2019 discarded]
  Data Loading
       ↓
  Preprocessing per study      (run_preprocessing)
  [offset → low counts removal → CLR]
       ↓
  ┌────────────┬─────────────────┐
  ↓            ↓
  Forward      Individual
  Selection    Filtering
  (Layer 3)    (Layer 2)
       ↓            ↓
  FS model     IF model
  (run_mint)   (run_mint)
  └─ choose lower BER ─┘
         ↓
     Final Model
  ┌──────────────────┐
  ↓                  ↓
Biological       Generalisation
Interpretation   (align_columns → test on HGMA)
```

Both strategies (Forward Selection and Individual Filtering) are run. The one with lower BER is the Final Model. Post-modelling splits into two parallel branches: Biological Interpretation (training data analysis) and Generalisation (external validation on HGMA).

**Training studies (5):** FengQ_2015, KarlssonFH_2013, LiJ_2014, YuJ_2015, SankaranarayananK_2015
**Testing studies (2):** PRJEB1786 (Swedish), PRJNA422434 (Chinese) — testing set only, never used in training

---

## 2. Data Sources

### Training data — CuratedMetagenomicData (CMD) v3.16.1

**Five** studies are used as training cohorts. **HMP_2019_t2d is confirmed DISCARDED** (per Research_Structure.jpg). Each is loaded with `curatedMetagenomicData()`:

| Study name (short) | CMD query string | Population | Notes |
|---|---|---|---|
| FengQ_2015 | `FengQ_2015.relative_abundance` | Chinese | Colorectal cancer cohort; T2D/Healthy subset |
| KarlssonFH_2013 | `KarlssonFH_2013.relative_abundance` | Swedish | Canonical T2D gut microbiome study |
| LiJ_2014 | `LiJ_2014.relative_abundance` | Chinese | Large reference metagenome cohort |
| YuJ_2015 | `YuJ_2015.relative_abundance` | Chinese | Colorectal cancer cohort; T2D/Healthy subset |
| SankaranarayananK_2015 | `SankaranarayananK_2015.relative_abundance` | Indian | Indian population |

CMD taxa names include a `"species: "` prefix (e.g. `"species: Bacteroides uniformis"`) that must be stripped before comparing with HGMA taxa names.

### Testing data — Human Gut Microbiome Atlas (HGMA)

Two external studies used for generalisability evaluation. Loaded directly from pre-built `.rds` files (no longer built via `draw_TestSet()`):

| BioProject | Population | File used |
|---|---|---|
| PRJEB1786 | Swedish | `PRJEB1786.rds` — pre-built by `HGMA_Processing.Rmd` |
| PRJNA422434 | Chinese | `PRJNA422434.rds` — pre-built by `HGMA_Processing.Rmd` |

HGMA taxa are indexed by numeric IDs and mapped to species names via `corresponding_taxa.csv`. Some HGMA species names use `/` to indicate multiple possible identities (e.g. `"Firmicutes bacterium CAG:65 / [Clostridium] sp. 2789STDY5608883"`). These are resolved by `align_columns()`.

### Reference Papers

Key papers for the statistical methods used in this pipeline:

| File | Content |
|---|---|
| `references/model_paper/MINT.pdf` | Rohart et al. 2017 — the MINT method paper |
| `references/model_paper/MINT additional file.pdf` | MINT supplementary material |
| `references/model_paper/Le-cao-2008-A-sparse-pls-for-variable-selection.pdf` | sPLS / variable selection foundation |
| `references/model_paper/PLSDA.pdf` | PLS-DA methodology reference |
| `references/model_paper/Rohart-2017-Mixomics-an-r-package-for-omics-fea.pdf` | mixOmics R package paper |

---

## 3. Script Map

All files below are relative to `ClaudeHelp/Gut_Microbiome_Project-Edward-s-Paper-/Pkg_paper/`.

| File | Role | Status |
|---|---|---|
| `Function_Helpers.R` | All reusable helper functions. Always `source()` this at the top of every Rmd. | Active |
| `Test_Framework.Rmd` | **THE active combined sandbox** — full pipeline end-to-end. Sections: Data Loading → Data Processing → Model Development (IF + FS) → Strategy Selection → Biological Interpretation (stub) → Generalisation (HGMA preprocessing, `align_columns()`, external validation). Replaces the two former sandbox files. | Active |
| `Test_Framework_old.Rmd` | Original pre-split sandbox (kept as backup/reference). Do not edit. | Reference only |
| `full_framework.Rmd` | THE formal pipeline file. All finalised code goes here once sandbox is confirmed. | Skeleton only |
| `HGMA_Processing.Rmd` | Builds HGMA test sets; standalone notebook for HGMA ingestion. | Reference |
| `Batch_effect_management.Rmd` | Demo/reference only. Not part of active pipeline. | Reference only |
| `Joyce/` | Joyce's intermediate working files. Not original, not final. Reference only. | Reference only |
| `references/pipeline_reference.md` | **This file** — pipeline documentation. Mirrored at `ClaudeHelp/references/pipeline_reference.md`. Always update both copies together. | Active |

---

## 4. Variable Glossary

All variables listed here appear in `Test_Framework.Rmd` (the combined active sandbox).

### Data Loading section

| Variable | Type | Description |
|---|---|---|
| `req_pkgs` | `character vector` | All R packages required by the pipeline. Passed to `install_if_missing()`. |
| `study_names` | `character vector` | Five CMD query strings in format `"<Author>_<Year>.relative_abundance"`. Used to download data from CMD. |
| `study_names_v2` | `character vector` | Shortened study names (e.g. `"FengQ_2015"`). Used as list names and study labels downstream. |
| `abundance_lst` | `named list of data.frames` | One entry per study. Each data.frame: rows = samples (sample IDs as row names), columns = taxa (species names). Values are raw counts. Named by `study_names_v2` after the name-assignment chunk. |
| `meta_lst` | `named list of data.frames` | One entry per study. Each data.frame: rows = samples, columns = all CMD metadata fields (disease, gender, body_site, age, BMI, etc.). Named by `study_names_v2`. |

### Data Processing section

| Variable | Type | Description |
|---|---|---|
| `keep_cols` | `character vector` | Extra metadata columns to retain alongside `disease`: `c("gender", "body_site")`. |
| `diseases` | `character vector` | Disease classes to discriminate: `c("Healthy", "T2D")`. All other disease labels are dropped during preprocessing. |
| `res` | `named nested list` | Output from `run_preprocessing()`. Keys = study short names. Each entry: `$abund` (CLR-transformed abundance, samples × taxa) and `$meta` (filtered metadata with disease/gender/body_site columns). This is the primary data structure passed to filtering functions. |

### Model Development section

The modelling section runs **both** data integration strategies in parallel, each producing its own set of named variables and output files distinguished by the strategy suffix (`_individual` or `_forward`).

#### Individual Filtering strategy variables

| Variable | Type | Description |
|---|---|---|
| `result_individual` | `named nested list` | Output from `individual_filtering()`. Same structure as `res` (`$abund` + `$meta` per study), containing only studies that passed the per-study BER filter. |
| `result_df_individual` | `data.frame (list result of bind_rows)` | All individually-filtered studies merged into one data frame. Sub-elements: `$abund` (combined CLR abundance) and `$meta` (combined metadata including `$study` column). |
| `p_selected_individual` | `data.frame` | Alias for `result_df_individual$abund`. The individually-filtered CLR training abundance matrix. Used in `align_columns()` in the Generalisation section. |
| `mint_results_individual` | `named list` | Output from `run_mint(..., label="individual")`. Contains `$model`, `$optimal.keepX`, `$performance`, `$selected.features.comp1`, `$selected.features.comp2`, `$features.all`. |
| `predict_individual` | `list` | `predict()` output on the individual-strategy training set. Access via `$class$centroids.dist[, comp]`. |
| `conf_mat_individual_c1` | `matrix` | Confusion matrix (truth × predicted) for individual-strategy training predictions at component 1. |
| `conf_mat_individual_c2` | `matrix` | Confusion matrix (truth × predicted) for individual-strategy training predictions at component 2. |

#### Forward Selection strategy variables

| Variable | Type | Description |
|---|---|---|
| `result_forward` | `named nested list` | Output from `forward_selection()`. Same structure as `res`, containing only the greedy-optimal study combination. |
| `result_df_forward` | `data.frame (list result of bind_rows)` | All forward-selected studies merged into one data frame. Sub-elements: `$abund` and `$meta` (including `$study` column). |
| `p_selected_forward` | `data.frame` | Alias for `result_df_forward$abund`. The forward-selected CLR training abundance matrix. Used in `align_columns()` in the Generalisation section. |
| `mint_results_forward` | `named list` | Output from `run_mint(..., label="forward")`. Same structure as `mint_results_individual`. |
| `predict_forward` | `list` | `predict()` output on the forward-strategy training set. Access via `$class$centroids.dist[, comp]`. |
| `conf_mat_forward_c1` | `matrix` | Confusion matrix (truth × predicted) for forward-strategy training predictions at component 1. |
| `conf_mat_forward_c2` | `matrix` | Confusion matrix (truth × predicted) for forward-strategy training predictions at component 2. |

#### Generalisation section variables (HGMA test sets)

| Variable | Type | Description |
|---|---|---|
| `test1` | `named list` | PRJEB1786 (Swedish) test set loaded from `readRDS("PRJEB1786.rds")`. After preprocessing: `$abund` = CLR-transformed abundance (samples × taxa), `$meta` = metadata with `$disease`. |
| `test2` | `named list` | PRJNA422434 (Chinese) test set loaded from `readRDS("PRJNA422434.rds")`. Note: `$meta$disease` labels are recoded from `"NGT"` → `"Healthy"` before preprocessing. |
| `offset1` | `numeric` | Minimum positive value in `test1$abund` — used as CLR offset for PRJEB1786 (HGMA uses relative abundances, not counts, so offset=1 is inappropriate). |
| `offset2` | `numeric` | Minimum positive value in `test2$abund` — CLR offset for PRJNA422434. |
| `processed1` / `processed2` | `list` | Output from `preprocess_single()` on each test set. `$data` = CLR-transformed abundance; `$meta` = aligned metadata. Written back to `test1$abund`/`test2$abund` after call. |
| `combined_obj1` | `list` | Output from `align_columns(selected_set$abund, test1$abund)`. Sub-elements: `$training` (CMD training abundance restricted to common taxa), `$testing` (PRJEB1786 abundance restricted to common taxa), `$change_log`, `$messages`. |
| `combined_obj2` | `list` | Same structure as `combined_obj1` but for PRJNA422434 (test2). |
| `mint_results_test1` | `named list` | Output from `run_mint(..., label="prjeb1786", model_mode="testing")` trained on `combined_obj1$training`. No xlsx saved. Used only to obtain `$model` for `predict()`. |
| `mint_results_test2` | `named list` | Same structure as `mint_results_test1` but trained on `combined_obj2$training` for PRJNA422434. |
| `predict_test1` / `predict_test2` | `list` | `predict()` output on the respective HGMA test set. Access via `$class$centroids.dist[, comp]`. |
| `cm_test1_c1` / `cm_test1_c2` | `matrix` | PRJEB1786 confusion matrices (truth × predicted) at comp 1 and comp 2 respectively. |
| `cm_test2_c1` / `cm_test2_c2` | `matrix` | PRJNA422434 confusion matrices (truth × predicted) at comp 1 and comp 2 respectively. |

#### Shared temporary variables (reused by both strategies)

| Variable | Type | Description |
|---|---|---|
| `X` | `data.frame` | Abundance matrix input for `run_mint()`. Assigned from `result_df_individual$abund` or `result_df_forward$abund` before each `run_mint()` call. |
| `Y` | `factor` | Disease outcome vector. Levels: `"Healthy"`, `"T2D"`. Reassigned before each `run_mint()` call. |
| `study` | `factor` | Study membership vector. MUST exist in the calling environment — `run_mint()` reads it directly. Reassigned before each `run_mint()` call. |

---

## 5. Function Reference

All functions live in `Function_Helpers.R`. Source it at the top of every notebook.

---

### `install_if_missing(packages)`

Installs any packages not yet installed, then loads all of them.

| Parameter | Type | Description |
|---|---|---|
| `packages` | `character vector` | Package names to install + load |

**Returns**: nothing (invisible)

---

### `filter_meta(meta, diseases, cols_to_keep)`

Filters a metadata data.frame to samples in target disease groups. Normalises raw CMD disease labels (e.g. `"T2D_prediabetes"` → `"T2D"`) via case-insensitive regex matching.

| Parameter | Type | Description |
|---|---|---|
| `meta` | `data.frame` | Metadata; must contain a `disease` column; row names = sample IDs |
| `diseases` | `character vector` | Disease labels to keep (e.g. `c("Healthy","T2D")`) |
| `cols_to_keep` | `character vector` or NULL | Additional columns to retain (e.g. `c("gender","body_site")`) |

**Returns**: filtered metadata data.frame with `disease` + requested columns; row names preserved.

---

### `low_count_removal(data, percent)`

Removes taxa whose total relative abundance across all samples is below `percent`%.

| Parameter | Type | Description |
|---|---|---|
| `data` | `data.frame` or `matrix` | Samples × taxa (must be numeric) |
| `percent` | `numeric` | Minimum % of total abundance to keep a taxon (default `0.01`) |

**Returns**: `list($data.filter, $keep.otu)` — filtered data frame and indices of kept taxa.

---

### `split_count_filtering(labelled_abund, diseases, percent, filtered_meta)`

Applies `low_count_removal` within each disease group separately, then takes the UNION of surviving taxa. Ensures taxa present in at least one class are not accidentally removed.

| Parameter | Type | Description |
|---|---|---|
| `labelled_abund` | `data.frame` | Samples × taxa abundance |
| `diseases` | `character vector` | Disease groups to split by |
| `percent` | `numeric` | Passed to `low_count_removal` |
| `filtered_meta` | `data.frame` | Filtered metadata; row names must match `labelled_abund` rows |

**Returns**: filtered abundance data.frame (samples × surviving taxa).

---

### `clr.transformation(abundance.filtered, offset_value)`

Applies the Centered Log-Ratio (CLR) transformation. Converts compositional count data into real-valued log-ratio coordinates suitable for multivariate analysis (PLS-DA, MINT).

| Parameter | Type | Description |
|---|---|---|
| `abundance.filtered` | `matrix` or `data.frame` | Samples × taxa (numeric) |
| `offset_value` | `numeric` | Added to all values before log to handle zeros. Use `1` for count data; use the second-smallest non-zero value for relative abundances (HGMA) |

**Returns**: CLR-transformed data.frame (same dimensions).

---

### `preprocess_single(meta, abund, diseases, keep_cols, percent, offset_value, CLR)`

Full single-study preprocessing: `filter_meta` → align abundance rows → `split_count_filtering` → `clr.transformation`.

| Parameter | Type | Description |
|---|---|---|
| `meta` | `data.frame` | Study metadata |
| `abund` | `data.frame` | Study abundance (samples × taxa) |
| `diseases` | `character vector` | Disease groups to keep |
| `keep_cols` | `character vector` | Extra metadata columns |
| `percent` | `numeric` | Low-count removal threshold |
| `offset_value` | `numeric` | CLR offset |
| `CLR` | `logical` | If FALSE, skip CLR and return filtered counts only (default TRUE) |

**Returns**: `list($data, $meta)` — CLR-transformed abundance + filtered metadata.

---

### `preprocess_multiple(abundance_list, meta_list, diseases, keep_cols, percent, offset_value, threshold)`

Applies `preprocess_single` to each study in the list, then calls `integrate_filter` to align taxa across studies.

| Parameter | Type | Description |
|---|---|---|
| `abundance_list` | `named list` | One abundance data.frame per study |
| `meta_list` | `named list` | One metadata data.frame per study (same names) |
| `diseases` | `character vector` | Disease groups to keep |
| `keep_cols` | `character vector` | Extra metadata columns |
| `percent` | `numeric` | Low-count removal threshold |
| `offset_value` | `numeric` | CLR offset |
| `threshold` | `integer` | Minimum number of studies a taxon must appear in (default NULL) |

**Returns**: named nested list from `integrate_filter`: `result[[study_name]]$abund` + `$meta`.

---

### `integrate_filter(results, threshold)`

Finds taxa present in at least `threshold` studies and filters each study to only those taxa, giving all studies a common (though not necessarily identical) feature space.

| Parameter | Type | Description |
|---|---|---|
| `results` | `named list` | Output from study-wise preprocessing (each entry has `$data` + `$meta`) |
| `threshold` | `integer` | Minimum study count for a taxon to survive |

**Returns**: named nested list: `final[[study_name]]$abund` + `$meta`.

---

### `run_preprocessing(abundance_list, meta_list, diseases, keep_cols, percent, offset_value, threshold)`

Top-level dispatcher. Routes to `preprocess_single` (1 study) or `preprocess_multiple` (many studies). This is the function called in the Rmd notebooks.

**Returns**: single study → `list($data, $meta)`; multiple studies → nested list from `integrate_filter`.

---

### `individual_filtering(data, err_thres, ncomp, validation, folds, nrepeat, dist)`

**Layer 2**. Fits a per-study PLS-DA model and removes studies where any per-class error rate exceeds `err_thres`.

| Parameter | Type | Description |
|---|---|---|
| `data` | `named list` or single study | Output from `run_preprocessing` |
| `err_thres` | `numeric` | Max tolerated per-class BER (default `0.55`; thesis specifies `0.7`; sandbox uses `0.5`) |
| `ncomp` | `integer` | PLS-DA components to fit (default `2`) |
| `validation` | `"Mfold"` | CV strategy |
| `folds` | `integer` | CV folds (default `5`) |
| `nrepeat` | `integer` | CV repeats (default `10`) |
| `dist` | `character` | Classification rule (default `"centroids.dist"`) |

**Returns**: filtered named list of studies that passed, or single logical (for single study input).

---

### `calculate_BER(single_study, ncomp, validation, folds, nrepeat, dist)`

Fits PLS-DA on one study and returns BER at `ncomp`. Used internally by `forward_selection` to rank seed studies.

**BER definition**: `BER = (1/K) × Σ(misclassified_k / n_k)` — weighted average error rate across K classes; robust to class imbalance.

**Returns**: numeric BER value.

---

### `eliminate_single_taxa(lst_df, study_col, n_taxa)`

Removes taxa that appear in only 1 study (or `n_taxa` studies). Taxa unique to one study cannot be integrated across cohorts in MINT.

| Parameter | Type | Description |
|---|---|---|
| `lst_df` | `data.frame` | Output of `bind_rows(result)` — combined data frame with `$abund` and `$meta$study` |
| `study_col` | `character` | Column name holding study label (default `"study"`) |
| `n_taxa` | `integer` | Remove taxa with exactly this many unique values within a study (default `1` = remove constant-across-study taxa) |

**Returns**: filtered abundance data.frame including the study column.

---

### `calculate_multipleDB_BER(lst, ncomp, study_col, ...)`

Computes BER for a candidate combination of studies using MINT sPLS-DA. Used iteratively inside `forward_selection`.

**Returns**: numeric BER scalar for this study combination.

---

### `forward_selection(data, ncomp, study_col, validation, folds, nrepeat, dist, n_taxa)`

**Layer 3**. Greedy study selection algorithm.

**Algorithm**:
- Round 0 (seed): evaluate all individual studies; select the one with the lowest BER.
- Round 1+: for each remaining study, temporarily add it to the current combination and compute BER via `calculate_multipleDB_BER`. Add the study that gives the lowest BER. Stop if no addition improves BER.

| Parameter | Type | Description |
|---|---|---|
| `data` | `named list` | Preprocessed study pool (output of `run_preprocessing`) |
| `ncomp` | `integer` | MINT components (default `2`) |
| `nrepeat` | `integer` | CV repeats for BER calculation (default `20`) |
| `dist` | `character` | Classification rule (default `"centroids.dist"`) |
| `n_taxa` | `integer` | Passed to `eliminate_single_taxa` (default `1`) |

**Returns**: named list of selected study objects (the optimal combination).

---

### `align_columns(p_test, p_selected)`

Reconciles taxa column names between CMD training data and HGMA testing data so both matrices can be directly compared.

**Steps**:
1. Strip `"species: "` prefix from both matrices (CMD adds this; HGMA does not).
2. For each HGMA column with `/` aliases, match against training column names:
   - Single match → rename test column.
   - Multiple matches → sum training columns into new combined column, rename test column.
   - No match → leave unchanged.
3. Subset both matrices to common taxa.

| Parameter | Type | Description |
|---|---|---|
| `p_test` | `data.frame` | HGMA test abundance (samples × taxa) |
| `p_selected` | `data.frame` | CMD training abundance (samples × taxa) |

**Returns**: `list($training, $testing, $change_log, $messages)`.

---

### `run_mint(X, Y, ncomp, dist, nrepeat_perf, fold_perf, label, model_mode)`

**Layers 4 & 5**. The main MINT sPLS-DA modelling function. Performs sequential hyperparameter tuning, fits the final model, and evaluates it.

**Important**: `study` must exist as a factor in the **calling environment** — `run_mint()` references it directly. Assign `study <- as.factor(result_df_individual$meta$study)` (or `result_df_forward`) before calling.

**Sequential tuning rationale**: tuning comp 1 first (with `ncomp=1`), then tuning comp 2 with `already.tested.X = optimal.comp1` is the mixOmics-recommended approach. It is superior to joint tuning (`tune(ncomp=2)`) because it prevents comp 2 from interfering with comp 1 optimisation.

| Parameter | Type | Description |
|---|---|---|
| `X` | `matrix`/`data.frame` | Combined training abundance (samples × taxa) |
| `Y` | `factor` | Disease outcome per sample |
| `ncomp` | `integer` | Number of MINT components (default `2`) |
| `dist` | `character` | Classification rule: `"centroids.dist"`, `"mahalanobis.dist"`, or `"max.dist"` |
| `nrepeat_perf` | `integer` | Repeats for final `perf()` evaluation (default `20`) |
| `fold_perf` | `integer` | Folds for final `perf()` evaluation (default `5`) |
| `label` | `character` | Strategy label used to prefix output filenames (default `"mint"`). Use `"individual"` or `"forward"` for training models; `"prjeb1786"` / `"prjna422434"` for test-set refits. |
| `model_mode` | `character` | `"training"` (default) saves the `.xlsx` output; `"testing"` skips the xlsx save. Use `"testing"` when refitting on the common-taxa intersection for external validation, where only the fitted model object is needed. |

**Returns**: named list:

| Field | Type | Description |
|---|---|---|
| `$model` | `mint.splsda` | Fitted MINT sPLS-DA model object |
| `$optimal.keepX` | `integer vector` | `c(keepX_comp1, keepX_comp2)` — number of taxa kept per component |
| `$performance` | `perf object` | Full `perf()` output. Access `$global.error$BER` for BER matrix, `$global.error$error.rate.class` for per-class errors |
| `$selected.features.comp1` | `character vector` | Taxa selected on component 1 (length = `optimal.keepX[1]`) |
| `$selected.features.comp2` | `character vector` | Taxa selected on component 2 (length = `optimal.keepX[2]`) |
| `$features.all` | `character vector` | De-duplicated union of comp 1 + comp 2 taxa (lowercased) |

**Side effects** (in addition to the return value):
- Prints a formatted taxa summary to console: count + ranked list of names for each component.
- Creates `output_files/` directory if it does not exist.
- **Only when `model_mode == "training"`**: saves `output_files/{label}_selected_taxa.xlsx` — 4 sheets: `Comp1_Taxa_Loadings`, `Comp2_Taxa_Loadings`, `Comp1_Error_Rate_per_Class`, `Comp2_Error_Rate_per_Class`. Requires `openxlsx` (listed in `req_pkgs`).
- When `model_mode == "testing"`, the xlsx save is skipped. External validation performance is instead reported via `predict()` + `get.confusion_matrix()` in the calling notebook.

---

### `plot_mint_results(results, Y, study, label)`

Visualises the fitted MINT model. Produces three plots saved as PNG (dpi=300):

| File saved | Content |
|---|---|
| `output_files/{label}_scoreplot.png` | 2D score plot: samples in MINT latent space, coloured by disease, shaped by study, with 95% ellipses |
| `output_files/{label}_loadings_comp1.png` | Top 20 taxa by loading magnitude on component 1; coloured by contributing disease class |
| `output_files/{label}_loadings_comp2.png` | Top 20 taxa by loading magnitude on component 2; coloured by contributing disease class |

Files are saved to `output_files/` relative to the working directory (created automatically if absent). The `label` prefix distinguishes output from the two strategies — e.g. `individual_scoreplot.png` vs `forward_scoreplot.png`.

**Implementation note**: `plotIndiv()` returns a list with a `$graph` ggplot element → saved via `ggsave()`. `plotLoadings()` in newer mixOmics versions returns a **data.frame** (not a ggplot) and renders the plot as a side effect — calling `print()` or `ggsave()` on this data.frame causes a `grid.draw` error. Fix: call `plotLoadings()` directly (auto-renders inline), then save separately using `png()`/`dev.off()`.

| Parameter | Type | Description |
|---|---|---|
| `results` | `list` | Output from `run_mint()` |
| `Y` | `factor` | Disease factor (same row order as the model's training data) |
| `study` | `factor` | Study factor (same row order) |
| `label` | `character` | Strategy label used to prefix output PNG filenames (default `"mint"`). Use `"individual"` or `"forward"` to match the `run_mint()` label. |

---

## 6. Key Data Structures

### The `res` / `result_individual` / `result_forward` nested list (most common structure)

```
res  (or result_individual, result_forward — same shape)
├── FengQ_2015
│   ├── $abund   # data.frame: samples × taxa (CLR-transformed, numeric)
│   └── $meta    # data.frame: samples × (disease, gender, body_site)
├── KarlssonFH_2013
│   ├── $abund
│   └── $meta
... (one entry per study)
```

Row names of `$abund` and `$meta` are sample IDs (character strings). They are always aligned (same order).

### The `result_df_individual` / `result_df_forward` combined frames (after `bind_rows`)

```
result_df_individual  (or result_df_forward — same shape)
├── $abund   # data.frame: (all_samples) × taxa — all studies stacked
└── $meta    # data.frame: (all_samples) × (disease, gender, body_site, study)
```

### CLR transformation and offset

The CLR transform is: `CLR_i = log(x_i / geometric_mean(x))`. Zeros are handled by adding `offset_value` before computing. For CMD count data, `offset_value = 1` is standard. For HGMA relative abundances, `offset_value` is the second-smallest non-zero value in the dataset (dataset-specific small constant).

### `study` factor and calling environment

`run_mint()`, `calculate_multipleDB_BER()`, and `tune()` all require `study` as a factor vector. In the Rmd notebooks, this must be assigned **before** calling `run_mint()`. In `Test_Framework(BioInterp).Rmd`, this is done twice — once per strategy:
```r
# Individual filtering strategy:
study <- as.factor(result_df_individual$meta$study)
mint_results_individual <- run_mint(X, Y, ncomp=2, dist="centroids.dist", label="individual")

# Forward selection strategy:
study <- as.factor(result_df_forward$meta$study)
mint_results_forward <- run_mint(X, Y, ncomp=2, dist="centroids.dist", label="forward")
```
This is a known quirk of the mixOmics API — `study` is not passed as a function argument to `run_mint()` but is read from the calling environment.

---

## 7. Known Issues and Pending Decisions

| # | Issue | Status |
|---|---|---|
| 1 | `err_thres` value: sandbox uses `0.5`, `individual_filtering` default is `0.55`, thesis specifies `0.7`. Needs reconciling before writing to `full_framework.Rmd`. | **OPEN** |
| 2 | ~~`Test_Framework(Reproducibility).Rmd` not yet completed~~ — **RESOLVED by merge**: `Test_Framework(BioInterp).Rmd` and `Test_Framework(Reproducibility).Rmd` have been merged into a single `Test_Framework.Rmd`. The Generalisation section (HGMA preprocessing via `preprocess_single()`, offset computation, NGT recoding) is in progress within the merged file. | **IN PROGRESS** |
| 3 | `full_framework.Rmd` not yet started. All finalised code must be copied here once `Test_Framework.Rmd` is confirmed working end-to-end. | **PENDING** |
| 4 | ~~HMP2019 study not yet decided~~ — **RESOLVED**: HMP_2019_t2d is confirmed DISCARDED (per Research_Structure.jpg). Training set is 5 studies. | **CLOSED** |
| 5 | `QinJ_2012` in `full_framework.Rmd` studies vector — this study is not in the 5-study training set and must be removed before finalising `full_framework.Rmd`. | **OPEN** |
| 6 | `draw_TestSet()` — **REMOVED** from `Function_Helpers.R`. HGMA test data is now loaded directly from `PRJEB1786.rds` / `PRJNA422434.rds`. | **CLOSED** |

---

*This document is auto-maintained. Update it whenever functions, variables, or pipeline decisions change.*
