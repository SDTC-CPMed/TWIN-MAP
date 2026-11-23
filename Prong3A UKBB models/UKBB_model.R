# ==============================================================================
# PROJECT: UK Biobank Proteomics Risk Prediction Analysis (Ridge AUC + Covariates)
# DESCRIPTION:
# One script for all three organs. Organ settings (disease files, signature, known
# marker, covariates, output folder) are in config/<organ>.json.
#   Rscript UKBB_model.R prostate      (or colorectal, lung)
# Optional environment variables:
#   UKBB_NCORES=15      parallel workers
#   PREFLIGHT_ONLY=1    stop after the signature -> Olink mapping table
#   UKBB_NBOOT, UKBB_OUT_DIR   fewer iterations / another output folder (for testing)
#
# Ridge regression (alpha = 0, cv.glmnet) over 1,000 bootstrap iterations (random
# stratified 70/30 train/test splits), comparing six model specifications (model_specs):
#   - Signature:            signature proteins only
#   - Covariates:           clinical covariates only
#   - Signature_Covariates: signature proteins + covariates
#   - known marker:         KLK3 (PSA) for prostate, CEACAM5 (CEA) for colorectal and lung
#   - Random, Random_Covariates: a random protein panel the size of the measured
#     signature, redrawn in every iteration (null comparator), without / with covariates
# Covariates:
#   prostate:   age, BMI, HbA1c, ethnicity, diabetes (sex is listed but is constant
#               among men, so drop_constant_factors() removes it)
#   colorectal: age, sex, BMI, smoking, alcohol, diabetes
#   lung:       age, sex, BMI, smoking, alcohol
#
# Controls selected via MatchIt nearest-neighbour without replacement, exact on sex,
# age caliper in absolute years (std.caliper = FALSE): NonCancer and progression
# 5:1 within +-2 years; healthy controls (HC) 1:1 within +-1 year (the HC pool is
# small and younger than the cases, so 5:1 left an age gap; unmatched cases dropped).
# Progression (Dis1andDis2) cases are incident Dis2 only (diagnosed after blood draw).
# Features passed jointly into glmnet; categorical covariates one-hot encoded
# via model.matrix (treatment contrasts, reference level dropped).
# Imputation: kNN on proteins in training set; training-mean substitution on test.
# Evaluation: AUC (pROC) on held-out 30% test set across observation windows
# 2, 4, 6, 8, 10, 12, 14 years before/at diagnosis.
# ==============================================================================

rm(list = ls())

# --- 0. PACKAGES --- #
packages <- c(
  "dplyr", "ggplot2", "patchwork", "stringr", "data.table", "glmnet",
  "MatchIt", "purrr", "pROC", "foreach", "doParallel",
  "tidyr", "tibble", "impute", "HGNChelper"
)
invisible(lapply(packages, library, character.only = TRUE))

# One BLAS/OpenMP thread per process: the forked workers would otherwise each start
# multi-threaded BLAS calls and oversubscribe the CPU (silent worker crashes/hangs)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
} else {
  message("RhpcBLASctl not installed -- install.packages('RhpcBLASctl') for a ",
          "runtime BLAS thread cap; falling back to env vars only, which may ",
          "not affect an already-initialized BLAS runtime.")
}

#### --- 1. CONFIGURATION --- ####
# Data and result folders: set in paths.env (see config.R)
if (!exists("earlyt_path")) local({ a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); d <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", a[1], fixed = TRUE))) else getwd()
  while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config.R"))) stop("config.R not found: run this script from inside the EARLY-T code folder")
  Sys.setenv(EARLYT_CODE = d); source(file.path(d, "config.R")) })

# Organ settings: config/<organ>.json next to this script
organ <- commandArgs(trailingOnly = TRUE)[1]
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(gsub("~+~", " ", script_file, fixed = TRUE)))   # Rscript writes spaces as ~+~
cfg_file <- file.path(script_dir, "config", paste0(organ, ".json"))
if (!file.exists(cfg_file)) stop("Usage: Rscript UKBB_model.R <prostate|colorectal|lung>")
cfg <- jsonlite::read_json(cfg_file, simplifyVector = TRUE)

input.dir  <- earlyt_path("EARLYT_UKBB_SERVER")
input.dir2 <- earlyt_path("EARLYT_UKBB_SERVER")
out.dir    <- Sys.getenv("UKBB_OUT_DIR",
                         unset = file.path(earlyt_path("EARLYT_UKBB_SERVER"), "UKBB_Results", cfg$output_folder))
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

dir.bootstrap   <- file.path(out.dir, "bootstrap_raw")
dir.predictions <- file.path(out.dir, "predictions")
dir.thresholds  <- file.path(out.dir, "thresholds")
dir.summary     <- file.path(out.dir, "summary")
dir.models      <- file.path(out.dir, "models")
dir.importance  <- file.path(out.dir, "importance")

invisible(lapply(
  list(dir.bootstrap, dir.predictions, dir.thresholds,
       dir.summary, dir.models, dir.importance),
  dir.create, showWarnings = FALSE, recursive = TRUE
))

dis1             <- cfg$dis1   # premalignant / precursor condition
dis2             <- cfg$dis2   # cancer
observation_years <- c(2, 4, 6, 8, 10, 12, 14)
n_bootstrap      <- as.integer(Sys.getenv("UKBB_NBOOT", unset = "1000"))
alpha            <- 0          # Ridge
match_ratio      <- 5
# Workers: set UKBB_NCORES inside Docker, where detectCores() sees the host's cores
ncores           <- as.integer(Sys.getenv("UKBB_NCORES",
                                          unset = as.character(min(8, max(1, parallel::detectCores() - 2)))))
set.seed(123)

# Protein sets
known_marker <- cfg$known_marker
signature    <- cfg$signature

# Covariate names as they appear in All_clinic after processing below
covars_num <- cfg$covariates_numeric
covars_fac <- cfg$covariates_factor
all_covars <- c(covars_num, covars_fac)

model_specs <- list(
  Signature            = list(proteins = signature,    use_covariates = FALSE),
  Covariates           = list(proteins = character(0), use_covariates = TRUE),
  Signature_Covariates = list(proteins = signature,    use_covariates = TRUE)
)
model_specs[[cfg$known_marker_model]] <- list(proteins = known_marker, use_covariates = FALSE)
# Null comparator: a new random panel of non-signature proteins in every iteration
model_specs$Random            <- list(proteins = character(0), use_covariates = FALSE,
                                      random = TRUE, random_n = length(signature))
model_specs$Random_Covariates <- list(proteins = character(0), use_covariates = TRUE,
                                      random = TRUE, random_n = length(signature))

#### --- 2. DATA LOADING --- ####

# Proteomics
if (!exists("UKBB_Proteomics")) {
  UKBB_Proteomics <- read.table(
    file        = file.path(input.dir, "Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt"),
    sep         = "\t",
    header      = TRUE,
    fill        = TRUE,
    row.names   = "PID",
    check.names = FALSE
  )
}

# Filter proteins: keep columns with < 25% missingness
UKBB_Proteomics <- UKBB_Proteomics[
  , colSums(is.na(UKBB_Proteomics)) <= 0.25 * nrow(UKBB_Proteomics),
  drop = FALSE
]
UKBB_Proteomics$eid <- as.character(rownames(UKBB_Proteomics))

# --- Map signature / marker gene symbols to Olink assay names --- #
# Single-cell signature genes may carry outdated symbols relative to the Olink
# assay names (e.g. ACPP -> ACP3, H2AFZ -> H2AZ1). Both sides are mapped to
# current HGNC approved symbols (HGNChelper) and matched on those. Olink assays
# covering several genes (names joined by "_") are matched on any member gene.
hgnc_map <- tryCatch(HGNChelper::getCurrentHumanMap(), error = function(e) {
  message("HGNChelper: could not download current map, using bundled map (",
          conditionMessage(e), ")")
  NULL
})
to_approved <- function(x) {
  chk <- if (is.null(hgnc_map)) {
    HGNChelper::checkGeneSymbols(x, unmapped.as.na = FALSE)
  } else {
    HGNChelper::checkGeneSymbols(x, map = hgnc_map, unmapped.as.na = FALSE)
  }
  lapply(strsplit(chk$Suggested.Symbol, " /// ", fixed = TRUE), trimws)
}
match_to_olink <- function(genes, olink_cols) {
  col_parts <- strsplit(olink_cols, "_", fixed = TRUE)
  part_vec  <- unlist(col_parts)
  part_col  <- rep(olink_cols, lengths(col_parts))
  part_appr <- to_approved(part_vec)
  gene_appr <- to_approved(genes)
  res <- vapply(seq_along(genes), function(i) {
    if (genes[i] %in% olink_cols) return(genes[i])
    cand <- unique(c(genes[i], gene_appr[[i]]))
    hit  <- unique(part_col[vapply(seq_along(part_vec), function(j)
      part_vec[j] %in% cand || any(part_appr[[j]] %in% cand), logical(1))])
    if (length(hit) > 1) {
      message(sprintf("  %s matches several Olink assays (%s) -- left unmatched",
                      genes[i], paste(hit, collapse = ", ")))
      return(NA_character_)
    }
    if (length(hit) == 1) hit else NA_character_
  }, character(1))
  setNames(res, genes)
}

olink_cols         <- setdiff(colnames(UKBB_Proteomics), "eid")
signature_input    <- signature
known_marker_input <- known_marker
sig_map    <- match_to_olink(signature_input,    olink_cols)
marker_map <- match_to_olink(known_marker_input, olink_cols)

signature_mapping <- tibble::tibble(
  set         = c(rep("signature", length(sig_map)), rep("marker", length(marker_map))),
  input_gene  = c(names(sig_map), names(marker_map)),
  olink_assay = c(unname(sig_map), unname(marker_map))
)
message("Signature/marker -> Olink assay mapping (NA = not assayed or removed by the 25% missingness filter):")
print(as.data.frame(signature_mapping))
write.csv(signature_mapping, file.path(out.dir, "signature_olink_mapping.csv"), row.names = FALSE)
message(sprintf("Signature genes measured: %d of %d", sum(!is.na(sig_map)), length(sig_map)))
if (Sys.getenv("PREFLIGHT_ONLY") == "1") {
  message("PREFLIGHT_ONLY=1: stopping after the mapping table")
  quit(save = "no")
}

signature                <- unique(stats::na.omit(unname(sig_map)))
known_marker             <- unique(stats::na.omit(unname(marker_map)))
exist_signature_proteins <- signature
# Random pool excludes every assay matched to a signature gene
all_other_proteins       <- setdiff(olink_cols, c(signature, signature_input))

# Point the model specs at the matched assay names; random panels are the same
# size as the MEASURED signature panel (not the input gene list)
for (m in names(model_specs)) {
  if (isTRUE(model_specs[[m]]$random)) {
    model_specs[[m]]$random_n <- length(exist_signature_proteins)
  } else if (identical(model_specs[[m]]$proteins, signature_input)) {
    model_specs[[m]]$proteins <- signature
  } else if (identical(model_specs[[m]]$proteins, known_marker_input)) {
    model_specs[[m]]$proteins <- known_marker
  }
}

# Core clinical data
All_clinic <- readRDS(file.path(input.dir, "ukb_clinic.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
  dplyr::mutate(
    dplyr::across(
      c(date_of_attending_assessment_centre_f53_0_0,
        date_of_death_f40000_0_0,
        date_of_death_f40000_1_0,
        date_lost_to_followup_f191_0_0),
      as.Date
    ),
    # Unordered factor: the raw field is ordered, which model.matrix would
    # encode as a polynomial contrast ("sex_f31_0_0.L") instead of a dummy
    sex_f31_0_0                    = factor(as.character(sex_f31_0_0)),
    ethnic_background_f21000_0_0 = factor(
      dplyr::case_when(
        trimws(as.character(ethnic_background_f21000_0_0)) %in% c(
          "British", "Irish", "White", "Any other white background"
        ) ~ "White",
        trimws(as.character(ethnic_background_f21000_0_0)) %in% c(
          "Indian", "Pakistani", "Bangladeshi", "Any other Asian background", "Chinese"
        ) ~ "Asian",
        trimws(as.character(ethnic_background_f21000_0_0)) %in% c(
          "Caribbean", "African", "Any other Black background", "Black or Black British"
        ) ~ "Black",
        trimws(as.character(ethnic_background_f21000_0_0)) %in% c(
          "White and Black Caribbean", "White and Black African",
          "White and Asian", "Any other mixed background", "Mixed"
        ) ~ "Mixed/Other",
        trimws(as.character(ethnic_background_f21000_0_0)) %in% c(
          "Other ethnic group"
        ) ~ "Mixed/Other",
        trimws(as.character(ethnic_background_f21000_0_0)) %in% c(
          "Prefer not to answer", "Do not know"
        ) ~ NA_character_
      ),
      levels = c("White", "Asian", "Black", "Mixed/Other")
    ),

    # Covariates with reference levels set here (will be relevelled again in model)
    smoking_status_f20116_0_0 = stats::relevel(
      factor(dplyr::na_if(trimws(as.character(smoking_status_f20116_0_0)), "Prefer not to answer")),
      ref = "Never"
    ),
    alcohol_intake_frequency_f1558_0_0 = stats::relevel(
      factor(dplyr::na_if(trimws(as.character(alcohol_intake_frequency_f1558_0_0)), "Prefer not to answer")),
      ref = "Never"
    ),
    diabetes_diagnosed_by_doctor_f2443_0_0 = stats::relevel(
      factor(dplyr::na_if(
        dplyr::na_if(trimws(as.character(diabetes_diagnosed_by_doctor_f2443_0_0)),
                     "Prefer not to answer"), "Do not know")),
      ref = "No"
    ),
    
    OSstatus.raw = ifelse(is.na(date_of_death_f40000_0_0) & is.na(date_of_death_f40000_1_0), 0, 1),
    OStime = ifelse(
      OSstatus.raw == 1,
      (date_of_death_f40000_0_0 - date_of_attending_assessment_centre_f53_0_0) / 365,
      ifelse(
        is.na(date_lost_to_followup_f191_0_0),
        (as.Date("2022-12-30") - date_of_attending_assessment_centre_f53_0_0) / 365,
        (pmin(date_lost_to_followup_f191_0_0, as.Date("2022-12-30")) -
           date_of_attending_assessment_centre_f53_0_0) / 365
      )
    ),
    age = as.numeric(age_at_recruitment_f21022_0_0)
  ) %>%
  dplyr::filter(
    !is.na(ethnic_background_f21000_0_0),
    !is.na(age_at_recruitment_f21022_0_0),
    !is.na(uk_biobank_assessment_centre_f54_0_0),
    !is.na(sex_f31_0_0)
  ) %>%
  dplyr::mutate(eid = as.character(eid))

# Disease-specific clinical data
load_disease_data <- function(dis_name) {
  readRDS(file.path(input.dir2, paste0("ukb_", dis_name, "_subset_clinic.rda"))) %>%
    dplyr::filter(eid %in% rownames(UKBB_Proteomics)) %>%
    dplyr::mutate(eid = as.character(eid), diagnosis_time = as.Date(diagnosis_time))
}

dis1_clinic <- load_disease_data(dis1)
dis2_clinic <- load_disease_data(dis2)

dis1_eid <- dis1_clinic %>% dplyr::pull(eid)
dis2_eid <- dis2_clinic %>% dplyr::pull(eid)

# Prevalent / incident splits
extract_cases <- function(clinic_df, prevalent = TRUE) {
  tmp <- clinic_df %>%
    dplyr::select(-dplyr::any_of("date_of_attending_assessment_centre_f53_0_0")) %>%
    dplyr::left_join(
      All_clinic %>% dplyr::select(eid, date_of_attending_assessment_centre_f53_0_0),
      by = "eid"
    )
  if (prevalent) dplyr::filter(tmp, diagnosis_time <= date_of_attending_assessment_centre_f53_0_0)
  else           dplyr::filter(tmp, diagnosis_time >  date_of_attending_assessment_centre_f53_0_0)
}

dis2_prev_clinic <- extract_cases(dis2_clinic, prevalent = TRUE)
dis2_inci_clinic <- extract_cases(dis2_clinic, prevalent = FALSE)
dis2_prev_eid    <- dis2_prev_clinic %>% dplyr::pull(eid)
dis2_inci_eid    <- dis2_inci_clinic %>% dplyr::pull(eid)

dis1_prev_clinic <- extract_cases(dis1_clinic, prevalent = TRUE)
dis1_prev_eid    <- dis1_prev_clinic %>% dplyr::pull(eid)

# Healthy Controls
HC_clinic <- readRDS(file.path(input.dir2, "ukb_healthy_controls.rda")) %>%
  dplyr::filter(eid %in% rownames(UKBB_Proteomics))
HC_clinic <- All_clinic %>% dplyr::filter(eid %in% as.character(HC_clinic$eid))
HC_eid    <- HC_clinic %>% dplyr::pull(eid)

All_cancer_clinic_eid = readRDS(paste0(input.dir2,'/ukb_AllCancers_subset_clinic.rda')) %>% 
  filter(eid %in% rownames(UKBB_Proteomics)) %>% pull(eid)

non_cancer_eid = All_clinic %>%
  filter(!eid %in% All_cancer_clinic_eid & !eid %in% dis2_eid) %>% 
  pull(eid)

NonCancer_clinic <- All_clinic %>% dplyr::filter(eid %in% non_cancer_eid)
NonCancer_eid    <- NonCancer_clinic %>% dplyr::pull(eid)

# Progression subgroups (Dis1 -> Dis2 vs Dis1 -> no Dis2)
All_cancer_clinic_2 <- All_clinic %>%
  dplyr::left_join(dis1_clinic %>% dplyr::select(eid, diagnosis_time.dis1 = diagnosis_time), by = "eid") %>%
  dplyr::left_join(dis2_clinic %>% dplyr::select(eid, diagnosis_time.dis2 = diagnosis_time), by = "eid") %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("diagnosis_time"), as.Date))

Dis1andDis2_clinic <- All_cancer_clinic_2 %>%
  dplyr::mutate(is_dis2 = ifelse(eid %in% dis2_eid, 1, 0)) %>%
  dplyr::filter(!is.na(diagnosis_time.dis1),
                is.na(diagnosis_time.dis2) | diagnosis_time.dis1 <= diagnosis_time.dis2,
                is_dis2 == 1,
                # Incident Dis2 only (diagnosed after blood draw), as in the HC/NonCancer contrasts
                !eid %in% dis2_prev_eid) %>%
  dplyr::rename(diagnosis_time = diagnosis_time.dis2)
Dis1andDis2_eid <- Dis1andDis2_clinic %>% dplyr::pull(eid)

Dis1noDis2_clinic <- All_cancer_clinic_2 %>%
  dplyr::mutate(is_dis2 = ifelse(eid %in% dis2_eid, 1, 0)) %>%
  dplyr::filter(!is.na(diagnosis_time.dis1), is_dis2 == 0) %>%
  dplyr::rename(diagnosis_time = diagnosis_time.dis1)
Dis1noDis2_eid <- Dis1noDis2_clinic %>% dplyr::pull(eid)


#### --- 3. HELPER: BUILD FEATURE MATRIX FOR GLMNET --- ####

# Proteins and numeric covariates as columns; factor covariates one-hot encoded
# (reference level dropped). na.pass keeps every row; NA dummies become 0 (reference).

build_feature_matrix <- function(df, proteins, use_covariates,
                                 covars_num, covars_fac) {
  n <- nrow(df)
  
  # --- Protein block ---
  if (length(proteins) > 0) {
    prot_cols <- intersect(proteins, colnames(df))
    prot_mat  <- as.matrix(df[, prot_cols, drop = FALSE])
    mode(prot_mat) <- "numeric"
    if (nrow(prot_mat) != n || ncol(prot_mat) != length(prot_cols)) {
      prot_mat <- matrix(as.numeric(prot_mat), nrow = n, ncol = length(prot_cols),
                         dimnames = list(NULL, prot_cols))
    }
    stopifnot(nrow(prot_mat) == n)
  } else {
    prot_mat <- NULL
  }
  
  if (use_covariates) {
    num_present <- intersect(covars_num, colnames(df))
    if (length(num_present) > 0) {
      num_mat <- as.matrix(df[, num_present, drop = FALSE])
      mode(num_mat) <- "numeric"
      stopifnot(nrow(num_mat) == n)
    } else {
      num_mat <- NULL
    }
    
    fac_present <- intersect(covars_fac, colnames(df))
    
    if (length(fac_present) > 0) {
      fac_df <- df[, fac_present, drop = FALSE]
      
      fac_df <- as.data.frame(
        lapply(fac_df, function(x) {
          x <- as.factor(x)
          if (anyNA(x)) levels(x) <- c(levels(x), NA)
          x
        }),
        stringsAsFactors = FALSE
      )
      
      fac_formula <- as.formula(paste("~", paste(fac_present, collapse = " + ")))
      fac_mat <- tryCatch(
        model.matrix(
          fac_formula,
          data = fac_df,
          na.action = na.pass
        ),
        error = function(e) {
          matrix(0, nrow = n, ncol = 0)
        }
      )
      
      fac_mat <- fac_mat[, !colnames(fac_mat) %in% "(Intercept)", drop = FALSE]
      fac_mat[is.na(fac_mat)] <- 0
      na_level_cols <- grepl("NA$", colnames(fac_mat))
      if (any(na_level_cols)) fac_mat <- fac_mat[, !na_level_cols, drop = FALSE]
      
      stopifnot(nrow(fac_mat) == n)
    } else {
      fac_mat <- NULL
    }
    
    cov_mat <- do.call(cbind, Filter(Negate(is.null), list(num_mat, fac_mat)))
  } else {
    cov_mat <- NULL
  }
  
  feat_mat <- do.call(cbind, Filter(Negate(is.null), list(prot_mat, cov_mat)))
  if (is.null(feat_mat) || ncol(feat_mat) == 0) return(NULL)
  stopifnot(nrow(feat_mat) == n)
  feat_mat
}


# --- 3b. HELPER: ASSIGN SUBGROUP LABEL (priority-ordered) --- #
assign_subgroup <- function(eids) {
  dplyr::case_when(
    eids %in% Dis1andDis2_eid  ~ "Dis1andDis2",
    eids %in% dis2_inci_eid    ~ "dis2_inci_only",
    eids %in% dis2_prev_eid    ~ "dis2_prev",
    eids %in% Dis1noDis2_eid   ~ "Dis1noDis2",
    eids %in% NonCancer_eid  ~ "NonCancer",
    eids %in% HC_eid           ~ "HC",
    TRUE                       ~ "Other"
  )
}


# --- 3c. HELPER: FIT / PREDICT --- #
fit_risk_model <- function(x_train, y_train, alpha, nfolds = 10) {
  
  # One predictor (the known-marker model): plain logistic regression
  if (ncol(x_train) == 1) {
    train_df <- data.frame(y = y_train, x = as.numeric(x_train[, 1]))
    fit <- tryCatch(glm(y ~ x, data = train_df, family = binomial),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    return(list(type = "glm", fit = fit, alpha = NA_real_, lambda = NA_real_))
  }
  
  cv_fit <- tryCatch(
    glmnet::cv.glmnet(
      x = x_train, y = y_train, family = "binomial",
      alpha = alpha, type.measure = "auc",
      standardize = TRUE, nfolds = nfolds
    ),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) return(NULL)
  
  list(type = "glmnet", fit = cv_fit, alpha = alpha, lambda = cv_fit$lambda.min)
}

predict_risk_model <- function(model_obj, newx) {
  if (is.null(model_obj)) return(rep(NA_real_, nrow(newx)))
  
  if (model_obj$type == "glm") {
    newdf <- data.frame(x = as.numeric(newx[, 1]))
    return(tryCatch(
      predict(model_obj$fit, newdata = newdf, type = "response"),
      error = function(e) rep(NA_real_, nrow(newx))
    ))
  }
  
  tryCatch(
    as.numeric(predict(model_obj$fit, newx = newx, s = model_obj$lambda, type = "response")),
    error = function(e) rep(NA_real_, nrow(newx))
  )
}

# Ridge coefficients for one fitted glmnet model, plus a standardized version
# (coef x training SD). glmnet reports coefficients on the original feature
# scale, so raw values aren't comparable across features (years vs 0/1 vs NPX).
# Returns NULL for single-feature glm fits, where importance is trivial.
extract_coefs <- function(model_fit, x_train) {
  if (is.null(model_fit) || !identical(model_fit$type, "glmnet")) return(NULL)
  cf <- as.matrix(stats::coef(model_fit$fit, s = "lambda.min"))
  cf <- cf[rownames(cf) != "(Intercept)", 1]
  feats <- intersect(names(cf), colnames(x_train))
  if (length(feats) == 0) return(NULL)
  sds <- apply(x_train[, feats, drop = FALSE], 2, stats::sd)
  tibble::tibble(
    feature  = feats,
    coef     = as.numeric(cf[feats]),
    sd_train = as.numeric(sds),
    std_coef = as.numeric(cf[feats]) * as.numeric(sds)
  )
}

# Drop factor covariates that are constant within an analysis set (e.g. sex for
# prostate): a constant column is collinear with glmnet's intercept.
drop_constant_factors <- function(df, fac_cols) {
  keep <- vapply(fac_cols, function(col) {
    dplyr::n_distinct(df[[col]], na.rm = TRUE) > 1
  }, logical(1))
  fac_cols[keep]
}

#### --- 4. MAIN ANALYSIS FUNCTION --- ####
analyze_auc_multiyear <- function(
    model_name,
    case_group    = "dis2_inci",
    control_group = "HC",
    n_obs_yrs_vec = observation_years,
    match_ratio   = 5,
    n_bootstrap   = 1000,
    alpha         = 0) {
  
  spec         <- model_specs[[model_name]]
  is_random    <- isTRUE(spec$random)
  random_n     <- if (is_random) spec$random_n else NA_integer_
  proteins     <- if (is_random) character(0) else intersect(spec$proteins, colnames(UKBB_Proteomics))
  use_cov      <- spec$use_covariates
  
  purrr::map_dfr(n_obs_yrs_vec, function(n_obs_yrs) {
    message(paste0("  Model: ", model_name, " | Case: ", case_group,
                   " | Control: ", control_group, " | Years: ", n_obs_yrs))
    
    # ---- 4.1 Case data ----
    case_eid_var    <- paste0(case_group, "_eid")
    case_clinic_var <- paste0(case_group, "_clinic")
    
    case_data <- All_clinic %>%
      dplyr::filter(eid %in% get(case_eid_var)) %>%
      dplyr::left_join(UKBB_Proteomics, by = "eid") %>%
      dplyr::left_join(
        get(case_clinic_var) %>% dplyr::select(eid, diagnosis_time),
        by = "eid"
      ) %>%
      dplyr::filter(
        abs(as.numeric(diagnosis_time - date_of_attending_assessment_centre_f53_0_0)) <=
          n_obs_yrs * 365
      ) %>%
      dplyr::mutate(
        binary_group = 1L,
        sex_num      = ifelse(sex_f31_0_0 == "Female", 0, 1)
      )
    
    if (nrow(case_data) < 20) {
      message(sprintf("    Skipping: only %d cases", nrow(case_data)))
      return(NULL)
    }
    
    # ---- 4.2 Control data ----
    control_clinic <- switch(control_group,
                             HC         = HC_clinic,
                             Dis1noDis2 = Dis1noDis2_clinic,
                             NonCancer  = NonCancer_clinic
    )
    if(is.null(control_clinic)) stop(paste("Unknown control group:", control_group))
    
    control_data <- All_clinic %>%
      dplyr::filter(eid %in% control_clinic$eid) %>%
      # Exclude controls who died before the observation window closes
      dplyr::filter(OSstatus.raw != 1 | (OStime * 365) > n_obs_yrs * 365) %>%
      dplyr::left_join(UKBB_Proteomics, by = "eid") %>%
      dplyr::mutate(
        binary_group = 0L,
        sex_num      = ifelse(sex_f31_0_0 == "Female", 0, 1)
      )
    
    # ---- 4.3 Age/sex matching ----
    pool <- dplyr::bind_rows(case_data, control_data) %>%
      dplyr::mutate(binary_group = as.integer(binary_group))
    
    # Age caliper in absolute years (not SD units). HC: 1:1 within +-1 y (see header);
    # NonCancer / Dis1noDis2: match_ratio within +-2 y.
    is_hc <- control_group == "HC"
    match_obj <- tryCatch(
      MatchIt::matchit(
        binary_group ~ age,
        data        = pool,
        method      = "nearest",
        ratio       = if (is_hc) 1 else match_ratio,
        replace     = FALSE,
        exact       = ~ sex_num,
        caliper     = c(age = if (is_hc) 1 else 2),
        std.caliper = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(match_obj)) {
      message("Matching failed -- skipping window")
      return(NULL)
    }
    
    full_data <- MatchIt::match.data(match_obj) %>%
      dplyr::mutate(binary_group = factor(binary_group, levels = c(0, 1)))
    
    if (nrow(full_data) < 30) return(NULL)

    # Balance diagnostics for this window (logged with every bootstrap row)
    n_case_matched <- sum(full_data$binary_group == 1)
    age_smd <- (mean(full_data$age[full_data$binary_group == 1]) -
                mean(full_data$age[full_data$binary_group == 0])) / sd(pool$age)
    
    covars_fac_run <- drop_constant_factors(full_data, intersect(covars_fac, colnames(full_data)))
    # ---- 4.4 Parallel bootstrap ----
    # No .combine: foreach's default returns a flat list of exactly n_bootstrap results.
    results_boot_raw <- foreach(
      i         = 1:n_bootstrap,
      .errorhandling = "pass",
      .packages = c("dplyr", "glmnet", "pROC", "tibble")
    ) %dopar% {
      
      set.seed(1000 + i)
      
      # Random models draw a new protein panel in every iteration
      proteins_iter <- if (is_random) sample(all_other_proteins, random_n) else proteins
      
      # 70/30 stratified split
      case_idx    <- which(full_data$binary_group == 1)
      control_idx <- which(full_data$binary_group == 0)
      train_idx   <- c(
        sample(case_idx,    round(0.7 * length(case_idx))),
        sample(control_idx, round(0.7 * length(control_idx)))
      )
      test_idx  <- setdiff(seq_len(nrow(full_data)), train_idx)
      
      train_data <- full_data[train_idx, , drop = FALSE]
      test_data  <- full_data[test_idx,  , drop = FALSE]
      
      if (sum(train_data$binary_group == 1) < 5 ||
          sum(test_data$binary_group  == 1) < 3)
        return(list(metrics = NULL, preds = NULL, reason = "guard_small_n"))
      
      # ---- kNN imputation on protein columns (training only) ----
      protein_cols_avail <- intersect(proteins_iter, colnames(full_data))
      
      if (length(protein_cols_avail) > 0) {
        prot_train_mat <- apply(
          as.matrix(train_data[, protein_cols_avail, drop = FALSE]), 2, as.numeric
        )
        if (is.null(dim(prot_train_mat)))
          prot_train_mat <- matrix(prot_train_mat, ncol = 1,
                                   dimnames = list(NULL, protein_cols_avail))
        
        imputed_train <- tryCatch(
          impute::impute.knn(prot_train_mat)$data,
          error = function(e) prot_train_mat
        )
        train_data[, protein_cols_avail] <- imputed_train
        
        # Mean-substitute test proteins using training means
        train_prot_means <- colMeans(imputed_train, na.rm = TRUE)
        for (col in protein_cols_avail) {
          test_data[[col]][is.na(test_data[[col]])] <- train_prot_means[[col]]
        }
      }
      
      # Mean-substitute missing numeric covariates using training means
      for (col in intersect(covars_num, colnames(train_data))) {
        train_mean <- mean(train_data[[col]], na.rm = TRUE)
        if (is.na(train_mean)) train_mean <- 0
        train_data[[col]][is.na(train_data[[col]])] <- train_mean
        test_data[[col]][is.na(test_data[[col]])]   <- train_mean
      }
      
      # Mode-substitute missing factor covariates using training mode
      for (col in intersect(covars_fac, colnames(train_data))) {
        train_vals <- train_data[[col]][!is.na(train_data[[col]])]
        if (length(train_vals) == 0) next
        train_mode <- names(sort(table(train_vals), decreasing = TRUE))[1]
        train_data[[col]][is.na(train_data[[col]])] <- train_mode
        test_data[[col]][is.na(test_data[[col]])]   <- train_mode
      }
      
      # ---- Build feature matrices ----
      x_train <- build_feature_matrix(
        train_data, proteins_iter, use_cov, covars_num, covars_fac_run
      )
      x_test  <- build_feature_matrix(
        test_data, proteins_iter, use_cov, covars_num, covars_fac_run
      )
      
      if (is.null(x_train) || is.null(x_test))
        return(list(metrics = NULL, preds = NULL, reason = "guard_null_feature_matrix"))
      
      # Align columns (model.matrix may differ between train/test if a level is absent)
      shared_cols <- intersect(colnames(x_train), colnames(x_test))
      if (length(shared_cols) == 0)
        return(list(metrics = NULL, preds = NULL, reason = "guard_no_shared_cols"))
      x_train <- x_train[, shared_cols, drop = FALSE]
      x_test  <- x_test[,  shared_cols, drop = FALSE]
      
      mode(x_train) <- "numeric"
      mode(x_test)  <- "numeric"
      
      y_train <- as.numeric(train_data$binary_group) - 1  # 0/1
      
      # ---- Ridge regression (logistic regression for one predictor) ----
      model_fit <- fit_risk_model(x_train, y_train, alpha = alpha, nfolds = 10)
      if (is.null(model_fit))
        return(list(metrics = NULL, preds = NULL, reason = "guard_model_fit_null"))
      pred_prob <- predict_risk_model(model_fit, x_test)
      
      y_test <- as.numeric(test_data$binary_group) - 1
      
      # ---- AUC ----
      auc_val <- tryCatch({
        roc_obj <- pROC::roc(
          response  = y_test,
          predictor = pred_prob,
          quiet     = TRUE,
          direction = "<"
        )
        as.numeric(pROC::auc(roc_obj))
      }, error = function(e) NA_real_)
      
      # ---- Per-eid predictions ----
      pred_df <- tibble::tibble(
        eid               = test_data$eid,
        pred_prob         = pred_prob,
        y_true            = y_test,
        Bootstrap         = i,
        Model             = model_name,
        Case              = case_group,
        Control_Group     = control_group,
        Observation_Years = n_obs_yrs
      )
      
      # Coefficients for fixed panels only (random panels change every iteration)
      coefs <- if (is_random) NULL else extract_coefs(model_fit, x_train)
      if (!is.null(coefs)) coefs$Bootstrap <- i
      
      list(
        metrics = tibble::tibble(
          Model              = model_name,
          Case               = case_group,
          Control_Group      = control_group,
          Observation_Years  = n_obs_yrs,
          Bootstrap          = i,
          N_case_window      = nrow(case_data),
          N_matched_total    = nrow(full_data),
          N_case_matched     = n_case_matched,
          N_control_matched  = nrow(full_data) - n_case_matched,
          Age_SMD            = age_smd,
          N_proteins_used    = length(protein_cols_avail),
          AUC                = auc_val
        ),
        preds = pred_df,
        coefs = coefs
      )
    }
    
    # Tally failed draws before dropping them: every guard in the loop returns a tagged
    # list(metrics = NULL, reason = ...), and an errored draw comes back as a condition object.
    status <- vapply(results_boot_raw, function(x) {
      if (inherits(x, "error") || inherits(x, "try-error")) "worker_error"
      else if (is.null(x)) "unexpected_null"
      else if (is.null(x$metrics)) paste0("skipped:", x$reason %||% "unknown")
      else "ok"
    }, character(1))
    n_fail <- sum(status != "ok")
    if (n_fail > 0) {
      fail_tab <- table(status[status != "ok"])
      message(sprintf(
        "    [%s | %s vs %s | %dyr] %d/%d bootstrap draws failed (%s)",
        model_name, case_group, control_group, n_obs_yrs, n_fail, n_bootstrap,
        paste(names(fail_tab), fail_tab, sep = "=", collapse = ", ")
      ))
    }
    # Keep the successful draws
    results_boot_raw <- Filter(function(x) !is.null(x) && !inherits(x, "error") &&
                                 !inherits(x, "try-error") && !is.null(x$metrics),
                               results_boot_raw)
    
    auc_rows  <- dplyr::bind_rows(lapply(results_boot_raw, `[[`, "metrics"))
    pred_rows <- dplyr::bind_rows(lapply(results_boot_raw, `[[`, "preds"))
    
    # ---- Per-bootstrap feature importance (fixed-panel models only) ----
    coef_rows <- dplyr::bind_rows(lapply(results_boot_raw, `[[`, "coefs"))
    if (nrow(coef_rows) > 0) {
      coef_rows <- coef_rows %>%
        dplyr::mutate(Model = model_name, Case = case_group,
                      Control_Group = control_group, Observation_Years = n_obs_yrs)
      dir.create(dir.importance, showWarnings = FALSE, recursive = TRUE)
      saveRDS(
        coef_rows,
        file.path(dir.importance, sprintf("coefs_boot_%s_%s_vs_%s_%dyr.rds",
                                          model_name, case_group, control_group, n_obs_yrs))
      )
    }
    
    # ---- Per-eid summary across bootstraps ----
    eid_summary <- pred_rows %>%
      dplyr::group_by(eid, Model, Case, Control_Group, Observation_Years) %>%
      dplyr::summarise(
        mean_pred_prob  = mean(pred_prob, na.rm = TRUE),
        sd_pred_prob    = sd(pred_prob,   na.rm = TRUE),
        n_appearances   = sum(!is.na(pred_prob)),
        y_true          = dplyr::first(y_true),
        .groups         = "drop"
      ) %>%
      dplyr::mutate(subgroup = assign_subgroup(eid))
    
    # Save per observation window
    saveRDS(eid_summary, file.path(dir.predictions, sprintf(
      "pred_prob_eid_summary_%s_%s_vs_%s_%dyr.rds",
      model_name, case_group, control_group, n_obs_yrs)))
    
    compute_threshold_stats <- function(eid_summary, obs_years = n_obs_yrs) {
      df <- eid_summary %>%
        dplyr::filter(Observation_Years == obs_years, !is.na(mean_pred_prob), !is.na(y_true))
      
      if (nrow(df) == 0 || length(unique(df$y_true)) < 2) {
        return(list(threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_,
                    accuracy = NA_real_, ppv = NA_real_, npv = NA_real_,
                    auc = NA_real_, n_eids = nrow(df)))
      }
      
      roc_obj <- pROC::roc(
        response  = df$y_true,
        predictor = df$mean_pred_prob,
        quiet     = TRUE,
        direction = "<"
      )
      
      # Youden's J: maximizes sensitivity + specificity - 1
      best_coords <- pROC::coords(
        roc_obj, x = "best", best.method = "youden",
        ret = c("threshold", "sensitivity", "specificity", "accuracy", "ppv", "npv")
      )
      
      # Ties at the Youden-optimal point return several rows: take the first
      if (is.matrix(best_coords) || is.data.frame(best_coords)) {
        if (nrow(best_coords) > 1) best_coords <- best_coords[1, , drop = FALSE]
        best_coords <- unlist(best_coords)
      }
      
      list(
        threshold   = as.numeric(best_coords["threshold"]),
        sensitivity = as.numeric(best_coords["sensitivity"]),
        specificity = as.numeric(best_coords["specificity"]),
        accuracy    = as.numeric(best_coords["accuracy"]),
        ppv         = as.numeric(best_coords["ppv"]),
        npv         = as.numeric(best_coords["npv"]),
        auc         = as.numeric(pROC::auc(roc_obj)),
        n_eids      = nrow(df)
      )
    }
    
    perf_stats <- compute_threshold_stats(
      eid_summary,
      obs_years = n_obs_yrs
    )
    threshold_row <- tibble::tibble(
      Model              = model_name,
      Case               = case_group,
      Control_Group      = control_group,
      Observation_Years  = n_obs_yrs,
      threshold          = perf_stats$threshold,
      sensitivity        = perf_stats$sensitivity,
      specificity        = perf_stats$specificity,
      accuracy           = perf_stats$accuracy,
      ppv                = perf_stats$ppv,
      npv                = perf_stats$npv,
      auc_from_mean_pred = perf_stats$auc,
      n_eids             = perf_stats$n_eids
    )
    
    saveRDS(
      threshold_row,
      file.path(dir.thresholds, sprintf("threshold_stats_%s_%s_vs_%s_%dyr.rds",
                                        model_name, case_group, control_group, n_obs_yrs))
    )
    
    auc_rows  # return for map_dfr
  })
}

#### --- 5. RUNS LIST --- ####
# Three contrasts, each run for all six models:
#   - incident cases vs healthy controls (HC)
#   - incident cases vs non-cancer controls (NonCancer)
#   - progression (Dis1andDis2) vs non-progression (Dis1noDis2)
contrasts <- list(c("dis2_inci",   "HC"),
                  c("dis2_inci",   "NonCancer"),
                  c("Dis1andDis2", "Dis1noDis2"))
# Final 14-year models are saved for every model except the random panels
runs_list <- unlist(lapply(contrasts, function(cc) lapply(names(model_specs), function(m)
  list(model_name = m, case_group = cc[1], control_group = cc[2],
       save_final_model = !isTRUE(model_specs[[m]]$random)))), recursive = FALSE)

#### --- 6. EXECUTE --- ####

registerDoParallel(cores = ncores)

start_time <- Sys.time()

for (i in seq_along(runs_list)) {
  run <- runs_list[[i]]
  message(sprintf("\n[Run %d/%d] Model: %s | Case: %s | Control: %s",
                  i, length(runs_list),
                  run$model_name, run$case_group, run$control_group))
  
  res <- tryCatch(
    analyze_auc_multiyear(
      model_name    = run$model_name,
      case_group    = run$case_group,
      control_group = run$control_group,
      n_obs_yrs_vec = observation_years,
      match_ratio   = match_ratio,
      n_bootstrap   = n_bootstrap,
      alpha         = alpha
    ),
    error = function(e) {
      message(sprintf("  !! [Run %d/%d] %s FAILED: %s -- skipping, continuing to next run",
                      i, length(runs_list), run$model_name, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) next
  
  rds_file <- file.path(dir.bootstrap, sprintf(
    "%s_%s_vs_%s_run%d_bt%d.rds",
    run$model_name, run$case_group, run$control_group, i, n_bootstrap))
  saveRDS(res, rds_file)
  
  if (isTRUE(run$save_final_model)) {
    message(sprintf("  Fitting final saved model for %s (14yr window)...", run$model_name))
    
    spec_final   <- model_specs[[run$model_name]]
    proteins_fin <- intersect(spec_final$proteins, colnames(UKBB_Proteomics))
    use_cov_fin  <- spec_final$use_covariates
    
    # Rebuild matched dataset for 14yr window
    case_data_fin <- All_clinic %>%
      dplyr::filter(eid %in% get(paste0(run$case_group, "_eid"))) %>%
      dplyr::left_join(UKBB_Proteomics, by = "eid") %>%
      dplyr::left_join(
        get(paste0(run$case_group, "_clinic")) %>% dplyr::select(eid, diagnosis_time),
        by = "eid"
      ) %>%
      dplyr::filter(
        abs(as.numeric(diagnosis_time - date_of_attending_assessment_centre_f53_0_0)) <= 14 * 365
      ) %>%
      dplyr::mutate(binary_group = 1L,
                    sex_num      = ifelse(sex_f31_0_0 == "Female", 0, 1))
    
    control_clinic_fin <- switch(run$control_group,
                                 HC         = HC_clinic,
                                 Dis1noDis2 = Dis1noDis2_clinic,
                                 NonCancer  = NonCancer_clinic
    )
    
    control_data_fin <- All_clinic %>%
      dplyr::filter(eid %in% control_clinic_fin$eid) %>%
      dplyr::filter(OSstatus.raw != 1 | (OStime * 365) > 14 * 365) %>%
      dplyr::left_join(UKBB_Proteomics, by = "eid") %>%
      dplyr::mutate(binary_group = 0L,
                    sex_num      = ifelse(sex_f31_0_0 == "Female", 0, 1))
    
    pool_fin <- dplyr::bind_rows(case_data_fin, control_data_fin) %>%
      dplyr::mutate(binary_group = as.integer(binary_group))
    
    # Same matching as in the bootstrap
    is_hc_fin <- run$control_group == "HC"
    match_fin <- tryCatch(
      MatchIt::matchit(
        binary_group ~ age,
        data        = pool_fin,
        method      = "nearest",
        ratio       = if (is_hc_fin) 1 else match_ratio,
        replace     = FALSE,
        exact       = ~ sex_num,
        caliper     = c(age = if (is_hc_fin) 1 else 2),
        std.caliper = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(match_fin)) {
      set.seed(42)
      full_fin <- MatchIt::match.data(match_fin) %>%
        dplyr::mutate(binary_group = factor(binary_group, levels = c(0, 1)))
      
      covars_fac_fin_run <- drop_constant_factors(full_fin, intersect(covars_fac, colnames(full_fin)))
      
      # ---- Factor levels from the full matched set ---- #
      factor_levels_full <- lapply(
        full_fin[, intersect(covars_fac_fin_run, colnames(full_fin)), drop = FALSE],
        # Existing level order keeps the reference categories (Never / No / White)
        function(x) if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x)))
      )
      
      # Put every factor covariate back onto those levels
      force_factor_levels <- function(df, lvls_list) {
        for (col in names(lvls_list)) {
          if (col %in% colnames(df)) {
            df[[col]] <- factor(as.character(df[[col]]), levels = lvls_list[[col]])
          }
        }
        df
      }
      
      # Fixed 70/30 split
      case_idx_fin    <- which(full_fin$binary_group == 1)
      control_idx_fin <- which(full_fin$binary_group == 0)
      train_idx_fin   <- c(
        sample(case_idx_fin,    round(0.7 * length(case_idx_fin))),
        sample(control_idx_fin, round(0.7 * length(control_idx_fin)))
      )
      train_fin <- full_fin[train_idx_fin, , drop = FALSE]
      
      test_idx_fin <- setdiff(seq_len(nrow(full_fin)), train_idx_fin)
      test_fin  <- full_fin[test_idx_fin,  , drop = FALSE]
      
      # ---- Impute proteins on train ---- #
      if (length(proteins_fin) > 0) {
        prot_cols_fin <- intersect(proteins_fin, colnames(train_fin))
        prot_mat_fin  <- apply(as.matrix(train_fin[, prot_cols_fin, drop = FALSE]), 2, as.numeric)
        if (is.null(dim(prot_mat_fin)))
          prot_mat_fin <- matrix(prot_mat_fin, ncol = 1, dimnames = list(NULL, prot_cols_fin))
        imputed_fin <- tryCatch(
          impute::impute.knn(prot_mat_fin)$data,
          error = function(e) prot_mat_fin
        )
        train_fin[, prot_cols_fin] <- imputed_fin
      }
      
      # Mean-substitute test proteins using training means
      train_prot_means_fin <- colMeans(imputed_fin, na.rm = TRUE)
      for (col in prot_cols_fin) {
        test_fin[[col]][is.na(test_fin[[col]])] <- train_prot_means_fin[[col]]
      }
      
      # ---- Impute numeric covariates ---- #
      for (col in intersect(covars_num, colnames(train_fin))) {
        m <- mean(train_fin[[col]], na.rm = TRUE)
        train_fin[[col]][is.na(train_fin[[col]])] <- ifelse(is.na(m), 0, m)
        test_fin[[col]][is.na(test_fin[[col]])]   <- ifelse(is.na(m), 0, m)
      }
      
      # ---- Impute factor covariates (mode); this can turn factors into character ---- #
      for (col in intersect(covars_fac, colnames(train_fin))) {
        vals <- train_fin[[col]][!is.na(train_fin[[col]])]
        if (length(vals) == 0) next
        mode_val <- names(sort(table(vals), decreasing = TRUE))[1]
        train_fin[[col]][is.na(train_fin[[col]])] <- mode_val
        test_fin[[col]][is.na(test_fin[[col]])]   <- mode_val
      }
      
      # ---- Restore the factor levels before building the matrices ---- #
      train_fin <- force_factor_levels(train_fin, factor_levels_full)
      test_fin  <- force_factor_levels(test_fin,  factor_levels_full)
      
      x_train_fin <- build_feature_matrix(train_fin, proteins_fin, use_cov_fin, covars_num, covars_fac_fin_run)
      y_train_fin <- as.numeric(train_fin$binary_group) - 1
      
      if (!is.null(x_train_fin)) {
        model_fin <- fit_risk_model(x_train_fin, y_train_fin, alpha = alpha, nfolds = 10)
        
        if (!is.null(model_fin)) {
          
          # ---- Evaluate the saved model on its true 30% holdout ---- #
          x_test_fin <- build_feature_matrix(test_fin, proteins_fin, use_cov_fin, covars_num, covars_fac_fin_run)
          x_test_fin <- x_test_fin[, colnames(x_train_fin), drop = FALSE]  # enforce same column order
          mode(x_test_fin) <- "numeric"
          y_test_fin <- as.numeric(test_fin$binary_group) - 1
          
          pred_test_fin <- predict_risk_model(model_fin, x_test_fin)
          
          # ---- Calibration and accuracy on the same 30% holdout ---- #
          eps          <- 1e-6
          pred_clipped <- pmin(pmax(pred_test_fin, eps), 1 - eps)
          logit_pred   <- stats::qlogis(pred_clipped)
          
          # Calibration slope & intercept: regress observed outcome on logit(pred).
          # Slope=1, Intercept=0 -> perfectly calibrated. Slope<1 -> predictions
          # too extreme (overconfident); Slope>1 -> too conservative.
          cal_fit <- tryCatch(
            glm(y_test_fin ~ logit_pred, family = binomial),
            error = function(e) NULL
          )
          calibration_intercept <- if (!is.null(cal_fit)) unname(coef(cal_fit)[1]) else NA_real_
          calibration_slope     <- if (!is.null(cal_fit)) unname(coef(cal_fit)[2]) else NA_real_
          
          # Calibration-in-the-large: slope fixed at 1, tests only whether mean
          # predicted risk matches mean observed risk (systematic over/under-prediction).
          cal_large_fit <- tryCatch(
            glm(y_test_fin ~ 1, family = binomial, offset = logit_pred),
            error = function(e) NULL
          )
          calibration_in_the_large <- if (!is.null(cal_large_fit)) unname(coef(cal_large_fit)[1]) else NA_real_
          
          # Brier score: mean squared error between predicted prob and outcome (0 = perfect)
          brier_score <- mean((pred_test_fin - y_test_fin)^2)
          
          # Log-loss: penalizes confident-but-wrong predictions more heavily than Brier/AUC
          log_loss <- -mean(y_test_fin * log(pred_clipped) + (1 - y_test_fin) * log(1 - pred_clipped))
          
          # Calibration table: observed vs mean predicted risk by decile (fewer bins if n < 50)
          n_bins <- if (nrow(test_fin) >= 50) 10 else max(2, floor(nrow(test_fin) / 10))
          calibration_table <- tibble::tibble(
            pred_prob = pred_test_fin,
            y_true    = y_test_fin
          ) %>%
            dplyr::mutate(bin = dplyr::ntile(pred_prob, n_bins)) %>%
            dplyr::group_by(bin) %>%
            dplyr::summarise(
              n              = dplyr::n(),
              mean_pred_prob = mean(pred_prob),
              observed_rate  = mean(y_true),
              .groups        = "drop"
            )
          
          roc_fin <- pROC::roc(
            response  = y_test_fin,
            predictor = pred_test_fin,
            quiet     = TRUE,
            direction = "<"
          )
          
          best_coords_fin <- pROC::coords(
            roc_fin, x = "best", best.method = "youden",
            ret = c("threshold", "sensitivity", "specificity", "accuracy", "ppv", "npv")
          )
          
          # Ties at the Youden-optimal point: take the first row (as above)
          if (is.matrix(best_coords_fin) || is.data.frame(best_coords_fin)) {
            if (nrow(best_coords_fin) > 1) best_coords_fin <- best_coords_fin[1, , drop = FALSE]
            best_coords_fin <- unlist(best_coords_fin)
          }
          
          performance_stats <- list(
            threshold                = as.numeric(best_coords_fin["threshold"]),
            sensitivity               = as.numeric(best_coords_fin["sensitivity"]),
            specificity               = as.numeric(best_coords_fin["specificity"]),
            accuracy                  = as.numeric(best_coords_fin["accuracy"]),
            ppv                       = as.numeric(best_coords_fin["ppv"]),
            npv                       = as.numeric(best_coords_fin["npv"]),
            auc                       = as.numeric(pROC::auc(roc_fin)),
            n_test                    = nrow(test_fin),
            n_test_cases              = sum(y_test_fin == 1),
            brier_score               = brier_score,
            log_loss                  = log_loss,
            calibration_intercept     = calibration_intercept,
            calibration_slope         = calibration_slope,
            calibration_in_the_large  = calibration_in_the_large
          )
          
          message(sprintf(
            "  Calibration -- Slope: %.3f | Intercept: %.3f | Brier: %.4f | Log-loss: %.4f",
            performance_stats$calibration_slope, performance_stats$calibration_intercept,
            performance_stats$brier_score, performance_stats$log_loss
          ))
          
          message(sprintf(
            "  Holdout performance -- AUC: %.3f | Sens: %.3f | Spec: %.3f | Acc: %.3f | Threshold: %.3f",
            performance_stats$auc, performance_stats$sensitivity,
            performance_stats$specificity, performance_stats$accuracy,
            performance_stats$threshold
          ))
          
          final_model_obj <- list(
            cv_fit            = model_fin$fit,
            lambda            = model_fin$lambda,
            alpha_used        = model_fin$alpha,
            proteins_used     = proteins_fin,
            use_covariates    = use_cov_fin,
            covars_num        = covars_num,
            covars_fac        = covars_fac_fin_run,
            factor_levels     = factor_levels_full,
            feature_colnames  = colnames(x_train_fin),
            train_means       = colMeans(x_train_fin, na.rm = TRUE),
            train_sds         = apply(x_train_fin, 2, stats::sd),
            performance_stats = performance_stats,
            calibration_table = calibration_table,
            model_name        = run$model_name,
            case_group        = run$case_group,
            window_years      = 14L
          )
          saveRDS(final_model_obj, file.path(dir.models, sprintf(
            "final_model_%s_%s_vs_%s_14yr.rds",
            run$model_name, run$case_group, run$control_group)))
          message(sprintf("  Saved final model: %s", run$model_name))
        }
      }
    }
  }
}

stopImplicitCluster()

end_time <- Sys.time()
run_time <- end_time - start_time
message(sprintf("\nTotal run time: %.1f %s", as.numeric(run_time), units(run_time)))


# Merge the per-run bootstrap results saved in bootstrap_raw/
bootstrap_files <- list.files(dir.bootstrap, pattern = sprintf("_bt%d\\.rds$", n_bootstrap), full.names = TRUE)
bootstrap_files <- bootstrap_files[!grepl("^AUC_results_all", basename(bootstrap_files))]
all_results <- dplyr::bind_rows(lapply(bootstrap_files, readRDS))
saveRDS(all_results, file.path(dir.bootstrap, sprintf("AUC_results_all_bt%d.rds", n_bootstrap)))
write.csv(all_results, file.path(dir.bootstrap, sprintf("AUC_results_all_bt%d.csv", n_bootstrap)),
          row.names = FALSE)

# --- 6b. MERGE THRESHOLD STATS ACROSS ALL WINDOWS & MODELS --- #
threshold_files <- list.files(dir.thresholds, pattern = "^threshold_stats_.*\\.rds$", full.names = TRUE)
threshold_stats_all <- dplyr::bind_rows(lapply(threshold_files, readRDS))
saveRDS(threshold_stats_all, file.path(dir.thresholds, "threshold_stats_all.rds"))
write.csv(threshold_stats_all, file.path(dir.thresholds, "threshold_stats_all.csv"), row.names = FALSE)

#### --- 7. SUMMARY --- ####
summary_df <- all_results %>%
  dplyr::group_by(Model, Case, Control_Group, Observation_Years) %>%
  dplyr::summarise(
    N_boot_valid   = sum(!is.na(AUC)),
    AUC_median     = median(AUC, na.rm = TRUE),
    AUC_lower      = quantile(AUC, 0.025, na.rm = TRUE),
    AUC_upper      = quantile(AUC, 0.975, na.rm = TRUE),
    N_case_window  = dplyr::first(N_case_window),
    N_case_matched    = dplyr::first(N_case_matched),
    N_control_matched = dplyr::first(N_control_matched),
    Age_SMD           = dplyr::first(Age_SMD),
    .groups        = "drop"
  ) %>%
  dplyr::mutate(
    contrast = paste0(Model, " - ", Case, " vs ", Control_Group),
    Model    = factor(Model,
                      levels = cfg$summary_model_order)
  )

saveRDS(summary_df, file.path(dir.summary, sprintf("AUC_summary_bt%d.rds", n_bootstrap)))
write.csv(summary_df, file.path(dir.summary, sprintf("AUC_summary_bt%d.csv", n_bootstrap)),
          row.names = FALSE)

#### --- 8. FEATURE IMPORTANCE --- ####
# Standardized ridge coefficients (coef x training SD) from every bootstrap fit.
# Ridge never zeroes a coefficient, so the useful questions are: how large, and
# is the sign consistent across bootstraps.
clean_feature_label <- function(x) {
  x <- gsub("_f[0-9]+_[0-9]+_[0-9]+", ": ", x)
  gsub("_", " ", sub(": $", "", x))
}

coef_files <- list.files(dir.importance, pattern = "^coefs_boot_.*\\.rds$", full.names = TRUE)

if (length(coef_files) > 0) {
  coefs_all <- dplyr::bind_rows(lapply(coef_files, readRDS)) %>%
    dplyr::group_by(Model, Case, Control_Group, Observation_Years, Bootstrap) %>%
    dplyr::mutate(boot_rank = rank(-abs(std_coef), ties.method = "first")) %>%
    dplyr::ungroup()
  
  importance_summary <- coefs_all %>%
    dplyr::group_by(Model, Case, Control_Group, Observation_Years, feature) %>%
    dplyr::summarise(
      n_boot         = dplyr::n(),
      coef_mean      = mean(coef),
      std_coef_mean  = mean(std_coef),
      std_coef_lower = quantile(std_coef, 0.025),
      std_coef_upper = quantile(std_coef, 0.975),
      sign_stability = mean(sign(std_coef) == sign(mean(std_coef))),
      median_rank    = median(boot_rank),
      pct_top5       = mean(boot_rank <= 5),
      .groups        = "drop"
    ) %>%
    dplyr::group_by(Model, Case, Control_Group, Observation_Years) %>%
    dplyr::mutate(rel_importance = abs(std_coef_mean) / sum(abs(std_coef_mean))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      feature_type  = ifelse(feature %in% colnames(UKBB_Proteomics), "Protein", "Covariate"),
      feature_label = clean_feature_label(feature)
    ) %>%
    dplyr::arrange(Model, Case, Control_Group, Observation_Years, dplyr::desc(abs(std_coef_mean)))
  
  for (mdl in unique(importance_summary$Model)) {
    imp_m <- dplyr::filter(importance_summary, Model == mdl)
    saveRDS(imp_m, file.path(dir.importance, sprintf("feature_importance_summary_%s.rds", mdl)))
    write.csv(imp_m, file.path(dir.importance, sprintf("feature_importance_summary_%s.csv", mdl)),
              row.names = FALSE)
  }
}

message(sprintf("\nDone. Outputs written to: %s", out.dir))

