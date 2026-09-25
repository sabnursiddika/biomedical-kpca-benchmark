# ============================================================
# LEAKAGE-FREE NESTED-CV BENCHMARK:
# Six Biomedical Classification Datasets
# Wisconsin + Pima + Colon + Gravier+ Prostate + Leukemia
#
# Representations:
#   1. Raw standardized features
#   2. PCA
#   3. Global RBF Kernel PCA
#   4. Locally scaled Kernel PCA
#
# Classifiers:
#   1. Linear SVM
#   2. RBF SVM
#
# Primary tuning metric: ROC-AUC
# Includes expanded tuning, PR-AUC, inner-CV threshold selection,
# high-dimensional biomedical datasets, and runtime tracking
#
# IMPORTANT:
# - Imputation and scaling are estimated from training data only.
# - PCA/KPCA is fitted on training data only.
# - Test observations are projected with training-derived quantities.
# - Hyperparameters are selected in inner CV.
# - Final performance is estimated in outer CV.
# ============================================================
# ------------------------------------------------------------
# REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "e1071",
  "pROC",
  "PRROC",
  "plsgenomics",
  "sda",
  "SIS",
  "datamicroarray"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nPlease install the required packages before running the benchmark."
    )
  )
}

# ------------------------------------------------------------
# PROJECT PATHS
# ------------------------------------------------------------

# Run this script from the root of the repository:
# biomedical-kpca-benchmark/

PROJECT_DIR <- getwd()

DATA_DIR <- file.path(
  PROJECT_DIR,
  "data"
)

RESULTS_DIR <- file.path(
  PROJECT_DIR,
  "results"
)

WISCONSIN_FILE <- file.path(
  DATA_DIR,
  "wdbc.data"
)

PIMA_FILE <- file.path(
  DATA_DIR,
  "pima.csv"
)

# Stop immediately if either external input file is missing.
missing_input_files <- c(WISCONSIN_FILE, PIMA_FILE)[
  !file.exists(c(WISCONSIN_FILE, PIMA_FILE))
]

if (length(missing_input_files) > 0) {
  stop(
    "The following data file(s) were not found in PROJECT_DATA_DIR:
",
    paste(missing_input_files, collapse = "
"),
    "

Expected filenames are exactly: wdbc.data and pima.csv"
  )
}

message("Wisconsin file: ", normalizePath(WISCONSIN_FILE))
message("Pima file: ", normalizePath(PIMA_FILE))

# ------------------------------------------------------------
# STEP 1: BENCHMARK SETTINGS
# ------------------------------------------------------------
# Use "debug" for a quick test run.
# Use "final" only when intentionally running the full benchmark.
RUN_MODE <- "final"

# Safety switch: sourcing the script will NOT start the long benchmark unless TRUE.
# Keep FALSE while editing/checking. Change to TRUE only when intentionally running.
RUN_BENCHMARK <- FALSE

if (RUN_MODE == "debug") {
  OUTER_FOLDS   <- 5
  OUTER_REPEATS <- 1
  INNER_FOLDS   <- 3
  
  NCOMP_GRID           <- c(2, 5)
  KPCA_GAMMA_MULT_GRID <- c(0.25, 1)
  LOCAL_K_GRID          <- c(5, 15)
  SVM_COST_GRID         <- c(0.5, 2)
  SVM_GAMMA_MULT_GRID   <- c(0.5, 2)
  MAX_TUNING_CONFIGS    <- 20
  
} else if (RUN_MODE == "final") {
  OUTER_FOLDS   <- 5
  OUTER_REPEATS <- 3
  INNER_FOLDS   <- 3
  
  # The script automatically caps components at p and n - 1.
  NCOMP_GRID <- c(
    2, 3, 5, 8, 10, 15, 20, 25, 30, 40, 50
  )
  
  # Expanded toward smaller gamma because 0.125
  # was frequently selected at the lower boundary.
  KPCA_GAMMA_MULT_GRID <- c(
    0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4
  )
  
  # Expanded upward because k = 50 was frequently selected.
  # The script should cap infeasible k values automatically.
  LOCAL_K_GRID <- c(
    3, 5, 7, 10, 15, 20, 30, 50, 75, 100
  )
  
  # Expanded in both directions because C = 0.125
  # was frequently selected at the lower boundary.
  SVM_COST_GRID <- 2^c(
    -7, -5, -3, -1, 1, 3, 5, 7
  )
  
  # Expanded toward smaller gamma because 0.25
  # was frequently selected at the lower boundary.
  SVM_GAMMA_MULT_GRID <- c(
    0.0625, 0.125, 0.25, 0.5, 1, 2, 4
  )
  
  # Fixed-budget reproducible search prevents the full Cartesian grid
  # from becoming unnecessarily large. Boundary configurations are retained.
  MAX_TUNING_CONFIGS <- 60

  
} else {
  stop('RUN_MODE must be "debug", or "final".')
}

# Keep results from each run mode in a separate folder.
OUTPUT_DIR <- file.path(
  RESULTS_DIR,
  paste0("benchmark_results_", RUN_MODE)
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Results folder: ", normalizePath(OUTPUT_DIR, mustWork = FALSE))

MASTER_SEED <- 20260727
POSITIVE_CLASS <- "1"
TARGET_SPECIFICITY <- 0.80
TARGET_SENSITIVITY <- 0.80

# Optional sensitivity analysis. Keep FALSE for the primary comparison.
USE_BALANCED_CLASS_WEIGHTS <- FALSE

# ------------------------------------------------------------
# STEP 3: FIX THE STUDY QUESTION BEFORE ANALYSIS
# ------------------------------------------------------------
PRIMARY_RESEARCH_QUESTION <- paste(
  "Under what data conditions do PCA, global RBF-KPCA, and locally",
  "scaled KPCA improve classification compared with unreduced predictors?"
)

SECONDARY_RESEARCH_QUESTIONS <- c(
  "Does local scaling outperform global scaling when local density is heterogeneous?",
  "Does dimensionality reduction help more when the p-to-n ratio is high?",
  "Does nonlinear dimensionality reduction improve a linear classifier?",
  "Are performance differences stable across resamples?",
  "What predictive and computational trade-offs arise across representations?"
)

# ------------------------------------------------------------
# 1. DATA LOADERS
# ------------------------------------------------------------

load_wisconsin <- function(path) {
  if (!file.exists(path)) {
    stop("Wisconsin file not found: ", path)
  }
  
  dat <- read.csv(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (ncol(dat) != 32) {
    stop(
      "Wisconsin file should contain 32 columns: ID, diagnosis, and 30 predictors. ",
      "Found ", ncol(dat), " columns."
    )
  }
  
  diagnosis <- dat[[2]]
  
  if (!all(diagnosis %in% c("M", "B"))) {
    stop("Wisconsin diagnosis column must contain only M and B.")
  }
  
  x <- dat[, 3:32, drop = FALSE]
  colnames(x) <- paste0("X", seq_len(ncol(x)))
  
  # Malignant = 1, benign = 0
  y <- factor(ifelse(diagnosis == "M", "1", "0"), levels = c("0", "1"))
  
  list(
    name = "Wisconsin",
    x = as.data.frame(x),
    y = y
  )
}

load_pima <- function(path) {
  if (!file.exists(path)) {
    stop("Pima file not found: ", path)
  }
  
  dat <- read.csv(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (ncol(dat) != 9) {
    stop(
      "Pima file should contain 9 columns. Found ",
      ncol(dat), " columns."
    )
  }
  
  colnames(dat) <- c(
    "Pregnancies",
    "Glucose",
    "BloodPressure",
    "SkinThickness",
    "Insulin",
    "BMI",
    "DiabetesPedigree",
    "Age",
    "Outcome"
  )
  
  numeric_columns <- names(dat)
  dat[numeric_columns] <- lapply(dat[numeric_columns], function(z) {
    as.numeric(as.character(z))
  })
  
  if (anyNA(dat$Outcome) || !all(dat$Outcome %in% c(0, 1))) {
    stop("Pima outcome column must contain only 0 and 1.")
  }
  
  # Zeros in these physiological variables are treated as missing.
  zero_as_missing <- c(
    "Glucose",
    "BloodPressure",
    "SkinThickness",
    "Insulin",
    "BMI"
  )
  
  for (variable in zero_as_missing) {
    dat[[variable]][dat[[variable]] == 0] <- NA_real_
  }
  
  x <- dat[, setdiff(names(dat), "Outcome"), drop = FALSE]
  y <- factor(as.character(dat$Outcome), levels = c("0", "1"))
  
  list(
    name = "Pima",
    x = as.data.frame(x),
    y = y
  )
}


# ------------------------------------------------------------
# COLON CANCER
# ------------------------------------------------------------

load_colon <- function() {
  data("Colon", package = "plsgenomics", envir = environment())
  
  x <- as.data.frame(Colon$X)
  
  # Original package coding: 1 = normal, 2 = tumor.
  # For this benchmark: 0 = normal, 1 = tumor.
  y <- factor(
    ifelse(Colon$Y == 2, "1", "0"),
    levels = c("0", "1")
  )
  
  list(
    name = "Colon",
    x = x,
    y = y
  )
}


# ------------------------------------------------------------
# PROSTATE CANCER
# ------------------------------------------------------------

load_prostate <- function() {
  data("singh2002", package = "sda", envir = environment())
  
  x <- as.data.frame(singh2002$x)
  
  # For this benchmark: 0 = healthy, 1 = cancer.
  y <- factor(
    ifelse(singh2002$y == "cancer", "1", "0"),
    levels = c("0", "1")
  )
  
  list(
    name = "Prostate",
    x = x,
    y = y
  )
}


# ------------------------------------------------------------
# LEUKEMIA
# ------------------------------------------------------------

load_leukemia <- function() {
  data("leukemia.train", package = "SIS", envir = environment())
  data("leukemia.test", package = "SIS", envir = environment())
  
  # Combine the historical train/test partitions because this study
  # creates new nested-CV partitions from the complete 72-patient set.
  leukemia_all <- rbind(
    leukemia.train,
    leukemia.test
  )
  
  x <- as.data.frame(
    leukemia_all[, 1:7129, drop = FALSE]
  )
  
  # SIS coding: 0 = ALL, 1 = AML.
  y <- factor(
    leukemia_all[, 7130],
    levels = c(0, 1),
    labels = c("0", "1")
  )
  
  list(
    name = "Leukemia",
    x = x,
    y = y
  )
}

# ------------------------------------------------------------
# GRAVIER BREAST CANCER PROGNOSIS
# ------------------------------------------------------------

load_gravier <- function() {
	
  data(
    "gravier",
    package = "datamicroarray",
    envir = environment()
  )

  x <- as.data.frame(gravier$x)

  # Gravier coding for this benchmark:
  # 0 = good prognosis
  # 1 = poor prognosis
  y <- factor(
    ifelse(
      gravier$y == "poor",
      "1",
      "0"
    ),
    levels = c("0", "1")
  )

  list(
    name = "Gravier",
    x = x,
    y = y
  )
}

# ------------------------------------------------------------
# 1B. DATASET VALIDATION SUMMARY
# ------------------------------------------------------------

print_dataset_summary <- function(dataset) {
  zero_variance_predictors <- sum(
    vapply(
      dataset$x,
      function(z) {
        s <- stats::sd(z, na.rm = TRUE)
        is.na(s) || s < 1e-12
      },
      logical(1)
    )
  )
  
  cat("
============================================================
")
  cat("Dataset:", dataset$name, "
")
  cat("Observations:", nrow(dataset$x), "
")
  cat("Predictors:", ncol(dataset$x), "
")
  cat("p/n ratio:", round(ncol(dataset$x) / nrow(dataset$x), 2), "
")
  cat("Outcome counts:
")
  print(table(dataset$y, useNA = "ifany"))
  cat("Raw missing predictor cells:", sum(is.na(dataset$x)), "
")
  cat("Exact duplicate predictor rows:", sum(duplicated(dataset$x)), "
")
  cat("Zero-variance predictors:", zero_variance_predictors, "
")
  cat("============================================================
")
}

# ------------------------------------------------------------
# 2. STRATIFIED FOLDS
# ------------------------------------------------------------

make_stratified_folds <- function(y, k = 5, seed = 1) {
  y <- factor(y)
  set.seed(seed)
  
  fold_id <- integer(length(y))
  
  for (class_level in levels(y)) {
    class_indices <- which(y == class_level)
    class_indices <- sample(class_indices)
    
    assignments <- rep(seq_len(k), length.out = length(class_indices))
    assignments <- sample(assignments)
    
    fold_id[class_indices] <- assignments
  }
  
  lapply(seq_len(k), function(fold_number) {
    which(fold_id == fold_number)
  })
}

# ------------------------------------------------------------
# 3. TRAINING-ONLY PREPROCESSING
# ------------------------------------------------------------

fit_preprocessor <- function(x_train) {
  x_train <- as.data.frame(x_train)
  
  medians <- vapply(x_train, function(z) {
    z <- as.numeric(z)
    
    if (all(is.na(z))) {
      stop("A predictor is entirely missing in the training data.")
    }
    
    median(z, na.rm = TRUE)
  }, numeric(1))
  
  x_imputed <- x_train
  
  for (j in seq_along(x_imputed)) {
    missing_index <- is.na(x_imputed[[j]])
    x_imputed[[j]][missing_index] <- medians[j]
  }
  
  means <- vapply(x_imputed, mean, numeric(1))
  sds   <- vapply(x_imputed, stats::sd, numeric(1))
  
  # Prevent division by zero for constant predictors.
  sds[is.na(sds) | sds < 1e-12] <- 1
  
  list(
    medians = medians,
    means = means,
    sds = sds,
    columns = names(x_train)
  )
}

apply_preprocessor <- function(x, preprocessor) {
  x <- as.data.frame(x)
  
  if (!identical(names(x), preprocessor$columns)) {
    stop("Predictor columns do not match the fitted preprocessor.")
  }
  
  for (j in seq_along(x)) {
    missing_index <- is.na(x[[j]])
    x[[j]][missing_index] <- preprocessor$medians[j]
  }
  
  x_matrix <- as.matrix(x)
  storage.mode(x_matrix) <- "double"
  
  x_matrix <- sweep(x_matrix, 2, preprocessor$means, FUN = "-")
  x_matrix <- sweep(x_matrix, 2, preprocessor$sds, FUN = "/")
  
  # Remove predictors that are still non-finite for any reason.
  if (any(!is.finite(x_matrix))) {
    stop("Non-finite values remain after preprocessing.")
  }
  
  x_matrix
}

# ------------------------------------------------------------
# 4. DISTANCES AND KERNELS
# ------------------------------------------------------------

squared_distance_matrix <- function(a, b = NULL) {
  a <- as.matrix(a)
  storage.mode(a) <- "double"
  
  if (is.null(b)) {
    b <- a
  } else {
    b <- as.matrix(b)
    storage.mode(b) <- "double"
  }
  
  a_sq <- rowSums(a^2)
  b_sq <- rowSums(b^2)
  
  d2 <- outer(a_sq, b_sq, FUN = "+") - 2 * tcrossprod(a, b)
  d2[d2 < 0 & d2 > -1e-10] <- 0
  
  if (any(d2 < 0)) {
    stop("Substantial negative squared distances were produced.")
  }
  
  d2
}

median_rbf_gamma <- function(x) {
  d2 <- squared_distance_matrix(x)
  positive_d2 <- d2[upper.tri(d2) & d2 > 0]
  
  if (length(positive_d2) == 0) {
    stop("Unable to estimate an RBF bandwidth from identical observations.")
  }
  
  1 / median(positive_d2)
}

kth_neighbor_scale_train <- function(d2_train, k) {
  n <- nrow(d2_train)
  
  if (k < 1 || k >= n) {
    stop("For training data, local k must satisfy 1 <= k < n.")
  }
  
  d_train <- sqrt(d2_train)
  diag(d_train) <- Inf
  
  apply(d_train, 1, function(row_distances) {
    sort(row_distances, partial = k)[k]
  })
}

kth_neighbor_scale_test <- function(d2_test_train, k) {
  n_train <- ncol(d2_test_train)
  
  if (k < 1 || k > n_train) {
    stop("For test data, local k must satisfy 1 <= k <= n_train.")
  }
  
  d_test_train <- sqrt(d2_test_train)
  
  apply(d_test_train, 1, function(row_distances) {
    sort(row_distances, partial = k)[k]
  })
}

# ------------------------------------------------------------
# 5. MANUAL KPCA WITH OUT-OF-SAMPLE PROJECTION
# ------------------------------------------------------------

center_training_kernel <- function(k_train) {
  row_means <- rowMeans(k_train)
  col_means <- colMeans(k_train)
  overall_mean <- mean(k_train)
  
  k_centered <- sweep(k_train, 1, row_means, FUN = "-")
  k_centered <- sweep(k_centered, 2, col_means, FUN = "-")
  k_centered <- k_centered + overall_mean
  
  list(
    centered = k_centered,
    train_col_means = col_means,
    train_overall_mean = overall_mean
  )
}

center_test_kernel <- function(
    k_test,
    train_col_means,
    train_overall_mean
) {
  test_row_means <- rowMeans(k_test)
  
  k_centered <- sweep(k_test, 1, test_row_means, FUN = "-")
  k_centered <- sweep(k_centered, 2, train_col_means, FUN = "-")
  k_centered <- k_centered + train_overall_mean
  
  k_centered
}

fit_kpca_from_kernel <- function(k_train, ncomp, eigen_tolerance = 1e-8) {
  centered_info <- center_training_kernel(k_train)
  k_centered <- centered_info$centered
  
  eigen_result <- eigen(k_centered, symmetric = TRUE)
  
  positive_indices <- which(eigen_result$values > eigen_tolerance)
  
  if (length(positive_indices) == 0) {
    stop("The centered kernel has no usable positive eigenvalues.")
  }
  
  usable_components <- min(ncomp, length(positive_indices))
  selected_indices <- positive_indices[seq_len(usable_components)]
  
  eigenvalues <- eigen_result$values[selected_indices]
  eigenvectors <- eigen_result$vectors[, selected_indices, drop = FALSE]
  
  alphas <- sweep(
    eigenvectors,
    2,
    sqrt(eigenvalues),
    FUN = "/"
  )
  
  train_scores <- k_centered %*% alphas
  colnames(train_scores) <- paste0("KPC", seq_len(ncol(train_scores)))
  
  negative_values <- eigen_result$values[eigen_result$values < -eigen_tolerance]
  negative_mass <- if (length(negative_values) == 0) {
    0
  } else {
    sum(abs(negative_values)) / sum(abs(eigen_result$values))
  }
  
  list(
    train_scores = train_scores,
    alphas = alphas,
    train_col_means = centered_info$train_col_means,
    train_overall_mean = centered_info$train_overall_mean,
    eigenvalues = eigenvalues,
    negative_eigenvalue_count = length(negative_values),
    negative_eigenvalue_mass = negative_mass
  )
}

project_kpca_test <- function(k_test, kpca_model) {
  k_test_centered <- center_test_kernel(
    k_test = k_test,
    train_col_means = kpca_model$train_col_means,
    train_overall_mean = kpca_model$train_overall_mean
  )
  
  test_scores <- k_test_centered %*% kpca_model$alphas
  colnames(test_scores) <- colnames(kpca_model$train_scores)
  
  test_scores
}

# ------------------------------------------------------------
# 6. REPRESENTATION FITTING
# ------------------------------------------------------------

fit_representation <- function(
    x_train,
    x_test,
    representation,
    ncomp = NULL,
    kpca_gamma_multiplier = NULL,
    local_k = NULL
) {
  x_train <- as.matrix(x_train)
  x_test  <- as.matrix(x_test)
  
  if (representation == "raw") {
    return(list(
      train = x_train,
      test = x_test,
      diagnostics = list()
    ))
  }
  
  if (representation == "pca") {
    maximum_components <- min(ncol(x_train), nrow(x_train) - 1)
    ncomp_use <- min(ncomp, maximum_components)
    
    pca_model <- prcomp(
      x_train,
      center = FALSE,
      scale. = FALSE,
      rank. = ncomp_use
    )
    
    train_scores <- predict(pca_model, newdata = x_train)[, seq_len(ncomp_use), drop = FALSE]
    test_scores  <- predict(pca_model, newdata = x_test)[, seq_len(ncomp_use), drop = FALSE]
    
    return(list(
      train = train_scores,
      test = test_scores,
      diagnostics = list()
    ))
  }
  
  d2_train <- squared_distance_matrix(x_train)
  d2_test_train <- squared_distance_matrix(x_test, x_train)
  
  if (representation == "rbf_kpca") {
    base_gamma <- median_rbf_gamma(x_train)
    gamma <- base_gamma * kpca_gamma_multiplier
    
    k_train <- exp(-gamma * d2_train)
    k_test  <- exp(-gamma * d2_test_train)
    
    kpca_model <- fit_kpca_from_kernel(
      k_train = k_train,
      ncomp = ncomp
    )
    
    test_scores <- project_kpca_test(k_test, kpca_model)
    
    return(list(
      train = kpca_model$train_scores,
      test = test_scores,
      diagnostics = list(
        kernel_gamma = gamma,
        negative_eigenvalue_count = kpca_model$negative_eigenvalue_count,
        negative_eigenvalue_mass = kpca_model$negative_eigenvalue_mass
      )
    ))
  }
  
  if (representation == "local_kpca") {
    maximum_k <- nrow(x_train) - 1
    k_use <- min(local_k, maximum_k)
    
    sigma_train <- kth_neighbor_scale_train(d2_train, k = k_use)
    sigma_test  <- kth_neighbor_scale_test(d2_test_train, k = k_use)
    
    denominator_train <- outer(sigma_train, sigma_train)
    denominator_test  <- outer(sigma_test, sigma_train)
    
    denominator_train[denominator_train < 1e-12] <- 1e-12
    denominator_test[denominator_test < 1e-12] <- 1e-12
    
    k_train <- exp(-d2_train / denominator_train)
    k_test  <- exp(-d2_test_train / denominator_test)
    
    kpca_model <- fit_kpca_from_kernel(
      k_train = k_train,
      ncomp = ncomp
    )
    
    test_scores <- project_kpca_test(k_test, kpca_model)
    
    return(list(
      train = kpca_model$train_scores,
      test = test_scores,
      diagnostics = list(
        local_k = k_use,
        negative_eigenvalue_count = kpca_model$negative_eigenvalue_count,
        negative_eigenvalue_mass = kpca_model$negative_eigenvalue_mass
      )
    ))
  }
  
  stop("Unknown representation: ", representation)
}

# ------------------------------------------------------------
# 7. SVM FITTING AND PROBABILITY PREDICTION
# Representation scores are standardized inside each SVM fit using only
# that fit's training data (e1071 stores the scaling for prediction).
# ------------------------------------------------------------

fit_predict_svm <- function(
    x_train,
    y_train,
    x_test,
    classifier,
    cost,
    svm_gamma_multiplier = NULL
) {
  x_train <- as.matrix(x_train)
  x_test  <- as.matrix(x_test)
  y_train <- factor(y_train, levels = c("0", "1"))
  
  class_weights <- NULL
  if (USE_BALANCED_CLASS_WEIGHTS) {
    class_counts <- table(y_train)
    class_weights <- length(y_train) / (length(class_counts) * class_counts)
    class_weights <- as.numeric(class_weights)
    names(class_weights) <- names(class_counts)
  }
  
  if (classifier == "linear") {
    svm_model <- e1071::svm(
      x = x_train,
      y = y_train,
      type = "C-classification",
      kernel = "linear",
      cost = cost,
      scale = TRUE,
      probability = TRUE,
      class.weights = class_weights
    )
  } else if (classifier == "radial") {
    gamma <- svm_gamma_multiplier / ncol(x_train)
    
    svm_model <- e1071::svm(
      x = x_train,
      y = y_train,
      type = "C-classification",
      kernel = "radial",
      cost = cost,
      gamma = gamma,
      scale = TRUE,
      probability = TRUE,
      class.weights = class_weights
    )
  } else {
    stop("Unknown classifier: ", classifier)
  }
  
  predicted_class <- predict(
    svm_model,
    newdata = x_test,
    probability = TRUE
  )
  
  probability_matrix <- attr(predicted_class, "probabilities")
  
  if (is.null(probability_matrix)) {
    stop("SVM probability estimates were not returned.")
  }
  
  if (!(POSITIVE_CLASS %in% colnames(probability_matrix))) {
    stop(
      "Positive-class probability column was not found. Available columns: ",
      paste(colnames(probability_matrix), collapse = ", ")
    )
  }
  
  positive_probability <- probability_matrix[, POSITIVE_CLASS]
  
  list(
    class = factor(predicted_class, levels = c("0", "1")),
    probability = as.numeric(positive_probability),
    model = svm_model
  )
}

# ------------------------------------------------------------
# 8. PERFORMANCE METRICS
# ------------------------------------------------------------

safe_divide <- function(numerator, denominator) {
  if (
    length(numerator) != 1 ||
    length(denominator) != 1 ||
    !is.finite(numerator) ||
    !is.finite(denominator) ||
    denominator == 0
  ) {
    return(NA_real_)
  }
  
  numerator / denominator
}

confusion_metrics_at_threshold <- function(
    y_true,
    probability,
    threshold
) {
  y_true <- factor(y_true, levels = c("0", "1"))
  y_num <- as.integer(as.character(y_true))
  predicted_num <- ifelse(probability >= threshold, 1L, 0L)
  
  tp <- sum(predicted_num == 1L & y_num == 1L)
  tn <- sum(predicted_num == 0L & y_num == 0L)
  fp <- sum(predicted_num == 1L & y_num == 0L)
  fn <- sum(predicted_num == 0L & y_num == 1L)
  
  sensitivity <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  
  # Zero-division convention for classification summaries:
  # if a classifier predicts no positive cases, precision and F1 are set to 0
  # instead of NA. This avoids dropping degenerate folds from method averages.
  precision <- if ((tp + fp) == 0) {
    0
  } else {
    tp / (tp + fp)
  }
  
  f1_denominator <- 2 * tp + fp + fn
  f1 <- if (f1_denominator == 0) {
    0
  } else {
    2 * tp / f1_denominator
  }
  
  # Convert confusion counts to double precision before multiplication.
  # This prevents integer overflow when MCC is evaluated repeatedly on
  # moderate/large validation sets during threshold selection.
  tp_num <- as.numeric(tp)
  tn_num <- as.numeric(tn)
  fp_num <- as.numeric(fp)
  fn_num <- as.numeric(fn)
  
  mcc_denominator <- sqrt(
    (tp_num + fp_num) *
      (tp_num + fn_num) *
      (tn_num + fp_num) *
      (tn_num + fn_num)
  )
  
  # Mathematically MCC is undefined if one predicted class is absent.
  # For benchmark aggregation we use the explicit convention MCC = 0 in
  # these degenerate folds and retain a Constant_Prediction flag so the
  # behavior is transparent rather than silently omitting the fold.
  constant_prediction <- length(unique(predicted_num)) < 2
  
  mcc <- if (!is.finite(mcc_denominator) || mcc_denominator == 0) {
    0
  } else {
    (tp_num * tn_num - fp_num * fn_num) / mcc_denominator
  }
  
  c(
    Accuracy = safe_divide(tp + tn, tp + tn + fp + fn),
    Balanced_Accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    MCC = mcc,
    Constant_Prediction = as.numeric(constant_prediction)
  )
}

roc_auc_score <- function(y_true, probability) {
  y_num <- as.integer(as.character(factor(y_true, levels = c("0", "1"))))
  
  roc_object <- pROC::roc(
    response = y_num,
    predictor = probability,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  as.numeric(pROC::auc(roc_object))
}

pr_auc_score <- function(y_true, probability) {
  y_num <- as.integer(as.character(factor(y_true, levels = c("0", "1"))))
  
  positive_scores <- probability[y_num == 1L]
  negative_scores <- probability[y_num == 0L]
  
  if (length(positive_scores) == 0 || length(negative_scores) == 0) {
    return(NA_real_)
  }
  
  as.numeric(
    PRROC::pr.curve(
      scores.class0 = positive_scores,
      scores.class1 = negative_scores,
      curve = FALSE
    )$auc.integral
  )
}

operating_point_metrics <- function(
    y_true,
    probability,
    target_specificity = TARGET_SPECIFICITY,
    target_sensitivity = TARGET_SENSITIVITY
) {
  y_true <- factor(y_true, levels = c("0", "1"))
  
  thresholds <- unique(c(Inf, sort(probability, decreasing = TRUE), -Inf))
  operating_table <- do.call(
    rbind,
    lapply(thresholds, function(threshold) {
      metrics <- confusion_metrics_at_threshold(
        y_true = y_true,
        probability = probability,
        threshold = threshold
      )
      
      data.frame(
        Threshold = threshold,
        Sensitivity = metrics["Sensitivity"],
        Specificity = metrics["Specificity"]
      )
    })
  )
  
  sensitivity_candidates <- operating_table$Sensitivity[
    operating_table$Specificity >= target_specificity
  ]
  
  specificity_candidates <- operating_table$Specificity[
    operating_table$Sensitivity >= target_sensitivity
  ]
  
  sensitivity_at_target_specificity <- if (
    length(sensitivity_candidates) > 0
  ) {
    max(sensitivity_candidates, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  specificity_at_target_sensitivity <- if (
    length(specificity_candidates) > 0
  ) {
    max(specificity_candidates, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  c(
    Sensitivity_at_80_Specificity = sensitivity_at_target_specificity,
    Specificity_at_80_Sensitivity = specificity_at_target_sensitivity
  )
}

select_threshold_balanced_accuracy <- function(y_true, probability) {
  candidate_thresholds <- sort(unique(c(0, probability, 1)))
  
  threshold_table <- do.call(
    rbind,
    lapply(candidate_thresholds, function(threshold) {
      metrics <- confusion_metrics_at_threshold(
        y_true = y_true,
        probability = probability,
        threshold = threshold
      )
      
      data.frame(
        Threshold = threshold,
        Balanced_Accuracy = metrics["Balanced_Accuracy"],
        Sensitivity = metrics["Sensitivity"],
        Specificity = metrics["Specificity"]
      )
    })
  )
  
  best_value <- max(threshold_table$Balanced_Accuracy, na.rm = TRUE)
  best_rows <- which(
    abs(threshold_table$Balanced_Accuracy - best_value) < 1e-12
  )
  
  # Tie-break toward the conventional 0.50 threshold.
  best_row <- best_rows[
    which.min(abs(threshold_table$Threshold[best_rows] - 0.50))
  ]
  
  threshold_table[best_row, , drop = FALSE]
}

classification_metrics <- function(
    y_true,
    probability,
    tuned_threshold = 0.5
) {
  y_true <- factor(y_true, levels = c("0", "1"))
  y_num <- as.integer(as.character(y_true))
  
  roc_object <- pROC::roc(
    response = y_num,
    predictor = probability,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  ranking_metrics <- c(
    ROC_AUC = as.numeric(pROC::auc(roc_object)),
    PR_AUC = pr_auc_score(y_true, probability)
  )
  
  default_metrics <- confusion_metrics_at_threshold(
    y_true = y_true,
    probability = probability,
    threshold = 0.5
  )
  names(default_metrics) <- paste0(names(default_metrics), "_Default")
  
  tuned_metrics <- confusion_metrics_at_threshold(
    y_true = y_true,
    probability = probability,
    threshold = tuned_threshold
  )
  names(tuned_metrics) <- paste0(names(tuned_metrics), "_Tuned")
  
  operating_metrics <- operating_point_metrics(
    y_true = y_true,
    probability = probability
  )
  
  c(
    ranking_metrics,
    Tuned_Threshold = tuned_threshold,
    default_metrics,
    tuned_metrics,
    operating_metrics
  )
}

# ------------------------------------------------------------
# 9. HYPERPARAMETER GRIDS
# ------------------------------------------------------------

build_parameter_grid <- function(
    representation,
    classifier,
    p,
    n_train
) {
  valid_ncomp <- unique(
    pmax(
      1,
      pmin(NCOMP_GRID, p, n_train - 1)
    )
  )
  
  valid_local_k <- unique(
    LOCAL_K_GRID[LOCAL_K_GRID < n_train]
  )
  
  if (length(valid_local_k) == 0) {
    valid_local_k <- max(1, min(3, n_train - 1))
  }
  
  base <- data.frame(
    representation = representation,
    classifier = classifier,
    stringsAsFactors = FALSE
  )
  
  if (representation == "raw") {
    rep_grid <- data.frame(
      ncomp = NA_integer_,
      kpca_gamma_multiplier = NA_real_,
      local_k = NA_integer_
    )
  } else if (representation == "pca") {
    rep_grid <- data.frame(
      ncomp = valid_ncomp,
      kpca_gamma_multiplier = NA_real_,
      local_k = NA_integer_
    )
  } else if (representation == "rbf_kpca") {
    rep_grid <- expand.grid(
      ncomp = valid_ncomp,
      kpca_gamma_multiplier = KPCA_GAMMA_MULT_GRID,
      local_k = NA_integer_,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else if (representation == "local_kpca") {
    rep_grid <- expand.grid(
      ncomp = valid_ncomp,
      kpca_gamma_multiplier = NA_real_,
      local_k = valid_local_k,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown representation in grid builder.")
  }
  
  if (classifier == "linear") {
    svm_grid <- data.frame(
      cost = SVM_COST_GRID,
      svm_gamma_multiplier = NA_real_
    )
  } else if (classifier == "radial") {
    svm_grid <- expand.grid(
      cost = SVM_COST_GRID,
      svm_gamma_multiplier = SVM_GAMMA_MULT_GRID,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unknown classifier in grid builder.")
  }
  
  full_grid <- merge(rep_grid, svm_grid, by = NULL)
  full_grid$representation <- representation
  full_grid$classifier <- classifier
  
  full_grid[
    ,
    c(
      "representation",
      "classifier",
      "ncomp",
      "kpca_gamma_multiplier",
      "local_k",
      "cost",
      "svm_gamma_multiplier"
    )
  ]
}

limit_tuning_grid <- function(
    grid,
    max_configs,
    seed
) {
  if (nrow(grid) <= max_configs) {
    return(grid)
  }
  
  numeric_parameter_names <- intersect(
    c(
      "ncomp",
      "kpca_gamma_multiplier",
      "local_k",
      "cost",
      "svm_gamma_multiplier"
    ),
    names(grid)
  )
  
  boundary_rows <- integer(0)
  
  for (parameter_name in numeric_parameter_names) {
    values <- grid[[parameter_name]]
    finite_rows <- which(is.finite(values))
    
    if (length(finite_rows) > 0) {
      minimum_value <- min(values[finite_rows])
      maximum_value <- max(values[finite_rows])
      
      boundary_rows <- union(
        boundary_rows,
        which(values == minimum_value | values == maximum_value)
      )
    }
  }
  
  # Retain a manageable subset of boundary rows if there are many.
  set.seed(seed)
  if (length(boundary_rows) > floor(max_configs / 2)) {
    boundary_rows <- sample(
      boundary_rows,
      size = floor(max_configs / 2),
      replace = FALSE
    )
  }
  
  remaining_rows <- setdiff(seq_len(nrow(grid)), boundary_rows)
  remaining_slots <- max_configs - length(boundary_rows)
  
  sampled_rows <- if (remaining_slots > 0) {
    sample(
      remaining_rows,
      size = min(remaining_slots, length(remaining_rows)),
      replace = FALSE
    )
  } else {
    integer(0)
  }
  
  selected_rows <- sort(unique(c(boundary_rows, sampled_rows)))
  grid[selected_rows, , drop = FALSE]
}

# ------------------------------------------------------------
# 10. EVALUATE ONE CONFIGURATION IN INNER CV
# ------------------------------------------------------------

evaluate_configuration_inner_cv <- function(
    config,
    x,
    y,
    inner_folds
) {
  fold_auc <- rep(NA_real_, length(inner_folds))
  fold_pr_auc <- rep(NA_real_, length(inner_folds))
  
  for (fold_number in seq_along(inner_folds)) {
    validation_index <- inner_folds[[fold_number]]
    training_index <- setdiff(seq_len(nrow(x)), validation_index)
    
    x_inner_train_raw <- x[training_index, , drop = FALSE]
    x_inner_valid_raw <- x[validation_index, , drop = FALSE]
    y_inner_train <- y[training_index]
    y_inner_valid <- y[validation_index]
    
    fold_result <- tryCatch({
      preprocessor <- fit_preprocessor(x_inner_train_raw)
      
      x_inner_train <- apply_preprocessor(
        x_inner_train_raw,
        preprocessor
      )
      
      x_inner_valid <- apply_preprocessor(
        x_inner_valid_raw,
        preprocessor
      )
      
      representation_result <- fit_representation(
        x_train = x_inner_train,
        x_test = x_inner_valid,
        representation = config$representation,
        ncomp = config$ncomp,
        kpca_gamma_multiplier = config$kpca_gamma_multiplier,
        local_k = config$local_k
      )
      
      svm_result <- fit_predict_svm(
        x_train = representation_result$train,
        y_train = y_inner_train,
        x_test = representation_result$test,
        classifier = config$classifier,
        cost = config$cost,
        svm_gamma_multiplier = config$svm_gamma_multiplier
      )
      
      c(
        ROC_AUC = roc_auc_score(
          y_true = y_inner_valid,
          probability = svm_result$probability
        ),
        PR_AUC = pr_auc_score(
          y_true = y_inner_valid,
          probability = svm_result$probability
        )
      )
    }, error = function(e) {
      c(ROC_AUC = NA_real_, PR_AUC = NA_real_)
    })
    
    fold_auc[fold_number] <- fold_result["ROC_AUC"]
    fold_pr_auc[fold_number] <- fold_result["PR_AUC"]
  }
  
  c(
    mean_auc = mean(fold_auc, na.rm = TRUE),
    sd_auc = stats::sd(fold_auc, na.rm = TRUE),
    mean_pr_auc = mean(fold_pr_auc, na.rm = TRUE),
    sd_pr_auc = stats::sd(fold_pr_auc, na.rm = TRUE),
    valid_folds = sum(!is.na(fold_auc))
  )
}

collect_inner_oof_probabilities <- function(
    config,
    x,
    y,
    inner_folds
) {
  oof_probability <- rep(NA_real_, nrow(x))
  
  for (fold_number in seq_along(inner_folds)) {
    validation_index <- inner_folds[[fold_number]]
    training_index <- setdiff(seq_len(nrow(x)), validation_index)
    
    x_train_raw <- x[training_index, , drop = FALSE]
    x_valid_raw <- x[validation_index, , drop = FALSE]
    y_train <- y[training_index]
    
    preprocessor <- fit_preprocessor(x_train_raw)
    x_train <- apply_preprocessor(x_train_raw, preprocessor)
    x_valid <- apply_preprocessor(x_valid_raw, preprocessor)
    
    representation_result <- fit_representation(
      x_train = x_train,
      x_test = x_valid,
      representation = config$representation,
      ncomp = config$ncomp,
      kpca_gamma_multiplier = config$kpca_gamma_multiplier,
      local_k = config$local_k
    )
    
    svm_result <- fit_predict_svm(
      x_train = representation_result$train,
      y_train = y_train,
      x_test = representation_result$test,
      classifier = config$classifier,
      cost = config$cost,
      svm_gamma_multiplier = config$svm_gamma_multiplier
    )
    
    oof_probability[validation_index] <- svm_result$probability
  }
  
  if (anyNA(oof_probability)) {
    stop("Inner out-of-fold probabilities contain missing values.")
  }
  
  data.frame(
    Truth = factor(y, levels = c("0", "1")),
    Probability = oof_probability
  )
}

# ------------------------------------------------------------
# 11. INNER-CV TUNING
# ------------------------------------------------------------

tune_method <- function(
    x_outer_train,
    y_outer_train,
    representation,
    classifier,
    inner_folds,
    grid_seed
) {
  # Build the grid using the smallest inner-training sample size.
  # This prevents ncomp/local_k from being silently capped differently
  # across inner folds, which is especially important for small p >> n data.
  minimum_inner_training_n <- min(
    vapply(
      inner_folds,
      function(validation_index) {
        nrow(x_outer_train) - length(validation_index)
      },
      integer(1)
    )
  )
  
  full_grid <- build_parameter_grid(
    representation = representation,
    classifier = classifier,
    p = ncol(x_outer_train),
    n_train = minimum_inner_training_n
  )
  
  grid <- limit_tuning_grid(
    grid = full_grid,
    max_configs = MAX_TUNING_CONFIGS,
    seed = grid_seed
  )
  
  tuning_rows <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    config <- grid[i, , drop = FALSE]
    
    inner_result <- evaluate_configuration_inner_cv(
      config = config,
      x = x_outer_train,
      y = y_outer_train,
      inner_folds = inner_folds
    )
    
    tuning_rows[[i]] <- cbind(
      config,
      mean_inner_auc = unname(inner_result["mean_auc"]),
      sd_inner_auc = unname(inner_result["sd_auc"]),
      mean_inner_pr_auc = unname(inner_result["mean_pr_auc"]),
      sd_inner_pr_auc = unname(inner_result["sd_pr_auc"]),
      valid_inner_folds = unname(inner_result["valid_folds"])
    )
  }
  
  tuning_table <- do.call(rbind, tuning_rows)
  
  valid_rows <- which(
    is.finite(tuning_table$mean_inner_auc) &
      tuning_table$valid_inner_folds == length(inner_folds)
  )
  
  if (length(valid_rows) == 0) {
    stop(
      "No valid tuning configuration for ",
      representation,
      " + ",
      classifier
    )
  }
  
  best_index <- valid_rows[
    which.max(tuning_table$mean_inner_auc[valid_rows])
  ]
  
  best <- tuning_table[best_index, , drop = FALSE]
  
  inner_oof <- collect_inner_oof_probabilities(
    config = best,
    x = x_outer_train,
    y = y_outer_train,
    inner_folds = inner_folds
  )
  
  threshold_result <- select_threshold_balanced_accuracy(
    y_true = inner_oof$Truth,
    probability = inner_oof$Probability
  )
  
  list(
    best = best,
    table = tuning_table,
    selected_threshold = threshold_result$Threshold,
    threshold_inner_balanced_accuracy = threshold_result$Balanced_Accuracy,
    threshold_inner_sensitivity = threshold_result$Sensitivity,
    threshold_inner_specificity = threshold_result$Specificity
  )
}

# ------------------------------------------------------------
# 12. OUTER-CV EVALUATION
# ------------------------------------------------------------

method_specs <- expand.grid(
  representation = c(
    "raw",
    "pca",
    "rbf_kpca",
    "local_kpca"
  ),
  classifier = c(
    "linear",
    "radial"
  ),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

run_nested_cv_dataset <- function(dataset) {
  x <- dataset$x
  y <- dataset$y
  
  repeat_checkpoint_file <- file.path(
    OUTPUT_DIR,
    paste0(
      "repeat_checkpoint_",
      RUN_MODE,
      "_",
      tolower(dataset$name),
      ".rds"
    )
  )
  
  if (file.exists(repeat_checkpoint_file)) {
    
    repeat_state <- readRDS(repeat_checkpoint_file)
    
    all_outer_results <- repeat_state$all_outer_results
    all_tuning_results <- repeat_state$all_tuning_results
    all_outer_predictions <- repeat_state$all_outer_predictions
    
    result_counter <- repeat_state$result_counter
    tuning_counter <- repeat_state$tuning_counter
    prediction_counter <- repeat_state$prediction_counter
    
    completed_repeat <- repeat_state$completed_repeat
    
    message(
      "[", dataset$name,
      "] repeat checkpoint found. Completed repeat ",
      completed_repeat, "/", OUTER_REPEATS
    )
    
  } else {
    
    all_outer_results <- list()
    all_tuning_results <- list()
    all_outer_predictions <- list()
    
    result_counter <- 1
    tuning_counter <- 1
    prediction_counter <- 1
    
    completed_repeat <- 0
  }
  
  repeat_sequence <- if (completed_repeat < OUTER_REPEATS) {
    seq.int(completed_repeat + 1, OUTER_REPEATS)
  } else {
    integer(0)
  }
  
  for (repeat_number in repeat_sequence) {
    outer_folds <- make_stratified_folds(
      y = y,
      k = OUTER_FOLDS,
      seed = MASTER_SEED + 1000 * repeat_number
    )
    
    for (outer_fold_number in seq_along(outer_folds)) {
      outer_test_index <- outer_folds[[outer_fold_number]]
      outer_train_index <- setdiff(
        seq_len(nrow(x)),
        outer_test_index
      )
      
      x_outer_train_raw <- x[outer_train_index, , drop = FALSE]
      x_outer_test_raw  <- x[outer_test_index, , drop = FALSE]
      y_outer_train <- y[outer_train_index]
      y_outer_test  <- y[outer_test_index]
      
      inner_folds <- make_stratified_folds(
        y = y_outer_train,
        k = INNER_FOLDS,
        seed = MASTER_SEED +
          100000 * repeat_number +
          1000 * outer_fold_number
      )
      
      for (method_number in seq_len(nrow(method_specs))) {
        representation <- method_specs$representation[method_number]
        classifier <- method_specs$classifier[method_number]
        
        message(
          "[", dataset$name, "] ",
          "repeat ", repeat_number, "/", OUTER_REPEATS,
          ", outer fold ", outer_fold_number, "/", OUTER_FOLDS,
          ", ", representation, " + ", classifier
        )
        
        method_start_time <- proc.time()[["elapsed"]]
        
        tuning_result <- tune_method(
          x_outer_train = x_outer_train_raw,
          y_outer_train = y_outer_train,
          representation = representation,
          classifier = classifier,
          inner_folds = inner_folds,
          grid_seed = MASTER_SEED +
            1000000 * repeat_number +
            10000 * outer_fold_number +
            100 * method_number
        )
        
        best <- tuning_result$best
        selected_threshold <- tuning_result$selected_threshold
        
        tuning_table <- tuning_result$table
        tuning_table$Dataset <- dataset$name
        tuning_table$Outer_Repeat <- repeat_number
        tuning_table$Outer_Fold <- outer_fold_number
        
        all_tuning_results[[tuning_counter]] <- tuning_table
        tuning_counter <- tuning_counter + 1
        
        svm_result <- NULL
        
        outer_result <- tryCatch({
          preprocessor <- fit_preprocessor(x_outer_train_raw)
          
          x_outer_train <- apply_preprocessor(
            x_outer_train_raw,
            preprocessor
          )
          
          x_outer_test <- apply_preprocessor(
            x_outer_test_raw,
            preprocessor
          )
          
          representation_result <- fit_representation(
            x_train = x_outer_train,
            x_test = x_outer_test,
            representation = representation,
            ncomp = best$ncomp,
            kpca_gamma_multiplier = best$kpca_gamma_multiplier,
            local_k = best$local_k
          )
          
          svm_result <- fit_predict_svm(
            x_train = representation_result$train,
            y_train = y_outer_train,
            x_test = representation_result$test,
            classifier = classifier,
            cost = best$cost,
            svm_gamma_multiplier = best$svm_gamma_multiplier
          )
          
          metrics <- classification_metrics(
            y_true = y_outer_test,
            probability = svm_result$probability,
            tuned_threshold = selected_threshold
          )
          
          diagnostic_count <- representation_result$diagnostics$negative_eigenvalue_count
          diagnostic_mass  <- representation_result$diagnostics$negative_eigenvalue_mass
          
          if (is.null(diagnostic_count)) diagnostic_count <- NA_real_
          if (is.null(diagnostic_mass)) diagnostic_mass <- NA_real_
          
          data.frame(
            Dataset = dataset$name,
            Outer_Repeat = repeat_number,
            Outer_Fold = outer_fold_number,
            Representation = representation,
            Classifier = classifier,
            Ncomp = best$ncomp,
            KPCA_Gamma_Multiplier = best$kpca_gamma_multiplier,
            Local_K = best$local_k,
            SVM_Cost = best$cost,
            SVM_Gamma_Multiplier = best$svm_gamma_multiplier,
            Inner_ROC_AUC = best$mean_inner_auc,
            Inner_PR_AUC = best$mean_inner_pr_auc,
            Selected_Threshold = selected_threshold,
            Inner_Threshold_Balanced_Accuracy = tuning_result$threshold_inner_balanced_accuracy,
            ROC_AUC = metrics["ROC_AUC"],
            PR_AUC = metrics["PR_AUC"],
            Accuracy_Default = metrics["Accuracy_Default"],
            Balanced_Accuracy_Default = metrics["Balanced_Accuracy_Default"],
            Sensitivity_Default = metrics["Sensitivity_Default"],
            Specificity_Default = metrics["Specificity_Default"],
            Precision_Default = metrics["Precision_Default"],
            F1_Default = metrics["F1_Default"],
            MCC_Default = metrics["MCC_Default"],
            Constant_Prediction_Default = metrics["Constant_Prediction_Default"],
            Accuracy_Tuned = metrics["Accuracy_Tuned"],
            Balanced_Accuracy_Tuned = metrics["Balanced_Accuracy_Tuned"],
            Sensitivity_Tuned = metrics["Sensitivity_Tuned"],
            Specificity_Tuned = metrics["Specificity_Tuned"],
            Precision_Tuned = metrics["Precision_Tuned"],
            F1_Tuned = metrics["F1_Tuned"],
            MCC_Tuned = metrics["MCC_Tuned"],
            Constant_Prediction_Tuned = metrics["Constant_Prediction_Tuned"],
            Sensitivity_at_80_Specificity = metrics["Sensitivity_at_80_Specificity"],
            Specificity_at_80_Sensitivity = metrics["Specificity_at_80_Sensitivity"],
            Negative_Eigenvalue_Count = diagnostic_count,
            Negative_Eigenvalue_Mass = diagnostic_mass,
            stringsAsFactors = FALSE
          )
        }, error = function(e) {
          warning(
            "Outer evaluation failed for ",
            dataset$name, ", ",
            representation, " + ", classifier,
            ": ", conditionMessage(e)
          )
          
          data.frame(
            Dataset = dataset$name,
            Outer_Repeat = repeat_number,
            Outer_Fold = outer_fold_number,
            Representation = representation,
            Classifier = classifier,
            Ncomp = best$ncomp,
            KPCA_Gamma_Multiplier = best$kpca_gamma_multiplier,
            Local_K = best$local_k,
            SVM_Cost = best$cost,
            SVM_Gamma_Multiplier = best$svm_gamma_multiplier,
            Inner_ROC_AUC = best$mean_inner_auc,
            Inner_PR_AUC = best$mean_inner_pr_auc,
            Selected_Threshold = selected_threshold,
            Inner_Threshold_Balanced_Accuracy = tuning_result$threshold_inner_balanced_accuracy,
            ROC_AUC = NA_real_,
            PR_AUC = NA_real_,
            Accuracy_Default = NA_real_,
            Balanced_Accuracy_Default = NA_real_,
            Sensitivity_Default = NA_real_,
            Specificity_Default = NA_real_,
            Precision_Default = NA_real_,
            F1_Default = NA_real_,
            MCC_Default = NA_real_,
            Constant_Prediction_Default = NA_real_,
            Accuracy_Tuned = NA_real_,
            Balanced_Accuracy_Tuned = NA_real_,
            Sensitivity_Tuned = NA_real_,
            Specificity_Tuned = NA_real_,
            Precision_Tuned = NA_real_,
            F1_Tuned = NA_real_,
            MCC_Tuned = NA_real_,
            Constant_Prediction_Tuned = NA_real_,
            Sensitivity_at_80_Specificity = NA_real_,
            Specificity_at_80_Sensitivity = NA_real_,
            Negative_Eigenvalue_Count = NA_real_,
            Negative_Eigenvalue_Mass = NA_real_,
            stringsAsFactors = FALSE
          )
        })
        
        elapsed_seconds <- proc.time()[["elapsed"]] - method_start_time
        outer_result$Elapsed_Seconds <- elapsed_seconds
        
        if (!is.null(svm_result) && all(is.finite(svm_result$probability))) {
          prediction_rows <- data.frame(
            Dataset = dataset$name,
            Outer_Repeat = repeat_number,
            Outer_Fold = outer_fold_number,
            Observation_Index = outer_test_index,
            Representation = representation,
            Classifier = classifier,
            Truth = as.character(y_outer_test),
            Probability = svm_result$probability,
            Selected_Threshold = selected_threshold,
            Predicted_Default = ifelse(svm_result$probability >= 0.5, "1", "0"),
            Predicted_Tuned = ifelse(
              svm_result$probability >= selected_threshold,
              "1",
              "0"
            ),
            stringsAsFactors = FALSE
          )
          
          all_outer_predictions[[prediction_counter]] <- prediction_rows
          prediction_counter <- prediction_counter + 1
        }
        
        all_outer_results[[result_counter]] <- outer_result
        result_counter <- result_counter + 1
      }
    }
    
    saveRDS(
      list(
        completed_repeat = repeat_number,
        all_outer_results = all_outer_results,
        all_tuning_results = all_tuning_results,
        all_outer_predictions = all_outer_predictions,
        result_counter = result_counter,
        tuning_counter = tuning_counter,
        prediction_counter = prediction_counter
      ),
      repeat_checkpoint_file
    )
    
    message(
      "[", dataset$name,
      "] repeat checkpoint saved after repeat ",
      repeat_number, "/", OUTER_REPEATS
    )
  }
  
  list(
    outer_results = do.call(rbind, all_outer_results),
    tuning_results = do.call(rbind, all_tuning_results),
    outer_predictions = if (length(all_outer_predictions) > 0) {
      do.call(rbind, all_outer_predictions)
    } else {
      data.frame()
    }
  )
}

# ------------------------------------------------------------
# 13. RESULT SUMMARIES
# ------------------------------------------------------------

summarize_results <- function(results) {
  metric_names <- c(
    "ROC_AUC",
    "PR_AUC",
    "Accuracy_Default",
    "Balanced_Accuracy_Default",
    "Sensitivity_Default",
    "Specificity_Default",
    "F1_Default",
    "MCC_Default",
    "Constant_Prediction_Default",
    "Accuracy_Tuned",
    "Balanced_Accuracy_Tuned",
    "Sensitivity_Tuned",
    "Specificity_Tuned",
    "F1_Tuned",
    "MCC_Tuned",
    "Constant_Prediction_Tuned",
    "Sensitivity_at_80_Specificity",
    "Specificity_at_80_Sensitivity"
  )
  
  groups <- unique(
    results[, c("Dataset", "Representation", "Classifier")]
  )
  
  summary_rows <- list()
  counter <- 1
  
  for (i in seq_len(nrow(groups))) {
    group <- groups[i, , drop = FALSE]
    
    subset_rows <- results[
      results$Dataset == group$Dataset &
        results$Representation == group$Representation &
        results$Classifier == group$Classifier,
      ,
      drop = FALSE
    ]
    
    elapsed_values <- if ("Elapsed_Seconds" %in% names(subset_rows)) {
      suppressWarnings(as.numeric(subset_rows$Elapsed_Seconds))
    } else {
      numeric(0)
    }
    elapsed_values <- elapsed_values[is.finite(elapsed_values)]
    
    row_output <- data.frame(
      Dataset = group$Dataset,
      Representation = group$Representation,
      Classifier = group$Classifier,
      Valid_Outer_Folds = sum(is.finite(subset_rows$ROC_AUC)),
      Mean_Elapsed_Seconds = if (length(elapsed_values) > 0) {
        mean(elapsed_values)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
    
    for (metric in metric_names) {
      values <- subset_rows[[metric]]
      valid_values <- values[is.finite(values)]
      
      row_output[[paste0(metric, "_Mean")]] <- if (
        length(valid_values) > 0
      ) {
        mean(valid_values)
      } else {
        NA_real_
      }
      
      row_output[[paste0(metric, "_SD")]] <- if (
        length(valid_values) > 1
      ) {
        stats::sd(valid_values)
      } else {
        NA_real_
      }
      
      row_output[[paste0(metric, "_SE")]] <- if (
        length(valid_values) > 1
      ) {
        stats::sd(valid_values) / sqrt(length(valid_values))
      } else {
        NA_real_
      }
    }
    
    summary_rows[[counter]] <- row_output
    counter <- counter + 1
  }
  
  summary_table <- do.call(rbind, summary_rows)
  summary_table[
    order(
      summary_table$Dataset,
      -summary_table$ROC_AUC_Mean
    ),
    ,
    drop = FALSE
  ]
}

# ------------------------------------------------------------
# 14. RUN ALL SIX DATASETS
# ------------------------------------------------------------

set.seed(MASTER_SEED)

research_question_lines <- c(
  "PRIMARY RESEARCH QUESTION",
  PRIMARY_RESEARCH_QUESTION,
  "",
  "SECONDARY RESEARCH QUESTIONS",
  paste0(
    seq_along(SECONDARY_RESEARCH_QUESTIONS),
    ". ",
    SECONDARY_RESEARCH_QUESTIONS
  )
)

writeLines(
  research_question_lines,
  con = file.path(OUTPUT_DIR, "research_questions.txt")
)

# Create each dataset exactly once.
wisconsin <- load_wisconsin(WISCONSIN_FILE)
pima      <- load_pima(PIMA_FILE)
colon     <- load_colon()
gravier   <- load_gravier()
prostate  <- load_prostate()
leukemia  <- load_leukemia()

# Keep the same object structure for all six datasets.
datasets <- list(
  wisconsin,
  pima,
  colon,
  gravier,
  prostate,
  leukemia
)
# Verify the final benchmark contains exactly six datasets.
expected_datasets <- c(
  "Wisconsin",
  "Pima",
  "Colon",
  "Gravier",
  "Prostate",
  "Leukemia"
)

actual_datasets <- vapply(
  datasets,
  function(d) d$name,
  character(1)
)

stopifnot(
  length(datasets) == 6,
  identical(actual_datasets, expected_datasets)
)
# Print validation summaries before any modeling begins.
invisible(lapply(datasets, print_dataset_summary))

# Automatic study-level dataset summary.
dataset_summary <- do.call(
  rbind,
  lapply(datasets, function(dat) {
    zero_variance_predictors <- sum(
      vapply(
        dat$x,
        function(z) {
          s <- stats::sd(z, na.rm = TRUE)
          is.na(s) || s < 1e-12
        },
        logical(1)
      )
    )
    
    data.frame(
      Dataset = dat$name,
      Observations = nrow(dat$x),
      Predictors = ncol(dat$x),
      P_to_N_Ratio = round(ncol(dat$x) / nrow(dat$x), 2),
      Negative_Class = sum(dat$y == "0"),
      Positive_Class = sum(dat$y == "1"),
      Raw_Missing_Cells = sum(is.na(dat$x)),
      Duplicate_Predictor_Rows = sum(duplicated(dat$x)),
      Zero_Variance_Predictors = zero_variance_predictors,
      stringsAsFactors = FALSE
    )
  })
)

print(dataset_summary)

write.csv(
  dataset_summary,
  file.path(OUTPUT_DIR, "dataset_summary.csv"),
  row.names = FALSE
)

message(
  "\nStarting ", RUN_MODE,
  " nested-CV benchmark for ", length(datasets), " datasets."
)
message(
  "Settings: outer folds = ", OUTER_FOLDS,
  ", outer repeats = ", OUTER_REPEATS,
  ", inner folds = ", INNER_FOLDS,
  ", max tuning configs = ", MAX_TUNING_CONFIGS
)

# The long benchmark is protected by RUN_BENCHMARK so an accidental
# click on Source does not start the entire experiment again.
if (RUN_BENCHMARK) {
  
  message("\nBenchmark execution is ENABLED.")
  
  # Run one dataset at a time and checkpoint immediately after each dataset.
  # If the run is interrupted, completed datasets remain recoverable.
  dataset_results <- vector("list", length(datasets))
  names(dataset_results) <- vapply(datasets, function(d) d$name, character(1))
  
 checkpoint_file <- file.path(
  OUTPUT_DIR,
  "dataset_results_checkpoint_6datasets_FINAL.rds"
)
  completed_datasets <- 0
  
  if (file.exists(checkpoint_file)) {
    
    saved_dataset_results <- readRDS(checkpoint_file)
    
    dataset_results[
      seq_along(saved_dataset_results)
    ] <- saved_dataset_results
    
    completed_datasets <- length(saved_dataset_results)
    
    message(
      "Loaded dataset checkpoint: ",
      completed_datasets, "/",
      length(datasets),
      " dataset(s) already complete."
    )
  }
  
  remaining_dataset_indices <- if (
    completed_datasets < length(datasets)
  ) {
    seq.int(
      completed_datasets + 1,
      length(datasets)
    )
  } else {
    integer(0)
  }
  
  for (dataset_index in remaining_dataset_indices) {
    current_dataset <- datasets[[dataset_index]]
    
    message(
      "\n===== Starting dataset ", dataset_index, "/", length(datasets),
      ": ", current_dataset$name, " ====="
    )
    
    dataset_results[[dataset_index]] <- run_nested_cv_dataset(current_dataset)
    
    saveRDS(
      dataset_results[seq_len(dataset_index)],
      checkpoint_file
    )
    
    repeat_checkpoint_file <- file.path(
      OUTPUT_DIR,
      paste0(
        "repeat_checkpoint_",
        RUN_MODE,
        "_",
        tolower(current_dataset$name),
        ".rds"
      )
    )
    
    if (file.exists(repeat_checkpoint_file)) {
      unlink(repeat_checkpoint_file)
    }
    
    message(
      "Checkpoint saved after ", current_dataset$name, ": ",
      normalizePath(checkpoint_file, mustWork = FALSE)
    )
  }
  
  outer_results <- do.call(
    rbind,
    lapply(dataset_results, `[[`, "outer_results")
  )
  
  tuning_results <- do.call(
    rbind,
    lapply(dataset_results, `[[`, "tuning_results")
  )
  
  outer_predictions <- do.call(
    rbind,
    lapply(dataset_results, `[[`, "outer_predictions")
  )
  
  summary_results <- summarize_results(outer_results)
  
  write.csv(
    outer_results,
    file.path(OUTPUT_DIR, "nested_cv_outer_results.csv"),
    row.names = FALSE
  )
  
  write.csv(
    tuning_results,
    file.path(OUTPUT_DIR, "nested_cv_tuning_results.csv"),
    row.names = FALSE
  )
  
  write.csv(
    outer_predictions,
    file.path(OUTPUT_DIR, "nested_cv_outer_predictions.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary_results,
    file.path(OUTPUT_DIR, "nested_cv_summary_6datasets_FINAL.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    list(
      settings = list(
        run_mode = RUN_MODE,
        outer_folds = OUTER_FOLDS,
        outer_repeats = OUTER_REPEATS,
        inner_folds = INNER_FOLDS,
        ncomp_grid = NCOMP_GRID,
        kpca_gamma_multiplier_grid = KPCA_GAMMA_MULT_GRID,
        local_k_grid = LOCAL_K_GRID,
        svm_cost_grid = SVM_COST_GRID,
        svm_gamma_multiplier_grid = SVM_GAMMA_MULT_GRID,
        max_tuning_configs = MAX_TUNING_CONFIGS,
        use_balanced_class_weights = USE_BALANCED_CLASS_WEIGHTS,
        target_specificity = TARGET_SPECIFICITY,
        target_sensitivity = TARGET_SENSITIVITY,
        seed = MASTER_SEED
      ),
      research_questions = list(
        primary = PRIMARY_RESEARCH_QUESTION,
        secondary = SECONDARY_RESEARCH_QUESTIONS
      ),
      dataset_summary = dataset_summary,
      outer_results = outer_results,
      tuning_results = tuning_results,
      outer_predictions = outer_predictions,
      summary_results = summary_results,
      session_info = utils::sessionInfo()
    ),
    file.path(
      OUTPUT_DIR,
      paste0("complete_benchmark_results_6dataserta_FINAL.rds")
    )
  )
  
  print(summary_results)
  

  
  # R warnings are printed here. Note: LIBSVM can also print solver messages
  # directly to the console; those messages may not appear in warnings().
  warnings()
  
} else {
  message(
    "\nSetup loaded successfully. RUN_BENCHMARK is FALSE, so no long benchmark was started."
  )
  message(
    "When ready, set RUN_BENCHMARK <- TRUE and Source the script intentionally."
  )
}