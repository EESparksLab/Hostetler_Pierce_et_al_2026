#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidymodels)
  library(readr)
  library(dplyr)
  library(janitor)
  library(rlang)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("\nUsage:\n")
  cat("  Rscript reviewer_evaluate.R <model_rds> <required_cols_txt> <labeled_csv>\n\n")
  cat("Example:\n")
  cat("  Rscript reviewer_evaluate.R model2_final_fit_workflow.rds model2_required_columns.txt reviewer_labeled.csv\n\n")
  cat("The labeled CSV must include the ground-truth column: smurf_cat\n\n")
  quit(status = 1)
}

model_path <- args[1]
cols_path  <- args[2]
input_csv  <- args[3]

# Load model + schema
if (!file.exists(model_path)) stop("Model file not found: ", model_path)
if (!file.exists(cols_path))  stop("Required columns file not found: ", cols_path)
if (!file.exists(input_csv))  stop("Input CSV not found: ", input_csv)

final_fit <- readRDS(model_path)
required_cols <- readLines(cols_path)

if (length(required_cols) == 0) stop("Required columns file is empty: ", cols_path)

# Read and clean data
df <- read_csv(input_csv, show_col_types = FALSE) %>%
  clean_names()

if (!("smurf_cat" %in% names(df))) {
  stop(
    "Evaluation requires a ground-truth column named 'smurf_cat' (after clean_names()).\n",
    "Add a column smurf_cat to your CSV to compute confusion matrix + metrics."
  )
}

# Ensure label is categorical
df <- df %>% mutate(smurf_cat = as.factor(smurf_cat))

missing <- setdiff(required_cols, names(df))
if (length(missing) > 0) {
  stop(
    paste0(
      "Missing required predictor columns (after clean_names()):\n  ",
      paste(missing, collapse = ", "),
      "\n\nTip: ensure your CSV header matches the schema file."
    )
  )
}

X <- df %>% select(all_of(required_cols))
y <- df %>% select(smurf_cat)

# Predict class
pred_class <- predict(final_fit, X, type = "class")
preds <- bind_cols(y, pred_class)

cat("\nConfusion matrix:\n")
print(conf_mat(preds, truth = smurf_cat, estimate = .pred_class))

cat("\nMetrics:\n")
print(metrics(preds, truth = smurf_cat, estimate = .pred_class))

# If probability predictions are available, compute ROC AUC for binary outcomes only
pred_prob <- tryCatch(
  predict(final_fit, X, type = "prob"),
  error = function(e) NULL
)

if (!is.null(pred_prob)) {
  classes <- levels(df$smurf_cat)
  if (length(classes) == 2) {
    # parsnip names prob columns like .pred_<level>
    pos_level <- classes[2]
    prob_col <- paste0(".pred_", pos_level)
    
    if (prob_col %in% names(pred_prob)) {
      roc_df <- bind_cols(y, pred_prob)
      cat("\nROC AUC (binary only):\n")
      print(roc_auc(roc_df, truth = smurf_cat, !!sym(prob_col)))
    }
  }
}

cat("\nEvaluation complete.\n")
cat("Model: ", model_path, "\n")
cat("Data:  ", input_csv, "\n\n")
