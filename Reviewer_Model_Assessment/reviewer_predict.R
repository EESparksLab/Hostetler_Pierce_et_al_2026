#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidymodels)
  library(readr)
  library(dplyr)
  library(janitor)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("\nUsage:\n")
  cat("  Rscript reviewer_predict.R <model_rds> <required_cols_txt> <input_csv> <output_csv>\n\n")
  cat("Example:\n")
  cat("  Rscript reviewer_predict.R model2_final_fit_workflow.rds model2_required_columns.txt reviewer_data.csv predictions.csv\n\n")
  quit(status = 1)
}

model_path <- args[1]
cols_path  <- args[2]
input_csv  <- args[3]
output_csv <- args[4]

# Load model + schema
if (!file.exists(model_path)) stop("Model file not found: ", model_path)
if (!file.exists(cols_path))  stop("Required columns file not found: ", cols_path)
if (!file.exists(input_csv))  stop("Input CSV not found: ", input_csv)

final_fit <- readRDS(model_path)
required_cols <- readLines(cols_path)

if (length(required_cols) == 0) stop("Required columns file is empty: ", cols_path)

# Read and clean input
new_data <- read_csv(input_csv, show_col_types = FALSE) %>%
  clean_names()

missing <- setdiff(required_cols, names(new_data))
extra   <- setdiff(names(new_data), required_cols)

if (length(missing) > 0) {
  stop(
    paste0(
      "Missing required columns (after clean_names()):\n  ",
      paste(missing, collapse = ", "),
      "\n\nTip: ensure your CSV header matches the schema file."
    )
  )
}

# Keep only required columns, in the exact schema order
X <- new_data %>% select(all_of(required_cols))

# Predict class
pred_class <- predict(final_fit, X, type = "class")

# Predict probabilities (if supported)
pred_prob <- tryCatch(
  predict(final_fit, X, type = "prob"),
  error = function(e) NULL
)

out <- bind_cols(X, pred_class)
if (!is.null(pred_prob)) out <- bind_cols(out, pred_prob)

write_csv(out, output_csv)

cat("\nPrediction complete.\n")
cat("Model:  ", model_path, "\n")
cat("Input:  ", input_csv, "\n")
cat("Output: ", output_csv, "\n")

if (length(extra) > 0) {
  cat("\nNote: the following extra columns were ignored:\n  ")
  cat(paste(extra, collapse = ", "))
  cat("\n")
}
cat("\n")
