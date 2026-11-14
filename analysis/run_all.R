#!/usr/bin/env Rscript
# ============================================================================
# Breast Cancer Survival Analysis — entry point
#
# Run from anywhere:
#     Rscript analysis/run_all.R
#
# The script locates the project root from its own path, sources the R/ modules
# in order, and writes figures/tables to output/.
# ============================================================================

# Locate the project root from this script's path (robust to working directory).
args      <- commandArgs(trailingOnly = FALSE)
file_arg  <- args[grep("^--file=", args)][1]
this_file <- sub("^--file=", "", file_arg)
.root     <- normalizePath(file.path(dirname(this_file), ".."))

setwd(.root)
cat("Project root:", .root, "\n")

source("R/00_setup.R")
source("R/01_load_clean.R")
source("R/02_eda.R")
source("R/03_inference.R")
source("R/04_ml_models.R")
source("R/05_evaluation.R")
source("R/06_interpret.R")

# Reproducibility snapshot
writeLines(capture.output(sessionInfo()), out_path("output/session_info.txt"))

cat("\n===========================================================\n")
cat("DONE. Figures in output/figures/, tables in output/tables/\n")
cat("Final C-index ranking:\n")
print(cindex_summary[, c("model", "c_index_mean", "c_index_sd")])
cat("===========================================================\n")
