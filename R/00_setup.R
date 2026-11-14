# 00_setup.R ---------------------------------------------------------------
# Project-wide configuration: package loading, paths, random seed, constants,
# and small helpers shared by every downstream module.
# This file is sourced by analysis/run_all.R AFTER `.root` is defined.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(glmnet)
  library(gbm)
  library(randomForestSRC)
  library(xgboost)
  library(mice)
  library(car)
  library(gtsummary)
  library(timeROC)
  library(pdp)
  library(shapviz)
  library(gridExtra)
  library(viridis)
  library(scales)
})

# ---- Reproducibility -----------------------------------------------------
SEED       <- 20240115          # master random seed
N_FOLDS    <- 5                 # cross-validation folds
N_REPEATS  <- 5                 # repeated CV repetitions (stability of C-index)
AUC_TIMES  <- c(12, 24, 36, 60) # months for time-dependent AUC
GBM_TREES  <- 300               # trees used inside the CV loop (fast, fixed)
RSF_TREES  <- 300               # trees used inside the CV loop
XGB_ROUNDS <- 150               # boosting rounds inside the CV loop

# ---- Styling -------------------------------------------------------------
theme_set(theme_minimal(base_size = 12))
OKBLUE   <- "#2166AC"
OKRED    <- "#B2182B"

# ---- Helpers -------------------------------------------------------------
# Resolve an output path relative to the project root.
out_path <- function(...) file.path(.root, ...)

save_fig <- function(plot, name, w = 7, h = 5) {
  dir.create(out_path("output/figures"), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_path("output/figures", name), plot = plot, width = w, height = h, dpi = 300)
  invisible(plot)
}

save_table <- function(df, name) {
  dir.create(out_path("output/tables"), showWarnings = FALSE, recursive = TRUE)
  write_csv(df, out_path("output/tables", name))
  invisible(df)
}

log_step <- function(msg) cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
