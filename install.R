# install.R ---------------------------------------------------------------
# One-shot dependency installer. Installs all required CRAN packages into the
# project-local `lib/` directory (loaded automatically by .Rprofile), keeping
# the project self-contained and reproducible without polluting the system.
#
# Usage:  Rscript install.R

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 900)

lib <- file.path(getwd(), "lib")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

pkgs <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "survival", "survminer",
  "glmnet", "gbm", "randomForestSRC", "xgboost", "mice", "car",
  "gtsummary", "timeROC", "pdp", "shapviz", "gridExtra",
  "viridis", "scales"
)

already <- rownames(installed.packages())
todo <- setdiff(pkgs, already)
cat("Installing", length(todo), "packages into", lib, "\n")

if (length(todo) > 0) install.packages(todo, lib = lib, dependencies = TRUE)

missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
  stop("Some packages failed to install: ", paste(missing, collapse = ", "))
}
cat("All dependencies installed.\n")
