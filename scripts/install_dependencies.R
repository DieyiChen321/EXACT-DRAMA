#!/usr/bin/env Rscript

cran <- c("data.table", "MASS", "remotes")
missing <- cran[!vapply(cran, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!requireNamespace("JOBS", quietly = TRUE)) {
  remotes::install_github("LidaWangPSU/JOBS", subdir = "JOBS", upgrade = "never")
}
