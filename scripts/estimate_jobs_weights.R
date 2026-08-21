#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript scripts/estimate_jobs_weights.R",
    "<chromosome:1-22> <beta_dir> <se_dir> <output_dir>"
  ))
}

chromosome <- as.integer(args[[1L]])
beta_dir <- args[[2L]]
se_dir <- args[[3L]]
output_dir <- args[[4L]]

if (is.na(chromosome) || chromosome < 1L || chromosome > 22L) {
  stop("chromosome must be an integer from 1 through 22.")
}
if (!requireNamespace("JOBS", quietly = TRUE)) {
  stop("Package 'JOBS' is required. Run Rscript scripts/install_dependencies.R first.")
}

beta_file <- file.path(beta_dir, sprintf("matched_eqtl_beta_chr%d.txt", chromosome))
se_file <- file.path(se_dir, sprintf("matched_eqtl_se_chr%d.txt", chromosome))
if (!file.exists(beta_file)) stop("Beta file not found: ", beta_file)
if (!file.exists(se_file)) stop("Standard-error file not found: ", se_file)

beta <- utils::read.delim(beta_file, row.names = 1L, check.names = FALSE)
se <- utils::read.delim(se_file, row.names = 1L, check.names = FALSE)
if (!identical(dim(beta), dim(se)) || !identical(rownames(beta), rownames(se))) {
  stop("beta and se inputs must have aligned rows and identical dimensions.")
}
if (!identical(names(beta), names(se))) stop("beta and se column names are not aligned.")
if (!"beta.est" %in% names(beta)) stop("beta input must contain a beta.est column.")
if (ncol(beta) < 16L) {
  stop("Expected at least 16 beta columns so columns 3:16 identify 14 cell types.")
}

# Align the effect direction used in the manuscript analysis.
beta$beta.est <- -beta$beta.est
weight <- JOBS::jobs.nnls.weights(beta, se)
weights <- as.data.frame(t(as.numeric(weight)), check.names = FALSE)

cell_type_names <- sub("_beta$", "", names(beta)[3:16])
if (ncol(weights) != length(cell_type_names)) {
  stop(
    "jobs.nnls.weights returned ", ncol(weights), " weights, but ",
    length(cell_type_names), " cell-type columns were identified."
  )
}
names(weights) <- cell_type_names

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, sprintf("weight_est_13par_chr%d.txt", chromosome))
utils::write.table(weights, output_file, sep = "\t", quote = FALSE, row.names = TRUE)
message("Wrote ", output_file)
