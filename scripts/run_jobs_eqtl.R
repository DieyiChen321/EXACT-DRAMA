#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript scripts/run_jobs_eqtl.R",
    "<beta.tsv> <se.tsv> <weights_dir> <output.tsv>"
  ))
}

beta_file <- args[[1L]]
se_file <- args[[2L]]
weights_dir <- args[[3L]]
output_file <- args[[4L]]

if (!requireNamespace("JOBS", quietly = TRUE)) {
  stop("Package 'JOBS' is required. Run Rscript scripts/install_dependencies.R first.")
}

weight_files <- file.path(weights_dir, sprintf("weight_est_13par_chr%d.txt", 1:22))
missing_files <- weight_files[!file.exists(weight_files)]
if (length(missing_files)) {
  stop("Missing weight files: ", paste(basename(missing_files), collapse = ", "))
}

weights_by_chr <- lapply(weight_files, function(path) {
  utils::read.delim(path, row.names = 1L, check.names = FALSE)
})
expected_columns <- names(weights_by_chr[[1L]])
if (!all(vapply(weights_by_chr, function(x) identical(names(x), expected_columns), logical(1L)))) {
  stop("Weight files do not have identical columns.")
}
weights <- colMeans(do.call(rbind, weights_by_chr), na.rm = TRUE)
if (any(!is.finite(weights))) stop("Mean weights contain non-finite values.")

beta <- utils::read.delim(beta_file, row.names = 1L, check.names = FALSE)
se <- utils::read.delim(se_file, row.names = 1L, check.names = FALSE)
if (!identical(dim(beta), dim(se)) || !identical(rownames(beta), rownames(se))) {
  stop("beta and se inputs must have aligned rows and identical dimensions.")
}
if (!"beta.est" %in% names(beta)) stop("beta input must contain a beta.est column.")

# Align the effect direction used in the manuscript analysis.
beta$beta.est <- -beta$beta.est
result <- JOBS::jobs.eqtls(beta, se, weight = as.numeric(weights), COR = TRUE)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.table(result, output_file, sep = "\t", quote = FALSE, row.names = TRUE)
