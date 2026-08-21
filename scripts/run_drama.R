#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: Rscript scripts/run_drama.R",
    "<beta.tsv> <se.tsv> <female_proportions.tsv> <chunk_id> <output.tsv>"
  ))
}

beta_file <- args[[1L]]
se_file <- args[[2L]]
proportion_file <- args[[3L]]
chunk_id <- as.integer(args[[4L]])
output_file <- args[[5L]]
if (is.na(chunk_id) || chunk_id < 1L) stop("chunk_id must be a positive integer.")

source(file.path("R", "drama.R"))
beta <- utils::read.delim(beta_file, check.names = FALSE)
se <- utils::read.delim(se_file, check.names = FALSE)
proportions <- utils::read.delim(proportion_file, check.names = FALSE)

required <- c("study", "female_proportion")
if (!all(required %in% names(proportions))) {
  stop("female_proportions.tsv must contain: study, female_proportion.")
}
study_names <- names(beta)[-1L]
if (!identical(study_names, names(se)[-1L])) stop("beta and se study columns differ.")
match_index <- match(study_names, proportions$study)
if (anyNA(match_index)) stop("A study is missing from female_proportions.tsv.")
female_proportion <- proportions$female_proportion[match_index]
if (any(!is.finite(female_proportion)) || any(female_proportion < 0 | female_proportion > 1)) {
  stop("Female proportions must be finite values from 0 to 1.")
}

chunk_count <- 10L
chunk_size <- ceiling(nrow(beta) / chunk_count)
start <- (chunk_id - 1L) * chunk_size + 1L
end <- min(chunk_id * chunk_size, nrow(beta))
if (start > nrow(beta)) stop("chunk_id is larger than the number of available chunks.")

design <- cbind(female_proportion, male_proportion = 1 - female_proportion)
result <- drama(beta[start:end, , drop = FALSE], se[start:end, , drop = FALSE], design)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.table(result, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
