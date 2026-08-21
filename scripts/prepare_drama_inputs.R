#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: Rscript scripts/prepare_drama_inputs.R",
    "<disease_id:1-12> <sex_stratified_dir> <combined_dir>",
    "<sex_ratio.tsv> <output_dir>"
  ))
}
if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")

disease_id <- as.integer(args[[1L]])
sex_dir <- args[[2L]]
combined_dir <- args[[3L]]
ratio_file <- args[[4L]]
output_dir <- args[[5L]]

diseases <- c("SLE", "CD", "PBC", "RA", "UC", "VIT", "AS", "CELIAC", "MS", "PSOA", "SJOGREN", "T1D")
file_stems <- c("sle", "cd", "pbc", "ra", "uc", "vit", "as", "celiac", "ms", "psoar", "sjogren", "t1d")
combined_counts <- c(6L, 3L, 2L, 2L, 3L, 2L, 1L, 1L, 3L, 3L, 1L, 1L)
female_neff <- round(c(6617.647648, 10855.9597, 2383.053774, 36193.66076, 14068.37708, 1685.528902,
  6103.398103, 9331.347388, 10604.79793, 3982.147533, 7747.355305, 10930.5235))
male_neff <- round(c(1445.743014, 8261.278507, 819.2759709, 16469.78557, 12985.49611, 1254.301943,
  6190.221469, 4442.554167, 4105.701087, 3593.98316, 1314.045168, 13393.54083))
if (is.na(disease_id) || disease_id < 1L || disease_id > length(diseases)) stop("disease_id must be 1 through 12.")

disease <- diseases[[disease_id]]
female_parts <- vector("list", 22L)
male_parts <- vector("list", 22L)
for (chr in 1:22) {
  female <- data.table::fread(file.path(sex_dir, sprintf("MAF_0.001_meta_female_chr%d_%s.txt", chr, disease)))
  male <- data.table::fread(file.path(sex_dir, sprintf("MAF_0.001_meta_male_chr%d_%s.txt", chr, disease)))
  key <- intersect(female$POS, male$POS)
  female_parts[[chr]] <- female[match(key, female$POS)]
  male_parts[[chr]] <- male[match(key, male$POS)]
}
female <- data.table::rbindlist(female_parts)
male <- data.table::rbindlist(male_parts)
if (!all(c("CHR", "POS", "ALLELE0", "ALLELE1", "BETA_meta", "SE_meta") %in% names(male)) ||
    !all(c("BETA_meta", "SE_meta") %in% names(female))) stop("Sex-stratified files lack required columns.")

female_std_beta <- (female$BETA_meta / female$SE_meta) / sqrt(female_neff[[disease_id]])
male_std_beta <- (male$BETA_meta / male$SE_meta) / sqrt(male_neff[[disease_id]])
snp <- paste(male$CHR, male$POS, male$ALLELE0, male$ALLELE1, "b38", sep = "_")
bp <- paste0(male$CHR, ":", male$POS)
beta <- data.frame(snp, female = female_std_beta, male = male_std_beta, check.names = FALSE)
se <- data.frame(snp, female = 1 / sqrt(female_neff[[disease_id]]), male = 1 / sqrt(male_neff[[disease_id]]), check.names = FALSE)

ratio_table <- data.table::fread(ratio_file)
if (!all(c("file", "ss", "female_sex_ratio") %in% names(ratio_table))) {
  stop("sex_ratio.tsv must contain file, ss, and female_sex_ratio columns.")
}
combined_ratio <- numeric(combined_counts[[disease_id]])
for (i in seq_len(combined_counts[[disease_id]])) {
  path <- file.path(combined_dir, disease, sprintf("MAF_0.001_hg_38_%s_eur%d_processed.txt", file_stems[[disease_id]], i))
  combined <- data.table::fread(path)
  if (!all(c("CHR", "POS", "beta", "se") %in% names(combined))) stop("Combined file lacks required columns: ", path)
  row <- match(path, ratio_table$file)
  if (is.na(row)) row <- match(basename(path), basename(ratio_table$file))
  if (is.na(row)) stop("No sample size/sex ratio entry for: ", path)
  z <- combined$beta / combined$se
  combined_bp <- paste0(combined$CHR, ":", combined$POS)
  beta[[paste0("combined_", i)]] <- (z / sqrt(ratio_table$ss[[row]]))[match(bp, combined_bp)]
  se[[paste0("combined_", i)]] <- (1 / sqrt(ratio_table$ss[[row]]))[match(bp, combined_bp)]
  combined_ratio[[i]] <- ratio_table$female_sex_ratio[[row]]
}

proportions <- data.frame(study = names(beta)[-1L], female_proportion = c(1, 0, combined_ratio))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(beta, file.path(output_dir, paste0(file_stems[[disease_id]], "_beta.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(se, file.path(output_dir, paste0(file_stems[[disease_id]], "_se.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(proportions, file.path(output_dir, paste0(file_stems[[disease_id]], "_female_proportions.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
