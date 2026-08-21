source(file.path("R", "drama.R"))

beta <- data.frame(
  snp = c("1_100_A_G_b38", "1_200_C_T_b38"),
  female = c(0.20, -0.10), male = c(0.10, -0.04), combined_1 = c(0.14, -0.06)
)
se <- data.frame(
  snp = beta$snp, female = c(0.05, 0.04), male = c(0.05, 0.04),
  combined_1 = c(0.03, 0.03)
)
design <- rbind(c(1, 0), c(0, 1), c(0.6, 0.4))
result <- drama(beta, se, design)

stopifnot(
  nrow(result) == 2L,
  identical(names(result), c("snp", "beta.female", "beta.male", "se.female", "se.male")),
  all(is.finite(as.matrix(result[-1L])))
)
message("DRAMA smoke test passed.")
