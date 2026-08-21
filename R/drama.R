drama <- function(beta, se, design) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Install it with install.packages('MASS').")
  }
  if (!is.data.frame(beta) || !is.data.frame(se)) {
    stop("beta and se must be data frames.")
  }
  if (!identical(dim(beta), dim(se))) {
    stop("beta and se must have identical dimensions.")
  }
  if (ncol(beta) < 3L || nrow(beta) == 0L) {
    stop("beta and se must contain a SNP column and at least two study columns.")
  }
  if (!identical(as.character(beta[[1L]]), as.character(se[[1L]]))) {
    stop("SNP identifiers in beta and se are not aligned.")
  }

  beta_values <- as.matrix(beta[-1L])
  se_values <- as.matrix(se[-1L])
  storage.mode(beta_values) <- "double"
  storage.mode(se_values) <- "double"
  design <- as.matrix(design)
  storage.mode(design) <- "double"

  if (nrow(design) != ncol(beta_values)) {
    stop("The design matrix must have one row per study column.")
  }
  if (ncol(design) != 2L) {
    stop("The DRAMA design matrix must have two columns: female and male proportions.")
  }
  if (any(!is.finite(design))) stop("The design matrix contains non-finite values.")

  keep <- rowSums(abs(beta_values), na.rm = TRUE) != 0
  keep <- keep & rowSums(!is.na(beta_values)) >= ncol(design)
  message(nrow(beta), " SNPs in total")
  message(sum(keep), " SNPs retained for estimation")

  beta_values <- beta_values[keep, , drop = FALSE]
  se_values <- se_values[keep, , drop = FALSE]

  estimates <- t(vapply(seq_len(nrow(beta_values)), function(i) {
    available <- which(
      is.finite(beta_values[i, ]) & is.finite(se_values[i, ]) & se_values[i, ] > 0
    )
    if (length(available) < ncol(design)) return(rep(NA_real_, 4L))

    x <- design[available, , drop = FALSE]
    precision <- diag(1 / se_values[i, available]^2, nrow = length(available))
    covariance <- MASS::ginv(crossprod(x, precision %*% x))
    coefficient <- covariance %*% crossprod(x, precision %*% beta_values[i, available])
    c(as.numeric(coefficient), sqrt(diag(covariance)))
  }, numeric(4L)))

  data.frame(
    snp = as.character(beta[[1L]][keep]),
    beta.female = estimates[, 1L],
    beta.male = estimates[, 2L],
    se.female = estimates[, 3L],
    se.male = estimates[, 4L],
    check.names = FALSE
  )
}
