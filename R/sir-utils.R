# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# Matrix helpers shared by the proposal machinery.

# Sampling Importance Resampling (SIR) ----------------------------------------
#
# Reference: Dosne et al. (2013), PAGE 22, Abstract 2907.
# Algorithm mirrors PsN's sir tool (psn_ref/sir, sir.pm, sir_userguide.tex).

# Helpers ----------------------------------------------------------------------

# Apply bounds elementwise across rows of a matrix; returns logical vector.
.sirInBounds <- function(mat, lower, upper) {
  apply(mat, 1L, function(x) all(x >= lower & x <= upper))
}

.sirCapCovCorrelation <- function(covMat, capCorrelation = 0.8) {
  checkmate::assertMatrix(covMat, mode = "numeric")
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  cov_mat <- (covMat + t(covMat)) / 2
  if (capCorrelation >= 1 || ncol(cov_mat) < 2L) {
    return(cov_mat)
  }

  param_names <- colnames(cov_mat)
  sd_vals <- sqrt(pmax(diag(cov_mat), 0))
  corr_mat <- suppressWarnings(cov2cor(cov_mat))
  corr_mat[!is.finite(corr_mat)] <- 0
  diag(corr_mat) <- 1

  lo <- lower.tri(corr_mat)
  corr_mat[lo] <- pmax(-capCorrelation, pmin(capCorrelation, corr_mat[lo]))
  corr_mat[upper.tri(corr_mat)] <- t(corr_mat)[upper.tri(corr_mat)]

  capped <- outer(sd_vals, sd_vals) * corr_mat
  dimnames(capped) <- list(param_names, param_names)
  capped
}

.sirEnsurePosDef <- function(covMat, minEigen = sqrt(.Machine$double.eps)) {
  checkmate::assertMatrix(covMat, mode = "numeric")
  dim_names <- dimnames(covMat)
  cov_mat <- (covMat + t(covMat)) / 2
  eig <- eigen(cov_mat, symmetric = TRUE)
  if (all(is.finite(eig$values)) && min(eig$values) >= minEigen) {
    dimnames(cov_mat) <- dim_names
    return(cov_mat)
  }

  values <- pmax(eig$values, minEigen)
  out <- eig$vectors %*% diag(values, nrow = length(values)) %*% t(eig$vectors)
  out <- (out + t(out)) / 2
  dimnames(out) <- dim_names
  out
}
