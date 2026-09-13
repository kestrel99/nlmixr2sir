# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# Low-level THETA and OMEGA samplers.

# Step 2 -----------------------------------------------------------------------

#' Sample THETA vectors from a (truncated) multivariate normal
#'
#' Draws `n` parameter vectors from a multivariate normal distribution and
#' discards any that violate the supplied bounds.  Sampling is repeated in
#' batches until `n` valid vectors are collected or the maximum attempt
#' budget is exhausted.
#'
#' @param mu Named numeric vector of THETA means (length p).
#' @param covMat p×p covariance matrix.
#' @param n Positive integer. Number of valid samples to collect.
#' @param lower Numeric vector (length p, or scalar recycled). Lower bounds.
#'   Default `-Inf` (no lower truncation).
#' @param upper Numeric vector (length p, or scalar recycled). Upper bounds.
#'   Default `Inf` (no upper truncation).
#' @return A named list:
#'   \describe{
#'     \item{`samples`}{Numeric matrix, at most `n` rows × p columns, with
#'       column names from `mu`.  Fewer than `n` rows if the budget was
#'       exhausted before `n` valid draws were found.}
#'     \item{`nRejected`}{Integer. Total number of draws that violated at
#'       least one bound.}
#'   }
#' @noRd
sirSampleTheta <- function(
  mu,
  covMat,
  n,
  lower = rep(-Inf, length(mu)),
  upper = rep(Inf, length(mu))
) {
  p <- length(mu)
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(covMat, mode = "numeric", nrows = p, ncols = p)
  checkmate::assertCount(n, positive = TRUE)
  if (length(lower) == 1L) {
    lower <- rep(lower, p)
  }
  if (length(upper) == 1L) {
    upper <- rep(upper, p)
  }
  checkmate::assertNumeric(lower, len = p, any.missing = FALSE)
  checkmate::assertNumeric(upper, len = p, any.missing = FALSE)

  n <- as.integer(n)
  max_attempts <- 10L * n
  param_names <- names(mu)

  collected <- matrix(NA_real_, nrow = n, ncol = p)
  n_filled <- 0L
  n_rejected <- 0L
  n_attempted <- 0L

  while (n_filled < n && n_attempted < max_attempts) {
    batch_n <- min(n - n_filled, max_attempts - n_attempted)
    draws <- mvtnorm::rmvnorm(batch_n, mean = mu, sigma = covMat)
    n_attempted <- n_attempted + batch_n

    ok <- .sirInBounds(draws, lower, upper)
    good <- draws[ok, , drop = FALSE]
    n_rejected <- n_rejected + sum(!ok)

    n_take <- min(nrow(good), n - n_filled)
    if (n_take > 0L) {
      idx <- seq(n_filled + 1L, n_filled + n_take)
      collected[idx, ] <- good[seq_len(n_take), ]
      n_filled <- n_filled + n_take
    }
  }

  if (n_filled < n) {
    cli::cli_warn(c(
      "Only {n_filled} of {n} requested samples were within bounds \\
       after {max_attempts} draw attempts.",
      "i" = "Consider widening the parameter bounds or reducing {.arg n}."
    ))
    collected <- collected[seq_len(n_filled), , drop = FALSE]
  }

  colnames(collected) <- param_names
  list(samples = collected, nRejected = n_rejected)
}

# Step 3 -----------------------------------------------------------------------

#' Sample OMEGA matrices from a multivariate normal, retaining only PD draws
#'
#' Vectorizes the lower triangle (including diagonal) of `omegaEst`, draws from
#' a multivariate normal with covariance `omegaCovMat`, reconstructs symmetric
#' matrices, and discards any that are not positive definite.
#'
#' @param omegaEst Square symmetric matrix. Current OMEGA estimate.
#' @param omegaCovMat Square numeric matrix of size `p × p`, where
#'   `p = length(lower.tri(omegaEst, diag = TRUE))`. Covariance of the
#'   lower-triangle elements (e.g., constructed from SEs in `fit$parFixedDf`).
#' @param n Positive integer. Number of valid (positive-definite) samples to
#'   collect.
#' @return A named list:
#'   \describe{
#'     \item{`samples`}{List of at most `n` positive-definite matrices with the
#'       same `dimnames` as `omegaEst`. Fewer than `n` entries if the budget was
#'       exhausted.}
#'     \item{`nRejected`}{Integer. Draws that failed the Cholesky check.}
#'   }
#' @noRd
sirSampleOmegaSigma <- function(omegaEst, omegaCovMat, n) {
  checkmate::assertMatrix(omegaEst, mode = "numeric")
  if (nrow(omegaEst) != ncol(omegaEst)) {
    cli::cli_abort("{.arg omegaEst} must be a square matrix.")
  }
  dim_omega <- nrow(omegaEst)
  lt_idx <- which(lower.tri(omegaEst, diag = TRUE))
  mu <- omegaEst[lt_idx]
  p <- length(mu)
  checkmate::assertMatrix(omegaCovMat, mode = "numeric", nrows = p, ncols = p)
  checkmate::assertCount(n, positive = TRUE)

  n <- as.integer(n)
  max_attempts <- 10L * n

  collected <- vector("list", n)
  n_filled <- 0L
  n_rejected <- 0L
  n_attempted <- 0L

  while (n_filled < n && n_attempted < max_attempts) {
    batch_n <- min(n - n_filled, max_attempts - n_attempted)
    draws <- mvtnorm::rmvnorm(batch_n, mean = mu, sigma = omegaCovMat)
    n_attempted <- n_attempted + batch_n

    for (i in seq_len(nrow(draws))) {
      mat <- matrix(0, dim_omega, dim_omega)
      mat[lt_idx] <- draws[i, ]
      mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
      dimnames(mat) <- dimnames(omegaEst)

      is_pd <- tryCatch(
        {
          chol(mat)
          TRUE
        },
        error = function(e) FALSE
      )
      if (is_pd) {
        n_filled <- n_filled + 1L
        collected[[n_filled]] <- mat
        if (n_filled == n) break
      } else {
        n_rejected <- n_rejected + 1L
      }
    }
  }

  if (n_filled < n) {
    cli::cli_warn(c(
      "Only {n_filled} of {n} OMEGA samples were positive definite \\
       after {max_attempts} draw attempts.",
      "i" = "Consider revising {.arg omegaCovMat} or reducing {.arg n}."
    ))
    collected <- collected[seq_len(n_filled)]
  }

  list(samples = collected, nRejected = n_rejected)
}
