# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# Importance weights, resampling, and per-iteration raw results.

# Step 5 -----------------------------------------------------------------------

#' Compute relative PDFs, importance ratios, and resampling probabilities
#'
#' Given a matrix of sampled parameter vectors, their delta-OFV values relative
#' to the original fit, and the proposal distribution (mean + covariance),
#' computes the importance ratio for each sample and normalises them into
#' resampling probabilities.
#'
#' The relative PDF is the multivariate normal density evaluated at each sample
#' *relative* to the density at `mu`, so it equals 1 at the proposal mean and
#' decreases for samples further away.  The likelihood ratio converts the OFV
#' difference back to a probability scale.  Samples with `NA` `dOFV` (failed
#' evaluations) receive `prob_resample = 0` and are excluded from resampling.
#'
#' @param samples Numeric matrix, one row per sample (same parameter space as
#'   `mu` / `covMat`).
#' @param mu Named numeric vector. Proposal mean (length = `ncol(samples)`).
#' @param covMat Positive-definite covariance matrix of the proposal (same
#'   dimension as `length(mu)`).
#' @param dOFV Numeric vector of length `nrow(samples)`.
#'   `dOFV[i] = OFV_sample[i] - OFV_original`.  May contain `NA`.
#' @return Data frame with one row per sample and columns:
#'   \describe{
#'     \item{`sample_id`}{Integer row index.}
#'     \item{`dOFV`}{As supplied.}
#'     \item{`likelihood_ratio`}{`exp(-0.5 * dOFV)`.}
#'     \item{`relPDF`}{Relative proposal density (1 at `mu`).}
#'     \item{`importance_ratio`}{`likRatio / relPDF`; `NA` when `dOFV` is `NA`.}
#'     \item{`prob_resample`}{Normalised resampling probability; 0 for `NA` rows.}
#'   }
#' @noRd
sirCalcWeights <- function(samples, mu, covMat, dOFV) {
  n <- nrow(samples)
  p <- ncol(samples)
  checkmate::assertMatrix(samples, mode = "numeric", min.rows = 1L)
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, len = p)
  checkmate::assertMatrix(covMat, mode = "numeric", nrows = p, ncols = p)
  checkmate::assertNumeric(dOFV, len = n)

  L <- tryCatch(
    chol(covMat),
    error = function(e) {
      cli::cli_abort("{.arg covMat} is not positive definite.")
    }
  )

  # log(relPDF_i) = -0.5 * ||L^{-T}(x_i - mu)||^2; equals 0 when x_i = mu.
  log_rel_pdf <- vapply(
    seq_len(n),
    function(i) {
      z <- backsolve(L, samples[i, ] - mu)
      -0.5 * sum(z^2)
    },
    numeric(1L)
  )

  log_lik_ratio <- -0.5 * dOFV
  log_ir <- log_lik_ratio - log_rel_pdf

  valid <- is.finite(log_ir)
  if (!any(valid)) {
    cli::cli_abort(c(
      "No finite SIR importance weights could be computed.",
      "i" = "All OFV evaluations may have failed, or the proposal density is numerically degenerate."
    ))
  }

  max_log_ir <- max(log_ir[valid])
  scaled <- pmax(exp(log_ir[valid] - max_log_ir), .Machine$double.xmin)
  scaled_sum <- sum(scaled)
  if (!is.finite(scaled_sum) || scaled_sum <= 0) {
    cli::cli_abort(c(
      "SIR importance weights could not be normalized.",
      "i" = "Check failed OFV evaluations, proposal covariance, and resampling diagnostics."
    ))
  }

  prob <- numeric(n)
  prob[valid] <- scaled / scaled_sum

  rel_pdf <- exp(log_rel_pdf)
  lik_ratio <- exp(log_lik_ratio)
  ir <- exp(log_ir)
  ir[!valid] <- NA_real_

  data.frame(
    sample_id = seq_len(n),
    dOFV = dOFV,
    likelihood_ratio = lik_ratio,
    relPDF = rel_pdf,
    importance_ratio = ir,
    prob_resample = prob
  )
}

# Step 6 -----------------------------------------------------------------------

#' Weighted resample of parameter vectors with optional appearance cap
#'
#' Draws `m` rows from `samples` using `weights$prob_resample` as
#' probabilities.  `capResampling = 1` follows PsN's default and samples
#' without replacement.  Values greater than one allow limited replacement so
#' no original sample appears more than `capResampling` times.
#'
#' @param samples Numeric matrix. Original sampled parameter vectors (one row
#'   per sample).
#' @param weights Data frame returned by `sirCalcWeights()`.  Must contain a
#'   `prob_resample` column of length `nrow(samples)`.
#' @param m Positive integer. Number of resampled vectors to return.
#' @param capResampling Positive number. Default `1` samples without
#'   replacement; values `> 1` allow limited replacement.
#' @return A named list:
#'   \describe{
#'     \item{`samples`}{Numeric matrix with `m` rows selected from the
#'       originals.}
#'     \item{`resampleCounts`}{Integer vector of length `nrow(samples)`.
#'       `resampleCounts[i]` is how many times original row `i` was selected.}
#'     \item{`sampleOrder`}{Character vector of length `nrow(samples)`.
#'       Semicolon-separated draw positions for each original row.}
#'     \item{`selectionOrder`}{Integer vector of selected row indices, in draw
#'       order.}
#'   }
#' @noRd
sirResample <- function(samples, weights, m, capResampling = 1) {
  n <- nrow(samples)
  checkmate::assertMatrix(samples, mode = "numeric", min.rows = 1L)
  checkmate::assertDataFrame(weights)
  checkmate::assertNumeric(
    weights$prob_resample,
    len = n,
    lower = 0,
    any.missing = FALSE
  )
  checkmate::assertCount(m, positive = TRUE)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)

  m <- as.integer(m)
  prob <- weights$prob_resample
  valid <- is.finite(prob) & prob > 0
  cap <- as.integer(floor(capResampling))

  if (!any(valid)) {
    cli::cli_abort(c(
      "No samples have non-zero SIR resampling probability.",
      "i" = "Check OFV failures and importance weight diagnostics."
    ))
  }
  if (m > sum(valid) * cap) {
    cli::cli_abort(c(
      "Cannot draw {m} SIR resamples with {.arg capResampling} = {cap}.",
      "i" = "Only {sum(valid)} samples have non-zero probability."
    ))
  }

  if (cap <= 1L) {
    idx <- sample.int(n, m, replace = FALSE, prob = prob)
  } else {
    expanded_idx <- rep(seq_len(n), each = cap)
    expanded_prob <- rep(prob, each = cap)
    keep <- expanded_prob > 0 & is.finite(expanded_prob)
    idx <- sample(
      expanded_idx[keep],
      m,
      replace = FALSE,
      prob = expanded_prob[keep]
    )
  }

  resampled <- samples[idx, , drop = FALSE]
  rownames(resampled) <- NULL
  resample_counts <- tabulate(idx, nbins = n)
  sample_order <- vapply(
    seq_len(n),
    function(i) {
      pos <- which(idx == i)
      if (length(pos) == 0L) "" else paste(pos, collapse = ";")
    },
    character(1L)
  )

  list(
    samples = resampled,
    resampleCounts = resample_counts,
    sampleOrder = sample_order,
    selectionOrder = idx
  )
}

.sirBuildRawResults <- function(
  paramMat,
  weights,
  dOFV,
  resampled,
  mu,
  capResampling = 1
) {
  checkmate::assertMatrix(paramMat, mode = "numeric", min.rows = 1L)
  checkmate::assertDataFrame(weights, nrows = nrow(paramMat))
  checkmate::assertNumeric(dOFV, len = nrow(paramMat))
  checkmate::assertList(resampled)
  checkmate::assertNumeric(mu, len = ncol(paramMat))
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)

  param_names <- colnames(paramMat)
  cap <- as.integer(floor(capResampling))
  raw_df <- data.frame(
    sample_id = seq_len(nrow(paramMat)),
    as.data.frame(paramMat, check.names = FALSE),
    check.names = FALSE
  )
  raw_df$dOFV <- dOFV
  raw_df$deltaofv <- dOFV
  raw_df$likelihood_ratio <- weights$likelihood_ratio
  raw_df$relPDF <- weights$relPDF
  raw_df$importance_ratio <- weights$importance_ratio
  raw_df$IR <- weights$importance_ratio
  raw_df$probability_resample <- weights$prob_resample
  raw_df$prob <- weights$prob_resample

  expanded <- lapply(seq_len(nrow(raw_df)), function(i) {
    out <- raw_df[rep(i, cap), , drop = FALSE]
    orders <- resampled$sampleOrder[[i]]
    orders <- if (identical(orders, "")) {
      character(0L)
    } else {
      strsplit(orders, ";", fixed = TRUE)[[1L]]
    }
    n_selected <- min(length(orders), cap)
    out$resamples <- 0L
    out$sample_order <- ""
    if (n_selected > 0L) {
      idx <- seq_len(n_selected)
      out$resamples[idx] <- 1L
      out$sample_order[idx] <- orders[idx]
    }
    out
  })
  expanded <- do.call(rbind, expanded)
  rownames(expanded) <- NULL

  mu_mat <- matrix(
    mu[param_names],
    nrow = 1L,
    dimnames = list(NULL, param_names)
  )
  mu_df <- data.frame(
    sample_id = 0L,
    as.data.frame(mu_mat, check.names = FALSE),
    check.names = FALSE
  )
  mu_df$dOFV <- 0
  mu_df$deltaofv <- 0
  mu_df$likelihood_ratio <- 1
  mu_df$relPDF <- 1
  mu_df$importance_ratio <- NA_real_
  mu_df$IR <- NA_real_
  mu_df$probability_resample <- 0
  mu_df$prob <- 0
  mu_df$resamples <- 0L
  mu_df$sample_order <- ""

  out <- rbind(mu_df, expanded)
  rownames(out) <- NULL
  out
}
