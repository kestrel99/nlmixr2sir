# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# Summaries, canonical raw results, and on-disk iteration artifacts.

.sirOriginalEstimates <- function(fit, paramNames) {
  vals <- setNames(rep(NA_real_, length(paramNames)), paramNames)
  ps <- .sirParamSpace(fit)
  shared <- intersect(paramNames, ps$sirName)
  vals[shared] <- ps$est[match(shared, ps$sirName)]
  vals
}

.sirCanonicalRawResults <- function(fit, fitName, resampledMat) {
  if (is.data.frame(resampledMat)) {
    resampledMat <- as.matrix(resampledMat)
  }
  checkmate::assertMatrix(resampledMat, mode = "numeric", min.rows = 1L)

  schema <- nlmixr2utils::rawResultsSchema(fit)
  param_cols <- c(schema$thetaCols, schema$omegaCols, schema$sigmaCols)
  na_se <- if (length(param_cols) > 0L) {
    stats::setNames(rep(NA_real_, length(param_cols)), param_cols)
  } else {
    numeric(0)
  }
  ps <- .sirParamSpace(fit)
  has_omega <- any(ps$kind %in% c("omegaDiag", "omegaOffdiag"))
  theta_cols <- intersect(schema$thetaCols, colnames(resampledMat))

  rows <- vector("list", nrow(resampledMat) + 1L)
  rows[[1L]] <- nlmixr2utils::rawResultsRow(
    fit = fit,
    source = "sir",
    hypothesis = "reference",
    sample = 0L,
    modelLabel = fitName,
    role = "reference",
    schema = schema
  )

  for (i in seq_len(nrow(resampledMat))) {
    row_vals <- resampledMat[i, , drop = TRUE]
    theta_vals <- if (length(theta_cols) > 0L) {
      stats::setNames(as.numeric(row_vals[theta_cols]), theta_cols)
    } else {
      NULL
    }
    omega_mat <- if (has_omega) {
      .sirReconstructOmega(ps, row_vals, fit$omega)
    } else {
      NULL
    }

    rows[[i + 1L]] <- nlmixr2utils::rawResultsRow(
      fit = fit,
      source = "sir",
      hypothesis = "resampled",
      sample = i,
      modelLabel = fitName,
      role = "sample",
      objf = NA_real_,
      minimizationSuccessful = NA_integer_,
      covarianceStepSuccessful = NA_integer_,
      estimateNearBoundary = NA_integer_,
      significantDigits = NA_real_,
      conditionNumber = NA_real_,
      theta = theta_vals,
      omega = omega_mat,
      se = na_se,
      schema = schema
    )
  }

  do.call(rbind, rows)
}

#' Summarize final SIR resampled parameter vectors
#'
#' `sirSummary()` computes empirical summary statistics from final SIR
#' resampled parameter vectors. It reports the reference estimate from the
#' input fit, empirical standard deviation and RSE, and percentile intervals
#' for each sampled parameter. The empirical covariance and correlation
#' matrices are attached as attributes.
#'
#' @param resampledMat Numeric matrix of final resampled parameter vectors, one
#'   row per vector and one column per parameter.
#' @param fit An nlmixr2 fit object used to obtain reference parameter
#'   estimates.
#' @return A data frame with columns `param`, `estimate`, `sd`, `rse`, `p2.5`,
#'   `p5`, `p25`, `p50`, `p75`, `p95`, and `p97.5`. Attributes `covMatrix`
#'   and `corMatrix` contain empirical covariance and correlation matrices.
#' @export
sirSummary <- function(resampledMat, fit) {
  if (is.data.frame(resampledMat)) {
    resampledMat <- as.matrix(resampledMat)
  }
  checkmate::assertMatrix(resampledMat, mode = "numeric", min.rows = 2L)
  checkmate::assertClass(fit, "nlmixr2FitCore")
  if (is.null(colnames(resampledMat))) {
    cli::cli_abort("{.arg resampledMat} must have parameter column names.")
  }

  param_names <- colnames(resampledMat)
  estimate <- .sirOriginalEstimates(fit, param_names)
  probs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  qmat <- t(apply(
    resampledMat,
    2L,
    stats::quantile,
    probs = probs,
    na.rm = TRUE,
    names = FALSE
  ))
  colnames(qmat) <- c("p2.5", "p5", "p25", "p50", "p75", "p95", "p97.5")
  sd_vals <- apply(resampledMat, 2L, stats::sd, na.rm = TRUE)
  rse <- ifelse(
    is.finite(estimate) & estimate != 0,
    sd_vals / abs(estimate) * 100,
    NA_real_
  )
  out <- data.frame(
    param = param_names,
    estimate = unname(estimate),
    sd = unname(sd_vals),
    rse = unname(rse),
    qmat,
    row.names = NULL,
    check.names = FALSE
  )
  structure(
    out,
    covMatrix = stats::cov(resampledMat),
    corMatrix = stats::cor(resampledMat)
  )
}

.sirSummarizeResamples <- function(resampledMat, fit) {
  sirSummary(resampledMat, fit)
}

.sirWriteIterationSummary <- function(iterSummary, directory) {
  out <- iterSummary
  out$requested_sample_resample_ratio <- out$nSamples / out$nResample
  out$actual_sample_resample_ratio <- out$nSuccessful / out$nResampled
  utils::write.csv(
    out,
    file.path(directory, "summary_iterations.csv"),
    row.names = FALSE
  )
  invisible(out)
}

.sirWriteRejectionSummary <- function(iterSummary, directory) {
  s <- iterSummary[nrow(iterSummary), , drop = FALSE]
  lines <- c(
    "SIR sample rejection summary",
    sprintf("Iteration: %s", s$iter),
    sprintf("Requested samples: %s", s$nSamples),
    sprintf("Attempted samples: %s", s$nAttempted),
    sprintf("Raw draw attempts: %s", s$nDrawAttempts),
    sprintf("Collected samples: %s", s$nCollected),
    sprintf("Successful OFV evaluations: %s", s$nSuccessful),
    sprintf("Failed OFV evaluations: %s", s$nFailed),
    sprintf("Inverse Box-Cox rejections: %s", s$inverseRejected),
    sprintf("Theta bound rejections: %s", s$thetaRejected),
    sprintf("Omega positive-definiteness rejections: %s", s$omegaRejected),
    sprintf("Sigma bound rejections: %s", s$sigmaRejected)
  )
  writeLines(lines, file.path(directory, "sample_rejection_summary.txt"))
  invisible(lines)
}
