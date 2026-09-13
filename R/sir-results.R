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
#' @return A data frame with columns `param`, `estimate`, `mean`, `sd`, `rse`,
#'   `rse_sd_scale`, and the percentiles `p2.5`, `p5`, `p10`, `p30`, `p50`,
#'   `p70`, `p90`, `p95`, `p97.5`. The percentile set is PsN's, derived from
#'   prediction intervals 0, 40, 80, 90 and 95.
#'
#'   `rse` is a **percentage**; PsN reports the same quantity as a fraction.
#'   `rse_sd_scale` is `rse / 2`, the RSE of a variance expressed on the
#'   standard-deviation scale, and is `NA` for parameters that are not
#'   variances. Note this differs from PsN, which halves everything that is
#'   not a NONMEM THETA: nlmixr2 parameterises residual error on the SD scale
#'   already, so halving `add.sd` would rescale a quantity that needs no
#'   rescaling.
#'
#'   Attributes `covMatrix`, `corMatrix` and `sdCorMatrix` hold the empirical
#'   covariance, correlation, and standard-deviation/correlation matrices.
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

  # PsN's percentile set, derived from prediction intervals 0, 40, 80, 90, 95:
  # 2.5, 5, 10, 30, 50, 70, 90, 95, 97.5.
  probs <- .sirPercentileProbs()
  qmat <- t(apply(
    resampledMat,
    2L,
    stats::quantile,
    probs = probs,
    na.rm = TRUE,
    names = FALSE
  ))
  colnames(qmat) <- paste0("p", .sirPercentileLabels())

  mean_vals <- colMeans(resampledMat, na.rm = TRUE)
  sd_vals <- apply(resampledMat, 2L, stats::sd, na.rm = TRUE)
  rse <- ifelse(
    is.finite(estimate) & estimate != 0,
    sd_vals / abs(estimate) * 100,
    NA_real_
  )

  # RSE expressed on the standard-deviation scale, for parameters estimated as
  # variances. PsN halves the RSE for everything that is not a NONMEM THETA.
  # That rule cannot be carried over literally: nlmixr2 parameterises residual
  # error on the SD scale already (add.sd is a standard deviation, not a
  # variance), so halving it would understate a quantity that needs no
  # rescaling. Only OMEGA elements are halved here.
  ps <- .sirParamSpace(fit)
  kind <- ps$kind[match(param_names, ps$sirName)]
  onVarianceScale <- !is.na(kind) & kind %in% c("omegaDiag", "omegaOffdiag")
  rse_sd_scale <- ifelse(onVarianceScale, rse / 2, NA_real_)

  out <- data.frame(
    param = param_names,
    estimate = unname(estimate),
    mean = unname(mean_vals),
    sd = unname(sd_vals),
    rse = unname(rse),
    rse_sd_scale = unname(rse_sd_scale),
    qmat,
    row.names = NULL,
    check.names = FALSE
  )
  structure(
    out,
    covMatrix = stats::cov(resampledMat),
    corMatrix = stats::cor(resampledMat),
    sdCorMatrix = .sirSdCorMatrix(resampledMat),
    rseUnits = "percent"
  )
}

# PsN reports percentiles derived from a set of prediction intervals rather
# than a flat list: interval 0 contributes the median, and each other interval
# pi contributes the pair (100 - pi)/2 and 100 - (100 - pi)/2.
.sirPercentileLabels <- function(predictionIntervals = c(0, 40, 80, 90, 95)) {
  out <- unlist(lapply(sort(predictionIntervals), function(pi) {
    if (pi == 0) {
      return(50)
    }
    c((100 - pi) / 2, 100 - (100 - pi) / 2)
  }))
  sort(unique(out))
}

.sirPercentileProbs <- function(predictionIntervals = c(0, 40, 80, 90, 95)) {
  .sirPercentileLabels(predictionIntervals) / 100
}

# Standard deviations on the diagonal, correlations off it -- PsN's sdcorr
# form, and what the RSE/correlation diagnostic plot reads.
.sirSdCorMatrix <- function(resampledMat) {
  cm <- stats::cov(resampledMat)
  out <- stats::cov2cor(cm)
  diag(out) <- sqrt(diag(cm))
  out
}

.sirSummarizeResamples <- function(resampledMat, fit) {
  sirSummary(resampledMat, fit)
}

# PsN's summary_iterations.csv column names, so a PsN-literate reader and any
# downstream tooling can read either file. The nlmixr2sir-specific columns
# (rejection counts, mean dOFV) are kept and appended after them.
.sirPsnIterationColumns <- function() {
  c(
    iter = "iteration",
    nSamples = "commandline.samples",
    nAttempted = "attempted.samples",
    nSuccessful = "successful.samples",
    nResample = "commandline.resamples",
    nResampled = "actual.resamples",
    requested_sample_resample_ratio = "requested.ratio",
    actual_sample_resample_ratio = "actual.ratio",
    nNegativeDOFV = "negative.dOFV",
    minDOFV = "minimum.sample.ofv"
  )
}

.sirWriteIterationSummary <- function(iterSummary, directory) {
  out <- iterSummary
  out$requested_sample_resample_ratio <- out$nSamples / out$nResample
  out$actual_sample_resample_ratio <- out$nSuccessful / out$nResampled

  map <- .sirPsnIterationColumns()
  psn <- intersect(names(map), names(out))
  rest <- setdiff(names(out), psn)
  out <- out[, c(psn, rest), drop = FALSE]
  names(out)[seq_along(psn)] <- unname(map[psn])

  utils::write.csv(
    out,
    file.path(directory, "summary_iterations.csv"),
    row.names = FALSE
  )
  invisible(out)
}

# PsN writes the empirical covariance of the final resampled vectors as
# <model>_sir.cov, and the sd/correlation form alongside it. Both were
# previously attached to the result as attributes only.
.sirWriteCovMatrices <- function(summary, directory, fitName = "sir") {
  cov_mat <- attr(summary, "covMatrix", exact = TRUE)
  if (is.null(cov_mat)) {
    return(invisible(NULL))
  }
  sdcor <- attr(summary, "sdCorMatrix", exact = TRUE)

  .write <- function(m, path) {
    df <- data.frame(NAME = rownames(m), m, check.names = FALSE)
    utils::write.table(
      df,
      path,
      row.names = FALSE,
      quote = FALSE,
      sep = "	"
    )
    path
  }

  written <- .write(cov_mat, file.path(directory, paste0(fitName, "_sir.cov")))
  if (!is.null(sdcor)) {
    written <- c(
      written,
      .write(sdcor, file.path(directory, paste0(fitName, "_sir.sdcorr")))
    )
  }
  invisible(written)
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
