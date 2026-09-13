# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# One SIR iteration, and the PsN-style attempted-sample adjustment.

sirRunIteration <- function(
  fit,
  mu,
  proposalCov,
  nSamples,
  nResample,
  iterNum,
  capResampling = 1,
  recenter = TRUE,
  boxcox = TRUE,
  directory = NULL,
  workers = NULL,
  rxThreads = NULL,
  boxcoxState = NULL,
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  omegaFallback = c("cov", "wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL,
  requestedSamples = nSamples,
  isLastIteration = FALSE
) {
  omegaFallback <- match.arg(omegaFallback)

  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(proposalCov, mode = "numeric")
  checkmate::assertCount(nSamples, positive = TRUE)
  checkmate::assertCount(nResample, positive = TRUE)
  checkmate::assertCount(iterNum, positive = TRUE)
  checkmate::assertCount(requestedSamples, positive = TRUE)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)
  checkmate::assertFlag(recenter)
  checkmate::assertFlag(boxcox)
  checkmate::assertFlag(isLastIteration)
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)
  checkmate::assertNumber(sigmaFallbackRse, lower = 0, finite = TRUE)
  if (!is.null(directory)) {
    checkmate::assertString(directory)
  }

  proposal <- .sirInitialProposal(
    fit = fit,
    mu = mu,
    proposalCov = proposalCov,
    thetaInflation = thetaInflation,
    omegaInflation = omegaInflation,
    sigmaInflation = sigmaInflation,
    capCorrelation = capCorrelation,
    omegaFallback = omegaFallback,
    sigmaFallbackRse = sigmaFallbackRse,
    omegaDf = omegaDf
  )

  ps <- proposal$paramSpace
  param_names <- proposal$paramNames
  theta_names <- ps$sirName[ps$kind == "theta"]
  sigma_names <- ps$sirName[ps$kind == "sigma"]
  bounds <- .sirParamBounds(ps)

  # ---- 1-3. Sample full proposal vectors and reject invalid draws ----
  bc_mu <- .sirBcTransformMu(proposal$mu, boxcoxState)
  sampled <- .sirSampleFullProposal(
    mu = bc_mu,
    covMat = proposal$covMat,
    n = nSamples,
    lower = bounds$lower,
    upper = bounds$upper,
    ps = ps,
    baseOmega = fit$omega,
    thetaNames = theta_names,
    sigmaNames = sigma_names,
    boxcoxState = boxcoxState
  )

  param_mat <- sampled$samples
  samples_for_weights <- sampled$samplesForPdf
  n_collected <- nrow(param_mat)
  if (n_collected == 0L) {
    cli::cli_abort(c(
      "No valid SIR samples were collected.",
      "i" = "Check proposal covariance, bounds, and omega positive-definiteness diagnostics."
    ))
  }

  # ---- 4. Evaluate OFV ----
  ofv_vals <- sirEvalOFV(
    fit,
    param_mat,
    workers = workers,
    rxThreads = rxThreads
  )
  dofv <- ofv_vals - fit$objf

  # ---- 5. Handle failures ----
  n_failed <- sum(is.na(dofv))
  n_success <- n_collected - n_failed
  if (n_success == 0L) {
    eval_errors <- attr(ofv_vals, "evalErrors")
    cli::cli_abort(c(
      "All SIR OFV evaluations failed.",
      "i" = "No valid importance weights can be computed.",
      if (length(eval_errors) > 0L) {
        c("x" = "First error: {eval_errors[[1L]]}")
      }
    ))
  }
  n_resample_adj <- .sirAdjustedResamples(
    requestedResamples = nResample,
    requestedSamples = requestedSamples,
    successfulCount = n_success
  )
  if (n_resample_adj != nResample) {
    cli::cli_warn(c(
      "{n_success}/{requestedSamples} requested SIR samples had usable OFV evaluations.",
      "i" = "Adjusting nResample to {n_resample_adj}."
    ))
  }
  if (n_resample_adj < 1L) {
    cli::cli_abort(c(
      "Turnout scaled the resample count below one.",
      "i" = "{n_success}/{requestedSamples} samples usable, {nResample} requested resamples.",
      "i" = "Increase {.arg nSamples} or {.arg nResample}."
    ))
  }

  # ---- 6. Compute weights in the same full proposal scale used for sampling ----
  weights <- sirCalcWeights(
    samples_for_weights,
    bc_mu,
    proposal$covMat,
    dOFV = dofv
  )

  # ---- 7. Resample ----
  resampled <- sirResample(
    param_mat,
    weights,
    m = n_resample_adj,
    capResampling = capResampling
  )

  # ---- 8. Recenter ----
  new_mu <- proposal$mu
  if (recenter) {
    valid_dofv <- ifelse(is.na(dofv), Inf, dofv)
    best_idx <- which.min(valid_dofv)
    if (isTRUE(valid_dofv[best_idx] < 0)) {
      new_mu <- param_mat[best_idx, ]
      cli::cli_inform(
        "  Iter {iterNum}: recentered mu (dOFV = {round(dofv[best_idx], 4)})."
      )
    }
  } else if (any(dofv < 0, na.rm = TRUE)) {
    cli::cli_warn(
      "At least one SIR sample had negative dOFV, but {.arg recenter} is FALSE."
    )
  }

  # ---- 9. Update full proposal in original scale ----
  updated <- sirUpdateProposal(
    resampled$samples[, param_names, drop = FALSE],
    boxcox = boxcox && !isLastIteration,
    capCorrelation = capCorrelation
  )
  new_cov <- updated$covMat
  new_bc_state <- updated$boxcoxParams # NULL when boxcox = FALSE

  # ---- 10. Raw results data frame ----
  raw_df <- .sirBuildRawResults(
    paramMat = param_mat,
    weights = weights,
    dOFV = dofv,
    resampled = resampled,
    mu = proposal$mu,
    capResampling = capResampling
  )

  # ---- Summary ----
  iter_summary <- data.frame(
    iter = iterNum,
    nSamples = requestedSamples,
    nAttempted = nSamples,
    nDrawAttempts = sampled$nAttempted,
    nCollected = n_collected,
    nSuccessful = n_success,
    nFailed = n_failed,
    nResample = nResample,
    nResampled = nrow(resampled$samples),
    minDOFV = if (all(is.na(dofv))) NA_real_ else min(dofv, na.rm = TRUE),
    meanDOFV = if (all(is.na(dofv))) NA_real_ else mean(dofv, na.rm = TRUE),
    nNegativeDOFV = sum(dofv < 0, na.rm = TRUE),
    thetaRejected = sampled$thetaRejected,
    omegaRejected = sampled$omegaRejected,
    sigmaRejected = sampled$sigmaRejected,
    inverseRejected = sampled$inverseRejected
  )

  list(
    resampledMat = resampled$samples,
    newMu = new_mu,
    newCov = new_cov,
    iterSummary = iter_summary,
    boxcoxState = new_bc_state,
    rawResults = raw_df
  )
}

# Port of PsN update_attempted_samples() (lib/tool/sir.pm). Compensates the
# next iteration's sample count for samples lost to failed OFV evaluation, so
# the requested count is what actually survives.
#
# Matches PsN exactly: triggers on loss only, at turnout <= 0.95 inclusive, and
# rounds half away from zero. The previous implementation used a strict `<`,
# `ceiling()`, and a max() clamp against the requested count, which gave 112
# where PsN's own oracle says 111.
.sirAdjustedAttemptedSamples <- function(
  requestedSamples,
  previousAttempted = NULL,
  previousSuccessful = NULL
) {
  checkmate::assertCount(requestedSamples, positive = TRUE)
  if (is.null(previousAttempted) || is.null(previousSuccessful)) {
    return(as.integer(requestedSamples))
  }
  checkmate::assertCount(previousAttempted, positive = TRUE)
  checkmate::assertCount(previousSuccessful, positive = TRUE)

  previousTurnout <- previousSuccessful / previousAttempted
  if (previousTurnout > 1) {
    cli::cli_abort(c(
      "More successful samples than attempted in the previous iteration.",
      "i" = "{previousSuccessful} successful of {previousAttempted} attempted."
    ))
  }
  if (previousTurnout <= 0.95) {
    return(.sirRound(requestedSamples / previousTurnout))
  }
  as.integer(requestedSamples)
}

# Port of PsN update_actual_resamples() (lib/tool/sir.pm). Scales the resample
# count by this iteration's turnout, on a gain *or* a loss of at least 5%.
#
# `turnout` is measured against the originally requested sample count, not the
# compensated attempted count -- PsN's oracle pins this: 109 successful of a
# requested 100 gives turnout 1.09 even though 109 were attempted.
.sirAdjustedResamples <- function(
  requestedResamples,
  requestedSamples,
  successfulCount
) {
  checkmate::assertCount(requestedResamples, positive = TRUE)
  checkmate::assertCount(requestedSamples, positive = TRUE)
  checkmate::assertCount(successfulCount)

  turnout <- successfulCount / requestedSamples
  if (abs(turnout - 1) >= 0.05) {
    return(.sirRound(requestedResamples * turnout))
  }
  as.integer(requestedResamples)
}
