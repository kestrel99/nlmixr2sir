# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# runSIR(): the top-level driver.

#' Run sampling importance resampling for an nlmixr2 fit
#'
#' `runSIR()` runs one or more sampling importance resampling iterations using
#' an nlmixr2 fit as the reference model. It follows the PsN SIR workflow where
#' practical for nlmixr2: sample parameter vectors from a proposal covariance,
#' evaluate them with population parameters fixed, compute importance weights,
#' resample, and use the empirical resampled covariance as the next proposal.
#'
#' @param fit An nlmixr2 fit object with a covariance matrix.
#' @param nSamples Integer vector. Requested number of samples per iteration.
#' @param nResample Integer vector. Requested number of resamples per
#'   iteration. Must have the same length as `nSamples`.
#' @param directory Output directory. If `NULL`, a numbered
#'   `<fitName>_sir_<N>` directory is created.
#' @param fitName Optional fit label used when creating an automatic output
#'   directory and canonical raw-results metadata. When `NULL` (default), the
#'   label is derived from the expression supplied to `fit`.
#' @param control A [runSIRControl()] object holding everything that tunes
#'   how the run behaves: inflation, caps, recentering, Box-Cox, parallelism,
#'   resume behaviour, and the OMEGA fallback.
#' @param ... Reserved for future PsN-compatible inputs. Passing a run
#'   setting here is an error; put it in `control` instead.
#' @return A data frame of final SIR summary statistics with class
#'   `c("nlmixr2SIR", "data.frame")`. Attributes include
#'   `iterationSummary`, `iterations`, `resampledMat`, `covMatrix`,
#'   `corMatrix`, and `outputDir`.
#' @examples
#' \dontrun{
#' sir <- runSIR(
#'   fit,
#'   nSamples = c(1000, 1000, 1000),
#'   nResample = c(200, 400, 500),
#'   control = runSIRControl(workers = 4, rxThreads = 2)
#' )
#' }
#' @seealso [runSIRControl()] for the run settings.
#' @export
runSIR <- function(
  fit,
  nSamples = c(1000, 1000, 1000, 2000, 2000),
  nResample = c(200, 400, 500, 1000, 1000),
  directory = NULL,
  fitName = NULL,
  control = runSIRControl(),
  ...
) {
  dots <- list(...)
  if (length(dots) > 0L) {
    cli::cli_abort(c(
      "Unsupported SIR argument(s): {.arg {names(dots)}}.",
      "i" = "Run settings now live in {.fn runSIRControl}."
    ))
  }
  checkmate::assertClass(control, "runSIRControl")

  thetaInflation <- control$thetaInflation
  omegaInflation <- control$omegaInflation
  sigmaInflation <- control$sigmaInflation
  capCorrelation <- control$capCorrelation
  capResampling <- control$capResampling
  recenter <- control$recenter
  boxcox <- control$boxcox
  workers <- control$workers
  rxThreads <- control$rxThreads
  recover <- control$recover
  addIterations <- control$addIterations
  omegaFallback <- control$omegaFallback
  sigmaFallbackRse <- control$sigmaFallbackRse
  omegaDf <- control$omegaDf

  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertIntegerish(
    nSamples,
    lower = 1,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertIntegerish(
    nResample,
    lower = 1,
    any.missing = FALSE,
    len = length(nSamples)
  )
  if (is.null(fitName)) {
    fitName <- nlmixr2utils::deriveFitName(substitute(fit))
  }

  nSamples <- as.integer(nSamples)
  nResample <- as.integer(nResample)
  ps <- .sirParamSpace(fit)
  run_dir <- nlmixr2utils::resolveRunDir(
    "sir",
    fitName,
    restart = !(recover || addIterations),
    outputDir = directory
  )
  output_dir <- run_dir$path
  if (identical(run_dir$mode, "overwrite") && dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  master_seed <- nlmixr2utils::withRunSeed(output_dir, prefix = "sir")
  saved_state <- if (recover || addIterations) {
    nlmixr2utils::readRunState(output_dir, "sir")
  } else {
    NULL
  }

  if (!is.null(saved_state) && isTRUE(addIterations)) {
    iter_offset <- saved_state$completedIterations
    iter_index <- seq_along(nSamples)
    iter_numbers <- iter_offset + iter_index
    mu <- saved_state$nextMu
    proposal_cov <- saved_state$nextCov
    boxcox_state <- saved_state$nextBoxcoxState
    iter_results <- saved_state$iterations
    iter_summary <- saved_state$iterationSummary
    prev_attempted <- tail(iter_summary$nAttempted, 1L)
    prev_successful <- tail(iter_summary$nSuccessful, 1L)
  } else if (!is.null(saved_state) && isTRUE(recover)) {
    completed <- saved_state$completedIterations
    if (completed >= length(nSamples) && !is.null(saved_state$result)) {
      cli::cli_inform("Recovered completed SIR run from {.path {output_dir}}.")
      return(saved_state$result)
    }
    iter_numbers <- seq.int(completed + 1L, length(nSamples))
    iter_index <- iter_numbers
    mu <- saved_state$nextMu
    proposal_cov <- saved_state$nextCov
    boxcox_state <- saved_state$nextBoxcoxState
    iter_results <- saved_state$iterations
    iter_summary <- saved_state$iterationSummary
    prev_attempted <- tail(iter_summary$nAttempted, 1L)
    prev_successful <- tail(iter_summary$nSuccessful, 1L)
  } else {
    iter_numbers <- seq_along(nSamples)
    iter_index <- iter_numbers
    initial <- .sirResolveInitialProposal(fit, ps, control)
    proposal_cov <- initial$covMat
    # The raw-results route derives its own centre and Box-Cox state from the
    # supplied vectors, the way PsN's iteration 0 does; every other route
    # centres on the fit's own estimates.
    mu <- initial$mu %||% .sirProposalMu(fit, ps)
    boxcox_state <- initial$boxcoxState
    if (!identical(initial$source, "cov")) {
      cli::cli_inform(
        "Initial SIR proposal built from {.arg {initial$source}}, not {.code fit$cov}."
      )
    }
    iter_results <- list()
    iter_summary <- data.frame()
    prev_attempted <- NULL
    prev_successful <- NULL
  }

  for (j in seq_along(iter_numbers)) {
    iter_num <- iter_numbers[[j]]
    schedule_idx <- iter_index[[j]]
    requested_samples <- nSamples[[schedule_idx]]
    attempted_samples <- .sirAdjustedAttemptedSamples(
      requested_samples,
      previousAttempted = prev_attempted,
      previousSuccessful = prev_successful
    )
    is_last <- j == length(iter_numbers)
    # Inflation widens the *initial* proposal only, as in PsN. From the second
    # iteration on the proposal is the previous empirical covariance, which
    # must not be re-inflated each time.
    inflate_now <- j == 1L && is.null(saved_state)

    cli::cli_inform(
      "Running SIR iteration {iter_num}: {attempted_samples} attempted samples, {nResample[[schedule_idx]]} requested resamples." # nolint: line_length_linter.
    )
    iter_res <- nlmixr2utils::withRunSeed(
      output_dir,
      key = paste0("sir-iteration-", iter_num),
      prefix = "sir",
      expr = sirRunIteration(
        fit = fit,
        mu = mu,
        proposalCov = proposal_cov,
        nSamples = attempted_samples,
        requestedSamples = requested_samples,
        nResample = nResample[[schedule_idx]],
        iterNum = iter_num,
        capResampling = capResampling,
        recenter = recenter,
        boxcox = boxcox,
        directory = output_dir,
        workers = workers,
        rxThreads = rxThreads,
        boxcoxState = boxcox_state,
        thetaInflation = if (inflate_now) thetaInflation else 1,
        omegaInflation = if (inflate_now) omegaInflation else 1,
        sigmaInflation = if (inflate_now) sigmaInflation else 1,
        capCorrelation = capCorrelation,
        omegaFallback = omegaFallback,
        sigmaFallbackRse = sigmaFallbackRse,
        omegaDf = omegaDf,
        isLastIteration = is_last
      )
    )

    iter_results[[as.character(iter_num)]] <- iter_res
    iter_summary <- rbind(iter_summary, iter_res$iterSummary)
    .sirWriteIterationSummary(iter_summary, output_dir)
    .sirWriteRejectionSummary(iter_summary, output_dir)

    mu <- iter_res$newMu
    proposal_cov <- iter_res$newCov
    boxcox_state <- iter_res$boxcoxState
    prev_attempted <- iter_res$iterSummary$nAttempted
    prev_successful <- iter_res$iterSummary$nSuccessful

    nlmixr2utils::writeRunState(
      output_dir,
      list(
        completedIterations = iter_num,
        nextMu = mu,
        nextCov = proposal_cov,
        nextBoxcoxState = boxcox_state,
        iterations = iter_results,
        iterationSummary = iter_summary,
        result = NULL
      ),
      "sir"
    )
  }

  final_iter <- iter_results[[length(iter_results)]]
  summary_df <- sirSummary(final_iter$resampledMat, fit)
  cov_mat <- attr(summary_df, "covMatrix")
  cor_mat <- attr(summary_df, "corMatrix")
  sdcor_mat <- attr(summary_df, "sdCorMatrix")
  utils::write.csv(
    summary_df,
    file.path(output_dir, "sir_results.csv"),
    row.names = FALSE
  )
  .sirWriteCovMatrices(summary_df, output_dir, fitName = fitName)
  raw_results <- .sirCanonicalRawResults(fit, fitName, final_iter$resampledMat)
  nlmixr2utils::writeRawResults(raw_results, output_dir)

  class(summary_df) <- c("nlmixr2SIR", "data.frame")
  attr(summary_df, "iterationSummary") <- iter_summary
  attr(summary_df, "iterations") <- iter_results
  attr(summary_df, "resampledMat") <- final_iter$resampledMat
  attr(summary_df, "covMatrix") <- cov_mat
  attr(summary_df, "corMatrix") <- cor_mat
  attr(summary_df, "sdCorMatrix") <- sdcor_mat
  attr(summary_df, "outputDir") <- output_dir
  attr(summary_df, "fitName") <- fitName
  attr(summary_df, "rawResults") <- raw_results
  attr(summary_df, "seed") <- master_seed
  attr(summary_df, "call") <- match.call()

  nlmixr2utils::writeRunState(
    output_dir,
    list(
      completedIterations = tail(iter_summary$iter, 1L),
      nextMu = mu,
      nextCov = proposal_cov,
      nextBoxcoxState = boxcox_state,
      iterations = iter_results,
      iterationSummary = iter_summary,
      result = summary_df
    ),
    "sir"
  )

  summary_df
}
