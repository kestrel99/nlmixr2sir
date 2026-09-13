#' Configure an `nlmixr2sir` run
#'
#' `runSIRControl()` constructs and validates the control object used by
#' [runSIR()]. It holds everything that tunes *how* SIR runs, leaving
#' [runSIR()] itself to take only the fit, the sampling schedule, and where
#' output goes.
#'
#' @param thetaInflation,omegaInflation,sigmaInflation Non-negative variance
#'   multipliers applied to the THETA, OMEGA, and residual-error blocks of the
#'   **initial** proposal. From the second iteration onwards the proposal is
#'   the previous iteration's empirical covariance and is not re-inflated.
#' @param capCorrelation Numeric in `[0, 1]`. Maximum absolute proposal
#'   correlation after covariance construction and updates.
#' @param capResampling Numeric greater than or equal to one. `1` resamples
#'   without replacement; larger values allow limited replacement.
#' @param recenter Logical. If `TRUE`, recenter the next proposal on the best
#'   sampled vector when any sample has negative dOFV.
#' @param boxcox Logical. If `TRUE`, use Box-Cox transformed empirical
#'   covariances for non-final iterations.
#' @param workers `NULL`, `"auto"`, `1`, or a positive integer. Controls
#'   parallel OFV evaluation through `future`. `NULL` leaves the current
#'   `future::plan()` unchanged, `1` forces sequential execution, a positive
#'   integer temporarily uses a multisession plan, and `"auto"` uses
#'   `future::availableCores(omit = 1L)`. SIR iterations are always sequential,
#'   because each iteration builds the next iteration's proposal.
#' @param rxThreads Integer, `"auto"`, or `NULL`; rxode2 OpenMP threads per
#'   worker. `NULL` (the default) uses the current `rxode2::getRxThreads()`
#'   value for every worker; `"auto"` divides the core count evenly across
#'   workers. Whenever `workers > 1`, `workers * rxThreads` must not exceed the
#'   machine's core count, since each worker is a separate process running its
#'   own rxode2 thread pool.
#' @param recover Logical. If `TRUE` and the output directory holds
#'   `sir_state.rds`, resume from the last completed iteration when possible.
#' @param addIterations Logical. If `TRUE`, append the supplied schedule after
#'   an existing completed state in the output directory.
#' @param omegaFallback How OMEGA uncertainty is obtained. `"cov"` (the
#'   default) takes it from `fit$cov`, including its correlations with THETA,
#'   and degrades to `"wishart"` automatically when `fit$cov` does not carry
#'   OMEGA. `"wishart"` always uses the Wishart-style approximation, which
#'   gives a block-diagonal proposal.
#' @param sigmaFallbackRse Percent relative standard error used for
#'   residual-error uncertainty when neither `fit$cov` nor `fit$parFixedDf`
#'   reports a standard error.
#' @param omegaDf Optional degrees of freedom for the Wishart-style OMEGA
#'   fallback; defaults to `nsub - 1`. Unused when OMEGA comes from `fit$cov`.
#'
#' @return An object of class `runSIRControl`.
#' @examples
#' runSIRControl(thetaInflation = 2, workers = 4, rxThreads = 2)
#' @export
runSIRControl <- function(
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  capResampling = 1,
  recenter = TRUE,
  boxcox = TRUE,
  workers = NULL,
  rxThreads = NULL,
  recover = TRUE,
  addIterations = FALSE,
  omegaFallback = c("cov", "wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL
) {
  omegaFallback <- match.arg(omegaFallback)

  checkmate::assertNumeric(
    thetaInflation,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertNumeric(
    omegaInflation,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertNumeric(
    sigmaInflation,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)
  checkmate::assertFlag(recenter)
  checkmate::assertFlag(boxcox)
  checkmate::assertFlag(recover)
  checkmate::assertFlag(addIterations)
  checkmate::assertNumber(sigmaFallbackRse, lower = 0, finite = TRUE)
  if (!is.null(omegaDf)) {
    checkmate::assertNumber(omegaDf, lower = 1, finite = TRUE)
  }
  nlmixr2utils::.validateWorkers(workers)

  structure(
    list(
      thetaInflation = thetaInflation,
      omegaInflation = omegaInflation,
      sigmaInflation = sigmaInflation,
      capCorrelation = capCorrelation,
      capResampling = capResampling,
      recenter = recenter,
      boxcox = boxcox,
      workers = workers,
      rxThreads = rxThreads,
      recover = recover,
      addIterations = addIterations,
      omegaFallback = omegaFallback,
      sigmaFallbackRse = sigmaFallbackRse,
      omegaDf = omegaDf
    ),
    class = "runSIRControl"
  )
}

#' @export
print.runSIRControl <- function(x, ...) {
  cli::cli_h2("nlmixr2sir control")
  cli::cli_dl(c(
    inflation = "theta {x$thetaInflation}, omega {x$omegaInflation}, sigma {x$sigmaInflation}",
    caps = "correlation {x$capCorrelation}, resampling {x$capResampling}",
    proposal = "recenter {x$recenter}, boxcox {x$boxcox}, omegaFallback {.val {x$omegaFallback}}",
    parallel = "workers {format(x$workers %||% 'current plan')}, rxThreads {format(x$rxThreads %||% 'rxode2 default')}",
    resume = "recover {x$recover}, addIterations {x$addIterations}"
  ))
  invisible(x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
