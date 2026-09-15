# B2: the candidate objective must be the fit's objective ---------------------
#
# sirEvalOFV() scores every candidate with a freshly built FOCEi call at
# maxOuterIterations = 0. That is the workflow analogue of NONMEM MAXEVAL=0,
# but it is not by itself proof that the resulting objective is the same
# surface that produced fit$objf. Importance sampling assumes one fixed target:
# a candidate-dependent difference between evaluators changes the shape of the
# target and so the retained distribution.
#
# The preflight settles it empirically, per run, by re-evaluating the fitted
# centre and comparing with the stored objective before any sampling happens.

test_that("the preflight accepts a fit whose centre reproduces its objective", {
  skip_on_cran()
  fit <- theoFit()
  expect_no_error(.sirCheckObjective(fit, workers = 1L))
})

test_that("the preflight returns the stored and reevaluated objectives", {
  skip_on_cran()
  fit <- theoFit()
  res <- .sirCheckObjective(fit, workers = 1L)
  expect_named(res, c("stored", "reevaluated", "absDiff", "relDiff", "stencil"))
  expect_equal(res$stored, fit$objf, tolerance = 1e-12)
  expect_lt(res$absDiff, 1e-3)
})

test_that("the preflight aborts when the centre does not reproduce the objective", {
  skip_on_cran()
  fit <- theoFit()
  # A tolerance tight enough that even the genuine numerical difference between
  # the stored and reevaluated objective fails it. This is the mismatch path:
  # the message must name both values so the user can judge the gap.
  expect_error(
    .sirCheckObjective(fit, workers = 1L, objfTolerance = 0),
    "objective"
  )
  err <- tryCatch(
    .sirCheckObjective(fit, workers = 1L, objfTolerance = 0),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, format(fit$objf, digits = 10), fixed = TRUE)
})

test_that("unsupported estimation methods are rejected before sampling", {
  skip_on_cran()
  # The check is on the recorded method, so it fires without needing a real fit
  # of that kind. saem is the case that matters: its objective comes from
  # Gaussian quadrature, and reevaluating it under FOCEi shifts the OFV by
  # ~2.7 units on theo_sd -- a different surface, not numerical noise.
  expect_error(.sirSupportedEstimation("saem"), "saem")
  expect_error(.sirSupportedEstimation("nlme"), "nlme")
  expect_error(.sirSupportedEstimation(NA_character_), "determine")
  expect_silent(.sirSupportedEstimation("focei"))
})

test_that("runSIR runs the objective preflight before sampling", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()
  expect_error(
    suppressMessages(runSIR(
      fit,
      nSamples = 16L,
      nResample = 8L,
      directory = dir,
      control = runSIRControl(
        recover = FALSE,
        workers = 1L,
        objfTolerance = 0
      )
    )),
    "objective"
  )
  # It aborted before doing any work, so there are no iteration artifacts.
  expect_false(file.exists(file.path(dir, "raw_results.csv")))
})

# B2: the stencil, and carrying the fit's likelihood settings ----------------

test_that("the evaluator carries the fit's likelihood-relevant controls", {
  skip_on_cran()
  fit <- theoFit()
  ec <- .sirEvalControl(fit)
  expect_s3_class(ec, "foceiControl")
  # Evaluation-only overrides are always ours.
  expect_equal(ec$maxOuterIterations, 0L)
  # foceiControl() normalises covMethod = "" to integer 0, meaning no step.
  expect_equal(as.integer(ec$covMethod), 0L)
  # Likelihood-defining settings come from the fit, not from the defaults.
  expect_identical(ec$interaction, fit$control$interaction)
  expect_identical(ec$addProp, fit$control$addProp)
})

test_that("the preflight probes the surface around the centre", {
  skip_on_cran()
  fit <- theoFit()
  res <- .sirCheckObjective(fit, workers = 1L)
  expect_false(is.null(res$stencil))
  # Two probes per parameter, less any clipped by a bound.
  expect_lte(res$stencil$nProbes, 2L * nrow(.sirParamSpace(fit)))
  expect_gt(res$stencil$nProbes, 0L)
  expect_equal(res$stencil$nFailed, 0L)
  # The fitted estimates are a local optimum, so no probe improves much on it.
  expect_gt(res$stencil$minDOFV, -1)
})

test_that("the stencil can be switched off", {
  skip_on_cran()
  res <- .sirCheckObjective(theoFit(), workers = 1L, stencil = FALSE)
  expect_null(res$stencil)
})

test_that("the preflight tolerance is absolute, not relative", {
  skip_on_cran()
  fit <- theoFit()
  # A relative rule would wave through a large absolute gap on a large
  # objective. The weights depend on differences in OFV, so only the absolute
  # scale is meaningful.
  expect_error(
    .sirCheckObjective(fit, workers = 1L, objfTolerance = 0, stencil = FALSE),
    "absolute"
  )
})
