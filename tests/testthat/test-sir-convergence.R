# The dOFV-versus-chi-square convergence diagnostic, and PsN's summary-
# statistics parity. Ported from PsN R-scripts/sir_default.R and
# lib/tool/sir.pm empirical_statistics().

test_that(".sirDofvQuantiles stops short of 1 so the reference stays finite", {
  q <- .sirDofvQuantiles(20L)
  expect_length(q, 19L)
  expect_equal(min(q), 0)
  expect_equal(max(q), 19 / 20)
  expect_true(all(is.finite(stats::qchisq(q, df = 5))))
})

test_that(".sirDofvQuantiles refuses a grid too small to be a curve", {
  expect_error(.sirDofvQuantiles(2L), "At least 3 resamples")
})

test_that(".sirDofvCurves gives proposal and SIR curves per iteration", {
  skip_on_cran()
  cv <- .sirDofvCurves(sirObj())
  expect_setequal(cv$type, c("reference", "proposal", "SIR"))
  # one reference curve, shared across iterations
  expect_equal(sum(cv$type == "reference"), sum(cv$type == "SIR"))
  expect_true(all(is.na(cv$iteration[cv$type == "reference"])))
  expect_false(anyNA(cv$iteration[cv$type != "reference"]))
})

test_that(".sirDofvCurves reference is a chi-square on the parameter count", {
  skip_on_cran()
  cv <- .sirDofvCurves(sirObj())
  ref <- cv[cv$type == "reference", ]
  expect_equal(
    ref$dOFV,
    stats::qchisq(ref$quantile, df = length(unique(sirObj()$param))),
    tolerance = 1e-10
  )
  expect_false(is.unsorted(ref$dOFV))
})

test_that(".sirProposalTooNarrow applies PsN's 25% rule", {
  # Hand-built curves: the proposal sits below the reference at 3 of 4
  # quantiles, which is above the quarter threshold.
  curves <- rbind(
    data.frame(
      iteration = NA_integer_,
      label = "reference",
      type = "reference",
      quantile = c(0.1, 0.2, 0.3, 0.4),
      dOFV = c(1, 2, 3, 4),
      stringsAsFactors = FALSE
    ),
    data.frame(
      iteration = 1L,
      label = "proposal 1",
      type = "proposal",
      quantile = c(0.1, 0.2, 0.3, 0.4),
      dOFV = c(0, 0, 0, 5),
      stringsAsFactors = FALSE
    )
  )
  chk <- .sirProposalTooNarrow(curves)
  expect_equal(chk$fraction, 0.75)
  expect_true(chk$warn)
})

test_that(".sirProposalTooNarrow stays quiet for a wide enough proposal", {
  curves <- rbind(
    data.frame(
      iteration = NA_integer_,
      label = "reference",
      type = "reference",
      quantile = c(0.1, 0.2, 0.3, 0.4),
      dOFV = c(1, 2, 3, 4),
      stringsAsFactors = FALSE
    ),
    data.frame(
      iteration = 1L,
      label = "proposal 1",
      type = "proposal",
      quantile = c(0.1, 0.2, 0.3, 0.4),
      dOFV = c(2, 3, 4, 5),
      stringsAsFactors = FALSE
    )
  )
  chk <- .sirProposalTooNarrow(curves)
  expect_equal(chk$fraction, 0)
  expect_false(chk$warn)
})

test_that("the too-narrow warning names inflation as the remedy", {
  curves <- rbind(
    data.frame(
      iteration = NA_integer_,
      label = "reference",
      type = "reference",
      quantile = c(0.1, 0.2),
      dOFV = c(1, 2),
      stringsAsFactors = FALSE
    ),
    data.frame(
      iteration = 1L,
      label = "proposal 1",
      type = "proposal",
      quantile = c(0.1, 0.2),
      dOFV = c(0, 0),
      stringsAsFactors = FALSE
    )
  )
  expect_warning(
    .sirWarnProposalTooNarrow(.sirProposalTooNarrow(curves)),
    "inflated proposal"
  )
})

test_that(".sirDofvNoise brackets the curve it is built from", {
  skip_on_cran()
  cv <- .sirDofvCurves(sirObj())
  quant <- cv$quantile[cv$type == "reference"]
  set.seed(1)
  nb <- .sirDofvNoise(sirObj(), 1L, quant, nReplicate = 50L)
  expect_equal(nrow(nb), length(quant))
  expect_true(all(nb$low <= nb$high))
})

test_that("plot(type = 'convergence') builds a faceted ggplot", {
  skip_on_cran()
  p <- suppressWarnings(
    plot(sirObj(), type = "convergence", nReplicate = 25L)
  )
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "dOFV")
  expect_equal(p$labels$x, "Quantile")
})

test_that("plot(type = 'convergence') can omit the noise band", {
  skip_on_cran()
  p <- suppressWarnings(plot(sirObj(), type = "convergence", noise = FALSE))
  expect_s3_class(p, "ggplot")
})

# Summary-statistics parity.

test_that(".sirPercentileLabels reproduces PsN's percentile set", {
  # From prediction intervals 0, 40, 80, 90, 95.
  expect_equal(
    .sirPercentileLabels(),
    c(2.5, 5, 10, 30, 50, 70, 90, 95, 97.5)
  )
})

test_that("sirSummary reports both mean and median", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  expect_true(all(c("mean", "p50") %in% names(s)))
  expect_equal(
    s$mean,
    unname(colMeans(iter1()$resampledMat)),
    tolerance = 1e-12
  )
})

test_that("rse_sd_scale halves the RSE of variance parameters only", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  ps <- .sirParamSpace(theoFit())
  kind <- ps$kind[match(s$param, ps$sirName)]

  isOmega <- kind %in% c("omegaDiag", "omegaOffdiag")
  expect_equal(s$rse_sd_scale[isOmega], s$rse[isOmega] / 2, tolerance = 1e-12)
  # THETA is not a variance, and nlmixr2's residual error is already on the SD
  # scale, so neither is rescaled.
  expect_true(all(is.na(s$rse_sd_scale[!isOmega])))
})

test_that("sirSummary attaches an sd/correlation matrix", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  sdcor <- attr(s, "sdCorMatrix")
  cm <- stats::cov(iter1()$resampledMat)
  expect_equal(diag(sdcor), sqrt(diag(cm)), tolerance = 1e-12)
  expect_equal(
    sdcor[lower.tri(sdcor)],
    stats::cov2cor(cm)[lower.tri(cm)],
    tolerance = 1e-12
  )
})

test_that("sirSummary records that rse is a percentage", {
  skip_on_cran()
  # PsN reports the same quantity as a fraction; the units are recorded so the
  # difference is not silent.
  expect_equal(
    attr(sirSummary(iter1()$resampledMat, theoFit()), "rseUnits"),
    "percent"
  )
})

# P2.2 intervals by iteration, P2.3 RSE/correlation, P2.5 on-disk parity.

test_that(".sirIterationIntervals covers both distributions per iteration", {
  skip_on_cran()
  iv <- .sirIterationIntervals(sirObj())
  expect_setequal(iv$type, c("proposal", "SIR"))
  expect_setequal(iv$param, colnames(attr(sirObj(), "resampledMat")))
  expect_true(all(iv$low <= iv$median))
  expect_true(all(iv$median <= iv$high))
})

test_that(".sirIterationIntervals honours the requested interval width", {
  skip_on_cran()
  narrow <- .sirIterationIntervals(sirObj(), ci = 50)
  wide <- .sirIterationIntervals(sirObj(), ci = 95)
  expect_true(all(
    (wide$high - wide$low) >= (narrow$high - narrow$low) - 1e-12
  ))
})

test_that("asymmetry is the ratio of the two half-widths", {
  skip_on_cran()
  iv <- .sirIterationIntervals(sirObj())
  expect_equal(
    iv$asymmetry,
    (iv$high - iv$median) / (iv$median - iv$low),
    tolerance = 1e-12
  )
})

test_that(".sirRseCorData keeps one triangle with RSE on the diagonal", {
  skip_on_cran()
  d <- .sirRseCorData(sirObj())
  n <- length(unique(sirObj()$param))
  expect_equal(nrow(d), n * (n + 1) / 2)
  expect_equal(sum(d$isDiagonal), n)
  # off-diagonal cells are correlations
  expect_true(all(abs(d$value[!d$isDiagonal]) <= 1 + 1e-12))
  # diagonal cells carry an asymmetry band, off-diagonal ones do not
  expect_false(anyNA(d$asymmetryBand[d$isDiagonal]))
  expect_true(all(is.na(d$asymmetryBand[!d$isDiagonal])))
})

test_that(".sirRseCorData can show the proposal instead of the posterior", {
  skip_on_cran()
  expect_equal(
    unique(.sirRseCorData(sirObj(), which = "proposal")$which),
    "proposal"
  )
  expect_equal(unique(.sirRseCorData(sirObj(), which = "SIR")$which), "SIR")
})

test_that("asymmetry bands follow PsN's breakpoints", {
  skip_on_cran()
  d <- .sirRseCorData(sirObj())
  diagRows <- d[d$isDiagonal, ]
  expected <- cut(
    diagRows$asymmetry,
    breaks = c(-Inf, 0.5, 1, 1.25, 2, Inf),
    labels = c("<0.5", "0.5-1", "1-1.25", "1.25-2", ">2"),
    right = FALSE
  )
  expect_equal(diagRows$asymmetryBand, expected)
})

test_that("plot(type = 'intervals') and plot(type = 'rsecor') build", {
  skip_on_cran()
  expect_s3_class(plot(sirObj(), type = "intervals"), "ggplot")
  expect_s3_class(plot(sirObj(), type = "rsecor"), "ggplot")
  expect_s3_class(
    plot(sirObj(), type = "rsecor", which = "proposal"),
    "ggplot"
  )
})

test_that("summary_iterations.csv leads with PsN's column names", {
  skip_on_cran()
  tmp <- tempfile("sir_cols_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  .sirWriteIterationSummary(attr(sirObj(), "iterationSummary"), tmp)
  got <- names(utils::read.csv(
    file.path(tmp, "summary_iterations.csv"),
    check.names = FALSE
  ))
  expect_identical(
    got[seq_len(10)],
    c(
      "iteration",
      "commandline.samples",
      "attempted.samples",
      "successful.samples",
      "commandline.resamples",
      "actual.resamples",
      "requested.ratio",
      "actual.ratio",
      "negative.dOFV",
      "minimum.sample.ofv"
    )
  )
})

test_that(".sirWriteCovMatrices exports the covariance and sd/correlation", {
  skip_on_cran()
  tmp <- tempfile("sir_cov_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  s <- sirSummary(iter1()$resampledMat, theoFit())
  .sirWriteCovMatrices(s, tmp, fitName = "demo")

  expect_true(file.exists(file.path(tmp, "demo_sir.cov")))
  expect_true(file.exists(file.path(tmp, "demo_sir.sdcorr")))
  back <- utils::read.delim(file.path(tmp, "demo_sir.cov"), check.names = FALSE)
  expect_equal(back$NAME, rownames(attr(s, "covMatrix")))
  expect_equal(
    as.matrix(back[, -1L]),
    unname(attr(s, "covMatrix")),
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})
