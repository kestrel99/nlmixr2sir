# Numeric parity with PsN, checked against oracle values taken from PsN's own
# unit tests. These exercise the pure arithmetic, so they need no model fit and
# no NONMEM -- which is what makes them worth porting.
#
# Source: PsN test/unit/tool/sir.t, lines 120-174.

test_that(".sirRound rounds half away from zero, as PsN's math.pm does", {
  expect_equal(.sirRound(20.5), 21L)
  expect_equal(.sirRound(21.5), 22L)
  expect_equal(.sirRound(20.4), 20L)
  expect_equal(.sirRound(20.6), 21L)
  expect_equal(.sirRound(-20.5), -21L)
  expect_equal(.sirRound(-20.4), -20L)
  # R's own round() is banker's rounding and disagrees on exact halves; this
  # is the reason .sirRound() exists.
  expect_equal(round(20.5), 20)
})

# PsN oracle: samples = 100 and resamples = 20 every iteration, with
# successful counts 90, 102, 109, 98. Attempted must come out 100, 111, 109,
# 100 and actual resamples 18, 20, 22, 20.
test_that(".sirAdjustedAttemptedSamples matches the PsN oracle sequence", {
  attempted <- integer(4)
  successful <- c(90L, 102L, 109L, 98L)

  attempted[1] <- .sirAdjustedAttemptedSamples(100L)
  for (i in 2:4) {
    attempted[i] <- .sirAdjustedAttemptedSamples(
      100L,
      previousAttempted = attempted[i - 1L],
      previousSuccessful = successful[i - 1L]
    )
  }
  expect_equal(attempted, c(100L, 111L, 109L, 100L))
})

test_that(".sirAdjustedResamples matches the PsN oracle sequence", {
  successful <- c(90L, 102L, 109L, 98L)
  actual <- vapply(
    successful,
    function(s) .sirAdjustedResamples(20L, 100L, s),
    integer(1L)
  )
  expect_equal(actual, c(18L, 20L, 22L, 20L))
})

test_that(".sirAdjustedAttemptedSamples compensates only for loss", {
  # 10% loss: 100 / 0.9 = 111.1 -> 111
  expect_equal(.sirAdjustedAttemptedSamples(100L, 100L, 90L), 111L)
  # No loss: unchanged
  expect_equal(.sirAdjustedAttemptedSamples(100L, 100L, 100L), 100L)
  # A gain in the previous iteration is impossible, and is a bug if seen
  expect_error(
    .sirAdjustedAttemptedSamples(100L, 100L, 101L),
    "More successful samples than attempted"
  )
})

test_that(".sirAdjustedAttemptedSamples triggers at exactly 0.95 turnout", {
  # PsN's condition is `previous_turnout <= 0.95`, inclusive. At exactly 0.95
  # the count is compensated; just above it, it is not.
  expect_equal(.sirAdjustedAttemptedSamples(100L, 100L, 95L), 105L)
  expect_equal(.sirAdjustedAttemptedSamples(100L, 100L, 96L), 100L)
})

test_that(".sirAdjustedResamples scales on gain as well as loss", {
  expect_equal(.sirAdjustedResamples(20L, 100L, 90L), 18L) # 10% loss
  expect_equal(.sirAdjustedResamples(20L, 100L, 110L), 22L) # 10% gain
  expect_equal(.sirAdjustedResamples(20L, 100L, 98L), 20L) # 2% loss, inert
  expect_equal(.sirAdjustedResamples(20L, 100L, 102L), 20L) # 2% gain, inert
})

test_that(".sirAdjustedResamples triggers at exactly 5% deviation", {
  # PsN's condition is `abs(turnout - 1) >= 0.05`, inclusive at the boundary.
  expect_equal(.sirAdjustedResamples(20L, 100L, 95L), 19L)
  expect_equal(.sirAdjustedResamples(20L, 100L, 105L), 21L)
  expect_equal(.sirAdjustedResamples(20L, 100L, 96L), 20L)
})

test_that(".sirAdjustedResamples measures turnout against requested samples", {
  # Not against the compensated attempted count. PsN's third oracle iteration
  # pins this: 109 successful of a requested 100 is a 9% gain, giving 22, even
  # though 109 samples were attempted that iteration.
  expect_equal(.sirAdjustedResamples(20L, 100L, 109L), 22L)
})

# Inflation. PsN setup_inflation() takes one value per class, or one per
# *diagonal* element of that class; an OMEGA off-diagonal is never given a
# factor directly but derives sqrt(infl_i) * sqrt(infl_j) from its diagonals.

test_that(".sirInflationVector recycles a scalar across the whole class", {
  skip_on_cran()
  v <- .sirInflationVector(.sirParamSpace(blockFit()), omegaInflation = 4)
  expect_equal(unname(v[c("eta.ka", "eta.cl:eta.ka", "eta.cl")]), c(4, 4, 4))
  expect_equal(unname(v[c("tka", "tcl", "tv", "add.sd")]), rep(1, 4))
})

test_that(".sirInflationVector derives off-diagonal inflation from diagonals", {
  skip_on_cran()
  v <- .sirInflationVector(
    .sirParamSpace(blockFit()),
    omegaInflation = c(4, 9)
  )
  expect_equal(v[["eta.ka"]], 4)
  expect_equal(v[["eta.cl"]], 9)
  expect_equal(v[["eta.cl:eta.ka"]], sqrt(4) * sqrt(9))
})

test_that(".sirInflationVector accepts one value per THETA", {
  skip_on_cran()
  v <- .sirInflationVector(
    .sirParamSpace(blockFit()),
    thetaInflation = c(1, 2, 3)
  )
  expect_equal(unname(v[c("tka", "tcl", "tv")]), c(1, 2, 3))
  expect_equal(v[["add.sd"]], 1)
})

test_that(".sirInflationVector rejects a length that is neither 1 nor n", {
  skip_on_cran()
  ps <- .sirParamSpace(blockFit())
  expect_error(
    .sirInflationVector(ps, omegaInflation = c(1, 2, 3)),
    "one value per"
  )
  expect_error(
    .sirInflationVector(ps, thetaInflation = c(1, 2)),
    "one value per"
  )
})

test_that("equal inflation factors leave the proposal correlation unchanged", {
  skip_on_cran()
  inflated <- sirGetProposalCov(
    blockFit(),
    omegaInflation = c(4, 4),
    capCorrelation = 1
  )
  base <- sirGetProposalCov(blockFit(), capCorrelation = 1)
  expect_equal(
    unname(cov2cor(inflated)),
    unname(cov2cor(base)),
    tolerance = 1e-12
  )
})

test_that("vector inflation reaches the assembled proposal covariance", {
  skip_on_cran()
  fit <- blockFit()
  mu <- .sirProposalMu(fit)
  pc <- sirGetProposalCov(fit)
  base <- .sirInitialProposal(fit, mu, pc, capCorrelation = 1)
  infl <- .sirInitialProposal(
    fit,
    mu,
    pc,
    omegaInflation = c(4, 9),
    capCorrelation = 1
  )
  expect_equal(
    infl$covMat["eta.ka", "eta.ka"],
    4 * base$covMat["eta.ka", "eta.ka"],
    tolerance = 1e-10
  )
  expect_equal(
    infl$covMat["eta.cl", "eta.cl"],
    9 * base$covMat["eta.cl", "eta.cl"],
    tolerance = 1e-10
  )
  expect_equal(
    infl$covMat["eta.cl:eta.ka", "eta.cl:eta.ka"],
    6 * base$covMat["eta.cl:eta.ka", "eta.cl:eta.ka"],
    tolerance = 1e-10
  )
})
