# Box-Cox domain and non-finite handling.
#
# Two independent robustness problems, both in the review as R5:
#   * inverse transformation can overflow to Inf without raising an R error,
#     and is.na(Inf) is FALSE; and
#   * the shift is chosen from the retained sample alone, so a recentred mean
#     below that sample's minimum can fall outside the transform domain.

test_that("inverse Box-Cox overflows to Inf rather than NA or an error", {
  # The premise of the next test, asserted separately so a change in overflow
  # behaviour is reported here rather than making the next test vacuous.
  # lambda = 0 inverts as exp(v) - delta, so a large transformed value
  # overflows. (A small non-zero lambda raises to a high power and can stay
  # finite: (0.05 * 1e5 + 1)^20 is only about 1e74.)
  v <- sirBoxCoxInverse(1000, lambda = 0, delta = 1e-6)
  expect_false(is.finite(v))
  expect_false(is.na(v))
})

test_that("an infinite inverse result is counted as an inverse failure", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)

  # One unbounded theta, so nothing can be rejected for being out of bounds and
  # there is no OMEGA to reject on: any rejection must be the inverse
  # transform. This matters because Inf > Inf is FALSE, so an infinite value
  # passes a bounds test against an infinite upper bound.
  nm <- ps$sirName[ps$kind == "theta"][[1L]]
  ps_one <- ps[ps$sirName == nm, , drop = FALSE]
  bc <- data.frame(param = nm, lambda = 0, delta = 1e-6)

  mu <- stats::setNames(1e5, nm)
  lower <- stats::setNames(-Inf, nm)
  upper <- stats::setNames(Inf, nm)

  out <- suppressWarnings(.sirSampleFullProposal(
    mu = mu,
    covMat = matrix(1, 1L, 1L, dimnames = list(nm, nm)),
    n = 4L,
    lower = lower,
    upper = upper,
    ps = ps_one,
    baseOmega = matrix(numeric(0), 0L, 0L),
    thetaNames = nm,
    sigmaNames = character(0),
    boxcoxState = bc
  ))

  expect_gt(out$inverseRejected, 0L)
  expect_equal(out$thetaRejected, 0L)
  expect_equal(out$omegaRejected, 0L)
  expect_equal(out$sigmaRejected, 0L)
  # Nothing infinite may reach the returned sample either.
  expect_true(all(is.finite(out$samples)))
})

test_that("the shift keeps a recentred mu inside the Box-Cox domain", {
  # The shift is chosen from the retained sample. With recenter = TRUE the next
  # centre is the best candidate, which need not be among the randomly retained
  # rows -- it can sit below their minimum. A delta chosen from the sample
  # alone then gives mu + delta <= 0, and .sirBcTransformMu() aborts on a run
  # that was proceeding normally.
  retained <- matrix(
    c(1, 2, 3, 4),
    ncol = 1L,
    dimnames = list(NULL, "tka")
  )
  centre <- c(tka = -5)

  updated <- sirUpdateProposal(retained, centre = centre)
  expect_silent(tr <- .sirBcTransformMu(centre, updated$boxcoxParams))
  expect_true(all(is.finite(tr)))
})
