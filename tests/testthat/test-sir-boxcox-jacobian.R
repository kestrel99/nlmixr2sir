# B3: the Box-Cox change-of-variables Jacobian --------------------------------
#
# With Box-Cox enabled, candidates are drawn as y = T(x) from a multivariate
# normal q_y, then mapped back to x. The density actually induced on the
# original parameter scale is
#
#   q_x(x) = q_y(T(x)) * |det J_T(x)|,
#
# so the importance weight for an original-scale likelihood target is
#
#   w(x) = L(x) / [ q_y(T(x)) * |det J_T(x)| ].
#
# Dropping the Jacobian (which is what PsN does) retains a sample from
# L(x)|det J_T(x)| instead: a distribution that changes when the same model is
# written in a different smooth parameterization. nlmixr2sir includes it, so
# the retained sample targets the original-scale normalized likelihood.
#
# Box-Cox is applied one coordinate at a time, so J_T is diagonal and
#
#   log|det J_T(x)| = sum_j (lambda_j - 1) * log(x_j + delta_j).

test_that(".sirBcLogJacobian matches the analytic Box-Cox derivative", {
  bc <- data.frame(
    param = c("a", "b"),
    lambda = c(0.3, 0),
    delta = c(2, 1),
    stringsAsFactors = FALSE
  )
  x <- matrix(
    c(1.0, 4.0, 2.5, 0.5),
    ncol = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b"))
  )

  # d/dx of (x+delta)^lambda/lambda is (x+delta)^(lambda-1); for lambda = 0 the
  # transform is log(x+delta), whose derivative is (x+delta)^(-1). Same formula.
  expected <- (0.3 - 1) * log(x[, "a"] + 2) + (0 - 1) * log(x[, "b"] + 1)

  expect_equal(.sirBcLogJacobian(x, bc), expected, tolerance = 1e-12)
})

test_that(".sirBcLogJacobian is zero for an identity transform and for NULL", {
  x <- matrix(c(1.5, 2.5), ncol = 1L, dimnames = list(NULL, "a"))
  # lambda = 1 is a pure shift: the Jacobian is 1, so its log is 0.
  bc <- data.frame(param = "a", lambda = 1, delta = 0.5, stringsAsFactors = FALSE)
  expect_equal(.sirBcLogJacobian(x, bc), c(0, 0), tolerance = 1e-12)
  expect_equal(.sirBcLogJacobian(x, NULL), c(0, 0), tolerance = 1e-12)
})

test_that("sirCalcWeights folds logJacobian into the relative density", {
  mu <- c(a = 0)
  cov1 <- matrix(1, 1L, 1L)
  x <- matrix(c(0, 0.5, -1.25), ncol = 1L, dimnames = list(NULL, "a"))
  lj <- c(0, 0.75, -0.4)

  plain <- sirCalcWeights(x, mu, cov1, dOFV = rep(0, 3L))
  withJ <- sirCalcWeights(x, mu, cov1, dOFV = rep(0, 3L), logJacobian = lj)

  expect_equal(log(withJ$relPDF), log(plain$relPDF) + lj, tolerance = 1e-12)
  # The centre is still the reference point: relPDF == 1 there.
  expect_equal(withJ$relPDF[1L], 1, tolerance = 1e-12)
})

test_that("a proposal equal to the target gives uniform weights under Box-Cox", {
  # The decisive check. Build a 1-D Box-Cox proposal, then set every dOFV so
  # that the likelihood is exactly proportional to the induced original-scale
  # density q_x. Every importance ratio is then identical, so every resampling
  # probability must be 1/n. This holds only if the Jacobian is included --
  # without it the weights tilt by exactly the Jacobian ratio.
  lambda <- 0.3
  delta <- 2
  bc <- data.frame(
    param = "a", lambda = lambda, delta = delta, stringsAsFactors = FALSE
  )
  tf <- function(x) ((x + delta)^lambda - 1) / lambda

  muX <- 1
  muY <- tf(muX)
  sdY <- 0.4

  xs <- c(1, 0.2, 2.75, 4.5, 0.6)
  ys <- tf(xs)

  xMat <- matrix(xs, ncol = 1L, dimnames = list(NULL, "a"))
  yMat <- matrix(ys, ncol = 1L, dimnames = list(NULL, "a"))

  logJ <- .sirBcLogJacobian(xMat, bc) -
    .sirBcLogJacobian(matrix(muX, ncol = 1L, dimnames = list(NULL, "a")), bc)

  # log q_x(x) - log q_x(muX), the induced original-scale log relative density.
  logQx <- stats::dnorm(ys, muY, sdY, log = TRUE) -
    stats::dnorm(muY, muY, sdY, log = TRUE) +
    logJ

  # Make the target equal the proposal: -0.5 * dOFV == log relative density.
  dofv <- -2 * logQx

  res <- sirCalcWeights(
    yMat,
    c(a = muY),
    matrix(sdY^2, 1L, 1L),
    dOFV = dofv,
    logJacobian = logJ
  )

  expect_equal(res$prob_resample, rep(1 / length(xs), length(xs)), tolerance = 1e-10)

  # And confirm the test is not vacuous: without the Jacobian the same inputs
  # give visibly non-uniform weights.
  bad <- sirCalcWeights(yMat, c(a = muY), matrix(sdY^2, 1L, 1L), dOFV = dofv)
  expect_gt(max(abs(bad$prob_resample - 1 / length(xs))), 1e-3)
})

test_that("sirRunIteration applies the Jacobian when Box-Cox is active", {
  skip_on_cran()
  fit <- theoFit()
  mu <- .sirProposalMu(fit)
  bc <- data.frame(
    param = names(mu),
    lambda = 0.5,
    delta = 5,
    stringsAsFactors = FALSE
  )
  set.seed(99)
  it <- suppressMessages(suppressWarnings(sirRunIteration(
    fit,
    mu = mu,
    proposalCov = sirGetProposalCov(fit),
    nSamples = 16L,
    nResample = 8L,
    iterNum = 2L,
    recenter = FALSE,
    boxcox = TRUE,
    directory = NULL,
    boxcoxState = bc
  )))

  raw <- it$rawResults[it$rawResults$role == "sample", , drop = FALSE]
  raw <- raw[!duplicated(raw$sample_id), , drop = FALSE]
  params <- colnames(it$resampledMat)
  xMat <- as.matrix(raw[, params, drop = FALSE])
  bcp <- bc[match(params, bc$param), ]

  # Rebuild both pieces of the density from the stored original-scale draws.
  tf <- function(v, lam, del) {
    xs <- v + del
    if (abs(lam) < 1e-10) log(xs) else (xs^lam - 1) / lam
  }
  yMat <- xMat
  for (j in seq_along(params)) {
    yMat[, j] <- tf(xMat[, j], bcp$lambda[j], bcp$delta[j])
  }
  muY <- vapply(
    seq_along(params),
    function(j) tf(mu[[params[j]]], bcp$lambda[j], bcp$delta[j]),
    numeric(1L)
  )

  covMat <- sirGetProposalCov(fit)[params, params, drop = FALSE]
  inv <- solve(covMat)
  logMvn <- apply(yMat, 1L, function(z) {
    d <- z - muY
    -0.5 * as.numeric(t(d) %*% inv %*% d)
  })

  logJ <- .sirBcLogJacobian(xMat, bc) -
    .sirBcLogJacobian(
      matrix(mu[params], nrow = 1L, dimnames = list(NULL, params)),
      bc
    )
  # The Jacobian must actually bite, or this test proves nothing.
  expect_gt(max(abs(logJ)), 1e-6)

  # relPDF is the ORIGINAL-scale relative density: normal part plus Jacobian.
  expect_equal(log(raw$relPDF), unname(logMvn + logJ), tolerance = 1e-8)

  # And it is not merely the normal part, which is what the old code reported.
  expect_gt(max(abs(log(raw$relPDF) - logMvn)), 1e-6)
})

test_that("the retained distribution is invariant to the Box-Cox reparameterization", {
  # The payoff test, and the justification for diverging from PsN here.
  #
  # Target a known skewed likelihood -- Gamma(3, 1), mean 3, second moment 12 --
  # and importance-sample it two ways: a normal proposal on the original scale,
  # and a normal proposal on a Box-Cox-transformed scale. Both are proposals for
  # the SAME target, so both weighted estimates must recover the same moments.
  # They do only when the Jacobian is included.
  skip_on_cran()
  set.seed(4242)
  n <- 40000L
  shape <- 3
  rate <- 1
  lambda <- 0.25
  delta <- 1e-6
  tf <- function(x) ((x + delta)^lambda - 1) / lambda
  itf <- function(y) (lambda * y + 1)^(1 / lambda) - delta

  nm <- list(NULL, "a")

  # (a) normal proposal on the original scale; no transform, no Jacobian.
  muX <- 3
  sdX <- 2.5
  x <- stats::rnorm(n, muX, sdX)
  dofvA <- ifelse(x > 0, -2 * stats::dgamma(x, shape, rate, log = TRUE), Inf)
  wA <- sirCalcWeights(
    matrix(x, ncol = 1L, dimnames = nm),
    c(a = muX),
    matrix(sdX^2, 1L, 1L),
    dOFV = dofvA
  )

  # (b) normal proposal on the Box-Cox scale.
  muY <- tf(3)
  sdY <- 0.9
  y <- stats::rnorm(n, muY, sdY)
  xb <- suppressWarnings(itf(y))
  ok <- is.finite(xb) & xb > 0
  dofvB <- ifelse(
    ok,
    -2 * stats::dgamma(pmax(xb, 1e-300), shape, rate, log = TRUE),
    Inf
  )
  bc <- data.frame(
    param = "a", lambda = lambda, delta = delta, stringsAsFactors = FALSE
  )
  logJ <- .sirBcLogJacobian(
    matrix(ifelse(ok, xb, 1), ncol = 1L, dimnames = nm),
    bc
  ) - .sirBcLogJacobian(matrix(3, ncol = 1L, dimnames = nm), bc)

  yMat <- matrix(y, ncol = 1L, dimnames = nm)
  wB <- sirCalcWeights(
    yMat, c(a = muY), matrix(sdY^2, 1L, 1L), dOFV = dofvB, logJacobian = logJ
  )
  wNoJ <- sirCalcWeights(
    yMat, c(a = muY), matrix(sdY^2, 1L, 1L), dOFV = dofvB
  )

  xv <- ifelse(ok, xb, 0)
  meanA <- sum(wA$prob_resample * x)
  meanB <- sum(wB$prob_resample * xv)
  meanNoJ <- sum(wNoJ$prob_resample * xv)
  m2B <- sum(wB$prob_resample * xv^2)

  # Both parameterizations recover the target. Tolerances are several Monte
  # Carlo standard errors wide at this sample size, with the seed fixed.
  expect_equal(meanA, shape / rate, tolerance = 0.05)
  expect_equal(meanB, shape / rate, tolerance = 0.05)
  expect_equal(m2B, shape * (shape + 1) / rate^2, tolerance = 0.2)

  # Without the Jacobian the same draws give a badly biased answer: ~2.25
  # against a true mean of 3. This is the defect, not a rounding difference.
  expect_lt(meanNoJ, 2.5)
})
