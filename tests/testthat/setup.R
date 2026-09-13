# Shared fixtures for the nlmixr2sir test suite.
#
# Fixtures are lazy and memoised: each is built on first use and the outcome --
# success or failure -- is cached. Lazy means nothing expensive is built on
# CRAN, because every test calls skip_on_cran() before touching a fixture.
# Memoised failure means a broken fixture is reported by the tests that need
# it, once each, instead of halting the rest of the file.

.sirLazy <- function(expr) {
  expr <- substitute(expr)
  env <- parent.frame()
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- tryCatch(
        list(ok = TRUE, value = eval(expr, env)),
        error = function(e) list(ok = FALSE, value = e)
      )
    }
    if (!cached$ok) {
      stop(cached$value)
    }
    cached$value
  }
}

# One-compartment model on theo_sd, single eta, with a covariance step.
theoOneCmt <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
    eta.ka ~ 0.6
    add.sd <- 0.7 # nolint: object_usage_linter.
  })
  model({
    ka <- exp(tka + eta.ka) # nolint: object_usage_linter.
    cl <- exp(tcl) # nolint: object_usage_linter.
    v <- exp(tv) # nolint: object_usage_linter.
    cp <- linCmt() # nolint: object_usage_linter.
    cp ~ add(add.sd)
  })
}

theoFit <- .sirLazy(suppressMessages(
  nlmixr2utils::nlmixr2(
    theoOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L, covMethod = "r")
  )
))

# Same fit with the covariance step suppressed; exercises the fallback paths.
theoFitNoCov <- .sirLazy(suppressMessages(
  nlmixr2utils::nlmixr2(
    theoOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L, covMethod = "")
  )
))

# Three-eta variant, used by the tests that need more than one omega element.
threeEtaOneCmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.00
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7 # nolint: object_usage_linter.
  })
  model({
    ka <- exp(tka + eta.ka) # nolint: object_usage_linter.
    cl <- exp(tcl + eta.cl) # nolint: object_usage_linter.
    v <- exp(tv + eta.v) # nolint: object_usage_linter.
    linCmt() ~ add(add.sd)
  })
}

threeEtaFit <- .sirLazy(suppressMessages(suppressWarnings(
  nlmixr2utils::nlmixr2(
    threeEtaOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L),
    table = list(npde = TRUE, cwres = TRUE)
  )
)))

# Correlated-eta variant. The only fixture with an OMEGA off-diagonal, so it
# is what pins the cov.<eta1>.<eta2> naming bridge.
blockOneCmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.00
    tv <- 3.45
    eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)
    add.sd <- 0.7 # nolint: object_usage_linter.
  })
  model({
    ka <- exp(tka + eta.ka) # nolint: object_usage_linter.
    cl <- exp(tcl + eta.cl) # nolint: object_usage_linter.
    v <- exp(tv) # nolint: object_usage_linter.
    linCmt() ~ add(add.sd)
  })
}

blockFit <- .sirLazy(suppressMessages(suppressWarnings(
  nlmixr2utils::nlmixr2(
    blockOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L)
  )
)))

# A single SIR iteration on theoFit(); tiny schedule for speed.
iter1 <- .sirLazy(local({
  fit <- theoFit()
  set.seed(42)
  suppressMessages(
    sirRunIteration(
      fit,
      mu = .sirProposalMu(fit),
      proposalCov = sirGetProposalCov(fit),
      nSamples = 8L,
      nResample = 4L,
      iterNum = 1L,
      recenter = TRUE,
      boxcox = TRUE,
      directory = NULL
    )
  )
}))

# A minimal nlmixr2SIR object for the S3 method tests.
sirObj <- .sirLazy(local({
  it <- iter1()
  out <- sirSummary(it$resampledMat, theoFit())
  class(out) <- c("nlmixr2SIR", "data.frame")
  attr(out, "iterationSummary") <- it$iterSummary
  attr(out, "iterations") <- list(it)
  attr(out, "resampledMat") <- it$resampledMat
  attr(out, "outputDir") <- tempdir()
  out
}))
