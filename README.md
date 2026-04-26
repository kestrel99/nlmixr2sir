# nlmixr2sir

`nlmixr2sir` provides Sampling Importance Resampling (SIR) for `nlmixr2`
population PK/PD models.

SIR is a simulation-based approach for parameter uncertainty estimation. Rather
than relying only on the asymptotic covariance matrix around the maximum
likelihood estimate, it repeatedly:

1. Samples parameter vectors from a proposal distribution centered on the MLE.
2. Evaluates each sampled vector by fixing parameters and recomputing the OFV.
3. Weights each sample by the ratio of target likelihood to proposal density.
4. Resamples according to those importance weights.
5. Updates the proposal from the empirical covariance of the resampled vectors.

The result is a set of empirical draws from the parameter uncertainty
distribution that can be summarized with nonparametric intervals, covariance
matrices, and diagnostic plots.

## What the Package Implements

`runSIR()` implements the SIR workflow for `nlmixr2` FOCE-I fits, following the
PsN `sir` algorithm where practical.

In particular, `nlmixr2sir`:

* Builds the initial proposal from `fit$cov`, with optional inflation and
  correlation capping.
* Samples THETA, OMEGA, and sigma-like parameters and enforces parameter-space
  constraints such as bounds and positive-definite OMEGA matrices.
* Re-evaluates sampled vectors by fixing parameters in the model and running
  an nlmixr2 FOCE-I equivalent of NONMEM `MAXEVAL=0`.
* Computes importance ratios, performs weighted resampling, and updates the
  proposal for the next iteration.
* Supports recentering, Box-Cox proposal updates, recovery from saved state,
  iteration summaries, and diagnostic plots.

The package is designed to work alongside `nlmixr2utils`, which provides the
shared worker-plan helpers and core infrastructure used across the split
`nlmixr2` extension packages.

## Installation

The package is not on CRAN. Install it from GitHub together with
`nlmixr2utils`.

Using `pak`:

```r
pak::pkg_install(c(
  "nlmixr2/nlmixr2utils",
  "nlmixr2/nlmixr2sir"
))
```

Using `remotes`:

```r
remotes::install_github("nlmixr2/nlmixr2utils")
remotes::install_github("nlmixr2/nlmixr2sir")
```

If you are working locally with the split repositories:

```r
devtools::install_local("../nlmixr2utils")
devtools::install_local("../nlmixr2sir")
```

## Basic Use

```r
library(nlmixr2data)
library(nlmixr2sir)

one_cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.00
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

fit <- nlmixr2utils::nlmixr2(
  one_cmt,
  data = nlmixr2data::theo_sd,
  est = "focei",
  control = list(print = 0L),
  table = list(npde = TRUE, cwres = TRUE)
)

sir <- runSIR(
  fit,
  nSamples = c(1000L, 1000L, 1000L, 2000L, 2000L),
  nResample = c(200L, 400L, 500L, 1000L, 1000L)
)

print(sir)
plot(sir, type = "parameters")
```

## Requirements and Practical Notes

`runSIR()` expects a fit with a successful covariance step because the initial
proposal is derived from `fit$cov`.

For practical use:

* Use a model with a successful covariance step before calling `runSIR()`.
* Use larger production schedules than toy examples; the default PsN-style
  schedule is usually a good starting point.
* Review the dOFV and resampling diagnostics to make sure the proposal is not
  too narrow or too wide.
* Use `workers` to parallelize OFV evaluation when runs are large enough to
  justify it.

## References

The SIR methodology implemented here is based primarily on:

* Dosne, A.-G., Bergstrand, M., & Karlsson, M.O. (2013). An automated sampling
  importance resampling procedure for estimating parameter uncertainty. *PAGE
  22*, Abstract 2907.
* Dosne, A.-G., Bergstrand, M., Harling, K., & Karlsson, M.O. (2016).
  Improving the estimation of parameter uncertainty distributions in nonlinear
  mixed effects models using sampling importance resampling. *Journal of
  Pharmacokinetics and Pharmacodynamics*, 43(6), 583-596.

The implementation also follows the original PsN SIR workflow:

* Lindbom, L., Ribbing, J., & Jonsson, E.N. (2004). Perl-speaks-NONMEM (PsN):
  a Perl module for NONMEM related programming. *Computer Methods and Programs
  in Biomedicine*, 75, 85-94.

For a fuller worked example, see the package vignette:
`vignette("runSIR", package = "nlmixr2sir")`.
