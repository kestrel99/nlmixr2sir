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

## The details

`runSIR()` implements the SIR workflow for `nlmixr2` Where practical, it follows
the same process as Perl-speaks-NONMEM.

`nlmixr2sir` builds the initial proposal from `fit$cov`, with optional inflation and
correlation capping. THETA, OMEGA, and sigma-like parameters are sampled, and 
parameter-space constraints such as bounds and positive-definite OMEGA matrices are
enforced.

Because `nlmixr2` does not currently report standard errors for OMEGA elements,
`nlmixr2sir` initializes OMEGA uncertainty separately from `fit$omega` instead
of reading it from `fit$cov` or `fit$parFixedDf`. The package takes the free
lower-triangular OMEGA elements, applies a Wishart-style fallback
(`omegaFallback = "wishart"`), and uses that to approximate an SE for each
sampled OMEGA element. By default the fallback uses `omegaDf = nSubjects - 1`,
so diagonal OMEGA elements use `sqrt(2 * omega^2 / df)` and off-diagonal
elements use `sqrt((omega[i, i] * omega[j, j] + omega[i, j]^2) / df)`.

Those fallback SEs define the initial OMEGA proposal block, with only the free
lower-triangular elements sampled directly. Each proposed vector is then
reconstructed into an OMEGA matrix, and non-positive-definite draws are
discarded. After each SIR iteration, the next proposal covariance is updated
from the empirical covariance of the retained samples, so later iterations are
not limited to the initial diagonal OMEGA approximation.

Sampled vectors are re-evaluated by fixing parameters in the model and running
Bayesian feedback (i.e. the model is evaluated aginst the data without any
estimation being performed). Importance ratios are computed, weighted resampling
is performed, and the proposal for the next iteration is updated.

The SIR tool supports recentering, Box-Cox proposal updates, recovery from saved state,
iteration summaries, and diagnostic plots.

The package is designed to work alongside `nlmixr2utils`, which provides the
shared worker-plan helpers and core infrastructure.

## Installation

The package is not on CRAN. Install it from GitHub together with
`nlmixr2utils`.

Using `pak`:

```r
pak::pkg_install(c(
  "kestrel99/nlmixr2utils",
  "kestrel99/nlmixr2sir"
))
```

Using `remotes`:

```r
remotes::install_github("kestrel99/nlmixr2utils")
remotes::install_github("kestrel99/nlmixr2sir")
```

## Basic Use

```r
library(nlmixr2)
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

fit <- nlmixr2(
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

## Acknowledgments

The SIR methodology used here is based primarily on the method of [Dosne 
et al](https://link.springer.com/article/10.1007/s10928-016-9487-8). The technical implementation is based heavily on the 
[PsN tool](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/sir_userguide.pdf).

## References 

* Dosne, A.-G., Bergstrand, M., Harling, K., & Karlsson, M.O. (2016).
  Improving the estimation of parameter uncertainty distributions in nonlinear
  mixed effects models using sampling importance resampling. *Journal of
  Pharmacokinetics and Pharmacodynamics*, 43(6), 583-596.

For a fuller worked example, see the package vignette:
`vignette("runSIR", package = "nlmixr2sir")`.
