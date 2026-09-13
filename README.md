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

`nlmixr2sir` builds the initial proposal from `fit$cov` by default, with
optional inflation and correlation capping, and can take it instead from
relative standard errors, a supplied covariance matrix, or the parameter
vectors in a raw-results file -- see *Requirements and Practical Notes*. THETA,
OMEGA, and residual-error parameters are sampled, and parameter-space
constraints such as bounds and positive-definite OMEGA matrices are enforced.

Since `nlmixr2est` 7, `foceiControl(covFull = TRUE)` is the default and `fit$cov`
covers THETA, residual error, and OMEGA together. `nlmixr2sir` uses that matrix
directly (`omegaFallback = "cov"`, the default), which means the initial
proposal carries the **correlations between THETA and OMEGA** rather than
treating the two as independent blocks.

When `fit$cov` does not carry OMEGA -- a fit run with `covFull = FALSE`, a
failed covariance step, or `covMethod = ""` -- `nlmixr2sir` falls back
automatically to a Wishart-style approximation, and `omegaFallback = "wishart"`
forces it. The fallback takes the free lower-triangular OMEGA elements and, with
`omegaDf = nSubjects - 1` by default, approximates diagonal SEs as
`sqrt(2 * omega^2 / df)` and off-diagonal SEs as
`sqrt((omega[i, i] * omega[j, j] + omega[i, j]^2) / df)`. That route gives a
block-diagonal proposal, so it discards the THETA-OMEGA correlations the
default route keeps. The route actually taken is reported in the run log.

Only the free lower-triangular elements are sampled directly. Each proposed
vector is reconstructed into an OMEGA matrix, and non-positive-definite draws
are discarded. After each SIR iteration, the next proposal covariance is
updated from the empirical covariance of the retained samples.

Sampled vectors are re-evaluated by fixing the population parameters and
recomputing the objective function against the data, without any estimation
being performed. Importance ratios are computed, weighted resampling is
performed, and the proposal for the next iteration is updated.

The SIR tool supports recentering, Box-Cox proposal updates, recovery from
saved state, extending a finished run with further iterations, iteration
summaries, and diagnostic plots.

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

Everything that tunes *how* the run behaves lives in `runSIRControl()`:

```r
sir <- runSIR(
  fit,
  nSamples = c(1000L, 1000L, 1000L, 2000L, 2000L),
  nResample = c(200L, 400L, 500L, 1000L, 1000L),
  control = runSIRControl(
    thetaInflation = 2,
    workers = 4,
    rxThreads = 2
  )
)
```

## Diagnostics

```r
plot(sir, type = "convergence")   # dOFV vs reference chi-square, per iteration
plot(sir, type = "intervals")     # proposal vs SIR interval, per parameter
plot(sir, type = "rsecor")        # RSE / correlation, with CI asymmetry
plot(sir, type = "parameters")    # resampled parameter distributions
```

`type = "convergence"` is the primary diagnostic. For each iteration it draws
the empirical dOFV quantile curve for the proposal and for the SIR posterior
against a reference chi-square on the number of estimated parameters.
Convergence reads as the SIR curve settling onto the reference. If the first
iteration's proposal falls below the reference for more than a quarter of the
quantiles, `runSIR()` warns: the proposal is too narrow, and resampling cannot
recover from that -- restart with inflation.

`type = "rsecor"` annotates each parameter's RSE with the confidence-interval
asymmetry ratio `(high - median) / (median - low)`. A symmetric
normal-approximation covariance reports one standard error per parameter and
cannot express that asymmetry, which is a large part of why SIR is run at all.

## Parity with PsN

| PsN option | `nlmixr2sir` | Status |
|---|---|---|
| `-samples` | `nSamples` | supported |
| `-resamples` | `nResample` | supported |
| covariance matrix from the fit | default | supported |
| `-rse_theta` / `-rse_omega` / `-rse_sigma` | `rseTheta` / `rseOmega` / `rseSigma` | supported |
| `-covmat_input=<file>` / `=identity` | `covmatInput` | supported |
| `-rawres_input` | `rawresInput` | supported |
| `-offset_rawres` | `offsetRawres` | supported |
| `-in_filter` | `inFilter` | supported |
| `-theta_inflation` etc., scalar or vector | `thetaInflation` etc. | supported |
| `-inflate_only_diagonal` semantics | always applied | supported |
| `-recenter` | `recenter` | supported |
| `-boxcox` | `boxcox` | supported |
| `-cap_resampling` | `capResampling` | supported |
| `-cap_correlation` | `capCorrelation` | supported |
| `-add_iterations` | `addIterations` | supported |
| sample / resample count adjustment | automatic | supported, oracle-tested |
| dOFV vs chi-square plot | `plot(type = "convergence")` | supported |
| CI-by-iteration plot | `plot(type = "intervals")` | supported |
| RSE / correlation plot | `plot(type = "rsecor")` | supported |
| `empirical_statistics()` output | `sirSummary()` | supported |
| `<model>_sir.cov` | `<fitName>_sir.cov` | supported |
| `-auto_rawres` | — | not implemented |
| `-print_iter` | — | not implemented |
| `-fast_posdef_checks` | — | not implemented |
| `rplots_level = 2` extras | — | not implemented |
| `-mceta`, `-copy_data`, `-problems_per_file`, `-nm_version` | — | not applicable (NONMEM execution) |

Numeric parity for the sample/resample adjustment, the inflation vector and
the RSE-to-variance conversion is checked against oracle values taken from
PsN's own unit tests.

Two deliberate differences. `sirSummary()` reports `rse` as a **percentage**
where PsN reports a fraction. And `rse_sd_scale` halves the RSE of OMEGA
elements only: PsN halves everything that is not a NONMEM THETA, which catches
`$SIGMA` because NONMEM parameterises residual error as a variance, whereas
nlmixr2 parameterises it on the standard-deviation scale already.

## Requirements and Practical Notes

`runSIR()` no longer requires a successful covariance step. When `fit$cov` is
unavailable, supply the proposal another way:

```r
# from relative standard errors
runSIR(fit, control = runSIRControl(rseTheta = 30))

# from a diagonal proposal, widened by inflation
runSIR(fit, control = runSIRControl(covmatInput = "identity",
                                    thetaInflation = 0.05))

# seeded from the parameter vectors in a raw-results file
runSIR(fit, control = runSIRControl(rawresInput = "raw_results.csv"))
```

For practical use:

* Use larger production schedules than toy examples; the default PsN-style
  schedule is usually a good starting point.
* Review the convergence diagnostic to make sure the proposal is not too
  narrow or too wide.
* Use `workers` and `rxThreads` to parallelize OFV evaluation when runs are
  large enough to justify it. Whenever `workers > 1`, `workers * rxThreads`
  must not exceed the machine's core count.
* `nlmixr2est::setCov(fit, "sir")` switches the fit's reported uncertainty to
  the SIR result after a run.

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

## Credit where it's due

`nlmixr2sir` is based on the [PsN implementation](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/sir_userguide.pdf) written by Lars Lindbom, 
Niclas Jonsson, Pontus Pihlgren, Mats Karlsson, Andrew Hooker, Kajsa Harling, 
Rikard Nordgren and Svetlana Freiberga.