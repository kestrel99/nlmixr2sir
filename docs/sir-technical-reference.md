# `nlmixr2sir` technical reference

This document describes the sampling importance resampling (SIR) procedure
implemented by `nlmixr2sir`. It is a user-facing specification of the current
source code, not a general recipe for parameter-uncertainty estimation. Source
references name functions rather than line numbers so that they remain useful
as the package evolves.

## Scope and statistical target

A SIR run produces an empirical sample from the parameter uncertainty
distribution of a fitted nonlinear mixed-effects model, without assuming that
distribution is multivariate normal. Given a fit with parameter estimate
$\widehat\psi$ and objective function value $\mathrm{OFV}(\widehat\psi)$, each
iteration:

1. draws $M$ candidate vectors $\psi_1,\ldots,\psi_M$ from a proposal
   distribution $g$;
2. evaluates each candidate's objective function with the population
   parameters fixed, giving
   $\Delta\mathrm{OFV}_i = \mathrm{OFV}(\psi_i) - \mathrm{OFV}(\widehat\psi)$;
3. forms importance ratios comparing the target density to $g$;
4. resamples $m < M$ vectors without replacement in proportion to those
   ratios; and
5. rebuilds $g$ from the empirical covariance of the resampled vectors,
   feeding the next iteration.

The retained vectors after the final iteration are the deliverable. They
support nonparametric intervals, an empirical covariance, and — unlike a
normal approximation — asymmetric intervals.

The method is due to [Dosne, Bergstrand, Harling & Karlsson
(2016)](https://doi.org/10.1007/s10928-016-9487-8), which introduced SIR for
this purpose and characterised its behaviour against the covariance step,
bootstrap, and log-likelihood profiling. The iterative, self-tuning form
implemented here — repeated rounds with the proposal rebuilt each time, sample
counts compensated for evaluation failures — follows [Dosne, Bergstrand &
Karlsson (2017)](https://doi.org/10.1007/s10928-017-9542-0). See
[Literature](#literature).

What the output is **not**: it is not a Bayesian posterior (there is no prior,
and the target is the likelihood surface around the maximum likelihood
estimate), it is not a bootstrap (the data are never resampled), and it does
not correct a misspecified model. A SIR interval is a statement about
parameter uncertainty under the fitted model, conditional on that model being
the right one.

## Public interface

| Function | Purpose |
| --- | --- |
| [`runSIR()`](../R/sir-run.R) | Run the iterative SIR workflow and write its artifacts |
| [`runSIRControl()`](../R/sir-control.R) | Construct and validate the run settings |
| [`sirSummary()`](../R/sir-results.R) | Empirical summary of a resampled matrix |
| `print()`, `plot()` | S3 methods on the returned `nlmixr2SIR` object |

`runSIR()` takes the fit, the sampling schedule, and where output goes;
everything that tunes *how* the run behaves lives in `runSIRControl()`.

```r
sir <- runSIR(
  fit,
  nSamples  = c(1000, 1000, 1000, 2000, 2000),
  nResample = c(200, 400, 500, 1000, 1000),
  control   = runSIRControl(thetaInflation = 2, workers = 4, rxThreads = 2)
)
```

## The parameter vector

Three naming conventions meet in this package, and deriving them independently
in each consumer is what once let the code paths drift apart.
[`.sirParamSpace()`](../R/sir-paramspace.R) is the single source of truth: one
row per estimated, non-fixed parameter, carrying every name that parameter is
known by.

| kind | SIR name | `fit$cov` rowname | raw-results column |
| --- | --- | --- | --- |
| `theta` | `tka` | `tka` | `tka` |
| `sigma` | `add.sd` | `add.sd` | `add.sd` |
| `omegaDiag` | `eta.ka` | `om.eta.ka` | `omega(eta.ka,eta.ka)` |
| `omegaOffdiag` | `eta.cl:eta.ka` | `cov.eta.cl.eta.ka` | `omega(eta.cl,eta.ka)` |

Rows are ordered THETA, then residual error, then the OMEGA lower triangle by
column and then row. That ordering reproduces `rownames(fit$cov)` exactly, so
the fitted covariance can be consumed without reordering; it also coincides
with PsN's `parameter_hash` order, which makes positional comparison against
PsN output meaningful.

Two details matter for correctness. OMEGA off-diagonals are ordered
`(neta1, neta2)`, not alphabetically, and eta names may contain periods — so
`cov.eta.cl.eta.ka` cannot be split back unambiguously. Every name is derived
from the eta index pair, never parsed. And `kind == "sigma"` is decided
structurally from `iniDf$err`, not from absence from `fit$cov`; the latter
silently changed meaning when `nlmixr2est` began reporting OMEGA in `fit$cov`.

## The initial proposal

Four sources are available, dispatched in this order by
[`.sirResolveInitialProposal()`](../R/sir-proposal-input.R). They are mutually
exclusive and validated at control construction, rather than resolved silently
by precedence.

### From the fitted covariance

The default. Since `nlmixr2est` 7, `foceiControl(covFull = TRUE)` is the
default and `fit$cov` spans THETA, residual error, and OMEGA jointly, so

$$
  g_1 = N\!\left(\widehat\psi,\ \widehat V\right),
  \qquad \widehat V = \texttt{fit\$cov},
$$

including the estimated THETA-OMEGA correlations. A block-diagonal
construction discards those correlations and is therefore needlessly wide in
exactly the directions SIR must then claw back over subsequent iterations.

### Wishart-style fallback

When `fit$cov` does not carry OMEGA — `covFull = FALSE`, a failed covariance
step, or `covMethod = ""` — or when `omegaFallback = "wishart"` forces it,
OMEGA uncertainty is approximated with $df = n_{\text{sub}} - 1$ by default:

$$
  \operatorname{Var}(\widehat\Omega_{jj}) = \frac{2\widehat\Omega_{jj}^2}{df},
  \qquad
  \operatorname{Var}(\widehat\Omega_{jk})
    = \frac{\widehat\Omega_{jj}\widehat\Omega_{kk}
            + \widehat\Omega_{jk}^2}{df}.
$$

This route gives a block-diagonal proposal. Which route was taken is reported
in the run log and on the returned object.

### From relative standard errors

`rseTheta`, `rseOmega` and `rseSigma` build a diagonal proposal with no
covariance step at all, porting PsN's `setup_variancevec_from_rse()`. Diagonal
variances are $(\mathrm{rse}\cdot\widehat\psi_j/100)^2$. OMEGA off-diagonals
use

$$
  N = \left(\frac{100}{\mathrm{rse}_j}\right)^2
    + \left(\frac{100}{\mathrm{rse}_k}\right)^2 + 1,
  \qquad
  \operatorname{Var}(\widehat\Omega_{jk})
    = \frac{\widehat\Omega_{jk}^2
            + \widehat\Omega_{jj}\widehat\Omega_{kk}}{N}.
$$

PsN's documentation describes this off-diagonal rule as choosing the variance
"so that the correlation from the final estimate is unchanged". Its
implementation is the Wishart-style expression above, which is not the same
thing. This package follows the implementation, since that is what PsN runs.

Each argument is a scalar for the whole class or one value per estimated
element of it. Following PsN, a scalar `rseTheta` fills in an unset `rseOmega`
and `rseSigma`, a vector `rseTheta` does not, and setting `rseOmega` without
`rseTheta` is an error. Inflation cannot be combined with this route: the RSE
already states the width.

### Supplied directly, or from previous parameter vectors

`covmatInput` accepts a matrix, a NONMEM-style `.cov` file, or `"identity"`;
the last together with inflation is the cheap "any diagonal proposal" route.
`rawresInput` seeds the first proposal from the parameter vectors in a
canonical raw-results file — PsN's iteration 0 — taking their empirical mean
and covariance, with `offsetRawres` and `inFilter` narrowing which rows are
used. Any canonical raw-results file works, including one written by
`nlmixr2boot`.

## Inflation, correlation capping, and positive-definiteness

Inflation multiplies each parameter's proposal *variance*, preserving
correlations, and is applied to the initial proposal only — from the second
iteration the proposal is the previous empirical covariance and re-inflating
it each round would compound. [`.sirInflationVector()`](../R/sir-proposal.R)
ports PsN's `setup_inflation()`: each argument is a scalar or one value per
*diagonal* element of its class, and an OMEGA off-diagonal is never given a
factor directly but derives

$$
  c_{jk} = \sqrt{c_{jj}}\,\sqrt{c_{kk}},
$$

which leaves the correlation unchanged when the two factors are equal. The
rescaling is implemented as
$\Sigma_{jk} \mapsto \Sigma_{jk}\sqrt{c_j c_k}$, algebraically identical to
rescaling standard deviations but requiring no `cov2cor()` — which would fail
on a variance of exactly zero, as the Wishart fallback produces for an OMEGA
element estimated at zero.

[`.sirCapCovCorrelation()`](../R/sir-utils.R) then clamps every off-diagonal
correlation to $\pm$ `capCorrelation` (default 0.8) while holding the
standard deviations fixed, and [`.sirEnsurePosDef()`](../R/sir-utils.R)
symmetrises and floors the eigenvalues at $\sqrt{\varepsilon}$. Capping a
correlation can itself destroy positive-definiteness, so the order matters.

## Sampling and rejection

[`.sirSampleFullProposal()`](../R/sir-proposal.R) draws from the multivariate
normal proposal in batches until $M$ valid vectors are collected or a budget
of `maxAttemptFactor * M` draws is exhausted. A draw is rejected if:

- back-transformation from the Box-Cox scale fails (`inverseRejected`);
- any THETA or residual-error parameter falls outside its `iniDf` bounds
  (`thetaRejected`, `sigmaRejected`); or
- the reconstructed OMEGA matrix is not positive-definite, tested by Cholesky
  (`omegaRejected`).

OMEGA elements are constrained by the positive-definiteness test rather than
by element-wise bounds, which keeps the four rejection counts separable — they
are reported per iteration and written to `sample_rejection_summary.txt`. A
run that rejects heavily in one category is diagnosable; a single pooled count
would not be.

## Objective function evaluation

[`sirEvalOFV()`](../R/sir-eval.R) sets each sampled vector into the model with
`rxode2::ini()` and evaluates the objective with
`foceiControl(maxOuterIterations = 0L)`, so no estimation occurs — the
population parameters are fixed at the proposed values and the inner problem
is solved. Evaluation is parallelised across `workers`, each worker using
`rxThreads` rxode2 threads; whenever `workers > 1`, `workers * rxThreads` must
not exceed the core count, since each worker is a separate process with its
own thread pool.

Failed evaluations return `NA` rather than aborting the run, but the
underlying error messages are retained and surfaced if *every* evaluation
fails. A configuration fault fails all samples identically, and reporting only
"all evaluations failed" hides the cause.

## Importance weights

For proposal $g = N(\mu,\Sigma)$ with Cholesky factor $L$,
[`sirCalcWeights()`](../R/sir-weights.R) computes a *relative* proposal
density, normalised so that $\psi_i = \mu$ gives exactly 1:

$$
  \log \mathrm{relPDF}_i
    = -\tfrac12 \left\| L^{-\mathsf T}(\psi_i - \mu) \right\|^2 .
$$

The likelihood ratio relative to the fit is
$\log \mathrm{LR}_i = -\tfrac12 \Delta\mathrm{OFV}_i$, so the importance ratio
and resampling probability are

$$
  \log \mathrm{IR}_i = \log \mathrm{LR}_i - \log \mathrm{relPDF}_i,
  \qquad
  p_i = \frac{\exp(\log \mathrm{IR}_i - \max_k \log \mathrm{IR}_k)}
             {\sum_j \exp(\log \mathrm{IR}_j - \max_k \log \mathrm{IR}_k)} .
$$

Everything is carried on the log scale and the maximum is subtracted before
exponentiating, because raw importance ratios overflow readily. Normalising
constants common to all samples cancel in $p_i$ and are never formed.

## Resampling

[`sirResample()`](../R/sir-weights.R) draws $m$ vectors **without
replacement** with probability proportional to $p_i$. Sampling without
replacement is what distinguishes SIR here from a naive importance sample: it
prevents a single high-weight vector from dominating the retained set, at the
cost of requiring $m < M$.

`capResampling` relaxes this. A value $c > 1$ expands each candidate into $c$
slots before drawing, so a vector may be selected up to $c$ times — limited
replacement, with the cap bounding how far any one vector can dominate.

## Proposal update and Box-Cox

[`sirUpdateProposal()`](../R/sir-boxcox.R) rebuilds the proposal from the
empirical covariance of the retained vectors. With `boxcox = TRUE` each
column is first transformed by

$$
  x^{(\lambda)} = \begin{cases}
    \dfrac{(x+\delta)^{\lambda} - 1}{\lambda}, & \lambda \neq 0 \\[2ex]
    \log(x+\delta), & \lambda = 0
  \end{cases}
$$

with $\delta = |\min x| + 10^{-6}$ guaranteeing positivity, and $\lambda$
chosen on $[-3,3]$ to maximise the correlation between the sorted transformed
values and normal scores — a normality-of-fit criterion rather than a
profile likelihood. The covariance is then taken on the transformed scale, so
the next iteration proposes in a space where the parameters are closer to
normal, and draws are back-transformed before evaluation.

The Box-Cox transform is deliberately **not** applied on the final iteration:
the last proposal is built on the original scale so the delivered vectors and
their covariance need no back-transformation.

If `recenter = TRUE` and any sample has $\Delta\mathrm{OFV} < 0$ — meaning a
proposed vector fits better than the reported maximum likelihood estimate, so
the fit was not fully converged — the proposal is recentred on that vector.
With `recenter = FALSE` the same condition warns instead.

## Sample count adjustment

Evaluation failures shrink the usable sample, so both counts are compensated
per iteration, porting PsN's `update_attempted_samples()` and
`update_actual_resamples()` in [`R/sir-iterate.R`](../R/sir-iterate.R). With
turnout $t$ = successful / requested:

$$
  M' = \begin{cases}
    \mathrm{round}(M / t_{\text{prev}}), & t_{\text{prev}} \le 0.95 \\
    M, & \text{otherwise}
  \end{cases}
  \qquad
  m' = \begin{cases}
    \mathrm{round}(m \cdot t), & |t - 1| \ge 0.05 \\
    m, & \text{otherwise}
  \end{cases}
$$

The attempted count compensates for loss only; the resample count scales on
gain *or* loss. Turnout for the resample adjustment is measured against the
originally requested sample count, not the compensated attempted count.

`round()` here is round-half-away-from-zero, PsN's convention, implemented as
[`.sirRound()`](../R/sir-utils.R). R's own `round()` is round-half-to-even and
disagrees on exact halves — `round(20.5)` is 21 in PsN and 20 in R. Both rules
are pinned to PsN's oracle sequence: requested 100 samples and 20 resamples
with successful counts 90, 102, 109, 98 must give attempted 100, 111, 109, 100
and resamples 18, 20, 22, 20.

## Diagnostics

### Convergence: dOFV against a reference chi-square

The primary SIR diagnostic, and the reason the method is trusted. For a
quantile grid $q$ stopping short of 1 so the reference stays finite,
[`.sirDofvCurves()`](../R/sir-convergence.R) draws three curves per iteration:

| curve | definition |
| --- | --- |
| reference | $\chi^2_{p}$ quantiles, $p$ = number of estimated parameters |
| proposal | empirical $\Delta\mathrm{OFV}$ quantiles over all evaluated samples |
| SIR | empirical $\Delta\mathrm{OFV}$ quantiles over the resampled subset |

Under the asymptotic theory the SIR curve should approach the reference.
Convergence reads as that curve settling onto it across iterations, with a
resampling-noise band on the last two iterations obtained by repeating the
weighted resampling and taking the 2.5th and 97.5th percentiles of the
resulting curves.

If the first iteration's proposal curve falls *below* the reference for more
than a quarter of the quantiles, the proposal is too narrow: the vectors SIR
would need were never drawn, and resampling cannot manufacture them. The run
warns and recommends restarting with inflation. This check is PsN's and is
most of the practical value of the plot.

### Intervals by iteration, and CI asymmetry

`plot(type = "intervals")` shows the proposal and SIR interval per parameter
per iteration; uncertainty that is still moving between the last two
iterations means the run has not settled.

`plot(type = "rsecor")` draws RSE% on the diagonal and correlations off it,
annotating the diagonal with

$$
  \text{asymmetry} = \frac{P_{\text{high}} - P_{\text{med}}}
                          {P_{\text{med}} - P_{\text{low}}},
$$

banded at 0.5, 1, 1.25 and 2. A symmetric normal-approximation covariance
reports one standard error per parameter and cannot express that a parameter's
upper interval half is twice its lower. Showing it is a large part of why SIR
is run at all.

## Summaries and artifacts

[`sirSummary()`](../R/sir-results.R) reports `estimate`, `mean`, `sd`, `rse`,
`rse_sd_scale`, and PsN's percentile set — 2.5, 5, 10, 30, 50, 70, 90, 95,
97.5, derived from prediction intervals 0, 40, 80, 90 and 95. Empirical
covariance, correlation, and standard-deviation/correlation matrices are
attached as attributes and written to disk.

Two documented differences from PsN. `rse` is a **percentage** where PsN
reports a fraction; the returned object records this in an `rseUnits`
attribute rather than leaving it implicit. And `rse_sd_scale` halves the RSE
of OMEGA elements only, where PsN halves everything that is not a NONMEM
THETA: that rule catches `$SIGMA` because NONMEM parameterises residual error
as a variance, whereas nlmixr2 parameterises it on the standard-deviation
scale already, so halving `add.sd` would rescale a quantity that needs no
rescaling.

A run writes `sir_results.csv`, `summary_iterations.csv` (leading with PsN's
column names), `<fitName>_sir.cov` and `.sdcorr`, a canonical
`raw_results.*` set, `sample_rejection_summary.txt`, and `sir_state.rds`.
`runSIR()` also registers the empirical covariance so that
`nlmixr2est::setCov(fit, "sir")` switches the fit's reported uncertainty to
the SIR result, skipping registration with a message if the parameters do not
match `fit$cov` or the covariance is not positive-definite.

## Persistence, resume and extension

State is written after every iteration. `recover = TRUE` resumes from the last
completed iteration, returning the stored result unchanged if the schedule was
already finished. `addIterations = TRUE` appends further iterations to a
completed run, carrying the existing iterations over rather than recomputing
them. Seeding is managed per iteration through `nlmixr2utils::withRunSeed()`,
so a resumed run reproduces the stream it would have had.

## Interpretation checklist and limitations

- **Read the convergence plot first.** If the iteration-1 proposal sits below
  the reference chi-square, the result is not usable; restart with inflation.
- **Check the sample-to-resample ratio.** Roughly 5:1 is the working default.
  Heavy rejection in any one category points at bounds, positive-definiteness,
  or a proposal on the wrong scale.
- **Check that the last two iterations agree.** If intervals are still moving,
  add iterations or increase sample counts.
- SIR characterises uncertainty *under the fitted model*. It cannot diagnose
  structural misspecification, and a tight SIR interval around a wrong model
  is still wrong.
- The target is the likelihood surface near $\widehat\psi$. With a
  poorly-identified parameter the surface may be flat or multimodal, and SIR
  will report that flatness faithfully rather than resolving it.
- Default schedules are expensive: the objective is evaluated once per sample,
  which at the default 7,000 samples is minutes to hours depending on the
  model. `workers` and `rxThreads` are the practical levers.

## Difference from PsN

The algorithm, the sample-count adjustment rules, the inflation semantics, the
RSE-to-variance conversion and the principal diagnostics follow PsN, and the
numeric cores are checked against oracle values taken from PsN's own unit
tests (`test/unit/tool/sir.t`). The implementations diverge in execution:
PsN generates and runs NONMEM control streams, whereas `nlmixr2sir` calls
`rxode2::ini()` and `nlmixr2est::nlmixr2()` directly and keys everything by
canonical named parameter schemas rather than by position.

Four PsN options are not implemented: `-auto_rawres`, `-print_iter`,
`-fast_posdef_checks`, and the `rplots_level = 2` extras (bin-exhaustion
diagnostics and inverse-Wishart degrees-of-freedom estimation). The
NONMEM-execution options — `-mceta`, `-copy_data`, `-problems_per_file`,
`-nm_version` and similar — have no analogue. A per-option parity matrix is
kept in the [README](../README.md).

The deliberate numeric differences are the two described under
[Summaries and artifacts](#summaries-and-artifacts), and the off-diagonal RSE
rule following PsN's code rather than its documentation, described under
[From relative standard errors](#from-relative-standard-errors).

This comparison is based on PsN's
[`tool::sir`](https://github.com/UUPharmacometrics/PsN/blob/master/lib/tool/sir.pm),
its diagnostic template
[`sir_default.R`](https://github.com/UUPharmacometrics/PsN/blob/master/R-scripts/sir_default.R),
and the
[SIR user guide](https://github.com/UUPharmacometrics/PsN/releases/download/v5.7.0/sir_userguide.pdf).

## Literature

1. Dosne A-G, Bergstrand M, Harling K, Karlsson MO. Improving the estimation
   of parameter uncertainty distributions in nonlinear mixed effects models
   using sampling importance resampling. *Journal of Pharmacokinetics and
   Pharmacodynamics*. 2016;43:583-596.
   [doi:10.1007/s10928-016-9487-8](https://doi.org/10.1007/s10928-016-9487-8).
   This is the original SIR paper and the primary reference for the method
   implemented here: the importance-ratio construction, the use of $\Delta$OFV
   against a reference chi-square as the convergence criterion, and the
   comparison against the covariance step, bootstrap, and log-likelihood
   profiling.

2. Dosne A-G, Bergstrand M, Karlsson MO. An automated sampling importance
   resampling procedure for estimating parameter uncertainty. *Journal of
   Pharmacokinetics and Pharmacodynamics*. 2017;44:509-520.
   [doi:10.1007/s10928-017-9542-0](https://doi.org/10.1007/s10928-017-9542-0).
   This develops the iterative, self-tuning procedure that this package
   implements — multiple rounds with the proposal rebuilt from each round's
   resampled vectors, and sample counts adjusted for evaluation failures.

3. Rubin DB. Using the SIR algorithm to simulate posterior distributions. In:
   *Bayesian Statistics 3*. Oxford University Press; 1988:395-402. The
   origin of sampling importance resampling as a general technique, on which
   the pharmacometric application above builds.

4. Box GEP, Cox DR. An analysis of transformations. *Journal of the Royal
   Statistical Society, Series B*. 1964;26:211-252.
   [doi:10.1111/j.2517-6161.1964.tb00553.x](https://doi.org/10.1111/j.2517-6161.1964.tb00553.x).
   The transformation used to make the inter-iteration proposal closer to
   normal.

5. Lindbom L, Ribbing J, Jonsson EN. Perl-speaks-NONMEM (PsN) — a Perl module
   for NONMEM related programming. *Computer Methods and Programs in
   Biomedicine*. 2004;75:85-94.
   [doi:10.1016/j.cmpb.2003.11.003](https://doi.org/10.1016/j.cmpb.2003.11.003).
   The reference implementation this package is checked against.

## Implementation index

- Orchestration, persistence, and return construction:
  [`runSIR()`](../R/sir-run.R)
- Control validation: [`runSIRControl()`](../R/sir-control.R)
- Parameter vector and naming bridge:
  [`.sirParamSpace()`, `.sirProposalMu()`](../R/sir-paramspace.R)
- Proposal construction, fallback uncertainty, inflation, and sampling:
  [`sirGetProposalCov()`, `.sirFallbackSe()`, `.sirInitialProposal()`,
  `.sirInflationVector()`, `.sirSampleFullProposal()`](../R/sir-proposal.R)
- Alternative proposal sources:
  [`.sirRseVariance()`, `.sirProposalFromCovmatInput()`,
  `.sirProposalFromRawResults()`,
  `.sirResolveInitialProposal()`](../R/sir-proposal-input.R)
- Objective function evaluation: [`sirEvalOFV()`](../R/sir-eval.R)
- Weights, resampling, and per-iteration raw results:
  [`sirCalcWeights()`, `sirResample()`,
  `.sirBuildRawResults()`](../R/sir-weights.R)
- Box-Cox and proposal update:
  [`sirBoxCox()`, `sirUpdateProposal()`](../R/sir-boxcox.R)
- One iteration, and the PsN sample-count adjustments:
  [`sirRunIteration()`, `.sirAdjustedAttemptedSamples()`,
  `.sirAdjustedResamples()`](../R/sir-iterate.R)
- Convergence diagnostic:
  [`.sirDofvCurves()`, `.sirDofvNoise()`,
  `.sirProposalTooNarrow()`](../R/sir-convergence.R)
- Interval and RSE/correlation diagnostics:
  [`.sirIterationIntervals()`, `.sirRseCorData()`](../R/sir-diagnostics.R)
- Summaries and on-disk artifacts:
  [`sirSummary()`, `.sirWriteIterationSummary()`,
  `.sirWriteCovMatrices()`](../R/sir-results.R)
- `setCov()` registration:
  [`.sirCovAsFitCov()`, `.sirRegisterCov()`](../R/sir-setcov.R)
- Matrix helpers: [`.sirCapCovCorrelation()`, `.sirEnsurePosDef()`,
  `.sirRound()`](../R/sir-utils.R)
- Print and plot methods: [`R/sir-methods.R`](../R/sir-methods.R)
