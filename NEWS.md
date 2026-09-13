# nlmixr2sir (development version)

* SIR now derives its parameter vector from a single internal description of the fit instead of re-deriving THETA, sigma, and OMEGA names independently in each function. This fixes `runSIR()` aborting with `Assertion on 'mu' failed: Contains missing values` against nlmixr2est 7, where `fit$cov` reports OMEGA alongside THETA and the old proposal mean came back `NA` for every OMEGA element.

* OFV evaluation no longer fails in a session where `nlmixr2sir` was attached but `nlmixr2` was not. `rxode2::ini()` evaluates the OMEGA line it is given as a `lotri({...})` call in the caller's environment, and `lotri` is an Imports rather than a Depends of both rxode2 and nlmixr2est, so it was never visible; `nlmixr2sir` now imports it. Every sample previously returned `NA` and `runSIR()` aborted with `All SIR OFV evaluations failed`.

* That abort now reports the first underlying evaluation error, instead of discarding it. A configuration problem fails every sample identically, which the old message hid.

* The package now requires nlmixr2est >= 7.0.0, nlmixr2utils >= 0.3, and rxode2 >= 5.0.0.

* `runSIR()` now takes its run settings through `control = runSIRControl()` rather than as flat arguments, cutting its signature from 21 arguments to 6. Passing a setting directly to `runSIR()` is an error naming `runSIRControl()`.

* `plot()` gains three diagnostics. `type = "convergence"` is the dOFV-versus-chi-square plot: per iteration, the empirical dOFV quantile curve for the proposal and for the SIR posterior against a reference chi-square on the number of estimated parameters, with a resampling-noise band on the last two iterations. Convergence reads as the SIR curve settling onto the reference, and a proposal that falls below the reference for more than a quarter of the quantiles now warns and recommends inflation. `type = "intervals"` compares the proposal and SIR interval per parameter per iteration. `type = "rsecor"` draws the RSE/correlation matrix with the diagonal annotated by the confidence-interval asymmetry ratio, which a symmetric normal approximation cannot show.

* `sirSummary()` now matches PsN's output: it adds `mean` alongside the median, reports PsN's percentile set (2.5, 5, 10, 30, 50, 70, 90, 95, 97.5, from prediction intervals 0/40/80/90/95) in place of the previous set, and adds `rse_sd_scale`. Note `p25` and `p75` are no longer reported. `rse` remains a percentage where PsN reports a fraction, and the returned object now records that in an `rseUnits` attribute.

* `runSIR()` writes `<fitName>_sir.cov` and `<fitName>_sir.sdcorr`, and `summary_iterations.csv` now leads with PsN's column names so either file can be read by PsN-literate tooling.

* `runSIR()` no longer requires a successful covariance step. `runSIRControl()` gains three alternative proposal sources, following PsN: `rseTheta`/`rseOmega`/`rseSigma` build a diagonal proposal from relative standard errors, `covmatInput` takes a covariance matrix, a NONMEM-style `.cov` file, or `"identity"`, and `rawresInput` seeds the first proposal from the parameter vectors in a canonical raw-results file (with `offsetRawres` and `inFilter`). These exist for the models whose covariance step fails, which is where SIR is most wanted.

* `runSIRControl()` inflation arguments accept a vector as well as a scalar: one value per estimated THETA, OMEGA diagonal, or residual-error parameter. OMEGA off-diagonals derive `sqrt(infl_i) * sqrt(infl_j)` from the two diagonals they connect, as PsN does.

* The sample and resample count adjustments now match PsN exactly, and are checked against the oracle values in PsN's own unit tests. The attempted-sample count previously used `ceiling()` where PsN uses round-half-away-from-zero, and the resample count used `floor()` and a strict `>` where PsN uses rounding and `>=`; both also clamped in ways PsN does not.

* `runSIR()` gains `rxThreads`, controlling rxode2 OpenMP threads per worker during OFV evaluation. `nlmixr2utils` 0.3 requires it whenever `workers > 1`, so parallel SIR runs previously aborted on most multicore machines.

* `runSIR()` now applies `sigmaInflation` to residual-error parameters. They were classified as THETA, which made the argument silently unreachable.

* `runSIR()` now takes OMEGA uncertainty from `fit$cov` by default, via the new `omegaFallback = "cov"`. The initial proposal therefore carries the correlations between THETA and OMEGA, which the previous block-diagonal construction discarded even when they were available. `omegaFallback = "wishart"` forces the old approximation, and it is still used automatically when `fit$cov` does not carry OMEGA.

# nlmixr2sir 0.2

* `runSIR()` now writes a final canonical shared `raw_results.*` artifact, uses the common `nlmixr2utils` run-state and seeding helpers, and no longer emits per-iteration raw-results CSV files.

# nlmixr2sir 0.1

* Initial package split from `nlmixr2extra`, providing `runSIR()`,
  `sirSummary()`, S3 print/plot methods, tests, and the SIR vignette as a
  standalone package depending on `nlmixr2utils`.
