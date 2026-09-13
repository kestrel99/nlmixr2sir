# nlmixr2sir (development version)

* SIR now derives its parameter vector from a single internal description of the fit instead of re-deriving THETA, sigma, and OMEGA names independently in each function. This fixes `runSIR()` aborting with `Assertion on 'mu' failed: Contains missing values` against nlmixr2est 7, where `fit$cov` reports OMEGA alongside THETA and the old proposal mean came back `NA` for every OMEGA element.

* OFV evaluation no longer fails in a session where `nlmixr2sir` was attached but `nlmixr2` was not. `rxode2::ini()` evaluates the OMEGA line it is given as a `lotri({...})` call in the caller's environment, and `lotri` is an Imports rather than a Depends of both rxode2 and nlmixr2est, so it was never visible; `nlmixr2sir` now imports it. Every sample previously returned `NA` and `runSIR()` aborted with `All SIR OFV evaluations failed`.

* That abort now reports the first underlying evaluation error, instead of discarding it. A configuration problem fails every sample identically, which the old message hid.

* The package now requires nlmixr2est >= 7.0.0, nlmixr2utils >= 0.3, and rxode2 >= 5.0.0.

* `runSIR()` now takes its run settings through `control = runSIRControl()` rather than as flat arguments, cutting its signature from 21 arguments to 6. Passing a setting directly to `runSIR()` is an error naming `runSIRControl()`.

* `runSIR()` gains `rxThreads`, controlling rxode2 OpenMP threads per worker during OFV evaluation. `nlmixr2utils` 0.3 requires it whenever `workers > 1`, so parallel SIR runs previously aborted on most multicore machines.

* `runSIR()` now applies `sigmaInflation` to residual-error parameters. They were classified as THETA, which made the argument silently unreachable.

* `runSIR()` now takes OMEGA uncertainty from `fit$cov` by default, via the new `omegaFallback = "cov"`. The initial proposal therefore carries the correlations between THETA and OMEGA, which the previous block-diagonal construction discarded even when they were available. `omegaFallback = "wishart"` forces the old approximation, and it is still used automatically when `fit$cov` does not carry OMEGA.

# nlmixr2sir 0.2

* `runSIR()` now writes a final canonical shared `raw_results.*` artifact, uses the common `nlmixr2utils` run-state and seeding helpers, and no longer emits per-iteration raw-results CSV files.

# nlmixr2sir 0.1

* Initial package split from `nlmixr2extra`, providing `runSIR()`,
  `sirSummary()`, S3 print/plot methods, tests, and the SIR vignette as a
  standalone package depending on `nlmixr2utils`.
