# nlmixr2sir 0.2

* `runSIR()` now writes a final canonical shared `raw_results.*` artifact, uses the common `nlmixr2utils` run-state and seeding helpers, and no longer emits per-iteration raw-results CSV files.

# nlmixr2sir 0.1

* Initial package split from `nlmixr2extra`, providing `runSIR()`,
  `sirSummary()`, S3 print/plot methods, tests, and the SIR vignette as a
  standalone package depending on `nlmixr2utils`.
