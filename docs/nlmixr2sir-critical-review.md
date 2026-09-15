# Follow-up critical review of `nlmixr2sir` 0.3

- **Review date:** 2026-09-14
- **Package source reviewed:** branch `main`, commit `612980e`, including the
  current uncommitted remediation work
- **Primary focus:** statistical correctness, numerical robustness, recovery
  safety, reproducibility, and the accuracy of
  [`sir-technical-reference.md`](sir-technical-reference.md)
- **Comparator:** local PsN 5.7.1 source at
  `C:/Users/justin/Documents/GitHub/PsN`
- **Comparator policy:** PsN is a useful reference implementation and source of
  numerical oracles, but PsN parity is **not** a requirement

## Executive summary

**Verdict: request changes before an inferential release.**

The remediation is substantial and several of the original high-severity
findings are now correctly resolved. In particular, the proposal-density
Cholesky solve is fixed, the Box-Cox Jacobian is included, recentering carries
the reference OFV, reference rows are explicitly excluded from empirical
summaries, OMEGA SD-scale RSE reporting is corrected, and useful
importance-weight diagnostics have been added.

This follow-up review does **not** require nlmixr2sir to reproduce PsN when
nlmixr2sir has a better justified design. A difference from PsN is a finding
only when it independently causes one of the following:

- the implemented statistical target differs from the documented target;
- results depend improperly on arbitrary parameter units;
- candidate OFVs may be evaluated on a different likelihood surface;
- saved output can be confused with another run;
- unrelated files can become eligible for recursive deletion;
- numerical failures are misclassified or allowed into evaluation; or
- documentation overstates what the code establishes.

Under that standard, four blocking issues remain:

1. empirical-rank and positive-definiteness decisions depend on the relative
   units of different parameters;
2. the objective preflight tests only one point and the candidate evaluator
   does not preserve the likelihood-relevant settings of the fitted model;
3. the default recovery path can claim an unrelated non-empty directory as an
   nlmixr2sir-owned directory; and
4. the recovery fingerprint omits the effective proposal contents and can fail
   open when identity cannot be verified.

The complete non-CRAN test suite passed during this review with zero failures
and 14 warnings. The remaining issues are counterexamples or state transitions
not represented in the current tests.

## Review standard and the role of PsN

The original report gave behavioral parity too much weight. This revision uses
the following hierarchy instead:

1. **Statistical correctness:** the implementation must sample the target it
   documents.
2. **Numerical robustness:** equivalent problems expressed in different units
   should not receive contradictory rank or regularization decisions.
3. **Scientific reproducibility:** recovered state must be demonstrably the
   state requested by the current call.
4. **Operational safety:** output handling must not put unrelated files at
   risk.
5. **Transparent design:** deliberate differences from PsN are acceptable when
   their consequences are documented and tested.

PsN remains valuable in two narrower roles:

- as an independently developed implementation against which shared numerical
  primitives can be checked; and
- as a source of workflow ideas and edge cases.

It is not the normative specification for nlmixr2sir. Exact agreement with
PsN's Box-Cox shift, lambda search, draw budget, rejection recovery, output
format, or NONMEM execution behavior is not required.

## Status of the original findings

| Original finding | Follow-up status | Assessment |
| --- | --- | --- |
| B1: incorrect Cholesky solve | **Resolved** | The upper-triangular Cholesky factor is now solved with the required transpose and exact oracles cover correlated proposals. |
| B2: objective semantics not preserved | **Partially resolved; blocking** | Unsupported estimation methods fail and the center is checked, but equality at one point does not establish equality of the candidate surface. |
| B3: Box-Cox Jacobian omitted | **Resolved** | The Jacobian is included consistently. This is a justified improvement over PsN for the documented original-scale target. |
| B4: unsafe/unfingerprinted recovery | **Partially resolved; blocking** | State versioning, manifests, and fingerprints exist, but directory claiming and incomplete proposal identity leave unsafe paths. |
| R1: Box-Cox differs from PsN | **Reframed** | Difference itself is not a defect. Only center-domain and non-finite-output failure modes remain findings. |
| R2: no-covariance fallback contract | **Resolved by contract** | The documentation now requires an alternative proposal source when `fit$cov` is absent. Automatic fallback is not required. |
| R3: reference OFV does not move | **Resolved** | Center and reference OFV now move and persist together. |
| R4: synthetic center contaminates summaries | **Resolved** | Explicit row roles are used and proposal consumers exclude the reference row. |
| R5: rank-deficient proposals are repaired | **Partially resolved; blocking replacement finding** | True rank deficiency now aborts, but the numerical-rank test itself is unit-dependent. |
| R6: off-diagonal `rse_sd_scale` | **Resolved** | The conversion is restricted to OMEGA diagonals. |
| S1: weight degeneracy diagnostics | **Resolved** | ESS and weight-concentration information are available. |
| S2: noise band uses different resampler | **Partially resolved** | It now uses `sirResample()` and the actual cap, but consumes an overflow-prone raw weight column. |
| S3: end-to-end PsN parity test | **No longer a requirement** | Internal mathematical invariants and package-specific end-to-end tests are the appropriate acceptance standard. |
| S4: incomplete provenance | **Partially resolved** | Much more is stored, but the effective proposal, cumulative schedule, and initial repair are incomplete. |
| S5: overstrong inferential language | **Largely resolved** | One important parameterization-invariance claim remains incorrect. |

## Findings at a glance

| ID | Severity | Finding | Primary consequence |
| --- | --- | --- | --- |
| B1 | Blocking | Rank and PD decisions use raw-coordinate eigenvalue ratios | Equivalent full-rank samples can be accepted or rejected solely because of parameter units |
| B2 | Blocking | Objective preflight establishes agreement only at the fitted center | Candidate likelihood ratios can still mix different objective surfaces |
| B3 | Blocking | Recovery can claim an unrelated non-empty directory | Unrelated files can later become eligible for recursive deletion |
| B4 | Blocking | Fingerprint omits effective proposal contents and fails open on unverifiable fields | A changed proposal or unidentified model/data can reuse stale state |
| R1 | Required | `addIterations` replaces rather than extends schedule provenance | A multi-iteration result can be returned under a shorter schedule identity |
| R2 | Required | Diagnostic noise bands use overflow-prone raw importance ratios | The highest-weight candidate can be discarded from the diagnostic |
| R3 | Required | Technical-reference claims do not match the implementation or estimand | Users can misunderstand rank repair, warnings, and parameterization dependence |
| R4 | Required | Pre-run feasibility and initial covariance-repair provenance are incomplete | Expensive doomed runs and incomplete audit trails remain possible |
| R5 | Required | Box-Cox center-domain and non-finite inverse handling are incomplete | A later center can be outside the transform domain and infinite candidates can be misclassified |

## Blocking findings

### B1. Numerical rank and covariance repair depend on parameter units

#### Evidence

[`.sirCheckProposalRank()`](../R/sir-utils.R) computes the empirical covariance
in the original parameter coordinates and defines numerical rank as the number
of eigenvalues greater than `rankTol * max(eigenvalue)`. This tests the
condition of a dimensional covariance matrix, not whether the centered sample
matrix contains independent directions.

The current tests multiply the complete covariance by one common scalar. That
proves invariance to a global change of scale, but pharmacometric parameters do
not generally share one unit. A clearance-like parameter, a log parameter, a
small covariance, and a residual standard deviation can differ by many orders
of magnitude.

The following full-rank sample was rejected during this review:

```r
set.seed(11)
x <- cbind(
  a = rnorm(100),
  b = 1e-6 * rnorm(100)
)
.sirCheckProposalRank(x)
```

The error reported numerical rank 1 for two parameters. Multiplying only the
second column by a unit-conversion constant therefore changes the algorithm's
conclusion from full rank to deficient rank.

[`.sirEnsurePosDef()`](../R/sir-utils.R) has the same coordinate problem at a
smaller tolerance. For example, `diag(c(1, 1e-14))` is changed to approximately
`diag(c(1, 1e-12))`, inflating the second marginal variance 100-fold because a
different coordinate happens to have variance 1.

#### Consequence

This can reject a statistically supported empirical proposal or alter a valid
initial proposal solely because the model uses different parameter units. It
is independent of PsN behavior and is blocking for a general-purpose
inferential implementation.

#### Required change

1. Center and standardize every non-constant column before determining rank,
   or determine rank from the correlation matrix/SVD of the standardized
   centered sample.
2. Treat an exactly constant column as unsupported and report its parameter
   name.
3. If numerical PD repair remains necessary, repair the correlation matrix and
   map it back with the original marginal standard deviations. Do not inflate
   marginal variances based on another parameter's units.
4. Add tests that independently rescale columns and require invariant
   rank/repair decisions under `X %*% D` and `D %*% Sigma %*% D`.

### B2. The objective preflight does not establish one fixed target surface

#### Evidence

[`.sirCheckObjective()`](../R/sir-preflight.R) is a useful fail-fast check. It
rejects unsupported estimation methods and reevaluates the fit at the stored
center. However, [`sirEvalOFV()`](../R/sir-eval.R) still constructs a new
`foceiControl()` with selected defaults rather than reconstructing the
likelihood-relevant controls from `fit$control`.

Potentially relevant fitted settings include, among others, interaction/FOCE
type, censoring behavior, likelihood adjustment, event behavior, ODE solver
controls and tolerances, and bad-solve objective handling. Optimization-only
settings need not be copied when `maxOuterIterations = 0`, but
likelihood-defining and likelihood-evaluation settings do.

The preflight evaluates one vector. Two objective implementations can agree at
that vector within tolerance while differing elsewhere. Importance sampling
requires all candidate likelihoods and the reference likelihood to belong to
one fixed surface, not merely to intersect near the fitted estimates.

The current acceptance rule also passes when either the absolute or relative
difference is below `objfTolerance`. Relative error in the absolute OFV is not
the relevant inferential scale: the OFV contains additive constants and grows
with the number of observations, while SIR weights depend on differences in
OFV.

#### Consequence

A run can pass preflight but compute candidate dOFVs using settings different
from those that produced `fit$objf`. That directly changes likelihood ratios
and importance weights.

#### Required change

1. Define which `fit$control` fields affect objective evaluation and carry
   those fields into the fixed-parameter evaluator.
2. Explicitly reject fitted configurations whose likelihood semantics cannot
   be reconstructed.
3. Evaluate a deterministic stencil: the center plus small, bounds-aware
   perturbations along each parameter direction. This is still not a formal
   proof, but it directly checks the local surface rather than one point.
4. Use an absolute tolerance on the dOFV scale, with separate numerical
   tolerances only where justified.
5. Attach and persist the preflight inputs and results.

### B3. Default recovery can claim an unrelated directory

#### Evidence

[`runSIR()`](../R/sir-run.R) calls `.sirAssertSafeToClear()` only when
`resolveRunDir()` returns overwrite mode. For an explicitly supplied existing
directory and the default `recover = TRUE`, `resolveRunDir()` returns resume
mode. If that directory contains no saved SIR state, `saved_state` is `NULL`,
but `.sirWriteManifest()` is still called.

This transition was reproduced with a temporary non-empty directory containing
`sentinel.txt`:

```text
resolve mode:           resume
owned before manifest:  FALSE
owned after manifest:   TRUE
sentinel still present: TRUE
```

[`.sirDirIsOwned()`](../R/sir-provenance.R) checks only whether a file named
`sir_manifest.dcf` exists. It does not parse the manifest or verify the
package, prefix, state version, or relationship to a saved state.

On a later fresh run, the existence of the newly written marker can authorize
recursive deletion of the same directory, including the unrelated sentinel
and any other pre-existing files.

#### Consequence

This is an operational data-loss path. PsN behavior is irrelevant; the package
must establish ownership before writing an ownership marker or deleting
contents.

#### Required change

1. Record whether the output directory existed and whether it was empty before
   any seed, state, or manifest write.
2. For recovery, require a valid manifest **and** compatible state before using
   a pre-existing non-empty directory.
3. Refuse to claim an existing non-empty unowned directory.
4. Parse and validate marker contents; filename existence alone is not proof of
   ownership.
5. Make manifest-write failure fatal because the marker participates in the
   deletion safety policy.
6. Test the two-call sequence: failed/empty recovery attempt followed by a
   fresh run, asserting that an unrelated sentinel can never be deleted.

### B4. Recovery fingerprints do not identify the effective proposal

#### Evidence

[`.sirRunFingerprint()`](../R/sir-provenance.R) includes model text, data,
parameter names, fitted estimates, objective value, estimation method,
schedule, and selected controls. It does not include:

- the contents of `fit$cov` used by the default proposal path;
- the validated numeric contents of `covmatInput`;
- the selected and validated contents of `rawresInput`;
- the resolved, inflated initial proposal covariance;
- parameter kinds, bounds, and fixed/free schema beyond names; or
- the initial Box-Cox/proposal state when recovery begins after an update.

For path-based input, the control digest stores the path string. During this
review, overwriting a covariance file at the same path produced an identical
control digest.

The comparison function also treats an `NA` digest as "unverifiable" and skips
the comparison. That is a fail-open policy: recovery can proceed precisely when
the package cannot establish that the model or data are the same.

Package versions are recorded but not enforced. A blanket requirement that all
dependency versions match would be unnecessarily strict, but the current
comment that the one-point objective preflight catches any meaningful version
effect is too strong. A version change can alter proposal construction,
Box-Cox estimation, random-number behavior, or result structure without
changing the center OFV.

#### Consequence

Changed scientific inputs can reuse a result or intermediate state generated
from a different proposal. This is a reproducibility and result-integrity
failure.

#### Required change

1. Fingerprint canonical numeric proposal inputs after parsing and validation,
   not only their paths.
2. Include the effective initial proposal covariance and parameter schema.
3. Make recovery fail closed when required identity fields cannot be digested.
   A new non-recovery run may still proceed with a warning if appropriate.
4. Enforce the nlmixr2sir state/algorithm version. Define and document a more
   selective compatibility policy for dependency versions.
5. Add recovery tests that mutate `fit$cov`, covariance-file contents,
   raw-results contents, bounds, and model/data identity one at a time.

## Required changes

### R1. `addIterations` does not maintain cumulative schedule identity

When `addIterations = TRUE`, the previous iterations and summaries are
retained, but the result's `schedule` attribute and saved fingerprint are built
only from the extension schedule supplied to the new call.

For example, after a two-iteration run is extended by one iteration, the saved
state can report three completed iterations while its fingerprint and result
schedule describe only the one-iteration extension. A later recovery request
with that one-iteration schedule can match and return the three-iteration
result immediately.

Store a cumulative schedule with actual iteration numbers and compute the saved
fingerprint from that cumulative schedule. Add a test that extends a run and
then recovers it, checking the schedule attribute, manifest, fingerprint,
completed count, and iteration summary together.

### R2. The convergence-noise calculation can discard the dominant weight

[`.sirDofvNoise()`](../R/sir-convergence.R) filters and renormalizes
`importance_ratio`. [`sirCalcWeights()`](../R/sir-weights.R) correctly computes
normalized probabilities in log space, but it also exposes
`importance_ratio = exp(log_ir)`, which can overflow.

The following case was reproduced:

```text
importance_ratio = 1, Inf
prob_resample     = 2.225074e-308, 1
```

The actual resampler correctly selects the second candidate. The diagnostic
filters it out because the raw ratio is not finite, so its noise band can
describe a different weighted population even though it now calls the correct
resampling function.

Use the already normalized `prob_resample`/`probability_resample` values
directly. Add an overflow test that requires the diagnostic to retain the row
with probability one.

### R3. The technical reference still contains material inaccuracies

#### Positive-definiteness floor

The reference says `.sirEnsurePosDef()` floors eigenvalues at
`sqrt(.Machine$double.eps)`. The code now uses `relTol * max(eigenvalue)`. Once
the unit-dependence in B1 is fixed, document the final standardized-coordinate
algorithm rather than either obsolete formula.

#### Narrow-proposal warning

The README and technical reference say `runSIR()` warns when the first proposal
is too narrow. The warning is currently emitted by the convergence plotting
path. Either run the check automatically before returning or state that
`plot(type = "convergence")` performs the check.

#### Parameterization invariance

The documents claim that including the Box-Cox Jacobian makes the retained
distribution invariant when the **model** is re-expressed in another smooth
parameterization. That is not correct for normalized likelihood under a flat
measure on the current nlmixr parameter scale.

If `x` is the nlmixr parameter and `z = T(x)` is only a proposal coordinate,
the induced proposal density is

$$
q_x(x) = q_z(T(x))\left|T'(x)\right|,
$$

and dividing by this density correctly targets `L(x) dx`. This makes the
answer invariant to the internal proposal transformation.

If the model itself is instead rewritten in parameter `y = h(x)` and the
estimand is again defined using a flat measure `dy`, mapping that target back
to `x` introduces `|h'(x)|`. Flat normalized likelihood is therefore not
invariant to arbitrary model reparameterization.

The documentation should say:

> The Jacobian makes Box-Cox an internal proposal transformation while
> preserving the normalized-likelihood target on nlmixr2's original parameter
> scale.

It should not claim invariance to redefining the model's parameter scale.

#### Comparability wording

Replace claims that every deliberate PsN difference is enumerated with a more
modest statement: PsN is used for selected workflow comparisons and numerical
oracles, while nlmixr2sir defines and tests its own statistical contract.

### R4. Feasibility checks and covariance-repair provenance are incomplete

`runSIR()` checks `nResample > nParameters` before objective evaluation. It
does not fully check whether the requested and attainable number of unique
candidate vectors can support the update:

- `nSamples` may itself be less than or equal to the parameter dimension;
- `capResampling` limits the number of distinct retained candidates;
- failed evaluations can reduce the usable candidate set; and
- dynamic resample reduction can lower the retained count to a value that
  cannot support full rank.

Validate deterministic impossibilities before expensive evaluation. After
evaluation, abort with an explicit feasibility message before attempting
resampling or covariance calculation when fewer than `p + 1` usable distinct
candidates remain.

The initial proposal path computes whether positive-definiteness repair was
applied, but iteration provenance records only repair of later empirical
updates. Store initial repair status, method, threshold, and magnitude
separately from per-iteration updates.

### R5. Box-Cox failure handling needs two independent robustness fixes

Exact replication of PsN's shift rule or lambda grid is not required. Two
current behaviors are nevertheless unsafe on their own merits.

First, the shift is chosen from the retained sample without guaranteeing that
the next proposal center lies inside the resulting Box-Cox domain. A recentered
best candidate need not be among the randomly retained rows. If it lies below
their minimum, the subsequent `.sirBcTransformMu()` call can abort because
`mu + delta <= 0`. Choose or adjust the shift using both the retained sample and
the center that will be transformed.

Second, inverse transformation can overflow to `Inf` without raising an R
error. `.sirSampleFullProposal()` classifies inverse failures using `is.na()`
rather than `!is.finite()`. An infinite unbounded parameter can therefore pass
the inverse-failure filter and reach bounds, OMEGA, or objective evaluation
logic under the wrong rejection category. Reject any row containing a
non-finite transformed value before all parameter-specific checks.

## Acceptable deliberate differences from PsN

The following are **not findings merely because PsN behaves differently**:

### Box-Cox Jacobian

Including the change-of-variables Jacobian is mathematically appropriate for
the package's stated target, normalized likelihood with respect to the current
nlmixr parameter scale. The analytic and Gamma-target tests are good evidence
for this choice. PsN's omission does not need to be copied.

### Box-Cox fitting heuristics

nlmixr2sir may use `stats::optimize()`, its own near-identity policy, and its own
shift convention. These are proposal-efficiency choices provided that:

- the transform is one-to-one over every sampled and centered value;
- the correct induced proposal density is used;
- all outputs are finite or rejected cleanly; and
- the choices are documented and tested.

### Draw-attempt budget and rejection strategy

The package does not need PsN's `2000 * nSamples` budget or PsN's
OMEGA/SIGMA-block adjustment. A smaller bounded budget with transparent
warnings can be preferable for interactive R use. The package must report
attempted and successful counts, avoid proceeding into an impossible
full-rank update, and make clear when fewer samples than requested were used.

### RSE conventions and summaries

Using percentage RSE and restricting SD-scale conversion to genuine variance
parameters are reasonable nlmixr2-native decisions. PsN's conventions need not
be reproduced.

### Execution and output differences

NONMEM-specific batching, copying, `mceta`, version selection, and file layouts
are not applicable requirements. nlmixr2-native persistence and diagnostics
should be judged by their own correctness and usability.

### End-to-end PsN equality

An exact paired PsN/NONMEM result is not an acceptance criterion. Estimation
engines, parameterizations, numerical solvers, and random-number streams make
exact equality neither necessary nor generally meaningful. More useful tests
are:

- analytic density and Jacobian invariants;
- exact oracles for genuinely shared mathematical primitives;
- deterministic recovery/state-transition tests;
- coordinate-rescaling invariance tests; and
- simulation-based calibration against targets with known moments or
  quantiles.

## Confirmed improvements

The remediation includes several strong changes worth preserving:

- [`sirCalcWeights()`](../R/sir-weights.R) now applies the correct triangular
  solve for an upper Cholesky factor.
- Exact correlated-density oracles cover the former Cholesky defect.
- Box-Cox weights include the forward-transform Jacobian and use stable
  log-space normalization.
- Recentering moves and persists both the proposal center and its reference
  OFV.
- Raw results distinguish proposal and reference rows explicitly, and
  empirical summaries use proposal rows only.
- `rse_sd_scale` is limited to OMEGA diagonal variances.
- ESS, maximum-weight, entropy/perplexity, and related diagnostics make severe
  importance-weight concentration visible.
- The diagnostic noise simulation now uses the actual limited-replacement
  resampler and effective cap; only its choice of weight column remains wrong.
- Missing fit covariance now has a clear contract: an alternative proposal
  source is required.
- State versioning, manifests, saved controls, schedules, proposal source, and
  reference-OFV history form a much better provenance foundation than the
  original implementation.

## Verification performed

The following checks were performed against the remediated tree:

1. `NOT_CRAN=true devtools::test(reporter = "summary")` completed with zero
   failures and 14 warnings.
2. The correlated Cholesky proposal-density implementation and its new oracle
   tests were inspected.
3. The Box-Cox Jacobian was checked against its analytic formula and integration
   path through iteration weighting.
4. A full-rank two-column sample with a `1e-6` unit ratio was rejected as rank
   one, reproducing B1.
5. A covariance with diagonal values `1` and `1e-14` was repaired to roughly
   `1` and `1e-12`, demonstrating coordinate-dependent marginal inflation.
6. A non-empty unowned directory entered resume mode and became "owned" after
   `.sirWriteManifest()`, reproducing B3 without deleting the sentinel.
7. Changing the contents of a path-based covariance input left the statistical
   control digest unchanged, reproducing part of B4.
8. An extreme but finite log-weight example produced raw importance ratios
   `1, Inf` while normalized probabilities remained finite and approximately
   `0, 1`, reproducing R2.
9. Recovery, `addIterations`, proposal-rank, provenance, Box-Cox, weight, and
   end-to-end tests were inspected for coverage of the counterexamples above.

The test warnings mostly arise from intentionally tiny fixtures: requested
resamples are clamped, effective sample sizes can be near one, and maximum
weights can approach 100%. These warnings are useful signals but make it harder
to distinguish expected test behavior from newly introduced warnings.

## Test improvements required with the fixes

Add focused tests for the following cases:

1. independently rescale each parameter column and require invariant rank and
   PD-repair decisions;
2. distinguish an exactly constant column from a small but independent column;
3. recover into a pre-existing non-empty unowned directory and verify no marker
   or state is written;
4. place a malformed or foreign `sir_manifest.dcf` in a directory and verify it
   does not establish ownership;
5. mutate covariance-file and raw-results-file contents without changing their
   paths and require recovery rejection;
6. mutate the fitted covariance and parameter bounds and require recovery
   rejection;
7. extend a run, then recover it and require a cumulative schedule everywhere;
8. feed `.sirDofvNoise()` a finite normalized probability paired with an
   infinite raw importance ratio;
9. choose a next center below the retained-sample minimum and require a valid
   Box-Cox domain; and
10. make Box-Cox inversion overflow and require an `inverseRejected` result
    before objective evaluation.

For new warning and error tests, snapshot the complete condition output so
changes to actionable diagnostics are reviewed deliberately.

## Technical-reference revision checklist

Before release, update [`sir-technical-reference.md`](sir-technical-reference.md)
to:

- define PsN as a comparator rather than a specification;
- state that the Jacobian preserves the target across internal proposal
  transformations, not arbitrary model reparameterizations;
- describe the final coordinate-standardized rank and PD policy;
- state exactly when the narrow-proposal warning is evaluated;
- describe recovery identity, directory ownership, and fail-closed behavior;
- explain cumulative schedules for added iterations;
- list the actual draw-attempt policy as an nlmixr2sir design choice without
  treating a difference from PsN as a defect;
- distinguish requested, attempted, successful, usable, and retained counts;
- describe initial and per-iteration covariance repair provenance; and
- avoid claims that every PsN difference has been exhaustively enumerated.

## Limitations of this follow-up review

- No simulation study was performed to establish frequentist interval coverage
  over a family of nonlinear mixed-effects models.
- No broad benchmark of proposal efficiency or parallel scaling was performed.
- The objective-surface concern is based on the evaluator construction and the
  insufficiency of a one-point equality test; a catalog of which nlmixr2 FOCEi
  controls produce material candidate-dependent differences remains to be
  built.
- The package is under active development, so source-function links are more
  stable than exact line numbers.
- Exact PsN parity was deliberately not used as an acceptance standard.

## Prioritized remediation plan

### Phase 0: correctness and safety

1. Make rank assessment and PD repair invariant to independent parameter-unit
   changes.
2. Prevent recovery from claiming a pre-existing unowned directory and validate
   manifest contents before deletion.
3. Fingerprint the effective proposal and fail closed on unverifiable recovery
   identity.
4. Preserve likelihood-relevant fitted controls and check the objective over a
   local parameter stencil.

No inferential release should be made before Phase 0 is complete.

### Phase 1: state and diagnostic integrity

1. Make `addIterations` provenance cumulative.
2. Use normalized probabilities in the convergence-noise calculation.
3. Add pre-run and post-evaluation feasibility checks.
4. Record initial covariance repair separately from later proposal updates.
5. Reject all non-finite inverse Box-Cox results and protect the next center's
   transformation domain.

### Phase 2: documentation and evidence

1. Correct the technical-reference statements listed above.
2. Add adversarial invariance and state-transition tests.
3. Reduce unexplained warnings in ordinary end-to-end fixtures and snapshot
   warnings that are intentional.
4. Keep PsN-derived oracles where both packages implement the same mathematical
   primitive, without implying that PsN defines the package's full behavior.

## Acceptance criteria for an inferential release

The package is ready for renewed review when all of the following are true:

- rank and PD decisions are unchanged by independent conversions of parameter
  units;
- candidate OFVs are evaluated with the likelihood semantics of the input fit
  and checked over more than the fitted center;
- no existing non-empty directory can become owned without a valid prior
  nlmixr2sir manifest and compatible state;
- recovery refuses changed or unverifiable model, data, proposal, parameter
  schema, schedule, or statistical controls;
- added iterations preserve a cumulative schedule and identity;
- diagnostic noise simulations use the same finite normalized probabilities as
  the actual resampler;
- impossible full-rank updates fail before avoidable expensive work;
- initial and later covariance repairs are fully recorded;
- Box-Cox centers and sampled vectors remain in a finite valid transform domain;
  and
- the README and technical reference describe the implemented estimand and
  behavior without claiming either unnecessary PsN equivalence or unsupported
  parameterization invariance.

The package need not reproduce PsN to satisfy these criteria. It needs to be
mathematically coherent, numerically invariant where it should be, safe to run,
and explicit about the target it computes.
