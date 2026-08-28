# Second-stage outer-iteration algorithm audit

## Status and authority

This document is the current decision and evidence authority for the outer
iteration implemented by `RunSecondStageLocalFitting`. The normative execution
description remains in [Second-stage local fitting](second-stage-local-fitting.md).
The four earlier audits are immutable historical records and are linked under
[Historical provenance](#historical-provenance).

The reviewed production baseline is `03cdb6ef` (`Unify
ConvergenceCertificate`). The audit consolidation and diagnostic cleanup do
not change `FitOptions`, command-line options, model persistence, convergence
thresholds, candidate selection, stop precedence, or the production
trajectory.

## Scope and canonical states

The review covers the second-stage outer loop from a validated accepted state
through proposal construction, candidate selection, post-processing,
convergence, final dependency polish, and persistence. The statistical
derivation of the local MDPDE estimator and the later group-fitting stage are
outside this audit.

Three states must remain distinct:

- `S(k)` is the previous validated accepted state.
- `F(S(k))` is the complete, undamped joint-offset-to-local-shape operator
  endpoint, including availability and solver evidence.
- `S(k+1)` is the candidate that is actually accepted after the geometric
  factor search, objective gates, joint polish, boundary reconciliation, and
  rescue.

Accepted movement is `T(S(k+1)) - T(S(k))`. The strict fixed-point residual is
`T(F(S(k))) - T(S(k))`. `T` uses three transformed coordinates: log peak,
log width, and shared physical offset normalized by peak.

Accepted movement samples active optimization DOFs. Every shape-active atom
contributes one log-peak and one log-width sample; every active
`(cluster, group_id)` shared offset contributes one group sample using its
largest absolute member change. Fixed and quarantined coordinates do not
dilute this population.

The operator residual instead samples the complete nominal-DOF population,
including fixed and quarantined shapes and shared offsets. Missing or
non-finite endpoint evidence makes the operator incomplete; it must never be
replaced by the previous state to manufacture a zero residual.

## End-to-end state machine

```text
validated S(k)
  -> complete undamped joint shared-offset endpoint
  -> complete undamped local-shape endpoint
  -> strict operator evidence F(S(k))
  -> geometric candidate factors: validity -> trust -> guard -> objective
  -> active-column joint polish
  -> boundary reconciliation and cooperative rescue
  -> complete-state global previous/best audit
  -> trust-radius and quarantine/probation transition
  -> assembled validated S(k+1)
  -> production convergence certificate
  -> stop policy selects a base final state
  -> final uncut dependency polish candidate
  -> converged path recertifies the persisted candidate
  -> persist Gaussian and peeling state
```

Validity establishes that a candidate can be represented. Trust limits the
step tested in the current iteration and updates the next radius. Guard tests
domain feasibility. Objective gates accept or reject candidates. None of
these responsibilities substitutes for fixed-point evidence.

## Authoritative production certificate

`ConvergenceCertificate::ProductionConverged()` is the only production stop
decision and requires all of the following:

```text
solver qualified
&& accepted active-DOF p99 < 1e-4
&& complete nominal-DOF operator
&& nominal fixed-point residual p99 < 1e-4
&& invariants clear
&& orthogonal blockers clear
```

The percentile predicate is coordinate-wise: the p99 for each of log peak,
log width, and shared offset must pass independently. Solver qualification
requires full, undamped, non-fallback active endpoints. Invariants include
complete shared-group activity and finite evidence. Orthogonal blockers cover
objective-domain changes, quarantine transitions, suspicious offset fallback,
and rejected clusters.

Maximum values and sparse-tail counts remain diagnostic measurements. The
historical `1e-3` maximum threshold is not part of
`ProductionConverged()` and does not define a separate production policy.

## Failure mode and safeguard coverage

| Failure mode | Accepted p99 | Strict operator p99 | Qualification | Invariants / blockers |
| --- | ---: | ---: | ---: | ---: |
| Trust clipping or objective backtracking makes the committed step small while the full endpoint remains material | Detects the small committed step | Blocks the false fixed point | Provides endpoint quality | Records the limiting state |
| Polish, reconciliation, or rescue moves the committed state after a small operator endpoint | Blocks convergence | Detects the small endpoint | Confirms the endpoint solve | Records post-processing blockers |
| Soft solver failure, damping, or fallback produces small numerical movement | Observes movement only | Observes residual only | Blocks convergence | Preserves failure classification |
| Fixed or quarantined coordinates hide an unavailable nominal endpoint | Excludes inactive DOFs by design | Fails closed on incomplete evidence | Reports restriction | Enforces population completeness |
| A shared group has mixed activity or non-finite member evidence | Cannot safely build a partial group sample | Cannot safely build a nominal group sample | Cannot establish a valid endpoint | Fails the invariant |
| Objective domain or quarantine changes during the iteration | May still be small | May still be small | May still pass | Orthogonal blocker prevents a premature stop |

No retained predicate is implied by the others. A zero-exposure corpus result
is empirical evidence, not a mathematical redundancy proof.

## Final dependency polish recertification

Final dependency polish is objective-accepted provisionally. On a `converged`
path, a changed polished state is persisted only when a new certificate built
at that state passes `StrictOperatorPassed()`: solver qualification, complete
nominal operator evidence, residual p99, and invariants must all pass. Failure,
incomplete evidence, or evaluation error retains the already converged base
state. Non-convergence stop reasons retain the existing objective-only polish
policy.

## Current diagnostic contract

The current Debug trajectory is schema 9 and serializes the production
certificate plus the measurements needed to explain it: active and nominal
populations, certificate bits, p99 and maximum values, unavailable and tail
counts, qualification, invariants, blockers, residual classification, and
candidate-path statistics. It has no historical policy vector, exposure flag,
accepted-only persistence record, or `production-maximum` decision.

Frozen-IRLS predicted-reduction and rho instrumentation is not part of a
normal or routine audit build. It is available only through the developer-only
`RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT` build option and never controls the
production trajectory. Its logs are consumed only by
`analyze_trust_model_experiment.py`; the production corpus analyzer does not
import or aggregate them.

## Evidence and retired experiments

The frozen baseline expands the checked manifest to 600 deterministic cases.
At `03cdb6ef`, all 600 cases completed with no safety regression, and the stop
distribution was:

| Stop reason | Cases |
| --- | ---: |
| `converged` | 42 |
| `audit-patience` | 372 |
| `all-rejected-backtracking-exhausted` | 163 |
| `maximum-iterations` | 23 |

The consolidation baseline records a manifest SHA-256 of
`2b0d74c249a4723575696d98f99df9d3c449a30837adf8117aaf744145d0a87a`,
case identity SHA-256 of
`63985f81c5b188cfc9992742cad1f4b9418cc36c88b12398e5775a324289bd98`,
and frozen truth SHA-256 of
`04e3cda3b49857b2d5e4f63e973b2392dfe3095360974db988770e7468edd628`.

Historical all-selected, active-proposal, cluster/maximum, and
production-maximum policies are retired. Accepted-only persistence is retained
only as a targeted negative unit scenario: small accepted movement cannot
declare convergence while the strict operator residual remains material.
Frozen-IRLS/rho remains a separate diagnostic experiment because the observed
coverage and action divergence did not justify a production controller.
Coarse-to-fine factor refinement remains rejected on objective and truth
outcomes. Fold-168 remains an optional quality regression, not convergence
policy evidence.

## Trajectory-neutral cleanup acceptance

Diagnostic cleanup is accepted only if paired runs complete 600/600 cases,
match the production semantic and normalized terminal-state digest for every
case, retain the stop distribution above, produce zero median and p90 deltas
for objective, transformed-truth RMSE, and accepted iterations, and report
zero safety regression. With the experimental flag disabled, a normal build
must not construct trust-model diagnostics. Routine audit elapsed-time median
and p90 must both decrease relative to the saved baseline.

### Cleanup recertification result (2026-08-28)

The paired run used AppleClang 21, `RelWithDebInfo`, the checked manifest and
truth, one estimator thread, sequential jobs, and an exported `03cdb6ef`
baseline executable. All blocking conditions passed:

| Gate | Result |
| --- | ---: |
| Baseline / candidate completed | `600/600` / `600/600` |
| Failed cases | `0` / `0` |
| Production semantic digest matches | `600/600` |
| Normalized terminal-state digest matches | `600/600` |
| Stop distribution match | exact |
| Objective delta median / p90 | `0 / 0` |
| Transformed-truth RMSE delta median / p90 | `0 / 0` |
| Accepted-iteration delta median / p90 | `0 / 0` |
| Safety regressions | `0` |
| Elapsed median, baseline → candidate | `0.101151 s → 0.088944 s` |
| Elapsed p90, baseline → candidate | `0.708533 s → 0.527565 s` |

Timing is excluded from both semantic digests. With the experiment flag off,
the trust-model data structures and calculations are not compiled. Schema 9
contains maximum and tail diagnostics but no rho, comparator, exposure, or
accepted-only persistence fields.

## Historical provenance

The following records preserve the original measurements and decisions in
chronological order:

1. [Convergence safeguard audit](audit-history/second-stage-convergence/convergence-safeguard-audit.md)
2. [Stationarity and active-coordinate population audit](audit-history/second-stage-convergence/stationarity-active-coordinate-audit.md)
3. [Counterfactual convergence continuation audit](audit-history/second-stage-convergence/counterfactual-convergence-continuation-audit.md)
4. [Convergence exposure and counterfactual outcome audit](audit-history/second-stage-convergence/convergence-exposure-counterfactual-outcome-audit.md)

The historical records are provenance, not current production specification.
