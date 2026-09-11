# Second-stage outer-iteration algorithm audit

## Status and authority

This document is the current decision and evidence authority for the outer
iteration implemented by `detail::RunSecondStageIterations`. The normative execution
description remains in [Second-stage local fitting](second-stage-local-fitting.md).
The earlier audits retain historical design reasoning and are linked under
[Historical provenance](#historical-provenance).

The current structural P0 implementation and its focused-test validation
are recorded in [Structural P0 result](#structural-p0-result-2026-09-11).

The earlier per-atom-offset review baseline is
`a4354e698e77398154009d231907ebf3ed4b1d52` (per-atom offsets). The audit consolidation and diagnostic cleanup do
not change `FitOptions`, command-line options, model persistence, convergence
thresholds, candidate selection, stop precedence, or the production
trajectory.

The [retained P1 changes](second-stage-p1-ablation.md) limit final polish to
converged stops with an independent radius; cooperative rescue and production
Grow retain their baseline policies.

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
log width, and per-atom physical offset normalized by peak.

Accepted movement samples active optimization DOFs. Every shape-active atom
contributes one log-peak and one log-width sample; every offset-active atom
contributes one absolute offset-to-peak-ratio change sample. Fixed and quarantined coordinates do not
dilute this population.

The operator residual instead samples the complete nominal-DOF population,
including fixed and quarantined shapes and per-atom offsets. Missing or
non-finite endpoint evidence makes the operator incomplete; it must never be
replaced by the previous state to manufacture a zero residual.

## Per-atom availability and conditioning

Each selected atom has independent shape and offset availability masks. Offset
availability is set only for a non-hard-failure joint-offset result whose model
passes validity checks. A soft failure may therefore supply an available endpoint
without being solver-qualified. Shape and offset availability are independent;
a missing shape sets both shape residual coordinates to infinity, and a missing
offset sets its own residual to infinity. `operator_complete` is the conjunction
of the masks for every nominal selected atom; transformed finiteness and solver
qualification are separate requirements. Inactive atoms are not removed from this
nominal population, and each coordinate has its own percentile population.

Joint-offset and joint-polish conditioning normalize each design-matrix column,
form the normalized Gram matrix, and use LDLT `min(D)/max(D)` as a conditioning
proxy before ridge. A ratio at or below `1e-8` triggers a ridge multiplier floor
of `10`. Empty/invalid columns, failed factorization, or nonpositive/nonfinite
pivots return the existing zero sentinel and require the guard. This ratio is
neither a singular-value condition number nor an effective-rank certificate.
Other existing ridge safeguards may also increase multipliers; diagnostics report
the actual minimum/maximum multipliers after all guards, not only the floor.

Schema-1 diagnostic records separately report conditioning, solve status, and
per-atom availability. Phases distinguish outer operator, candidate polish,
boundary reconciliation, final dependency polish, and final recertification.
Offset status codes follow `JointOffsetSolveStatus` (0 converged; 1 system build
failed; 2 empty; 3 initial solve failed; 4 IRLS solve failed; 5 objective
deteriorated; 6 iteration limit). Hard-failure classification is recorded
separately. Joint-polish solve records report `solved`/`failed`; this does not
imply objective acceptance or final persistence approval.

## End-to-end state machine

```text
validated S(k)
  -> complete undamped joint per-atom offset endpoint
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
  -> only converged and enabled: final uncut dependency polish at radius 1.0
  -> strict operator persistence safety check; otherwise retain chosen base
  -> persist Gaussian and peeling state
```

Validity establishes that a candidate can be represented. Trust limits the
step tested in the current iteration and updates the next radius. Guard tests
domain feasibility. Objective gates accept or reject candidates. None of
these responsibilities substitutes for fixed-point evidence.

Guard is feasibility-only: guard-only factor reduction does not request radius
shrink. Radius growth requires boundary utilization of at least `0.8` and actual
objective reduction larger than `1e-8 + 1e-3 * abs(previous)`. It does not use
predicted reduction or rho. Objective-induced shrink and radius bounds retain
their existing precedence and values.

## Authoritative production certificate

`ConvergenceCertificate::ProductionConverged()` is the only production stop
decision and requires all of the following:

```text
solver qualified
&& accepted active-DOF p99 < 1e-4
&& complete nominal-DOF operator
&& nominal fixed-point residual p99 < 1e-4
&& orthogonal blockers clear
```

The percentile predicate is coordinate-wise: the p99 for each of log peak,
log width, and per-atom offset must pass independently. Solver qualification
requires full, undamped, non-fallback active endpoints. Operator completeness and non-finite residuals fail
closed. Orthogonal blockers cover
objective-domain changes, quarantine transitions, suspicious offset fallback,
and rejected clusters.

Maximum values remain diagnostic measurements and do not define a separate
production policy.

## Failure mode and safeguard coverage

| Failure mode | Accepted p99 | Strict operator p99 | Qualification | Invariants / blockers |
| --- | ---: | ---: | ---: | ---: |
| Trust clipping or objective backtracking makes the committed step small while the full endpoint remains material | Detects the small committed step | Blocks the false fixed point | Provides endpoint quality | Records the limiting state |
| Polish, reconciliation, or rescue moves the committed state after a small operator endpoint | Blocks convergence | Detects the small endpoint | Confirms the endpoint solve | Records post-processing blockers |
| Soft solver failure, damping, or fallback produces small numerical movement | Observes movement only | Observes residual only | Blocks convergence | Preserves failure classification |
| Fixed or quarantined coordinates hide an unavailable nominal endpoint | Excludes inactive DOFs by design | Fails closed on incomplete evidence | Reports restriction | Enforces population completeness |
| One atom has unavailable or non-finite offset evidence | Includes its offset only when active | Unavailable evidence makes the operator incomplete; non-finite residual fails the percentile test | Availability alone does not qualify a solver | Each nominal atom retains its own coordinate |
| Objective domain or quarantine changes during the iteration | May still be small | May still be small | May still pass | Orthogonal blocker prevents a premature stop |

No retained predicate is implied by the others. Absence of observed failures
is not a mathematical redundancy proof.

## Final dependency polish recertification

Final dependency polish is objective-accepted provisionally. On a `converged`
path, a changed polished state is persisted only when a new certificate built
at that state passes `StrictOperatorPassed()`: solver qualification, complete
nominal operator evidence, and residual p99 must all pass. Failure,
incomplete evidence, or evaluation error retains the already converged base
state.

Non-convergence stop reasons persist the existing selected base state directly.
They do not run final dependency polish or operator recertification. The strict
candidate check is the only final-polish persistence policy.

## Current diagnostic contract

The current Debug trajectory is schema 10 and serializes the production
certificate plus its active and nominal populations, p99 and maximum values,
operator completeness, and four orthogonal blockers. Earlier trajectory
schemas are not accepted by the current analyzer.

Frozen-IRLS predicted-reduction and rho instrumentation is not part of a
normal or routine audit build. It is available only through the developer-only
`RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT` build option and never controls the
production trajectory. Its logs are consumed only by
`analyze_trust_model_experiment.py`. The experiment remains diagnostic-only.

## Structural P0 result (2026-09-11)

The structural baseline is `c5417717154896e53d31a1d43910cde15e3f80f6`.
No numerical policy or test expectation was relaxed.
See [P0 structural refactoring](second-stage-p0-structure.md) for ownership,
the flowchart and the per-path validation map.

The implementation extracts phase/trust-model observation, separates p99
certificate evidence from maximum/population diagnostics, centralizes typed
candidate evaluation and publishes builder-owned selection through a consuming
transaction. Existing tests only adapt to internal interface changes; no tests,
test cases or test targets were added.

The trust-enabled build passed 108/108 existing estimator defense tests,
including quiet/trace neutrality. Repository lint and whitespace checks passed.
Fold-168 was not run. Boundary/correction and fallback coverage relies on
existing focused tests; those tests do not establish coverage of every branch.

Historical all-selected, active-proposal, cluster/maximum, and
production-maximum policies are retired. Accepted-only persistence remains a
targeted negative unit scenario: small accepted movement cannot declare
convergence while the strict operator residual remains material. Geometric
factor search remains the production policy; coarse-to-fine refinement is not
part of that search. Fold-168 remains an optional quality regression.

## Historical provenance

The following records preserve independent design reasoning and evidence in
chronological order:

1. [Convergence safeguard audit](audit-history/second-stage-convergence/convergence-safeguard-audit.md)
2. [Stationarity and active-coordinate population audit](audit-history/second-stage-convergence/stationarity-active-coordinate-audit.md)
3. [Counterfactual convergence continuation audit](audit-history/second-stage-convergence/counterfactual-convergence-continuation-audit.md)

The historical records are provenance, not current production specification.
