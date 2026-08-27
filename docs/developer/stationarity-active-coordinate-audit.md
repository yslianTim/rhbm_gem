# Stationarity semantics and active-coordinate population audit

> Current baseline (2026-08-27): active shared-DOF p99 is the production
> convergence population, maximum change is diagnostic-only, and solver
> qualification remains a counterfactual comparator.

## Current resolution

Production samples shape-active atoms once for log peak and log width and
samples each active `(cluster, group_id)` offset DOF once, using the maximum
absolute member change. Fixed and quarantined blocks are excluded. Mixed
shared-group activity and non-finite member changes fail convergence. This
removes fixed-zero dilution and group-size weighting without changing the p99
threshold.

Production convergence qualification remains
`ClusterHealth::production_convergence_qualified`. Solver qualification is
retained in Debug/audit builds as `solver-qualified`; it is
not promoted without continuation evidence of material benefit.

## Purpose and constraints

This is the second-round convergence safeguard audit for
`RunSecondStageLocalFitting`. It follows the
[first-round audit](convergence-safeguard-audit.md) and compares two shadow
definitions against the production convergence gate:

- solver qualification for truly active parameter blocks;
- active-member and independent shared-offset-DOF change populations.

The audit is observation-only. It does not change the production conjunction,
thresholds, stop reason, trajectory, or final model. Debug-disabled runs do not
build the second-round summaries.

At the time of this audit, the implementation still computed production p99
and maximum from every selected atom even though the algorithm description
specified active-coordinate semantics. That discrepancy is now resolved.

## Qualification semantics

Activity and qualification are separate axes:

| Activity | Meaning | Included in solver-qualified conjunction |
| --- | --- | --- |
| Active | The block remains eligible to move this attempt. | Yes |
| Fixed | A current failure or fallback removed the block from the active parameterization. | No; reported as restricted state |
| Quarantined | A stage-local quarantine mask fixed the block. | No; reported as restricted state |

| Solver qualification | Shape atom | Shared-offset DOF |
| --- | --- | --- |
| Solver qualified | Local status `SUCCESS`, full guard factor, no fallback | Joint status `Converged`, full group factor, no fallback |
| Soft unqualified | Finite usable endpoint with a non-success status or damping | IRLS objective deterioration, maximum iterations, or damping |
| Hard failure | Missing/invalid refit evidence for an active block | System build, empty system, initial solve, or IRLS solve failure |

Proposal usability is deliberately not equivalent to solver qualification.
Soft endpoints retain their existing update and objective-gate behavior but
fail the `solver-qualified` comparator. Fixed and quarantined blocks are
excluded from its active conjunction. An empty active set passes the solver
qualification conjunction vacuously and remains explicitly labelled
`all-fixed/restricted`; this label does not add a production stop condition.

The production predicate remains
`ClusterHealth::production_convergence_qualified`. The Debug snapshot also
records solver qualification, shape/offset qualification counts, restricted
state, and `production=true, solver=false` disagreements.

## Historical population comparison

Accepted and raw transformed changes are summarized three ways:

| Population | log peak / log width | offset / peak |
| --- | --- | --- |
| Historical production | Every selected atom | Every selected atom |
| Active member | Shape-active atoms | Offset-active atoms |
| Active DOF | Shape-active atoms | One sample per `(cluster, group_id)` offset column |

The shared-DOF offset sample is the maximum absolute member
`delta(offset/peak)`. This prevents group size from becoming an implicit
statistical weight while retaining the most sensitive weak-peak member. A
non-finite member makes the group sample non-finite. Mixed active/fixed members
inside one shared-offset group are an invariant violation: the DOF summary is
forced to fail instead of choosing a partial group silently.

The historical comparison applied p99 `1e-4` and maximum `1e-3` to all three
populations. Current production applies only p99 `1e-4` to the active-DOF
population. A zero-sized coordinate population still passes vacuously and is
reported as restricted/all-fixed rather than full convergence.

## Failure Mode x Safeguard matrix

| Failure mode | Current stationarity | Strict stationarity | Active member | Shared DOF |
| --- | --- | --- | --- | --- |
| Active IRLS reaches its iteration limit with a finite endpoint | Can pass | Unique blocker | - | - |
| Active local refit has a finite non-success status | Can depend on cluster rollup | Unique blocker | - | - |
| Guard-damped material active block | Can be represented only by cluster bool | Direct block-level blocker | Complementary | Complementary |
| Fixed/quarantined zeros dilute p99 | No protection | Restricted-state label | Unique exposure | Unique exposure |
| Large shared group dominates offset p99 weighting | No protection | - | Still member-weighted | Removes size weighting |
| Extreme weak-peak member inside one group | - | - | Maximum catches it | Within-group maximum preserves it |
| Mixed activity inside one shared group | Can be hidden in rollup | Rejects strict shadow | Partial members remain visible | Forced non-finite failure |
| Invalid/non-finite active transformed coordinate | - | Hard/invariant evidence | Fails | Fails whole DOF sample |

## Reproduced synthetic evidence

- Active `IrlsMaximumIterationsReached`: production qualification is true while
  solver qualification is false; the offset group is classified soft
  unqualified.
- All shape and offset blocks quarantined: solver qualification passes
  vacuously and reports both `restricted` and `all-fixed`.
- All local refit statuses are classified explicitly. Only `SUCCESS` with a
  full-step qualification passes; a successful but damped refit remains soft
  nonstationary.
- One 100-member zero-change offset group plus one singleton at `5e-4`:
  member p99 passes while the two-DOF p99 fails. The difference is group-size
  weighting, not a threshold change.
- A member at `2e-3` remains the one-group shared-DOF maximum. A non-finite
  member makes both shared-DOF p99 and maximum fail.
- A mixed active/fixed shared group produces a solver-qualification violation
  and a failing shared-DOF change summary.

These cases establish mathematical and semantic independence. They do not by
themselves establish the frequency or quality impact on production datasets.

## Debug record and aggregation

Every accepted Debug attempt now emits schema `5` on the existing
`Convergence safeguard audit:` record. It includes:

- production active-DOF and legacy all-selected accepted/raw summaries;
- production, legacy-population, legacy-maximum, and solver-qualified
  predicate vectors;
- activity/qualification counts for shape and offset blocks;
- shared-offset group counts and min/p50/p99/max group sizes;
- shape-active, offset-active, and quarantine ratios;
- orthogonal-clear and four policy stop candidates;
- isolated legacy-population, maximum-gate, and solver-qualification exposures.

The former active-member p99/maximum/median track is no longer emitted or
aggregated. Active atom membership and group-size diagnostics remain because
they define the production DOF population and exposure strata.

Aggregate a captured Debug log with:

```bash
python3 resources/tools/developer/analyze_convergence_audit.py run.log
python3 resources/tools/developer/analyze_convergence_audit.py \
    run.log --format json --output audit.json
```

The analyzer reports blocker/unique-blocker counts, pairwise truth tables,
implication counterexamples, p99-to-maximum results split at `N=91`, and
strata for active/quarantine ratio, group size, solver status, and proposal
path. A production stop candidate rejected by a comparator is an exposure; a
mismatch alone does not claim downstream quality loss.

## Historical second-round decisions

| Mechanism | Decision | Evidence needed before production change |
| --- | --- | --- |
| Current stationarity rollup | Semantic redesign candidate | Count actual `production=true, solver=false` stop exposures on representative trajectories. |
| All-selected population | Semantic redesign needed | The implementation contradicts the documented active-coordinate semantics and admits fixed-zero dilution. |
| Active-member offset population | Retain as comparison | It preserves per-atom transformed geometry but remains group-size weighted. |
| Shared-DOF offset population | Preferred redesign candidate | It removes group-size weighting and preserves within-group extremes in targeted tests; confirm on real trajectories. |
| Maximum gate | No removal in this round | Recompute unique catches after active population sizes are observed; only `N<=91` cases without unique coverage are conditional ablation candidates. |

Zero observed unique catches remain empirical evidence only. A safeguard is not
called redundant without a mathematical implication or a subsequent
trajectory-changing ablation.

The subsequent production change adopted the preferred shared-DOF population
and removed the maximum gate. The active-member comparison was then retired.
Solver qualification remains diagnostic because the available counterfactual
evidence did not justify a production change.

## Verification and external evidence

The historical focused tests covered the stationarity status matrix, quarantine exclusion,
unequal group weighting, extreme/non-finite members, mixed activity, Debug
schema, serial/parallel trace equality, and exact Info/Debug final model values.
The developer analyzer now has a schema-5 fixture using isolated comparators.

Fold-168 inputs remain external to the repository. The current refresh uses
the locally available hash-matching inputs as a negative control without
changing fitting options or stopping behavior.

The trajectory-changing follow-up is the build-gated
[counterfactual convergence continuation audit](counterfactual-convergence-continuation-audit.md).
It suppresses only an exposed production convergence stop and evaluates policy
checkpoints on isolated finalization workspaces.

Current refresh verification on 2026-08-27:

- Audit-enabled CTest passes 21/21 and audit-disabled CTest passes 19/19.
- Fold-168 stops after seven accepted iterations with `audit-patience` in both
  builds; the audit report contains seven schema-5 records, no convergence
  trigger, and zero comparator exposures.
- Audit-enabled and audit-disabled fold-168 `actual.json` files are
  byte-identical; repository lint passes.
- The 600-case exposure corpus was not rerun.

Historical verification from the earlier 2026-08-27 audit:

- all seven second-round C++ audit tests and the schema-2 analyzer fixture test
  pass;
- the second-stage defense suite now passes 117 of 124 tests, adding seven
  passing audit tests without changing the same seven pre-existing summary-log
  spelling failures recorded by the first round;
- complete `ctest` passes 15 of 16 groups; the only failing group contains
  exactly those seven pre-existing assertions;
- repository lint and Python byte-code validation pass;
- fold-168 runner self-tests pass, but the external simulation is not enabled
  because `RHBM_GEM_FOLD_168_MODEL` and `RHBM_GEM_FOLD_168_MAP` are unset.
