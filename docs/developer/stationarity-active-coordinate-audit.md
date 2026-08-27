# Stationarity semantics and active-coordinate population audit

## Purpose and constraints

This is the second-round convergence safeguard audit for
`RunSecondStageLocalFitting`. It follows the
[first-round audit](convergence-safeguard-audit.md) and compares two shadow
definitions against the production convergence gate:

- convergence-qualified stationarity for truly active parameter blocks;
- active-member and independent shared-offset-DOF change populations.

The audit is observation-only. It does not change the production conjunction,
thresholds, stop reason, trajectory, or final model. Debug-disabled runs do not
build the second-round summaries.

The current Notion algorithm page says that quarantined fixed zeros cannot hide
active-variable change. The implementation still computes the production p99
and maximum from every selected atom. The statement therefore describes the
intended semantics, not the current stopping expression. This document records
that discrepancy rather than silently treating either source as authoritative.

## Stationarity semantics

Activity and qualification are separate axes:

| Activity | Meaning | Included in strict conjunction |
| --- | --- | --- |
| Active | The block remains eligible to move this attempt. | Yes |
| Fixed | A current failure or fallback removed the block from the active parameterization. | No; reported as restricted state |
| Quarantined | A stage-local quarantine mask fixed the block. | No; reported as restricted state |

| Qualification | Shape atom | Shared-offset DOF |
| --- | --- | --- |
| Qualified stationary | Local status `SUCCESS`, full guard factor, no fallback | Joint status `Converged`, full group factor, no fallback |
| Soft nonstationary | Finite usable endpoint with a non-success status or damping | IRLS objective deterioration, maximum iterations, or damping |
| Hard failure | Missing/invalid refit evidence for an active block | System build, empty system, initial solve, or IRLS solve failure |

`eligible for update` is deliberately not equivalent to `eligible for
convergence`. Soft endpoints retain their existing update and objective-gate
behavior but fail the strict shadow predicate. Fixed and quarantined blocks are
excluded from the active conjunction. An empty active set passes vacuously and
is explicitly labelled `all-fixed/restricted`, not full convergence.

The production predicate remains
`ClusterHealth::is_active_block_stationarity_eligible`. The Debug snapshot also
records the strict predicate, shape/offset qualification counts, restricted
state, and `current=true, strict=false` disagreements.

## Population definitions

Accepted and raw transformed changes are summarized three ways:

| Population | log peak / log width | offset / peak |
| --- | --- | --- |
| Production | Every selected atom | Every selected atom |
| Active member | Shape-active atoms | Offset-active atoms |
| Active DOF | Shape-active atoms | One sample per `(cluster, group_id)` offset column |

The shared-DOF offset sample is the maximum absolute member
`delta(offset/peak)`. This prevents group size from becoming an implicit
statistical weight while retaining the most sensitive weak-peak member. A
non-finite member makes the group sample non-finite. Mixed active/fixed members
inside one shared-offset group are an invariant violation: the DOF summary is
forced to fail instead of choosing a partial group silently.

The same p99 `1e-4` and maximum `1e-3` thresholds are applied to all three
populations. A zero-sized coordinate population has summary value zero and
passes vacuously; its population size and restricted-state label remain in the
record.

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

- Active `IrlsMaximumIterationsReached`: current stationarity is true while
  strict stationarity is false; the offset group is classified soft
  nonstationary.
- All shape and offset blocks quarantined: strict stationarity passes
  vacuously and reports both `restricted` and `all-fixed`.
- All local refit statuses are classified explicitly. Only `SUCCESS` with a
  full-step qualification passes; a successful but damped refit remains soft
  nonstationary.
- One 100-member zero-change offset group plus one singleton at `5e-4`:
  member p99 passes while the two-DOF p99 fails. The difference is group-size
  weighting, not a threshold change.
- A member at `2e-3` remains the one-group shared-DOF maximum. A non-finite
  member makes both shared-DOF p99 and maximum fail.
- A mixed active/fixed shared group produces a strict-stationarity violation
  and a failing shared-DOF change summary.

These cases establish mathematical and semantic independence. They do not by
themselves establish the frequency or quality impact on production datasets.

## Debug record and aggregation

Every accepted Debug attempt now emits schema `2` on the existing
`Convergence safeguard audit:` record. In addition to the first-round fields,
it includes:

- active-member and active-DOF populations and accepted/raw summaries;
- current, member-strict, and shared-DOF-strict predicate vectors;
- activity/qualification counts for shape and offset blocks;
- shared-offset group counts and min/p50/p99/max group sizes;
- shape-active, offset-member-active, and quarantine ratios;
- orthogonal-clear, production, member-shadow, and DOF-shadow stop candidates;
- stationarity and population premature-convergence exposure flags.

Aggregate a captured Debug log with:

```bash
python3 resources/tools/developer/analyze_convergence_audit.py run.log
python3 resources/tools/developer/analyze_convergence_audit.py \
    run.log --format json --output audit.json
```

The analyzer reports blocker/unique-blocker counts, pairwise truth tables,
implication counterexamples, p99-to-maximum results split at `N=91`, and
strata for active/quarantine ratio, group size, solver status, and proposal
path. A production stop candidate rejected by either shadow policy is called a
**premature-convergence exposure**. Because this audit never continues beyond
the production stop, that term does not claim observed downstream quality loss.

## Second-round decisions

| Mechanism | Decision | Evidence needed before production change |
| --- | --- | --- |
| Current stationarity rollup | Semantic redesign candidate | Count actual `current=true, strict=false` stop exposures on representative trajectories. |
| All-selected population | Semantic redesign needed | The implementation contradicts the documented active-coordinate semantics and admits fixed-zero dilution. |
| Active-member offset population | Retain as comparison | It preserves per-atom transformed geometry but remains group-size weighted. |
| Shared-DOF offset population | Preferred redesign candidate | It removes group-size weighting and preserves within-group extremes in targeted tests; confirm on real trajectories. |
| Maximum gate | No removal in this round | Recompute unique catches after active population sizes are observed; only `N<=91` cases without unique coverage are conditional ablation candidates. |

Zero observed unique catches remain empirical evidence only. A safeguard is not
called redundant without a mathematical implication or a subsequent
trajectory-changing ablation.

## Verification and external evidence

The focused tests cover the stationarity status matrix, quarantine exclusion,
unequal group weighting, extreme/non-finite members, mixed activity, Debug
schema, serial/parallel trace equality, and exact Info/Debug final model values.
The developer analyzer has a fixed schema-2 fixture test.

Fold-168 remains optional because the model and map are external to the
repository. When hash-matching inputs are available, its Debug log can be
aggregated without changing fitting options or stopping behavior.

The trajectory-changing follow-up is the build-gated
[counterfactual convergence continuation audit](counterfactual-convergence-continuation-audit.md).
It suppresses only an exposed production convergence stop and evaluates policy
checkpoints on isolated finalization workspaces.

Verification on 2026-08-27:

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
