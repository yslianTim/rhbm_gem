# Objective exhaustion quarantine policy audit

Baseline: `d511c88f` (clean working tree).
This is a numerical policy ablation: removing a freeze/recovery observation can
change later activity, proposals, convergence blockers, stopping and final output.
It is not a numerical-equivalence refactor.

## Evidence boundary

`QuarantineState::UpdateAfterIteration` is the sole production builder of the
quarantine failure map. Solver hard failures are inserted first. Both accepted
and rejected candidate terminal lists are then scanned, retaining only
`GuardInfeasible` and `InvalidCandidate`. Each terminal entry is filtered
independently, so another guard/invalid observation in the same list is retained.
The existing first-insertion precedence and shape/offset/cluster mapping remain.

`CandidateSelection` still produces `ObjectiveExhausted` when objective-rejected
trials exhaust the search; its priority over invalid trials is unchanged.
`Diagnosis` still prints `objective-exhausted`. Neither the enum nor the terminal
list is removed: search termination and quarantine evidence have different roles.

## State and recovery effects

The filtered map feeds both consecutive-failure tracking and affecting-failure
checks for domain retry. Objective exhaustion alone cannot create a tracking
entry or freeze a target, and an absent observation clears an Active target's
previous streak. Five consecutive matching retained failures still freeze.
`IterationProcess` also resets audit patience while Active failure tracking is
pending. Removing objective-only entries can therefore allow patience to advance
before any target would have reached the freeze threshold; that downstream
effect is part of this evidence removal, with the patience policy unchanged.

A Frozen target is not released merely because objective exhaustion was filtered.
A new domain revision must allow retry, all required coordinates must be active,
and recovery must have a finally accepted material change or a complete,
guard-safe, solver-qualified nonmaterial unrestricted endpoint. Other affecting
failures still block recovery. Failed retries retain Frozen state.

The state-update helper, domain revision rules, thresholds, retry radius/ridge,
mask overlap handling, transaction publication and final activity are unchanged.
Quarantine state is local to the second-stage execution; no persisted-state
migration or diagnostic schema change is required.

## Retained numerical policy and limitations

No edits are made to factor search, immediate guard isolation, exhausted-key
radius handling, candidate gates, rescue, global previous/best audit, all-rejected
stop or convergence tests. Their inputs may nevertheless change in later
iterations because objective exhaustion no longer creates quarantine evidence.
Nonimproving clusters can remain active and incur repeated proposal cost while
other clusters progress. Existing regression results alone do not establish
real-data quality or performance neutrality.

## Verification

Before/after Debug builds and the full existing CTest suite use all four
combinations of audit trace and trust-model experiment. No test files, cases,
assertions, tolerances or numerical baselines are changed. Build-local copies of
existing stdout-capture calls retain the original fixture assertions while saving
logs for trajectory, quarantine and final-state comparison.

All four before builds and all four after builds passed the full 16-group CTest
suite. Artifacts, commands, per-configuration results and exact comparison rules
are recorded in `build/quarantine-objective-exhausted/report.md`. The initial
after OFF/OFF documentation-link check ran before this newly linked file existed;
its failure log is preserved, and the complete suite passed on rerun after the
file was written. No test assertions were changed.

The trace-enabled, trust-model-disabled capture run passed all 107 enabled
existing defense cases on each side (108 case declarations remain in the source;
the shadow case is conditional). It produced 29 captures and 26,301 selected
production/diagnostic records per side. Twenty-three captures have identical
selected records; six frozen-background captures show the intended policy effect.
All 22 terminal records and 145 terminal-atom records are identical, including
stop reason, iteration counts, objective and final parameters. Phase event IDs
and their order are unchanged. There are 573 objective-exhaustion diagnostic
lines on each side.

In each of the six changed captures, objective exhaustion at attempts 6–10
formerly froze one atom at attempt 10; the new run keeps quarantine at zero.
At attempt 11, domain retry formerly applied the existing minimum-radius/10x-ridge
path. Its absence changes proposal and counterfactual operator diagnostics at
attempts 11–100. Baseline, local-search/polish, assembly, boundary and
final-selection objective records remain identical; production/unrestricted
proposal objectives differ. Both runs stop at 100 accepted iterations with
`maximum-iterations` and identical final atom records. These are observed
trajectory differences, not suppressed comparison failures.

Comparison preserves numerical values without rounding. Timing fields are
excluded and unordered solver/compatibility worker records are compared as
multisets. Terminal records are compared from their diagnostic marker because a
progress-table prefix on the same line contains the intentionally changed atom
quarantine count. Raw differences, event-order checks and field-specific results
are retained in `trace-comparison.json`, `behavior-comparison.json` and the diff
files. The captures exercise objective-driven freeze removal; they do not
establish recovery behavior for every mixed guard/invalid failure combination or
measure production performance.
