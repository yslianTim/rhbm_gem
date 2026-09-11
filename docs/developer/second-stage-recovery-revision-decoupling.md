# Objective and Frozen-recovery revision separation

Baseline: clean `e6296fdf604779f7259c99326de1258158df492a`.
This is a structural change retaining the existing recovery trigger policy.
It does not reduce background-triggered retries or change numerical acceptance.

## Ownership and dependency audit

`IterationState` retains `objective_domain_revision` and independently initializes
`frozen_recovery_revision` to 1. Recovery never copies, derives from, or compares
against objective revision. Objective revision remains the phase observer's
objective-context identifier. Global-best evaluation, history resets, convergence
and `objective_domain_changed` retain their existing inputs and behavior.

Only the recovery revision flows through `QuarantineState::BeginIteration`,
`CandidateTransactionBuilder::Finish`, `QuarantineState::UpdateAfterIteration`,
and `UpdateQuarantineFailureState`. Each target's `last_recovery_revision`
(default 0) records freezing, retry initiation, and failed retry. The old
`last_domain_revision` field is removed. The transaction still computes and
commits the same next-iteration quarantine state; no audited model is rewritten.

| Event | Objective revision | Frozen-recovery revision |
| --- | --- | --- |
| Initialization | 1 | 1 |
| Queued partition applied at the next attempt | +1 | +1 |
| Same partition, different `response_by_atom` | +1 | +1 |
| Partition applied together with background change | +1 once | +1 once |
| Same partition and exactly equal background response | unchanged | unchanged |
| Topology merely queued | unchanged | unchanged |

Objective revision's increment in `ResetIterationStateForPartition` stays there.
Recovery advances in the partition-application branch of
`BeginFrozenBackgroundIteration`, after reset succeeds, rather than inside the
objective reset function. The unchanged-partition/background-change branch also
advances recovery explicitly. Its exact response comparison and early return
are preserved. Objective reevaluation alone is not a recovery event.

Each Frozen target retries once per newer recovery revision. Failed retries
remain Frozen, and successful release still requires no affecting failure,
active required coordinates, and an accepted material update or qualified safe
nonmaterial endpoint. Failure streaks, reason priority, overlapping masks,
minimum radius, ridge, audit patience, stopping, and convergence are unchanged.
Background updates can still permit a retry every attempt. Reducing that policy
requires a separate numerical experiment.

## Interfaces and diagnostics

No public API, options, diagnostic fields, or schema change. Existing
`domain-retry`/probation labels and the existing test names remain for
compatibility; they now refer to recovery revision scheduling. Existing tests
only rename revision variables and the target-state field, preserving their
assertions and numerical tolerances. No new test case, framework, or baseline
is introduced. Historical audit documents are not rewritten.

## Existing validation matrix

All eight independent builds and complete CTest runs pass. Matrix wall times
are execution records, not performance comparisons; validation ran concurrently.

| Audit trace | Trust model | Before | After |
| --- | --- | ---: | ---: |
| OFF | OFF | 16/16 (95.92 s) | 16/16 (94.57 s) |
| ON | OFF | 16/16 (247.72 s) | 16/16 (246.12 s) |
| OFF | ON | 16/16 (116.42 s) | 16/16 (113.22 s) |
| ON | ON | 16/16 (264.78 s) | 16/16 (263.53 s) |

The implementation and adapted-test hashes in `source-verification.json`
match the validated source. `comparison.json` and `trace-comparison.json`
record behavior and diagnostics equality. Recovery no longer reads objective
revision; its two explicit event increments are outside objective reset and
reevaluation helpers. No test tolerance, numerical constant, diagnostic schema,
or user API was changed. The captured behavior is preserved within the coverage
limits below; this is not a new general quality or performance certification.

## Behavioral comparison and coverage

Both isolated trace-ON/trust-OFF capture executables pass all 861 existing cases
in 76 suites. The complete ordered comparison is exact, excluding only phase
timing values:

| Captured evidence | Identical records |
| --- | ---: |
| Revision, quarantine targets/reasons/lifecycles/counters and active masks | 9426 |
| Assembled-state trajectories | 3142 |
| Candidate eligibility, decisions and selected keys | 12598 |
| Finalization and persisted Gaussian/MSE/peeling evidence | 236 |
| Selected history/phase/convergence/terminal diagnostics in 29 stdout captures | 26369 |

Finalization evidence includes 75 base states with stop reasons and accepted
iteration counts, 75 persistence states, 75 post-persistence model records,
five polish summaries and six component summaries. It introduces no new
numerical baseline. Existing unordered solver/compatibility log records are
compared as multisets; other selected diagnostics retain order. No first
numerical or recovery-event divergence is observed.

Across 3,142 production attempts, the capture observes 75 initializations,
2,466 changed-background events, and 601 unchanged-background events. Every
record's objective and recovery revisions remain numerically equal under the
retained trigger policy, although their ownership and inputs are separate.
There are eight entries into Frozen, 190 retry attempts, 190 failed retries,
and no successful production releases. Frozen targets occur in 265 begin
records; target identities, failure reasons, streaks, last attempted revision,
retry lists and shape/offset/hard-failure masks match the baseline exactly.

The existing `PersistentQuarantineReasonRequiresStableReasonAndReleasesOnDomainRetry`
case retains its stable-reason, release and repeated failed-retry assertions
with only internal revision names adapted. Other existing endpoint and mask
assertions remain unchanged. No newly generated recovery scenario is used.

Coverage limits: the full executable capture does not queue or apply a changed
partition, including simultaneous partition/background change, and does not
exercise successful release through the full production loop. These are not
claimed as end-to-end dynamic coverage. Source inspection verifies that the
partition branch advances recovery once and returns before the background-only
increment; queuing and objective reset/reevaluation do not advance recovery.
The existing direct release assertions supplement, but do not replace, the
missing full-loop release coverage. No independent objective-only update event
exists in current production triggers; its lack of a recovery dependency is
verified by source wiring, not by inventing a numerical scenario.

## Reproduction

Workspace-local evidence is retained under `build/recovery-revision/`:
`before-source/` is an archive of the baseline; `verify.py before` and
`verify.py after` use independent Ninja/Debug builds with testing ON, Python
bindings OFF, SYSTEM dependencies, and all four audit-trace/trust-model switch
combinations. Existing dependency caches are reused. The after build reads the
working source, with only the independent counter wiring and interface renames
described above, without numerical policy changes.

`capture.py before|after` creates isolated capture copies from the trace-ON,
trust-OFF builds and runs the full existing C++ test executable. It records
revision events, target lifecycle/reason/retry state, transition counts and
active masks, candidate events, assembled-state trajectories, finalization
states, and post-persistence Gaussian/MSE/peeling data. Captures do not provide
decision values or modify production counters. Baseline capture labels the
original objective revision as its legacy recovery input; the after capture
reads the independent counter. Only the baseline field-name adapter differs.
The scripts and raw captures are build-local evidence, not versioned test
infrastructure; retain them when transferring the detailed comparison.
