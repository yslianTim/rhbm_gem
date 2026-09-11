# Shadow global-best-only rejection and complete-state ablation

Baseline: clean `5174961e2373f97738ec8201666b643531c1bbf1`.
Production stays ON. The independent OFF variant bypasses only complete-state
best rejection; it does not remove cooperative-component best rejection or
reuse conclusions from previous combined experiments.

## Isolation and gate paths

The OFF archive changes one argument in `EvaluateFinalSelectionAudit`:

```cpp
GlobalCandidateReference{affected_sample_ref_list, inputs.objective_domain,
    nullptr, &previous_audit_objective, inputs.performance_counters}
```

The baseline passes `best_audit_objective` instead of `nullptr`. Both initial
complete-state audit and greedy salvage reevaluation use this entry point.
The original best pointer remains available for shadow observation. Previous
availability, finite checks and previous non-regression are preserved.
Cooperative endpoint/correction/backtracking still use their original best
reference, member checks and strict-improvement requirement. Local and ordinary
boundary previous gates, salvage ordering, best updates/reevaluation, patience,
non-converged final-state selection, recovery, stopping and final polish remain
unchanged. A passing OFF audit can naturally avoid salvage; no alternate salvage
policy is introduced. Direct gate helper tests retain their original assertions.

## Shadow definition and attribution

Build-local observers read the already computed `EvaluateObjectiveDelta` result
inside `EvaluateCombinedObjective`. They neither recompute the candidate
objective nor add to production counters. Separate comparisons use the existing
`IsAuditObjectiveAcceptableForProgress`, `IsObjectiveDeteriorated`, and
`kObjectiveProgressTolerance`, with the original best pointer even in OFF:

```text
global_best_only_rejection = passed_previous && failed_global_best
```

Each record retains the case, run, attempt, objective revision, evaluation
sequence, source, keys, candidate/previous/original-best values, previous and
best outcomes, actual gate outcome, and candidate Gaussian parameters. Commit
records separately retain final selected keys and actual committed parameters.
A gate pass alone is not counted as committed retention: keys must survive and
their parameters must match. Exact whole-state retention is also distinguished.

Sources distinguish complete-state audit after ordinary boundary processing,
after cooperative processing, and their salvage reevaluations. Cooperative
endpoint, correction and backtracking are separate sources. Ordinary boundary
calls without a best reference are not pooled into best-gate denominators.
Comparisons count gate evaluations reaching the combined objective function,
not proposals already stopped by earlier member/guard checks. For cooperative
records, passing the previous/best comparisons does not imply passing the
subsequent strict-improvement requirement.

Candidate/previous unavailable or nonfinite, absent best, nonfinite best,
previous failure, finite best-only rejection, and passing both comparisons are
separate classifications. Raw flags remain available even when a classification
is unavailable. Only finite candidate/previous/best comparisons contribute to
the finite-comparable denominator; only its previous-passing subset contributes
to the independent finite-best rejection rate. Attempt and run counts avoid
confusing repeated salvage evaluations with independently affected executions.

## Existing validation matrix

All eight independent builds and complete CTest runs pass. The shadow captures
also pass the full existing executable. Production source and tests remain
unchanged; there is no new case, assertion adjustment, tolerance relaxation,
API option, or diagnostic schema. Wall times below record validation execution
only; concurrent builds/tests are not performance evidence.

| Audit trace | Trust model | ON | Complete-state OFF |
| --- | --- | ---: | ---: |
| OFF | OFF | 16/16 (94.40 s) | 16/16 (93.83 s) |
| ON | OFF | 16/16 (243.30 s) | 16/16 (242.72 s) |
| OFF | ON | 16/16 (112.72 s) | 16/16 (112.78 s) |
| ON | ON | 16/16 (262.27 s) | 16/16 (262.16 s) |

## Observed exposure and result

Both isolated shadow executables pass all 861 existing cases in 76 suites.
ON and OFF have the same counts and byte-identical gate/commit records:

| Gate source | Evaluations | Finite comparable and previous-passing | Best-only rejection |
| --- | ---: | ---: | ---: |
| Complete-state after ordinary boundary | 306 | 306 | 0 |
| Complete-state after cooperative processing | 6 | 6 | 0 |
| Complete-state salvage reevaluation | 0 | 0 | Not exposed |
| Cooperative endpoint | 6 | 6 | 0 |
| Cooperative correction | 6 | 6 | 0 |
| Cooperative backtracking reaching global comparison | 0 | 0 | Not exposed |

Complete-state comparisons cover 312 distinct attempts in eight runs of four
existing cases. Cooperative comparisons cover six attempts in three runs of
`BoundaryJointCorrectionMatchesSerialParallelAndIntensityScaling`; its endpoint
and correction records are not twelve independent attempts. Complete-state
exposure by case is 90, 90, 87, and 45 evaluations respectively for
`BoundaryComponentReconciliationMatchesSerialAndParallelSelection`,
`BoundaryComponentReconciliationIsIntensityScaleInvariant`,
`BoundaryJointCorrectionMatchesSerialParallelAndIntensityScaling`, and
`BoundaryComponentReconciliationBacktracksAndPreservesRemoteCluster`.

The finite best-only rates are 0/312 for complete-state and 0/12 for cooperative
global comparisons. Both also have zero raw `passed_previous && failed_best`
events, zero affected attempts/runs, zero best-only candidates released by OFF,
and zero such candidates retained at commit. The additional 582 ordinary
boundary evaluations have no best reference and are excluded from those rates.
No unavailable/nonfinite or previous-failing comparison reaches the measured
best-bearing gates in this corpus. Earlier rejection paths are not counted as
global evaluations. Direct gate helper assertions still execute in CTest, but
are not production best-gate exposure.

All 4,048 gate/commit records (906 evaluations and 3,142 commits) match exactly.
The shadow outcome is checked against the actual variant policy: complete-state
OFF uses previous alone; all other calls retain previous and best. The separate
behavior comparison also matches exactly:

| Evidence | Identical records |
| --- | ---: |
| Assembled-state trajectories | 3142 |
| Candidate eligibility/decision/selected-key events | 12598 |
| Recovery/quarantine revisions, targets and masks | 9426 |
| Finalization and post-persistence evidence | 236 |
| Selected history/phase/convergence/terminal diagnostics in 29 captures | 26369 |

Finalization evidence contains 75 base states with stop reason and accepted
iteration count, 75 persisted states, 75 model/MSE/peeling records, five final
polish summaries and six component summaries. Existing captured robust terminal
objectives match, as do persisted parameters, raw-response MSE and peeling
responses. No first production divergence is found because no best-only
rejection is bypassed. The earlier source-identical control without the best
shadow observer matches all candidate, trajectory, recovery and non-timing
finalization records as well. No regression is hidden by relaxed assertions.

## Decision and coverage limits

**Evidence is insufficient to downgrade the gate; retain production ON.** The
experiment measures no independent best rejection in the existing data, but
cannot evaluate the terminal consequences of releasing one. Neither its commit
tracking nor complete-state salvage/best-failure path is dynamically exercised
by a best-only event. This is not evidence that previous checks logically imply
the best check: tolerance-based previous non-regression need not imply
non-regression against a different historical best reference.

Cooperative best rejection remains enabled in OFF. Even a future positive
complete-state-only ablation would not establish that all global-best hard
gates can be demoted. Best update/reevaluation, audit patience, non-converged
final-state selection, diagnostics and final polish remain unchanged here.
No new corpus, quality threshold, test framework or numerical baseline is
introduced; existing serial/parallel, quiet/Debug, intensity scaling, hard
failure, remote-cluster and final-state assertions are retained. This is a
bounded existing-data experiment, not a general quality or performance
certification. No production gate is deleted or downgraded.

## Reproduction and evidence handling

Artifacts are workspace-local under `build/global-best-only/`:

- `before-source/` and `off-source/` are baseline archives with the one OFF
  argument change above. `source-isolation.json` verifies the difference.
- `verify.py before|off` configures independent Ninja/Debug builds, testing ON,
  Python bindings OFF, SYSTEM dependencies, and all OFF/ON combinations of
  audit trace and trust-model experiment. Existing dependency caches are reused.
- `capture.py before|off`, `best-capture.hpp`, and `best-wrapper.cpp` create
  isolated trace-ON/trust-OFF observer libraries and run the full existing C++
  test executable. `best.tsv` retains gate/commit evidence; the other captures
  retain recovery, candidate, trajectory, final model/MSE/peeling and diagnostics.
- `analyze-best.py` computes per-source denominators, affected attempts/runs and
  committed retention; `compare.py` and `compare-logs.py` compare outcomes.
  Timing fields are excluded; existing unordered solver/compatibility logs are
  compared as multisets. No numerical tolerances or golden data are changed.

The baseline's source identity against the earlier recovery-revision capture
is checked by hashes and the intervening Git diff. This supplies a control with
no best-gate shadow observer, in addition to this experiment's fresh
uninstrumented ON CTest matrix. The control does contain the earlier read-only
trajectory observers; it is not described as a wholly uninstrumented numerical
capture. Scripts and raw logs remain build-local evidence, not new versioned
test infrastructure. Preserve them when transferring the detailed evidence.
