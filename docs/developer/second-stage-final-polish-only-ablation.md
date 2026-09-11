# Converged-only final dependency polish ablation

Baseline: clean `bf8840b0ea9cbc6f199518783896fb09c8462104` (rescue/component
acceptance integration). This experiment isolates final dependency polish; it
does not reuse rescue/Grow or rescue-only conclusions. Production remains ON.

## Isolation and observation

The ON archive is unchanged. The OFF archive changes only the early-return
condition in `FinalizeSecondStageState`:

```cpp
if (true || stop_reason != SecondStageStopReason::Converged ||
    !options.enable_second_stage_dependency_polish)
```

OFF therefore persists the same selected base state without final dependency
polish or its recertification. Option defaults, zero-round validation, and the
direct `RunFinalDependencyPolish` entry point remain intact. Outer selection,
local/boundary polish, rescue, quarantine, global-best and stop policy, and the
persistence implementation are identical. This is an experiment in isolated
source trees, not a production policy or API change.

Separate build-local capture copies record candidate eligibility/decisions,
assembled states, finalization base and persist states, component acceptance
before global salvage, surviving component count, certificate status, actual
application, and phase timings. A second capture uses the existing response-MSE
helper after finalization and serializes selected Gaussian parameters and every
peeling sample response. Values use 17 significant digits. Observations do not
supply decision values or increment production counters. Instrumented libraries
are separate from the uninstrumented CTest builds and timing executables.

The direct polish test is counted separately from production finalization.
Accepting a component is not evidence that a converged production run passed
global audit, recertification, and persistence. No test case, framework,
assertion tolerance, or numerical baseline is added or changed.

## Existing validation matrix

All builds and all 16 existing CTest groups pass in each variant/configuration.
No assertions are adapted or relaxed. Matrix wall times below are validation
records only; concurrent builds/tests make them unsuitable for speedup claims.

| Audit trace | Trust model | ON CTest | OFF CTest |
| --- | --- | ---: | ---: |
| OFF | OFF | 16/16 (95.21 s) | 16/16 (95.30 s) |
| ON | OFF | 16/16 (252.83 s) | 16/16 (249.76 s) |
| OFF | ON | 16/16 (115.92 s) | 16/16 (113.21 s) |
| ON | ON | 16/16 (263.53 s) | 16/16 (262.19 s) |

## Exposure and numerical result

The isolated trace-OFF capture passes all 107 enabled defense cases. Its 2,759
assembled-state records and 11,066 candidate events match exactly; none of its
62 production finalizations is converged. Broader trace-ON capture passes all
861 existing cases in 76 suites in both variants and finds five converged
finalizations. These are repetitions across executions, not five unique cases.

| Production observation in the full executable | ON | OFF |
| --- | ---: | ---: |
| Finalizations reaching the state-selection/persistence path | 75 | 75 |
| Converged and configured to enable polish | 5 | 5 |
| Final polish calls | 5 | 0 |
| Components / attempted component solves | 1 / 1 | 0 / 0 |
| Components accepted before global salvage | 0 | 0 |
| Components surviving global audit/salvage | 0 | 0 |
| Recertification evaluations / passes | 0 / 0 | 0 / 0 |
| Actual polish applications | 0 | 0 |

Four converged executions have no dependency component. The remaining execution,
`LocalFittingResultRanksUseSelectedAtomsWithin2AInSecondStage`, attempts one
component. A build-local replay of that unchanged case records the correction
status `no-material-change`; no candidate reaches objective acceptance. All five
ON runs retain their base state, with residual safety `not-evaluated`.

The separate direct test `FinalDependencyPolishImprovesUncutComponent` accepts
one component in both variants. It is excluded from the production counts above:
it does not establish converged-path global/certificate/application coverage.

The full executable's 3,142 assembled-state records and 12,598 candidate events
match exactly. All 75 base records (including stop reason and accepted-iteration
count), 75 persistence records, and 75 post-persistence model/MSE/peeling records
match exactly. Stop counts are five converged, 25 audit-patience, 21 all-rejected
backtracking-exhausted, and 24 maximum-iterations. No numerical first divergence
is found. The first expected capture difference is the skipped polish call in
`RunPotentialFittingWorkflowProducesSingleGroupResultAfterLocalStages`.

| Converged case | Executions | Final response MSE, identical ON/OFF | ON robust objective before = after |
| --- | ---: | ---: | ---: |
| `RunPotentialFittingWorkflowProducesSingleGroupResultAfterLocalStages` | 2 | 7.2434412588721273e-11 | 1.8327931361410855e-6 |
| `LocalFittingResultRanksUseSelectedAtomsWithin2AInSecondStage` | 1 | 0.00012514463324786443 | 2.2010234999378541e-9 |
| `RunSecondStageIterationsImprovesBadFiniteEntryScale` | 1 | 7.2434412588721273e-11 | 1.8327931361410855e-6 |
| `RunSecondStageIterationsHandlesNearPerfectEntryScale` | 1 | 1.3712204091621674e-31 | 2.97138574302629e-18 |

The robust values are ON's existing polish before/after evaluation, not a new
OFF objective evaluation. OFF writes the identical base parameters through the
unchanged persistence path. Existing captured robust terminal diagnostics also
match: 29 stdout captures contain 26,369 identical selected diagnostic records,
including history, phase, convergence, and terminal evidence. Full finalization
captures intentionally differ by five polish-result and five component-summary
records; `final-comparison.json` separately verifies their identical base,
persist, and model projections. Timings are not compared as numerical evidence.
Existing downstream workflow, serial/parallel, scaling, quiet/Debug, failure,
remote-cluster, and final-state assertions are retained.

## Isolated timing

After all concurrent experiment builds and tests completed, the four converged
cases ran once per variant as warmup and then five times per variant in
alternating order (ON/OFF, then OFF/ON). Each repetition passes all four existing
cases. Uninstrumented trace-OFF/trust-OFF executables measure total process wall
time; separate minimal instrumented copies measure the polish call interval
without capture writes inside that interval. Their base/persist/model records
also match ON/OFF in every repetition.

| Measurement | Median | Range |
| --- | ---: | ---: |
| ON total process | 1.593441 s | 1.591390–1.599061 s |
| OFF total process | 1.583815 s | 1.578048–1.591772 s |
| ON polish calls combined per repetition | 10.195999 ms | 10.092793–10.284579 ms |
| ON one attempted component (rank case) | 9.689040 ms | 9.579960–9.755120 ms |
| Recertification | Not invoked | No cost sample |

The observed total median difference is about 9.6 ms in this four-case batch,
consistent with skipping roughly 10.2 ms of non-applied polish work. The four
no-component calls each take about 0.13 ms. OFF executes neither measured phase.
No evaluator call occurs in ON either: the tiny timer interval around its
certificate condition is checking/clock overhead, not recertification cost.
These measurements cover a no-material-change solve and empty-component work,
not accepted-patch polishing or recertification, and cannot estimate general
speedup or the cost of fully deleting the mechanism. Process totals include
startup; phase and process measurements are separate runs, not additive values.

## Decision and limits

**Evidence is insufficient to support deletion; retain production ON.** The
existing data covers entry and one unsuccessful solve, but no component accepted
and applied on a converged production path, no final global salvage of an
accepted patch, and no recertification. Passing OFF and observing identical
outputs in these cases cannot establish redundancy. The direct solver test's
improvement is likewise insufficient to prove final production benefit.

This experiment does not add a corpus, real-data benchmark, new tolerance, or
replacement policy. Results do not certify general quality or performance.
No unexplained numerical regression is observed in the captured data.

## Reproduction artifacts

Build-local artifacts are under `build/final-polish-only/` (ignored by Git):

- `before-source/` and `off-source/`: archives of the baseline, with the single
  OFF condition above; `verify.py before` and `verify.py off` configure, build,
  and run the four CTest combinations.
- Builds use Ninja, Debug, `BUILD_TESTING=ON`, `BUILD_PYTHON_BINDINGS=OFF`,
  `RHBM_GEM_DEP_PROVIDER=SYSTEM`, and all OFF/ON pairs of
  `RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE` and
  `RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT`. Existing dependency source caches
  are reused; each variant/configuration has an independent build directory.
- `capture.py before|off` runs the existing defense cases with trace OFF;
  `capture-trace.py before|off` runs the existing test executable with trace ON.
  Both use trust OFF and isolated instrumented libraries. Generated copies
  retain all original test assertions.
- `analyze.py`, `compare-logs.py`, `comparison.json`, and
  `trace-comparison.json`, and `final-comparison.json` retain exact-state/event and selected diagnostic
  comparisons. Timing fields are excluded; existing unordered solver and
  compatibility log records are compared as multisets.
- `timing.py` performs one warmup per variant, then five alternating ON/OFF
  runs without concurrent experiment work. It uses the existing
  converged workflow, selected-atom rank, bad-entry-scale, and near-perfect-entry-scale cases.
  `timing.json`, `phase-timing.json`, and individual logs retain every duration
  and outcome; `timing-build.py` prepares the separate phase capture executables.

The scripts and raw logs are workspace-local evidence, not new versioned test
infrastructure. On a fresh checkout, recreate both archives from the baseline,
apply the condition above only to OFF, and run the stated build matrix; preserve
these local artifacts when transferring the detailed capture evidence.
