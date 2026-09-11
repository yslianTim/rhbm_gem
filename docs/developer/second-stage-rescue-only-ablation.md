# Rescue-only ablation and retained component acceptance integration

Baseline: clean `ddc165052a2bce5828829cfd04d43385c3671a50`.
The independent experiment does not reproduce the historical combined
rescue/Grow regression. It finds a measurable quality/objective/cost tradeoff,
so rescue is retained rather than deleted. After reviewing that tradeoff, the
user explicitly requested applying the structural integration. The integration
keeps rescue-ON numerical policy and ordering; it does not adopt rescue-OFF.

## Isolation and existing-test results

The isolated OFF source changes only the condition invoking cooperative
component search. Finite rejected proposals remain collected. Ordinary boundary
endpoint/correction/backtracking, global previous/best audit and greedy salvage,
guards, quarantine, radii, convergence and final polish are unchanged. Grow
remains absent in both variants. Timing therefore measures skipping rescue
execution, not the total cost of deleting its storage.

Both ON and OFF pass all 16 existing CTest groups in each of the four audit-trace
and trust-model switch combinations. In particular, intensity scaling,
serial/parallel, hard failure, healthy remote cluster and final-state assertions
remain intact. No test case, framework, assertion tolerance or numerical baseline
is added or weakened. The historical combined experiment cannot establish that
rescue is necessary, nor that it is safe to delete.

## Exposure and final outcomes

Build-local copies of existing fixtures and two production translation units
record per-attempt/per-key finite proposals, component decisions and selection
after final salvage. The instrumented libraries are isolated from the production
CTest matrix. Both trace-ON/trust-OFF capture runs pass all 107 enabled defense
cases (108 declarations remain, including one conditional trust-model case).

| Observation | Rescue ON | Rescue OFF |
| --- | ---: | ---: |
| Attempt records | 2759 | 2756 |
| Finite rejected proposal key observations | 1144 | 1147 |
| Cooperative component attempts | 15 | 0 |
| Accepted cooperative keys | 6 | 0 |
| Cooperative keys retained after final salvage | 6 | 0 |

These are repeated observations across attempts/runs, not unique atoms or cases.
All actual rescue attempts occur in the existing
`BoundaryJointCorrectionMatchesSerialParallelAndIntensityScaling` fixture.
Each of its serial, parallel and intensity-scaled runs has five rescue attempts,
two accepted/retained keys, and three rejected rescue attempts. Finite rejected
proposals alone are insufficient: a component requires at least two keys linked
by shared boundary samples.

Across the full suite, first selection-state divergence is recorded at event
2611: attempt 8 of that fixture's first run. Rescue subsequently accepts again
at attempt 18. ON stops after 32 accepted attempts and selects best iteration 29;
OFF stops after 31 and selects best iteration 28. Both stop for `audit-patience`.

The existing response-MSE helper is evaluated after persistence in build-local
copies of the same fixture, retaining the original assertions:

| Run | ON final response MSE | OFF final response MSE |
| --- | ---: | ---: |
| Serial | 0.0010075571347304411 | 0.001009194886990455 |
| Parallel | 0.0010075571347304411 | 0.001009194886990455 |
| Intensity x100 | 10.075571347304086 | 10.0919488699049 |

ON improves response MSE by about 0.162%. The unscaled difference is 1.638e-6,
exceeding the existing objective-progress tolerance formula applied to this
metric (`1e-8 + 1e-3 * abs(reference)`, about 1.018e-6). This is an analysis
criterion using existing constants, not a new assertion or relaxed tolerance.

However, the terminal audit objective is approximately 1.681656e-4 for ON and
1.678245e-4 for OFF: OFF improves it by about 0.203% and uses one fewer attempt.
The Debug replay confirms both retain objective domain 1, identical initial
state/domain, and no domain changes. The audit-objective difference is also
larger than the existing progress tolerance at this scale. Robust audit objective
and raw response MSE are different metrics; neither improvement cancels the
other. Serial/parallel and x100-intensity replays reproduce the tradeoff.

## Decision and integrated production behavior

The ablation establishes neither a safety regression without rescue nor a
uniform advantage for either numerical policy. Rescue remains enabled under the
explicit decision to retain its capability and integrate its implementation.
The response-MSE/audit-objective/cost tradeoff above remains part of the evidence;
structural integration is not evidence that rescue dominates its removal.

The builder now keeps one provisional record per key, including candidate
status, cooperative patch, diagnostic, radius flags and rejection-event order.
Normal and cooperative components share endpoint/correction/backtracking and a
single result-application entry. Acceptance scope is explicit rather than read
from `diagnostic.is_rescue_attempt`. Final accepted/rejected lists are produced
only after global salvage; the special rejected-to-accepted promote function and
diagnostic transfers are removed. Rejection-event order remains deterministic.

The ordinary component sweep and its existing global audit/salvage still precede
the cooperative sweep and its conditional global audit/salvage. Component
construction, halo, member tolerances, strict/global-best gates, incomplete
baseline behavior and healthy remote updates remain unchanged. History observer
baseline, rollback, provenance and publication semantics, radius/quarantine
policy, stop conditions and final polish are retained. Public API and diagnostic
schemas are unchanged; rescue names remain diagnostic source labels.

## Integration verification

The exact integration initially explored during the ablation is now applied.
It matched rescue-ON for all 2,759 assembled-state trajectories, per-key candidate
and retention events, 29 stdout captures (26,369 selected records), and another
2,577 Debug/history/phase records in the rescue-active fixture. Final MSE and
Gaussian parameters match ON exactly in serial, parallel and x100 intensity runs.
No numerical tolerance is relaxed. The sole existing-test edit adapts the
internal boundary-reconciliation call signature; no case or assertion is added.

The applied source passes all 16 existing CTest groups in all four audit-trace /
trust-model combinations. The first two results come from the original checks
of these exact source bytes; the remaining two were completed when integration
was explicitly selected. Current documentation checks also pass. Per-build logs,
source hashes and the completed matrix are recorded in the build report.

## Reproduction and limits

Artifacts, commands, ON/OFF source archives, captured events, persisted quality
measurements and matrix logs are under `build/rescue-only/`; see `report.md` there.
The `final-*` artifacts and saved integration patch contain the same source now
applied; earlier `withdrawn-integration-*` names are preserved as historical
artifact labels. The current completion status is in the build report.
Comparisons exclude timing fields and unordered worker log order; no numerical
values are rounded to hide differences. The 29 ordinary stdout captures alone
show no ON/OFF differences because they do not capture the rescue-active fixture;
the per-key/trajectory instrumentation and targeted Debug replay expose it.

Production changes are limited to component evaluation/application and
provisional result ownership, plus the necessary existing-test interface edit. No 600-case corpus,
fold-168 numerical dataset, new benchmark or general real-data quality/performance
claim is made. This experiment does not establish every mixed-failure or
post-rescue rollback path's coverage.
