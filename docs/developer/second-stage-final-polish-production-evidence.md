# Final-polish production application evidence

## Baseline and decision rule

Experiment baseline: clean `fb926a5be58f0ade7ffd855bba44cc84789e79f0`.
Production remains unchanged after the experiment stopped at its baseline gate. The isolated OFF variant
only bypasses the converged finalization polish branch; it retains validation,
outer policies, persistence and the direct polish solver for comparison.

The agreed deletion gate requires all three specified workloads to complete with
zero actual polish applications at both four and one threads, followed by output
and regression validation of the removal. A failed or incomplete run is not a
zero-application result. A successful synthetic case alone would not override
this workload-based decision. Neither zero applications nor failure to find a
synthetic example establishes that the mechanism is unreachable or redundant.

## Inputs and execution

Exact input SHA-256 identities and commands are retained in
`build/final-polish-production/baseline.json` and each run's `command.json`.

| Case | Model | Map | Mode |
| --- | --- | --- | --- |
| 1 | `fold_test_model_0.cif` | `sim_map_gaus_grid0.10_charge1_bw0.50.map` | simulation, width 0.50 |
| 2 | `6z6u.cif` | `sim_map_gaus_charge1_grid0.30_6Z6U_bw0.50.map` | simulation, width 0.50 |
| 3 | `6z6u.cif` | `emd_11103_unsharpened.map` | experimental, default normalization |

Hydrogens are excluded, all remaining atoms are selected, and sampling uses the
existing deterministic default. Every run gets a new independent database.
The source files are not modified. Case 1 hashes match the existing external
fold-168 benchmark; its unchanged quality and iteration gates apply. Cases 1/2
are simulation evidence; only case 3 supplies experimental-map evidence.

## Bounded synthetic search

The search extends the existing rank/workflow fixture through
`BuildPotentialModelTestData` and `RunPotentialFittingWorkflow`. Its 128 fixed
inputs are the Cartesian product of spots O/N/C/CA, blurring widths
0.25/0.50/0.75/1.00, noise standard deviations 0/0.001/0.01/0.05 and seeds 42/43.
Spots vary the existing physical neighbor arrangements and atom counts. Noise
and width alter the generated input, not an intermediate fitting state. The
factory's `gaus_true` metadata is recorded but is not an independent control of
the generated potential.

No convergence, partition, accepted patch or certificate is injected. The
production default gates and stopping conditions remain intact. Build-local
cases check finite persisted amplitude and nonempty peeling and retain detailed
state/output observations. No successful applied regression is added unless
actual full-path application is observed and repeatable.

All 128 inputs completed successfully: 15 naturally converged and entered
polish, 46 stopped at audit patience, 48 at all-rejected and 19 at maximum
iterations. Each of the 15 polish calls attempted one component; none accepted
one, recertified a candidate or applied a polished state. A focused replay of
these same inputs supplies correction-status observations; it is not an
additional input search. Fourteen corrections returned NoMaterialChange; the
remaining candidate (Input010) failed the global improvement/unavailable
objective check before member checks and recertification. The result is **no applied case established within the
128-input search**, not proof of impossibility.

## Observation and build isolation

Eight ON/OFF configurations cover audit trace OFF/ON crossed with trust-model
experiment OFF/ON. Baseline Ninja objects are rebuilt from the current root;
isolated libraries and executables reuse unchanged objects and OFF recompiles
only IterationProcess. Their own rpaths resolve their own policy libraries.
Exact compile/link commands are retained. No public strategy switch is added.

Observer libraries are separate from these uninstrumented builds. They capture
selected/rejected keys, states, selection objectives, finalization base,
background response cache, native/common objectives, operator certificates,
actual application and post-persistence Gaussian/peeling/MSE. The common
objective uses the initial objective domain and rebuilds background from the
terminal Gaussian state. It is not a replacement for comparing persisted
background and peeling. Measurements do not feed algorithm decisions.

## Existing test validation

All eight complete existing CTest configurations passed. This does not include
passing the external fold-168 numerical dataset: its runner test is one of the
16 CTest groups, whereas the actual dataset gate failed as described below.

| Variant | Audit trace | Trust-model experiment | Existing CTest |
| --- | --- | --- | --- |
| ON | OFF | OFF | 16/16 |
| ON | OFF | ON | 16/16 |
| ON | ON | OFF | 16/16 |
| ON | ON | ON | 16/16 |
| OFF | OFF | OFF | 16/16 |
| OFF | OFF | ON | 16/16 |
| OFF | ON | OFF | 16/16 |
| OFF | ON | ON | 16/16 |

The final documentation core-contract check passed (1/1), `git diff --check`
passed, and source verification confirmed unchanged production source/tests.
Actual dynamic-library loading was checked for the isolated OFF executable.

## Workload results and final disposition

Case 1's unmodified trace-OFF/four-thread run completed successfully in
325.97 seconds under concurrent validation load (not a benchmark). It reached
100 accepted iterations. The existing fold-168 quality metrics passed, but its
unchanged accepted-iteration gate failed: expected <=25, observed 100. This is a
baseline regression, not an effect of disabling polish. No tolerance or baseline
was changed. `case1-on-j4-trace-off/fold-report.json` retains the exact result.

The instrumented ON case 1 run reproduced every database table exactly, including
Gaussian, grouping and peeling data; only row order was canonicalized. No database
columns were excluded. The dependent workload matrix was stopped at this failed baseline gate.
Production final polish is retained; no deletion or API change was made.

| Workload execution | Status | Final-polish conclusion |
| --- | --- | --- |
| Case 1 ON, four threads, trace OFF | Completed, uninstrumented and observer outputs identical | MaximumIterations; zero final-polish calls |
| Case 1 OFF four-thread capture, ON/OFF one-thread captures | Interrupted at baseline gate | Incomplete; no application conclusion |
| Cases 2/3 ON, four threads | Interrupted at baseline gate | Incomplete; no application conclusion |
| Remaining ON/OFF, thread and trace workload combinations | Not run | No evidence |
| Warmups and three-repetition performance matrix | Not run | No speedup claim |
| Production removal and post-removal verification | Not performed | Deletion conditions unmet |

Case 1's common objective was 0.037123982262462565, native final audit
0.66114841232488764 and selected-response MSE 7.2427410257401235. These are
baseline observations, not ON/OFF improvements. Its final nominal operator p99
vector was `(0.0063121282, 0.0084212314, 0.0010307997)`; operator completeness
was true, solver qualification and production convergence were false.

The regular terminal summary reports `final_uses_polish=true` for this run,
but that field includes outer accepted-state polish provenance. It does **not**
mean final dependency polish was applied. The explicit observer records zero
calls and the MaximumIterations branch bypasses final dependency polish. Do not
use the summary provenance flag as the application counter.

The unmodified fold-168 metrics were amplitude RMSE 0.01071979435, width RMSE
0.00034471244, offset RMSE 1.23698189587 and maximum absolute offset
0.90980697564. Their existing quality gates passed; only the iteration-count
gate failed. `fold-gate.py` invokes the existing runner's input validation,
metric calculations and gate functions against the retained database, changing
only its lookup key to this experiment's saved key. The executed command is
retained separately from the runner's reference command template.

This is a blocked experiment, not a completed zero-application decision for all
three datasets. Resolving the baseline regression is separate algorithm work;
neither accepting a larger iteration limit nor deleting final polish resolves
that failure under this plan.

## Reproduction and evidence preservation

The ignored `build/final-polish-production/` directory retains:

- `prepare.py`, `baseline.json`, ON/OFF source archives and source verification;
- `matrix.py`, per-configuration compiler/linker commands, CTest logs/results;
- `capture-build.py`, observer source copies and instrumented executables;
- `search.py`, all 128 input specifications, raw state/persistence records,
  focused replay and candidate-rejection capture;
- `run-data.py`, actual commands, completed output databases and explicitly
  interrupted run records; `compare-db.py` and observer equivalence results;
- `fold-gate.py`, unchanged-gate results, `summarize.py`, `summary.json`,
  `verify-evidence.py`, object hashes and the artifact manifest.

Reproduction order is prepare, matrix build/CTest, observer builds, bounded
search, uninstrumented case 1, observer case 1 and database equivalence, then
fold-gate evaluation. `run-all-data.py` refuses to start the dependent matrix
when the recorded baseline gate fails. Existing run outputs are never overwritten;
interrupted runs must be preserved and a fresh output directory used on resumption.
The build-local archives/tools must accompany this report to preserve raw evidence.

Production `src`, `include` and `tests` are unchanged from the baseline. No
successful applied fixture was found, so no new tracked test was added. The
result does not support deleting the mechanism and does not certify successful
recertification, salvage of an accepted final patch or applied persistence.
