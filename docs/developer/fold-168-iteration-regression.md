# fold-168 iteration-count regression investigation

Investigation date: 2026-09-12. Production source was fixed at
`6587542638d570685a11ef281aebd015d5eb8c06`.
This investigation added only this report to the production tree; it did not modify
production algorithms, the runner, tests, quality baselines, or tolerances.
Raw evidence, isolated sources, builds, and per-iteration data are preserved in
[`build/fold-168-regression-investigation/`](../../build/fold-168-regression-investigation/).
That directory is not tracked by Git. Reproduction requires retaining it and the external inputs.
The investigation completed 21 full dataset runs, covering the official runner,
historical versions, instrumentation, and single-change interventions.
`run-index.json` summarizes their results. Failed builds are retained separately
and are not counted as dataset runs.

## Conclusions and strength of evidence

The current fold-168 run accepts 100 iterations and stops with `maximum-iterations`.
All quality and atom-cutoff gates pass; only `accepted_iterations <= 25` fails.
This result predates the recent objective deduplication and observation-facade
refactors: the rebuilt `fb926a5b` and current version produce exactly the same
168 persisted local-potential records, including peeling blobs, after excluding
the storage key.

A historical passing endpoint under the current contract was found, and the
boundary was narrowed to adjacent commits: `f50a742c` accepts 11 iterations;
its child `49d3516a` accepts 100, with best iteration 27.
`49d3516a` changes both quarantine and cluster-history policy. Restoring only the
historical-best acceptance gate for local candidates and ordinary boundary members,
while retaining its new quarantine policy, all stopping conditions, and tolerances,
returns the run to 11 iterations with every original gate passing.
This is intervention evidence that removing the member-best gate sends this dataset
into a trajectory without net progress. The passing result was not manufactured
by weakening the certificate or stopping earlier.

The mechanism that sustains the run to 100 is a permitted small objective increase
followed by improvement relative to the previous state. Those improvements repeatedly
reset patience without improving the historical best. The unrestricted operator
residual and solver qualification both fail for this dataset, so the certificate
correctly refuses to declare fixed-point convergence. Best iteration 27 does not
imply that the run should stop at iteration 28.

Disabling ordinary joint correction alone still yields 100 accepted iterations,
although the cycle changes. The experiment therefore rules out correction as the
sole cause. This report does not recommend restoring the entire old algorithm,
deleting rescue, or relaxing the 25-iteration gate.

## Fixed conditions and reproduction

| Item | Fixed value |
| --- | --- |
| CIF | `/Users/yslian/Documents/simulation/fold_test_model_0.cif` |
| CIF SHA-256 | `156d35aa326f0d4408d726a999329d2ffede775489aeaa5d99a2cc9b9f663cab` |
| map | `/Users/yslian/Documents/simulation/sim_map_gaus_grid0.10_charge1_bw0.50.map` |
| map SHA-256 | `5e0dbb13fc3a76f8a944e6e2b18393d1896fafc2ec9020457cca8e8a421f120e` |
| Host | macOS 26.6.2 / arm64 |
| Compiler | Apple clang 21.0.0 (`clang-2100.1.1.101`) |
| CMake | 4.4.3 / Ninja |
| Build | Debug, SYSTEM dependencies, UMAP OFF, Python OFF |
| Diagnostic switches | SECOND_STAGE_AUDIT_TRACE OFF, TRUST_MODEL_EXPERIMENT OFF |
| Execution | Fresh SQLite database, separate run directory, `-j 4`; explicitly identified serial comparisons use `-j 1` |

Installed package versions are recorded in `dependency-versions.txt`: Eigen 5.0.1,
Boost 1.92.0, GSL 2.8, and libomp 23.1.0. ROOT is 6.40.04. The binary actually links
to system SQLite 3.51.0, not the Homebrew SQLite 3.53.4 installed on the same host;
see `linked-dependencies.json`.
The authoritative resolved paths and complete build options are in each build's
`*-build/CMakeCache.txt`, `compile_commands.json`, and `build.json`.
`current-build-provenance.json` and `current-CMakeCache.txt` describe the current
official executable. Each isolated build has an `executable.sha256`;
`binary-hashes.json` records SHA-256 hashes of both the CLI and `librhbm_gem.dylib`.

Run the current official runner from the repository root:

```sh
python3 tests/integration/fold_168_regression.py \
  --executable build/objective-evidence-tests/bin/RHBM-GEM \
  --model /Users/yslian/Documents/simulation/fold_test_model_0.cif \
  --map /Users/yslian/Documents/simulation/sim_map_gaus_grid0.10_charge1_bw0.50.map \
  --baseline tests/benchmarks/fold_168_simulation_baseline.json \
  --output-dir build/fold-168-regression-investigation/65875426-official
```

The CLI exits 0 and the runner exits 1. The latter is an iteration-gate failure,
not a crash. The official `run.log`, `actual.json`, and `report.json` are preserved.
Because the official runner removes its temporary database, the same command was
rerun as `65875426-persisted-j4` to retain the database without overwriting the
official evidence.

Isolated versions use `manage.py` to archive Git sources, build, and run, for example:

```sh
python3 build/fold-168-regression-investigation/manage.py build f50a742c
python3 build/fold-168-regression-investigation/manage.py run f50a742c new-run-id
```

Each run's `command.json` records the exact full command and input hashes.
`result.json` records the exit status, elapsed seconds, and parsing or gate failures.
The database for a new run ID must not already exist.
CLI commands for versions using the current contract are equivalent to:

```sh
RHBM-GEM potential_analysis -d NEW.sqlite \
  -a /Users/yslian/Documents/simulation/fold_test_model_0.cif \
  -m /Users/yslian/Documents/simulation/sim_map_gaus_grid0.10_charge1_bw0.50.map \
  -k fold_gaus_charge1 -j 4 -v 3 --simulation true -r 0.50 \
  --exclude-hydrogen true
```

`fb926a5b-legacy-command` deliberately retains `-k final_polish_case` and
`--only-backbone false` from the historical record.
`PotentialAnalysisRequest::only_backbone` already defaults to false, and the storage
key does not participate in numerical solving. The rebuilt historical command and
the current canonical command produce identical persisted results after excluding
the key. Their textual command differences must not be mistaken for the cause of
the iteration regression.

## Historical baseline and version boundary

`<=25` is a runner constant. The JSON does not store a measured baseline iteration
count. `ecd6bddd` introduced a limit of 10; `28eb8170` changed it to 25 when adding
the 10-residue cluster cutoff. That commit does not include verifiable per-iteration
measurements. This investigation reproduced 22 iterations with all historical gates
passing. That establishes that the old version can finish within 25 iterations,
but does not by itself validate the requirement under the current contract.

The contract evolved as follows:

- `28eb8170` uses a residue cutoff: 20 residues, limit 10, 3 clusters, and at most
  8 residues in a cluster.
- `b2c0bc5c` changes to an atom-count cutoff. The current gate requires 168 atoms,
  limit 100, at least 2 clusters, and at most 100 atoms per cluster.
- `1227c62b` separates signal and tail ranges and removes `--fit-max 1.0` from the
  official CLI command. The baseline JSON retains historical command metadata;
  the runner explicitly permits this difference. Older versions in this investigation
  use their valid historical commands with `--fit-max 1.0`, rather than being forced
  to accept the new command.
- The database gained results separated by fitting stage and later lost the third
  stage. Historical readers use the corresponding runner version. `adapt.py` only
  reads data and classifies gates; it changes no data, values, or limits.

| Version / run | accepted | best | stop reason | Quality / cutoff |
| --- | ---: | ---: | --- | --- |
| `28eb8170-historical` | 22 | 19 | audit-patience | Historical contract passes; current atom cutoff is not applicable |
| `before-patience-reference` (`dea92b67`) | 7 | 4 | audit-patience | Quality passes; historical residue contract |
| `0b495899-historical` | 7 | 4 | audit-patience | Quality passes; historical residue contract |
| `d55d7950-adjacent` | 7 | 4 | audit-patience | Quality passes; historical residue contract |
| `eb7a6daf-bisect` | 12 | 5 | all-rejected-backtracking-exhausted | All current gates pass |
| `7b54b65d-bisect` | 8 | 5 | audit-patience | All current gates pass |
| `53b4a4d7-bisect` | 11 | 11 | all-rejected-backtracking-exhausted | All current gates pass |
| `f50a742c-bisect` | 11 | 11 | all-rejected-backtracking-exhausted | All current gates pass |
| `49d3516a-bisect` | 100 | 27 | maximum-iterations | Quality / cutoff pass; iteration count fails |
| `ddc16505-bisect` | 100 | 27 | maximum-iterations | Quality / cutoff pass; iteration count fails |
| `fb926a5b-legacy-command` | 100 | 27 | maximum-iterations | Quality / cutoff pass; iteration count fails |
| `65875426-official` | 100 | 27 | maximum-iterations | Quality / cutoff pass; iteration count fails |

`0b495899` changes the improvement reference used by patience from best to previous.
Both adjacent versions reproduce 7 iterations with identical quality. This change
is therefore a pre-existing condition that allows the current cycle to persist,
not independently the first failing fold commit. Bisection followed first-parent
history and ended with a direct `f50a742c -> 49d3516a` comparison, rather than an
inference from commit titles or dates.

Unmodified `28eb8170` fails to compile with the current C++20 toolchain because
`path.u8string()` is incompatible with `std::string`. The original log and manifest
are preserved as `original-build.*`. The 22-iteration result uses only the
`-fno-char8_t` compatibility flag, without changing numerical source code.
Historical schema parsing failures remain in `result.json`; successful adaptation
is recorded separately in `historical-actual.json` and `adapted-report.json`.
Build or parsing failures were never classified as iteration regressions.
An initial detail logger used `std::osyncstream`, which is unavailable in the local
libc++. Those failed builds are preserved as `initial-failed-build.*`; replacing it
with a logger-only mutex fixed the build without changing solver code.
A passing version was found, so the investigation did not expand to rebuilding
`ecd6bddd` or running other large matrices.

## First divergence: local search at attempt 6

The adjacent sources were rerun as `before-detail-j4` (`f50a742c`) and
`after-detail-j4` (`49d3516a`). Only the first 12 attempts record previous, proposal,
and assembled parameters alongside existing local-evaluation evidence.
This is lightweight phase evidence for the narrowed path; it does not replay
solvers or objectives. Raw logs, `states.json`, and `local-events.json` retain
unrounded values.

The three Gaussian parameters of the previous, proposal, and assembled states
match for the first 5 attempts. At attempt 6, previous and proposal Gaussian
parameters still match exactly, but assembled parameters diverge for the first time.
Proposals begin to diverge at attempt 7 as their input states differ.
At attempt 12, selection rejects everything in the old version, leaving 11 accepted
iterations; the new version continues accepting candidates.
Matching records by key, candidate objective, and previous objective identifies
two local-search acceptance decisions that change at attempt 6:

| Cluster (internal zero-based atom indices) | previous | candidate | Best in the same overlay |
| --- | ---: | ---: | ---: |
| `[150,151,156–167]` | 0.022762934892239843 | 0.022779995116691074 | 0.02274849495238948 |
| `[0–20,23–29]` | 0.14939603980489818 | 0.14954517906583784 | 0.1493551576270224 |

For the first cluster, candidate−previous is `1.7060224451231204e-5`, below the
previous bound of `2.2772934892239844e-5`. However, candidate−best is
`3.150016430159344e-5`, above the best bound of `2.2758494952389483e-5`.
The second cluster is also within the previous bound but outside the best bound.
Candidate, previous, and best values match exactly on both sides. The old version
reports `accepted=0, previous_rejected=0, best_rejected=1`; the new version reports
`accepted=1, previous_rejected=0, best_rejected=0`.
See `matched-local-decision-differences.json`.

This locates the earliest divergence at removal of the local-search best gate,
rather than the initial solve, proposal, boundary correction, or floating-point
error emerging after iteration 27. The single-gate intervention retains the new
quarantine policy and still returns to 11 iterations, further ruling out the claim
that the quarantine change in the same commit is the sole cause. This does not
establish that quarantine is irrelevant for other inputs or combinations of changes.

## Why iteration continues after 27 and reaches 100

`current-diag-j4/iterations.json` retains model parameters, candidate sources, and
lightweight decision evidence for all 100 iterations. `iterations.csv` is the
per-iteration scalar table. `current-complete-j4` additionally records each key's
radius before and after the iteration, plus boundary previous, endpoint, correction
reference, correction, and accepted objectives, all read from existing values.
Instrumentation is isolated in `current-diag-source`; it reads existing evidence
without adding objective evaluations or solves.

Actual progress around iterations 25 and 27 is as follows, all in domain 2:

| accepted | candidate objective | best after | Patience before→after | Interpretation |
| ---: | ---: | ---: | --- | --- |
| 25 | 0.6611559895198077 | 0.6611486047610075 | 0→0 | Improves previous, but not best |
| 26 | 0.6611496661641321 | 0.6611486047610075 | 0→0 | Improves previous, but not best |
| 27 | 0.6611484123248916 | 0.6611484123248916 | 0→0 | Best genuinely improves again |
| 28 | 0.6611790738490916 | 0.6611484123248916 | 0→1 | Permitted increase; no reset |
| 29 | 0.6611591734906526 | 0.6611484123248916 | 1→0 | The decrease after that increase resets patience again |

Throughout these iterations, the operator is complete, the solver is not qualified,
and StrictOperatorPassed=false. Iteration 25 is only the runner's external gate;
the solver has not met its own stopping conditions at that point. The actual
improvement at iteration 27 also shows why the 25-iteration limit cannot serve as
a retrospective proof that stopping was justified.

Across the 73 iterations from 28 through 100:

| State | Observation |
| --- | --- |
| Objective domain / recovery revision | Both remain 2; no pending partition |
| Best updates | 0; best remains iteration 27's `0.6611484123248916` |
| Strict improvement relative to previous | 58 iterations |
| Patience after | 0 in 58 iterations, 1 in 15; never reaches 3 |
| Radius / quarantine / domain reset predicates | All false |
| Operator completeness | Always true |
| Solver qualification / strict operator passed | Always false |
| Rejected-cluster blocker | True in 30 iterations, false in 43 |
| Suspicious fallback / quarantine transition | Always false |

The complete radius records show that the last radius change occurs at iteration 19.
No key changes radius during iterations 28–100. Of the final 8 keys, 7 have radius
0.0625; the key starting at index 127 has radius 1. This is stronger evidence than
checking only the rejected-radius patience predicate.
`complete-trace-verification.json` records all changes and verifies exact equality
of the existing trace fields across 100 iterations and of all 168 persisted records.

Late iterations exhibit an almost exact five-iteration cycle. The maximum absolute
model-parameter difference between iterations 100 and 95 is
`1.4210854715202004e-14`. All objectives below belong to the same domain:

| accepted | candidate − historical best | Improves previous | patience after | boundary source |
| ---: | ---: | ---: | ---: | --- |
| 63 | 3.067838386217048e-5 | No | 1 | ordinary correction |
| 64 | 1.0767835740610465e-5 | Yes | 0 | ordinary correction; rescue exhausted |
| 65 | 7.60084385997839e-6 | Yes | 0 | ordinary correction; rescue exhausted |
| 66 | 1.2947598958534812e-6 | Yes | 0 | ordinary endpoint; rescue backtracking |
| 67 | 1.4730452635447477e-8 | Yes | 0 | ordinary backtracking |

At iteration 28, the objective rises from best to `0.6611790738490916`, then falls
again. The current progress tolerance is `1e-8 + 1e-3 * abs(reference)`, approximately
`6.61e-4` at this objective scale. The cycle's increase is approximately `3.07e-5`,
so it can be accepted within the existing non-regression bound. The strict tolerance
is `1e-10 + 1e-8 * abs(reference)`, so the descending part of the cycle qualifies as
strict improvement relative to previous and repeatedly resets patience.
Strict improvement of ordinary correction relative to its endpoint does not imply
improvement relative to historical best. Comparisons against these different
references cannot substitute for one another.

At iteration 28, boundary previous=`0.6611484123248876`,
endpoint / correction reference=`0.6611930309741361`, and
correction / accepted=`0.6611790738490916`. Correction strictly improves the endpoint
but does not improve previous; the combined previous bound permits that increase.
These are valid decisions against different references. There is no evidence here
of an incorrectly computed objective.

At iteration 100, accepted p99 for amplitude, width, and offset is
`[2.2949786309283168e-5, 1.166403857657765e-5, 1.1450518225355281e-5]`.
Although these values are below `1e-4`, nominal operator p99 is
`[0.006312128227289602, 0.008421231387339852, 0.001030799682985109]`.
Solver-qualified is also false. Small accepted movement does not establish that
the original operator has reached a fixed point.

The per-iteration trace can still show `stop=0` at iteration 100 because
maximum-iterations is assigned after the outer loop finishes. The final summary
is authoritative for the stopping reason; the logger's position does not create
a contradiction. Background is built each iteration, but this dataset has no
unselected atoms. Revision and reset evidence both show that late background or
partition changes are not the source of patience resets.

## Counterfactual experiments and parallel / observation checks

| Isolated change or comparison | Result | Supported conclusion |
| --- | --- | --- |
| Skip only ordinary joint correction in the current version | 100 iterations, best 16, maximum-iterations; quality / cutoff pass | Correction is not the sole cause; deleting this path is not a repair |
| Restore only the local / ordinary member-best gate in `49d3516a` | 11 iterations, best 11, all-rejected-backtracking-exhausted; all gates pass | Removing the member-best gate has a causal role in this version's regression |
| Restored gate with `-j 1` vs `-j 4` | Both accept 11; actual results and all 168 persisted local-potential records match exactly | The intervention result does not depend on worker count in this comparison |
| Current lightweight instrumentation with `-j 1` vs `-j 4` | All 100 iterations of trace match exactly | Worker-count differences do not cause this reproduction |

The gate-restoration source is `49-restore-member-best-source`. It changes only
local objective acceptance in `CandidateEvaluation.cpp`: it reevaluates the
historical-best cluster parameters in the same candidate neighbor overlay and
retains unavailable handling and the original best-deterioration gate.
Quarantine, patience, certificate, global audit, rescue, radius, background, and
the iteration limit remain unchanged. This historical intervention is not a
production repair submission.

The uninstrumented current persistent run, current diagnostic j4 run, and rebuilt
fb legacy run have exactly the same 168 `model_atom_local_potential` records after
excluding `key_tag`. The SHA-256 of the concatenated peeling blobs is
`62dcdce7dba4e6e3b2f27b6bef069693c238f1c15592363ccbfc16367fa75273`.
Atoms, quality, summary, and cutoff in the current official `actual.json` also match
the diagnostic version exactly. Complete comparisons are preserved in
`observation-equivalence.json` and `additional-equivalence.json`, including current
j1/j4, intervention j1/j4, and instrumented versus uninstrumented adjacent commits.
Actual results and persisted local-potential records also match exactly between
`49d3516a` and the current version. Wall time and nonnumerical log ordering are excluded.

The original quality requirements remain unchanged: 168 serial IDs, finite and
valid parameters, every quality metric at most 105% of its reference, and the
applicable cutoff contract. The current four metrics are:

| metric | reference | Current actual |
| --- | ---: | ---: |
| amplitude RMSE | 0.0538959080312383 | 0.010719794351592081 |
| width RMSE | 0.0013150213867548584 | 0.0003447124400131276 |
| offset RMSE | 1.234298853386113 | 1.236981895870946 |
| maximum absolute offset | 0.9086740100107469 | 0.9098069756400031 |

## Minimal repair proposal and subsequent validation

The proposed minimal repair candidate is to restore the necessary member-best
reference for acceptance decisions. Changing patience or selecting iteration 27
earlier is not presented as a validated repair. The historical single-change
intervention removes the regression, but porting it to current ownership and
observer interfaces still requires separate implementation and validation.

1. Restore the original best non-regression gate, evaluated in the same overlay
   and domain, in the shared local-objective and ordinary-boundary-member path in
   `CandidateEvaluation.cpp`. Retain strict local-polish prerequisites and previous
   gates. Do not apply it to cooperative members or infer strict improvement against
   one reference from improvement against another.
2. Pass the decision reference through `LocalCandidateReference` and the candidate
   transaction. Keep the minimal cluster-best parameters and domain lifecycle in
   `IterationProcess` / transaction state. Initialization, partition and background
   resets, candidate staging, and rollback must be explicit. Production acceptance
   must not read `ClusterHistoryObserver`, which can be disabled, fail, or depend on
   Debug conditions. The observer remains observational. Share the mathematical
   reference-evaluation code without moving the entire history observer back into
   the main flow.
3. Reevaluating a historical reference must preserve the candidate's other cluster
   and neighbor parameters and replace only the target key's best patch. Do not use
   a stored scalar best across different backgrounds. Continue sharing correction's
   precomputed raw objective, and count added reference-evaluation work accurately.

The expected effect is to reject the first candidate identified here that exceeds
the member-best bound and follow existing backtracking, fallback, or all-rejected
paths, changing the acceptance trajectory at its source. It must not label an
unconverged state as converged. Risks include additional objective-reference work,
more conservative progress in some cooperative scenarios, and reintroduced state
requiring correct domain and transaction management. Quality cannot be assumed
unchanged for every dataset.

Subsequent validation must include the current fold runner's original gates,
persisted Gaussian and peeling data, j1/j4 comparisons, existing defense tests under
all four observer-switch combinations, and unavailable, domain-reset, rollback,
rescue, and intensity-scaling paths. First verify that the initial divergence in
this report is removed, then require both <=25 iterations and all quality gates.
Iteration count alone is insufficient. If the port cannot satisfy both, retain it
as a candidate and investigate search or trust adaptation; do not conceal the
failure by raising the limit.

## Limitations

No tests were added or modified, and no large ablation matrix was rerun.
Fold is one concrete input pair, not a guarantee of general numerical equivalence
or quality. Expensive phase replay was not enabled for all 100 iterations.
The investigation uses lightweight observation of existing evidence and focused
records of early candidate-path differences between adjacent commits.
It does not exhaustively explain every inner solver's qualification failure.
The observed operator residual is far above the threshold, which is sufficient
to reject the argument that small accepted movement alone should imply convergence.
