# Second-stage decision audit

This is passive recording of evidence already produced by the algorithm.
Trust shadow strategies, full phase snapshots/replay, per-cluster historical
re-evaluation, and the six alternate coupling thresholds are retired. Production
trust radius, objective gates, salvage, quarantine, `best_audit_state`, convergence
and final operator recertification retain their decisions and ordering.

## Build and migration

The only extra-audit option is `RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT`, default OFF.
Library and private-header tests receive the same private definition through
`rhbm_gem_apply_observation_definitions`. ON works with `BUILD_TESTING=OFF`.
A defined retired option, including an OFF value left in a cache, is a configure
error. Remove old command arguments and clear their cache entries as instructed
by the error; no compatibility alias remains.

```sh
cmake -S . -B build/audit -DBUILD_PYTHON_BINDINGS=OFF \
  -DRHBM_GEM_ENABLE_SECOND_STAGE_AUDIT=ON
cmake --build build/audit -j
python3 resources/tools/developer/second_stage_audit.py run.log --output-dir audit-report
```

Extra recording also requires non-quiet execution and Debug verbosity. This gate
does not replace the existing raw Debug condition for parallel solver scheduling.
OFF, quiet and lower verbosity skip payload allocation, collection, formatting and
extra counters. Basic progress, warnings, final summary, atom-cutoff text and
elapsed time use scalar data independently of the audit.

## Observer and diagnostic contract

Observation is read-only with no decision return. It cannot request solver,
operator, objective, history or model-snapshot work. Gate evidence comes from the
original short-circuit branch, including which previous/best checks ran.
Session allocation, bounded collection and writer exceptions are contained and
disable extra recording for that run. No observation exception reaches the
production failure/fallback path. A missing terminal record therefore means an
incomplete diagnostic run, not evidence that fitting failed or succeeded.

Candidate evaluation owns the previous/best progress gate; its actual checked
references and reason are defined in `CandidateEvidence.hpp`. Objective evaluation
owns scoring, tolerance validation and scalar comparisons. Boundary correction
shares `BoundaryCandidateReference` and requires a non-null previous objective
reference; ordinary boundary evaluation can receive a missing previous objective.
An explicitly supplied unavailable delta remains unavailable and is never recomputed.

The numerical caller starts correction/final-polish observation with the existing
factor and trial/round number, then passes only numerical inputs to the candidate
reference. Session operations capture scores before best-state updates, retries,
convergence and final certification. Joint observers own rejection record fields
and reference labels. `Audit()` is read-only, and writable joint records are private;
logging and tests consume recorded evidence without changing it. Basic progress
diagnostics and performance counters retain their existing responsibilities.

The Logger line prefix is `Second-stage audit: schema=1, payload=` followed by
strict JSON. Start records contain version, thread count and fitting settings once.
Iteration records contain attempt/accepted counts, objective/recovery/background/
partition revisions, stage counts, committed selection, actual scores, quarantine
transitions and the outer convergence certificate. The operator reference is
`iteration_previous`; all-rejected branches that did not assess the certificate
say `not_evaluated`. Background revision identifies the frozen scoring environment;
a partition reset also advances that environment's revision.

Ordinary and after-rescue selection audits independently record `executed`, result,
reason, evaluation count, removed clusters and previous/candidate/best objectives.
Results are `skipped`, `passed`, `rejected`, `unavailable` or `empty_after_salvage`.
An optional final objective is never used to infer whether either audit executed.
Objectives have explicit scope/reference; missing scores are null with a reason.
Nonfinite numbers serialize as null with a reason, never NaN/Infinity or zero.

Normal successful trials contribute counts only. All abnormal categories share
five detail slots per attempt; finalization has its own five slots. Worker buffers
are bounded before merging, in stage/first-atom/key-size/trial order, independently
of completion order. All category totals remain, even for omitted details.
Component events identify their first member key (not an overlapping halo).
Each event contains an internal selected-atom key identity, factor/radius where
available, rejection/guard/solver or lifecycle reason and the evaluated references.
No full key/model snapshots or retained atom-by-atom output is stored.

Terminal records identify stop reason, final state source, best iteration, the
available score in the final state's frozen environment, and elapsed milliseconds.
Final polish reports attempted, objective accepted, operator certified and applied
separately. Certification status distinguishes not evaluated, failure and error;
its operator certificate does not fabricate an accepted-movement assessment.
The final logger performs no objective recomputation. `work_counters` contains,
in order, full state materializations, cache hits, cache misses, objective
samples recomputed, objective samples reused, and symbolic analyses.

## Analyzer

`second_stage_audit.py LOG --output-dir DIR` accepts one run in schema 1 only and
writes `audit.json` and `report.md`. It rejects invalid numbers, count mismatches,
more than five details, impossible applied-polish claims and mixed runs. Partial
logs produce an explicitly incomplete report. It has no replay, shadow or legacy
format branches. Old logs require the tools at their original Git revision; the
[retired guide and tools](https://github.com/yslianTim/rhbm_gem/blob/8bb7bf3a878c3c94fdbd6a2063bbc8c1b5552bec/docs/developer/second-stage-phase-audit.md)
remain historical references.

## Verification tiers

| Change | Required verification |
| --- | --- |
| This retirement | ON/OFF build and existing second-stage tests; small comparison to the pre-removal two-switch-OFF baseline; parser and documentation contracts |
| Text, fields or analyzer | Corresponding output/parser/documentation tests |
| Collection, lifetime or parallel recording | Small ON/OFF neutrality, bounded counts/order and allocation/collection/writer isolation |
| Solver, objective, acceptance, convergence, background or partition policy | Independent long numerical regression plus small ON/OFF neutrality |

`SecondStageNumericalProbe` captures actual commits, radii, quarantine, background,
Gaussian uncertainty/parameters and peeling samples in testing builds. Independent
instrumentation counts solver, operator, objective and model-snapshot entry calls;
new ON/OFF counts must match, not merely symbolic-analysis counts. Its code and
fault injection are absent from non-testing builds. Compare the three probe logs
with `python3 tests/integration/second_stage_neutrality_test.py BASELINE OFF ON`. Use the same toolchain,
threads and log level for baseline/OFF/ON comparisons. Tests also cover bounded
parallel merges, actual gate ordering, skipped/unavailable/empty audit states and
polish acceptance with failed certification.

For responsibility-only refactors, capture fresh OFF and ON baselines at the
pre-change commit. In addition to the existing comparator, require exact work
counts before/after (the comparator checks work only between OFF and ON), and
compare pre/post ON audit JSON after excluding version and elapsed time. Keep all
decision, score, stage, trial, ordering and count fields in that comparison.

External fold-168 regression is opt-in.
