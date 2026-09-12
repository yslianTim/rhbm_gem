# Shared component assembly and discrete shrink levels

## Baseline and intent

Baseline: `ac9ca72ce0e5093abb8e2c97966c41644b494053` (2026-09-12).
This is a structural refactor in two stages, with the existing numerical policies
retained. No tests, fixtures, assertions, public options or ablation switches are
added or changed. Build-local capture tools reuse existing defense cases.

## Component infrastructure

`CouplingGraph` shares participant merging and root-to-key collection between
boundary and uncut dependency components. Each entry point still controls its
participant source and filtering. Boundary components require two accepted keys;
final components may contain one multi-atom key. Existing DSU union order, root
iteration order, final component sorting, physical halo and sample filtering stay
with their original callers.

The internal `ComponentAssembly` module provides ordered application of borrowed
patches and base-state assembly with an optional excluded component position.
Null patch entries represent unselected components. Outer patch application stays
inside `CandidateTransactionBuilder`; final polish maintains its own accepted
patch list. There is no new mutable selection owner.

`AuditAndSalvageComponents` runs one initial evaluation and then policy-selected
removals until the caller's gate passes or no removal remains. The removal callback
returns the resulting audit, so a precomputed final-polish removal objective is
reused without another evaluation.

| Policy | Outer | Final polish |
| --- | --- | --- |
| Evaluation | Existing overlay/delta path and previous/best gates | Existing complete-state snapshot/audit path |
| Removal candidates | Units that do not independently strictly improve previous | Every currently accepted component, tried one at a time |
| Order | Worst objective first; missing evidence scores infinity; ties by lexical keys | Best available removal objective each round; ties retain first traversal position |
| Stop | Passing audit, or existing exhausted fallback | Strict improvement over base, or no improving removal |
| Updates | Builder rejection order, provenance, diagnostics and history | Final component accepted flags and assembled state |

Ordinary/cooperative sweep ordering, correction/backtracking, missing-baseline
handling, remote-cluster fallback, quarantine/recovery revision, global-best gates,
convergence and final strict operator recertification are unchanged. Sharing the
loop does not replace greedy salvage with atomic fixed-order acceptance.

## Discrete radius storage

`TrustRegionStateSet` now stores `std::map<ClusterKey, unsigned int>` levels. Its
existing methods and `double GetRadius()` interface are unchanged.

| Level | Radius |
| --- | --- |
| 0 | 1.0 |
| 1 | 0.5 |
| 2 | 0.25 |
| 3 | 0.125 |
| 4 | 0.0625 |

All five values are exactly representable in binary. Incrementing the level is
identical to the previous reachable multiply-by-half sequence. Level 4 saturates;
`ResetToMinimum` sets 4 directly. Reconciliation retains surviving levels, drops
removed keys and starts new keys at 0. Accepted shrink precedes retryable rejection
shrink, including the original repeated-key behavior and exhausted-key exclusion.
Missing-key exceptions and changed/saturated output ordering are retained.
Diagnostics still report actual radii. Final polish still uses its independent
fixed radius of 1.0 on every round.

## Validation

Artifacts and build-local reproduction scripts are in `build/component-assembly/`.
Both production builds use Debug, existing system/cached dependencies and trust-model
experiment OFF. Trace OFF and ON use the existing configured build directories
`build/rescue-only/final-trace-off-trust-off` and
`build/rescue-only/final-trace-on-trust-off`.

- Baseline: full trace-OFF CTest, 16/16 groups passed.
- Stage 1: full trace-OFF CTest, 16/16 groups passed, before changing radius storage.
- Final: full trace-OFF CTest, 16/16 groups passed.
- Final: full trace-ON CTest, 16/16 groups passed.
- Documentation contract check: passed after the documentation changes.

Baseline, stage 1 and final each pass all 107 enabled existing defense cases in
the isolated capture runs. Both stage 1 and final match the baseline in all 11,066
selection/rescue events, 2,759 assembled states and 26,369 selected production/diagnostic
records across 29 stdout captures; `comparison.json` records both comparisons.
The comparison excludes elapsed-time fields and normalizes unordered worker-log
order; it does not round numerical values. Existing serial/parallel/intensity-scale
quality captures match exactly in persisted Gaussian parameters and response MSE.
The final persisted-model/peeling capture also matches the baseline exactly.

Reproduction: rebuild each configured directory with `cmake --build <directory> -j 6`,
then run `ctest --test-dir <directory> --output-on-failure -j 4`. The existing
capture tools adapted under `build/component-assembly/` preserve original test
assertions and keep instrumented copies/libraries outside tracked sources;
`compare.py` compares baseline, stage1 and final captures.

## Limits

The existing suites cover component grouping, boundary/cooperative search,
remote-cluster isolation, serial/parallel behavior, direct final-polish improvement,
trust shrink/saturation and Frozen retry. They do not establish exhaustive coverage
of all salvage tie, missing-evidence and post-rescue rollback combinations, or a
production final-polish applied-path frequency. Persisted peeling is compared on
the captured existing persistence fixture, not every defense fixture. Repeated
radius keys and unchanged exception semantics are also reviewed directly in the
unchanged update ordering; no new edge-case tests were introduced.

No paired 600-case corpus, fold-168 numerical dataset, real-data quality benchmark
or performance claim is made. The existing fold-168 runner CTest validates the
runner, not that numerical dataset. The trust-model-ON matrix is outside this run.
