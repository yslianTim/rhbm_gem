# Joint analysis consolidation acceptance

Baseline: `c5619e0c` on `develop`. Changes are committed in S1–S5 order,
with each stage validated before the next one begins.

## Baseline

Release, EIGEN, system dependencies, ROOT/UMAP and Python bindings enabled:
`tests_all` built and all 25 default CTests passed. Dedicated build directory:
`/private/tmp/joint-consolidation-build`. Existing Debug/Release presets are not
used as test builds because they have testing disabled.

## S1 — provenance and invalidation

Reports derive estimator/peeling provenance from their actual stage population,
including mixed/unknown sources. Removing a diagnostic snapshot preserves the
published stage. UMAP uses the rows actually admitted to the embedding.

Analysis editor updates invalidate downstream state synchronously. Endpoint
replacement invalidates the affected joint runs (including halo-dependent
peeling) and groups; covariance changes invalidate evidence and group results;
raw replacement clears peeling without changing fixed points or parameter-domain
group results. Group membership/configuration changes clear group results.
Complete endpoint mappings are validated before publication. Standalone First
writeback retains Second/group results while invalidating sample-derived peeling.

Tests cover no-snapshot joint reports/UMAP, mixed sources, covariance changes under
the same run ID, cross-atom invalidation, same-size raw reordering, and invalid
batch mappings. An explicitly reattached mismatching diagnostic snapshot remains
a persistence error; an ordinary stage update now clears it automatically.

S1 validation: `tests_all` built; all 25 default CTests passed with ROOT, UMAP
and Python enabled. Frozen Joint numerical expectations and tolerances are
unchanged. `git diff --check` passed.

## S2 — canonical in-memory results

Each stage stores one published point, separate OLS/MDPDE uncertainty and native
fit diagnostics. Legacy result DTOs are value projections. A non-final seed or
invalid native outcome may be retained only while no final point is published.
Transient fitting matrices can be cleared without losing durable diagnostics.

Joint peeling stores paired responses only; compact samples, counts and ratios
are projections from raw samples and paired coverage. Independent compact-sample
or neighbor-count writes to a Joint stage are rejected before mutation.

Atom records own posterior results. Stored group summaries contain statistics
and member IDs, with member posterior DTOs assembled by the public view. Legacy
v1 JSON is still emitted as a compatibility projection until S5. Its duplicate
posterior fields are checked when decoding rather than retained in memory.

S2 validation: `tests_all` built and all 25 default CTests passed (ROOT/UMAP and
Python enabled). Added tests verify getter coherence after replacement, durable
diagnostics after transient clearing, rejection of independent Joint peeling
writes, copy consistency, and detached group projections. Existing persistence,
two-stage and frozen Joint regressions passed without tolerance changes.

## S3 — bounded postprocessing and serialization

Metadata updates modify only the diagnostic metadata. The workflow releases the
solver result before postprocessing and moves its captured snapshot into analysis
after the last uncertainty read. Saving encodes the joint document once.

Postprocessing defaults to recorded targets; an explicit output index span can
request other contributors, and inputs without selection metadata retain all-atom
behavior. Predictions and full-component covariance factorization retain every
contributor. Only target own-contribution maps and requested marginal covariance
blocks are materialized. Uncertainty reuses the immutable problem partition and
indexes support memberships by tile in the original atom/support order.

The opt-in `joint_postprocessing_benchmark` target and
`tests/integration/joint_postprocessing_benchmark.py` measure three independent
serial processes per full/halo/multi-component case, separately for complete
workflow+save and fixed-endpoint postprocessing+save. Synthetic maps use 0.16 A
spacing so component Jacobians span multiple 8192-row tiles. Process peak RSS
includes fixture construction; phase times exclude it. These bounded fixtures do
not establish a maximum supported problem size or general speedup guarantee.

Measured medians (seconds, MiB; baseline S2 versus S3 on the same host):

| Case | Total before / after | Peak RSS before / after |
| --- | --- | --- |
| full-workflow | 0.1440 / 0.1434 | 50.03 / 50.23 |
| full-post | 0.1026 / 0.1005 | 49.97 / 48.69 |
| halo-workflow | 0.1784 / 0.1706 | 50.20 / 48.77 |
| halo-post | 0.1157 / 0.1131 | 50.31 / 48.94 |
| multi-workflow | 0.3773 / 0.3544 | 89.33 / 89.09 |
| multi-post | 0.3037 / 0.2785 | 91.02 / 86.98 |

All 18 paired runs retained exactly equal endpoint A/C/B, objective, runtime
convergence and target peeling/covariance values. Raw runs, phase times and
diagnostic hashes are recorded in [the measurement report](joint-analysis-consolidation-benchmark.json).
Saving and halo peeling improved on these fixtures; uncertainty phase time was
slightly higher despite removal of repeated scans. Full-selection peak RSS was
essentially unchanged. No broader performance claim is inferred.

S3 validation: `tests_all` built; all 25 default CTests passed. Target-only and
explicit all-contributor postprocessing agree exactly on targets. Dense covariance
reference, coverage/negative contribution tests and frozen Joint regression pass.
CLI/persistence tests now explicitly require absent halo diagnostics. Numerical
tolerances are unchanged; `git diff --check` passed.

## S4 — workflow responsibilities

`PotentialFittingWorkflow.cpp` owns orchestration and `StageSummary.cpp` owns
reporting. Group training/inference declarations live in
`GroupPotentialFitting.hpp`, separate from uncertainty.

The explicit-workset First executor has two fixed modes: prepared-sample batch
training/fitting with the existing parallel behavior, and map-sampled contributor
initialization with per-atom exception isolation. Both use the same formal
First-atom fitting operation. Batch mode does not construct unused Joint
initialization diagnostics. The standalone wrapper still works on a model copy
and writes successful target First results back, preserving selection, Second,
group results and halo history; raw replacement invalidates peeling as in S1.

An additional test runs prepared-sample First on an unselected halo workset and
checks that the selected atom remains untouched, exactly one formal fit occurs,
and no sampling occurs. Existing exact First/direct-fit and standalone failure
isolation tests remain the numerical acceptance criteria.

S4 validation: `tests_all` built; all 25 default CTests passed, including
sampling, two-stage, standalone, CLI and Python binding/pipeline regression.
The explicit workset observer test passes; `git diff --check` passed.

## S5 — canonical persistence

SQLite v19 stores analysis only in neutral document v2. Final points, raw geometry,
paired Joint peeling, evidence and atom posteriors are written once. Native groups
store their statistics and IDs; Joint summaries contain no member posterior copy.
Alpha, OLS and its uncertainty, MDPDE uncertainty, native status/variance,
iteration and refinement diagnostics survive save/load and transient-state clearing.
Native failed-fit nonfinite diagnostics use explicit JSON tags; they remain
unavailable final points and UMAP keeps its existing exclusion behavior.
Legacy native uncertainty and diagnostics absent from the old format remain
unrecorded, rather than acquiring invented covariance or convergence information.

The writer no longer creates or writes the three legacy analysis tables, or the
old sample BLOB representation. v17/v18 adapters are restricted to load/migration.
First write converts every model key in a single transaction, drops replaced
tables, and saves the requested model. Conflicting duplicate points, samples,
paired responses or posterior projections are rejected. Conversion and subsequent
write failures restore both schema and all records. The checked-in SQL fixture
was generated by the actual v18 implementation at `df138da5`.

`SetJointResult` attaches diagnostic snapshots to existing stages; it no longer
relies on persistence to synthesize endpoints. `ApplyJointResult(result, run_id)`
is the explicit atomic endpoint-and-snapshot entry for C++ callers importing a
result. Saving checks snapshot point, contributor role, component and convergence
against canonical Second. Normal workflow already publishes the stages first.

S5 validation: `tests_all` built for all three configurations. EIGEN with
ROOT/UMAP ON and OFF each passed all 25 default CTests, including Python bindings
and pipeline examples. SPQR passed all 8 selected Joint/estimator/HRL numerical,
frozen regression, physical and CLI groups. Repository lint and installed C++
consumer execution passed. New acceptance tests cover exact persistent-document
round-trip (native diagnostics, Joint, mixed, unavailable and no snapshot), both
legacy versions remaining byte-identical on read, multi-key migration, conflicting
local/group/posterior representations, conversion rollback and database-write
rollback. The CLI smoke also asserts v2 and absence of the legacy tables/Joint
compact and member-posterior mirrors. No frozen expectations or numerical
tolerances were changed; `git diff --check` passed.

Final build directories: `/private/tmp/joint-consolidation-build` (EIGEN ON),
`/private/tmp/joint-consolidation-off` (EIGEN OFF), and
`/private/tmp/joint-consolidation-spqr` (SPQR OFF). Final run logs use
`/private/tmp/s5-*-complete-*`; these are local verification artifacts, not required
runtime files. The benchmark report above remains the S2-to-S3 comparison.
