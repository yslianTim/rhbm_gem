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
