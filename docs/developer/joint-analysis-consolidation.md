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
