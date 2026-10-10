# Joint sparse factorization backend

This document describes the shared sparse numerical backend used by the
current FixedNeighbor Joint estimator and its local `LegacyCompact` profile
solver. It is not an OperatorPcg or rank-prototype guide. Historical operator,
EIGEN comparison, and bounded-rank results are indexed in the
[canonical historical evidence](joint-component-evidence.md).

## Fixed production backend

The current Joint sparse factorization backend is EIGEN. It is the sole
production implementation; there is no sparse-backend selector, runtime
fallback, or automatic routing. The current build does not search for, link,
or install SuiteSparseQR, SPQR, CHOLMOD, or SuiteSparse.

```bash
cmake -S . -B build/joint-eigen -DCMAKE_BUILD_TYPE=Release \
  -DRHBM_GEM_ENABLE_UMAP=OFF
```

## Current numerical contract

`LinearSolve` performs the sparse solve with EIGEN `SparseQR` and
`COLAMDOrdering`; the compact reduction and derivative work feed the existing
local profile search. `StructuralPartition`, `CompactSvd` and
`TiledDerivative` provide the current partitioning and local numerical work.
Each width state performs its qualified EIGEN solve independently. No numeric
factor, active set or coefficient solution is carried across FixedNeighbor
sweeps.

The EIGEN implementation shares the estimator's observation domain,
parameterization, objective, threshold policy, endpoint certification and
`RuntimeConvergence` contract. Current configuration fingerprints and benchmark
metadata record `sparse_backend=EIGEN`.

FixedNeighbor's production numerical settings remain independent of backend
selection: core 12, Forward order, 30 maximum sweeps, one trusted accepted
local update per block visit, and local `LegacyCompact` profile search. The
existing KKT, width-gradient and eta-confirmation thresholds are unchanged.

## Qualification evidence

The one-time current FixedNeighbor qualification is retained as compact
promotion evidence in
[`joint-fixed-neighbor-backend-qualification-r1/`](figures/joint-fixed-neighbor-backend-qualification-r1/).
Small chain/cube, observable/nuisance, and partial-selection cases passed
with numerical parity. The completed 128/256 resource probes showed a
reproducible SPQR wall-time advantage, while EIGEN consistently used less
peak RSS; both backends timed out on the bounded 512 probes. The original
qualification selected SPQR under a wall-time-first policy. Production later
changed its policy to prioritize memory footprint, dependency simplicity,
installation portability and maintenance surface over that measured
wall-time benefit, so EIGEN is now the sole backend.

The qualification measurements and historical `chosen_backend=SPQR` decision
remain unchanged in the canonical evidence. They are retained as historical
readability and provenance, not as a current build or regression selector.
Historical saved-result metadata can still decode SPQR provenance; current v3
provenance requires `sparse_backend=EIGEN`.

## Historical rank and operator work

The former bounded SPQR rank certificate, OperatorPcg policy, preconditioner
and Schwarz experiments were research infrastructure. Their source drivers
and current route hooks are removed. They must not be reintroduced as hidden
fallbacks or current benchmark options. Compact summaries and figures remain
available for provenance in the
[historical evidence](joint-component-evidence.md) and
[retired OperatorPcg record](joint-operator-search.md).

This cleanup does not remove rank or identifiability evidence that is part of
the current endpoint assessment. It removes only the operator-only
rank/preconditioner infrastructure; the current estimator uses the EIGEN
sparse implementation.

## Current verification

Use the active benchmark wrapper for current FixedNeighbor measurements:

```bash
python3 tests/integration/joint_benchmark.py \
  --profile solve --case chain-8 --build-dir build/qualification \
  --output build/joint-solve.json
```

The benchmark wrapper has no sparse-backend, rank, operator, preconditioner or
Schwarz modes. Focused C++ tests cover the EIGEN linear solve, compact
assessment, prepared blocks, endpoint certification and FixedNeighbor
regressions. The historical EIGEN/SPQR comparison is not regenerated as a
current test.
