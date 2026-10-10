# Joint sparse factorization backend

This document describes the shared sparse numerical backend used by the
current FixedNeighbor Joint estimator and its local Profile LM search
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
local update per block visit, and local Profile LM search. The
existing KKT, width-gradient and eta-confirmation thresholds are unchanged.

## Qualification evidence

The one-time EIGEN/SPQR qualification and its superseded wall-time-first choice
are summarized in the [historical decision record](joint-operator-search.md).
Both backends passed bounded numerical parity; SPQR was faster on completed
128/256 probes, EIGEN consistently used less peak RSS, and both bounded 512
probes timed out. Production later prioritized memory footprint, dependency
simplicity and maintenance, making EIGEN the sole current backend. This is a
historical decision, not a current build or regression selector. Historical
saved-result metadata can still decode SPQR provenance; current v3 provenance
requires `sparse_backend=EIGEN`.

## Historical rank and operator work

The former bounded SPQR rank certificate, OperatorPcg policy, preconditioner
and Schwarz experiments were research infrastructure. Their source drivers,
current route hooks and machine-readable campaign outputs are removed. They
must not be reintroduced as hidden fallbacks or current benchmark options;
their decisions are summarized in the
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
