# Joint sparse factorization backend

This document describes the shared sparse numerical backend used by the
current FixedNeighbor Joint estimator and its local `LegacyCompact` profile
solver. It is not an OperatorPcg or rank-prototype guide. Historical operator,
EIGEN comparison, and bounded-rank results are indexed in the
[canonical historical evidence](joint-component-evidence.md).

## Fixed production backend

The current Joint sparse factorization backend is SPQR. SuiteSparseQR 4.x,
CHOLMOD and their transitive dependencies are required at configure time;
missing dependencies fail configuration. There is no Joint sparse-backend
selector, runtime fallback, or automatic routing.

```bash
cmake -S . -B build/joint-spqr -DCMAKE_BUILD_TYPE=Release \
  -DRHBM_GEM_ENABLE_UMAP=OFF
```

SPQR is GPL-2.0-or-later (alternate licenses are available from its author).
See the third-party notices. Installed shared builds need the linked runtime
libraries. Installed static builds resolve dependency targets through the
package configuration.

## Current numerical contract

The backend supplies the sparse factorization, compact reduction and derivative
work needed by the existing local profile search. `SparseFactor`,
`StructuralPartition`, `CompactSvd`, `TiledDerivative` and `LinearSolve` are
retained where they have current callers. Symbolic reuse is scoped to matching
matrix structure and policy; each new width state still performs fresh numeric
factorization and solve work. No numeric factor, active set or coefficient
solution is carried across FixedNeighbor sweeps.

The SPQR implementation shares the estimator's observation domain,
parameterization, objective, threshold policy, endpoint certification and
`RuntimeConvergence` contract. Configuration fingerprints and benchmark
metadata record `sparse_backend=SPQR` as scientific provenance, not as a
runtime choice.

FixedNeighbor's production numerical settings remain independent of backend
selection: core 12, Forward order, 30 maximum sweeps, one trusted accepted
local update per block visit, and local `LegacyCompact` profile search. The
existing KKT, width-gradient and eta-confirmation thresholds are unchanged.

## Qualification evidence

The one-time current FixedNeighbor qualification is retained as compact
promotion evidence in
[`joint-fixed-neighbor-backend-qualification-r1/`](figures/joint-fixed-neighbor-backend-qualification-r1/).
Small chain/cube, observable/nuisance, and partial-selection cases passed
with numerical parity. SPQR also showed a reproducible current-production
wall-time advantage on the bounded chain/cube resource probes, while using
more RSS; both backends timed out on the bounded 512 probes. That result,
along with the correctness and completion checks, supports the SPQR choice
despite its SuiteSparse/CHOLMOD deployment cost.

The EIGEN comparison is historical promotion evidence only. EIGEN is not a
current selectable build or regression backend.

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
rank/preconditioner infrastructure; shared sparse factorization remains active
because the current estimator still uses it.

## Current verification

Use the active benchmark wrapper for current FixedNeighbor measurements:

```bash
python3 tests/integration/joint_benchmark.py \
  --profile solve --case chain-8 --build-dir build/qualification \
  --output build/joint-solve.json
```

The benchmark wrapper has no sparse-backend, rank, operator, preconditioner or
Schwarz modes. Focused C++ tests cover the SPQR factorization, compact
assessment, prepared blocks, endpoint certification and FixedNeighbor
regressions. The EIGEN comparison is not regenerated as a current test.
