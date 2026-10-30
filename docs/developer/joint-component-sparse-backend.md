# Joint local sparse solver backend

The production FixedNeighbor local A/C solve uses Eigen `SparseQR` with
`COLAMDOrdering`. EIGEN is the sole backend; there is no sparse-backend
selector, runtime fallback, or automatic routing. The build does not search for
or link SuiteSparseQR, SPQR, or CHOLMOD.

```sh
cmake -S . -B build/joint-eigen -DCMAKE_BUILD_TYPE=Release \
  -DRHBM_GEM_ENABLE_UMAP=OFF
```

Each Profile LM evaluation solves the constrained local A/C problem at fixed
log widths, preserving `A >= 0`, signed `C`, and the original objective. Each
block visit accepts at most one trusted update. Numeric factors and active
sets are not carried across sweeps. Global stationarity and endpoint replay
remain separate streaming checks and do not invoke a second global A/C solver.

## Local numerical utilities

`LinearSolve`, `StructuralPartition`, `CompactSvd`, and `TiledDerivative`
provide the local sparse solve, block structure, rank checks, and profile
derivative calculations. Local conditioning uses bounded FixedNeighbor block
workspace after the accepted endpoint. It reports conditional block evidence;
it does not certify global rank. The complete derivative retains its
nonzero-residual coefficient correction.

## Memory and historical backend results

Global search and endpoint verification traverse support entries and do not
form a global dense design. The former global rank factors and projection
coefficients have been removed. Current stage-specific memory protections cover
construction, FixedNeighbor preparation, streaming evaluation, endpoint
verification and local conditioning.

The historical EIGEN/SPQR comparison found bounded parity on completed
128/256-atom probes, faster SPQR timings in those cases, and lower EIGEN peak
RSS. Both 512-atom probes timed out. The later production decision selected
EIGEN for memory footprint, dependency simplicity, and maintenance. These are
historical measurements, not current backend options or full-scale guarantees.
Historical operator, preconditioner, and global rank results are summarized in
[Joint component evidence](joint-component-evidence.md).

The production local solver and endpoint qualification checks are covered by
the focused `rhbm_tests_joint_component_numerical` group. Historical backend
comparison results and their limits are summarized in the
[Joint evidence record](joint-component-evidence.md).
