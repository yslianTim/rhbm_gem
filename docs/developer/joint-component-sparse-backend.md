# Joint sparse factorization backends

`RHBM_GEM_JOINT_SPARSE_BACKEND` selects the linear algebra implementation at
configure time. `EIGEN` is the default and needs no SuiteSparse installation.
`SPQR` requires an installed SPQR 4.x CMake package, CHOLMOD and their transitive
dependencies. Both SYSTEM and FETCH dependency modes use that installed package;
missing dependencies cause configuration to fail, without downloading SPQR or
changing the selected backend.

```sh
cmake -S . -B build/joint-spqr -DCMAKE_BUILD_TYPE=Release \
  -DRHBM_GEM_JOINT_SPARSE_BACKEND=SPQR
```

SPQR is GPL-2.0-or-later (alternate licenses are available from its author).
See the third-party notices. Installed shared builds need the linked runtime
libraries. Installed static builds resolve the dependency targets through the
package configuration; an upstream OpenMP-enabled SuiteSparse package also
requires a C compiler. AppleClang dependency discovery uses the same Homebrew
libomp hint as the project's OpenMP support.

## Numerical behavior

The primary SPQR solver factors the original weighted, column-normalized free
matrix with COLAMD ordering and the existing absolute pivot threshold. Each
component search owns one symbolic workspace. Reuse requires identical scope,
linear policy, actual pivot threshold, free-column identities and actual
compressed storage pattern.
Every numeric evaluation refactorizes its values. Face or storage-pattern changes
replace the symbolic analysis; no unbounded face cache or warm start is used.
The implementation follows the expert symbolic/numeric, orthogonal transform
and triangular solve interfaces in the [SPQR user guide](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/SPQR/Doc/spqr_user_guide.tex).

The reference independently rebuilds the original design and uses SPQR without
tolerance-based direction truncation to obtain a compact orthogonal reduction.
Column permutations are undone before the existing mixed-constrained active-set
and SVD solves. SVD rank thresholds retain the original observation dimensions
and numerical policy. Reference factors never come from the primary workspace.
Scalar support replay, KKT and gradient checks still govern trial trust.

A scoped free-design handle allows the derivative to reuse the current primary
factor. It checks matrix values, column identities and generation. Numerical
refactorization expires previous handles. Derivative columns are all C columns
plus strictly positive A columns; a different canonical face gets its own
factor. Evaluating supplied coefficients never replaces those coefficients.

Projection and residual correction use orthogonal transforms and triangular
solves, with at most 16 right-hand sides per batch. No dense observation-sized Q,
pseudo-inverse, normal-equation matrix or full Jacobian is constructed. The
free-design SVD checks the compact factor without computing unused U/V. Compact
free-design and reference solves use the [internal compact SVD policy](joint-component-compact-svd.md),
including a counted Jacobi retry near the original rank threshold. When projection nearly cancels
a raw width column (relative norm at most `64 * epsilon / 1e-10`), derivative
preparation retains the existing tiled free-design arithmetic to protect the
normalized-spectrum precision contract. This counted cancellation reduction
adds work and is included in measurements; it does not change a rank or
convergence threshold. The Jacobian reduction, compact LM
and endpoint spectra retain their existing algorithms and thresholds. Compact
factors and SVD retain quadratic storage and potentially cubic work. Sparse
fill-in and RHS batches add workspace; this is not a fixed memory bound.

The EIGEN route retains the existing row reduction, SparseQR, tiled reference
and tiled derivative. Public estimator APIs, selection/initialization behavior,
result schema and search budgets are shared by both routes. The configuration
fingerprint records the backend, SPQR, CHOLMOD, SuiteSparse_config and the BLAS
link information published by the dependency package. Benchmark receipts also
retain executable/library hashes and linked-library version information.

## Verification and bounded measurement

The [acceptance record](joint-component-sparse-acceptance.md) contains measured
costs, numerical comparisons, all timing samples and current scaling limits.

Use the default, extended and offline joint tests with both builds. The sparse
unit controls cover weighted and zero-weight rows, near-degenerate SVD rank,
cache invalidation, expired handles, canonical faces and derivative reuse.
Frozen inputs, numerical tolerances and nonconvergence limitations are unchanged.

`joint_sparse_benchmark` is test-only. `prepare MODEL MAP WIDTHS` runs the production
initializer without fitting; `initial MODEL MAP WIDTHS OUTPUT` evaluates the
largest component at the saved widths, including an independent reference,
derivative preparation and a repeated primary evaluation. `fixture DATASET CASE
OUTPUT` uses a frozen fixture. Stage records are flushed before expensive work,
so a resource stop does not masquerade as a numerical result.
The matrix-preparation counter covers basis construction and CSC index conversion;
it is not a disjoint total for all weighting and column assembly. Factor fill is
the library's reported upper bound on `nnz(R)`, not a measured byte allocation.
Phase wall times include their nested counters and must not be added to them.

```sh
python3 tests/integration/joint_sparse_validation.py \
  --baseline build/joint-sparse-baseline \
  --candidate build/joint-sparse-spqr \
  --work-dir build/joint-sparse-measurements
```

The baseline driver is compiled from the same benchmark source with
`SPARSE_BASELINE_DRIVER` against the frozen baseline library/test support. That
flag excludes candidate-only instrumentation and initialization entry points.
The runner uses identical widths and input hashes, one numerical worker, three
fresh processes for completed cases, and a single attempt after a resource stop.
Fixed-state workloads are 128, heterogeneous-168, 512 and the supplied 6Z6U
largest component. Full command workloads are 128 and 512, including SQLite and
JSON/CSV export. The analysis/export pair shares 600 seconds and a sampled
4 GiB process-tree threshold. The full measurement campaign is bounded to
80 minutes; skipped work remains explicit. The watchdog is not a hard OS memory
limit. Do not run builds or other benchmarks concurrently with measurements.

On macOS, `--profile-only` samples a completed campaign's Single 512 candidate
timeout at 120, 300 and 540 seconds in a separate diagnostic run. It retains the
same 600-second / 4-GiB cap and verifies executable and input fingerprints.
Those samples are excluded from timing repetitions.
