# Joint compact SVD

Current SVD and fixed-state measurements use the
[Joint benchmark guide](joint-benchmark.md). Historical results and
source-pinned provenance are indexed in the
[canonical historical evidence](joint-component-evidence.md).

Derivative preparation checks the column-normalized free design with a
singular-values-only decomposition. Projection coefficients and the full
nonzero-residual correction still use the existing QR and triangular solves.
The independent reference computes its own weighted reduction and retains the
singular-vector factors needed by its compact free-face least-squares solve.

The private compact helper uses JacobiSVD when the compact dimension is below
16 and Eigen BDCSVD otherwise. The current SPQR-backed path uses this dispatch;
there is no user-facing SVD option or additional dependency. The direct dense
oracle, block diagnostic, endpoint spectra, local-correction SVD and LM are
unchanged.

Rank and solve use the original caller's relative threshold, including the
original observation dimensions and any absolute-threshold override. Compact
dimensions do not replace those dimensions. Eigen rank-boundary and
zero-spectrum semantics are retained.

After an unsuccessful or nonfinite BDCSVD, or when a singular value lies within
64 * epsilon * max(compact.rows, compact.cols) * sigma_max of the caller's
absolute threshold, the same matrix is recomputed once with JacobiSVD. This
selects arithmetic, not a different rank threshold, and is not an error-bound
certificate. A failed retry remains invalid. This behavior does not relax rank
or convergence gates.

The existing tiled derivative has a separate near-cancellation fallback. When
the projected width column is within 64 * epsilon / 1e-10 relative roundoff of
cancellation, the derivative accumulation follows the independent reference's
row order and fused multiply-adds before normalization, then recomputes the
projection with tiled QR. This protects the normalized-spectrum precision
contract without changing a threshold. Permanent numerical tests exercise the
fallback.

The implementation uses Eigen 5 template options and its default BDCSVD
switching size. Unsafe floating-point optimization is not used; see the
[Eigen BDCSVD documentation](https://libeigen.gitlab.io/eigen/docs-5.0/classEigen_1_1BDCSVD.html).

## Counters and timing

Current internal work counters include reference QR, per-face compact
reduction, derivative compact extraction/reduction, reference and free-design
SVD, reference solve, Jacobi retry and cancellation reduction. Counters cover
the current SPQR path, including early returns and exceptions. An SVD count is a
helper invocation; BDCSVD attempts and Jacobi retries have separate counts.

derivative_inclusive_seconds covers the whole preparation. Reference phase
wall time includes reference QR, compact reduction, SVD and solve. Retry time
overlaps its SVD/solve counters; cancellation time overlaps the additional
tiled reduction and SVD. These inclusive and nested times must not be added.
SVD counters time decomposition; the reference-solve counter times application
of the factors. OS peak RSS is not an allocation-level attribution.

## Current measurement

Use the current FixedNeighbor solve profile in the SPQR build:

    python3 tests/integration/joint_benchmark.py \
      --profile solve --case chain-8 \
      --build-dir build/qualification --output build/joint-solve.json

The one-time EIGEN comparison is retained in the compact
[backend qualification evidence](figures/joint-fixed-neighbor-backend-qualification-r1/).
Historical campaign results, limitations, and source provenance are indexed in
the [canonical historical evidence](joint-component-evidence.md).
