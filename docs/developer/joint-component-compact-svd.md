# Joint compact SVD

Current SVD and fixed-state measurements use the
[Joint benchmark guide](joint-benchmark.md). Historical results and
source-pinned provenance are indexed in the
[canonical historical evidence](joint-component-evidence.md).

Profile LM retains the full derivative, including the nonzero-residual
correction, through the existing QR and triangular solves. Numerical
qualification builds only the projected-width derivative needed for its rank
assessment; it does not form a correction Jacobian or run a correction SVD or
step solve.
The independent reference computes its own weighted reduction and retains the
singular-vector factors needed by its compact free-face least-squares solve.
Derived Target Identifiability reuses Design Rank and Projected-width Rank and
does not run a full `(A, C, log B)` Jacobian SVD. The independent uncertainty
path may run a compact full-rank factorization for covariance only after its
separate `target-covariance` memory preflight; it does not certify target
identifiability.

The private compact helper uses JacobiSVD when the compact dimension is below
16 and Eigen BDCSVD otherwise. The current EIGEN-backed path uses this dispatch;
there is no user-facing SVD option or additional dependency. The direct dense `DesignSpectrum()` oracle remains available to offline and test
controls; production endpoint qualification no longer calls it. The one
free-design SVD in projected-width preparation uses the original observation
row count and full A/C column count for its relative threshold, even when the
free face has fewer columns. It applies the original strict `singular_value >
threshold` rule. Singular values within one representable step of the threshold
are marked sensitive and cannot qualify. Profile LM's derivative is unchanged;
the production correction-qualification SVD and solve have been removed.

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
the current EIGEN path, including early returns and exceptions. An SVD count is a
helper invocation; BDCSVD attempts and Jacobi retries have separate counts.

derivative_inclusive_seconds covers the whole preparation. Reference phase
wall time includes reference QR, compact reduction, SVD and solve. Retry time
overlaps its SVD/solve counters; cancellation time overlaps the additional
tiled reduction and SVD. These inclusive and nested times must not be added.
The projected-width preparation counts one free-design TiledQR and one SVD rank
assessment; that same `R` factor supplies free-design rank evidence and the
triangular solve for projection coefficients. The separate projected-width
QR/SVD remains necessary. SVD counters time decomposition; the reference-solve
counter times application of factors. OS peak RSS is not an allocation-level
attribution.

## Projected-width rank and weak-direction diagnostics

Projected-width rank, threshold sensitivity and `identified` use singular
values only. The right-vector matrix is used only to produce up to three
weak-direction diagnostics. That helper is now separate from the rank and
qualification assignments, but current runtime JSON and offline audit tools
consume those directions, so production still requests Eigen's thin right
vectors. Removing that request from production would change the present
diagnostic output contract.

`CompactSvdVectors::None` selects a values-only Eigen mode; `Right` selects a
thin right-vector mode. The projected-width research control compares their
singular values, rank, original threshold, and threshold-sensitivity decision,
including threshold-adjacent inputs. Its reported right-vector byte count is
the explicit `P × min(N, P)` matrix size. It is not an observed function-level
RSS reduction: process peak RSS is a whole-process high-water mark, and the
production diagnostic still retains the right-vector mode. No production
memory benefit is claimed from this separation.

## Current measurement

Use the current FixedNeighbor solve profile in the EIGEN-only build:

    python3 tests/integration/joint_benchmark.py \
      --profile solve --case chain-8 \
      --build-dir build/qualification --output build/joint-solve.json

The historical EIGEN/SPQR parity and resource conclusion is summarized in
[joint-operator-search.md](joint-operator-search.md). Its campaign outputs are
retired; source and raw results remain available from Git history.
