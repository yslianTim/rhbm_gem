# Joint sparse factorization backends

This document owns current EIGEN/SPQR backend behavior, including the bounded
SPQR rank-certificate contract. Current measurements use the
[Joint benchmark guide](joint-benchmark.md); retired campaign outcomes and
source-pinned provenance are indexed in the
[canonical historical evidence](joint-component-evidence.md).

RHBM_GEM_JOINT_SPARSE_BACKEND selects the linear algebra implementation at
configure time. EIGEN is the default and needs no SuiteSparse installation.
SPQR requires an installed SPQR 4.x CMake package, CHOLMOD and their transitive
dependencies. Both SYSTEM and FETCH dependency modes use that installed package;
missing dependencies cause configuration to fail without downloading SPQR or
changing the selected backend.

    cmake -S . -B build/joint-spqr -DCMAKE_BUILD_TYPE=Release \
      -DRHBM_GEM_JOINT_SPARSE_BACKEND=SPQR

SPQR is GPL-2.0-or-later (alternate licenses are available from its author).
See the third-party notices. Installed shared builds need the linked runtime
libraries. Installed static builds resolve dependency targets through the
package configuration; an upstream OpenMP-enabled SuiteSparse package also
requires a C compiler. AppleClang dependency discovery uses the same Homebrew
libomp hint as the project's OpenMP support.

## Numerical behavior

The primary SPQR solver factors the original weighted, column-normalized free
matrix with COLAMD ordering and the existing absolute pivot threshold. Each
component search owns one symbolic workspace. Reuse requires identical scope,
linear policy, actual pivot threshold, free-column identities and actual
compressed storage pattern. Every numeric evaluation refactorizes its values.
Face or storage-pattern changes replace the symbolic analysis; no unbounded
face cache or warm start is used. The implementation follows the expert
symbolic/numeric, orthogonal transform and triangular solve interfaces in the
[SPQR user guide](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/SPQR/Doc/spqr_user_guide.tex).

The reference independently rebuilds the original design and uses SPQR without
tolerance-based direction truncation to obtain a compact orthogonal reduction.
Column permutations are undone before the existing mixed-constrained active-set
and SVD solves. SVD rank thresholds retain the original observation dimensions
and numerical policy. Reference factors never come from the primary workspace.
Scalar support replay, KKT and gradient checks still govern trial trust.

A scoped free-design handle allows derivative preparation to reuse the current
primary factor. It checks matrix values, column identities and generation.
Numerical refactorization expires previous handles. Derivative columns are all C
columns plus strictly positive A columns; a different canonical face gets its
own factor. Evaluating supplied coefficients never replaces those coefficients.
Projection and residual correction use orthogonal transforms and triangular
solves, with at most 16 right-hand sides per batch. No dense observation-sized Q,
pseudo-inverse, normal-equation matrix or full Jacobian is constructed.

Free-design SVD checks the compact factor without computing unused U/V.
Compact free-design and reference solves use the
[compact SVD policy](joint-component-compact-svd.md), including a counted
Jacobi retry near the original rank threshold. When projection nearly cancels a
raw width column (relative norm at most 64 * epsilon / 1e-10), derivative
preparation retains the existing tiled free-design arithmetic to protect the
normalized-spectrum precision contract. This adds work but changes no rank or
convergence threshold. Compact factors retain quadratic storage and potentially
cubic work; sparse fill-in and right-hand-side batches also use workspace.

The EIGEN route retains the existing row reduction, SparseQR, tiled reference,
and tiled derivative. Public estimator APIs, selection/initialization behavior,
result schema and search budgets are shared by both routes. Configuration
fingerprints record the backend and relevant SPQR, CHOLMOD, SuiteSparse_config
and BLAS link information. Benchmark receipts also retain executable/library
hashes and linked-library version information.

## Bounded SPQR rank certificate

The bounded SPQR rank route is internal and is selected through the explicit
Operator-PCG policy or the rank-only benchmark profile. The production search
default remains `LegacyCompact`. `EvaluateFreeDesignRank` accepts the normalized
sparse design, an optional immutable factor, the existing `RankRequest`, and a
`RankBudget`. Its result is `FullRank`, `Deficient`, or `Unavailable`, with rank
bounds, threshold/spectral bounds, certificate path, charged work, workspace
bound, and elapsed time. Without a factor, structural checks and a
sufficient-only local-support witness can still decide rank; if neither
certifies, the result is unavailable. Dense `EvaluateRank` remains the oracle
and is never a hidden fallback.

For the public SPQR representation

    Z P = Q [R; 0] + E

the certificate proceeds in this order:

1. Sparse column, row and Frobenius norms enclose the largest singular value
   and therefore the original rank threshold. An absolute override includes
   the dense backend's division/multiplication rounding; unresolved subnormal
   cutoffs are unavailable.
2. Exact zero/duplicate columns and deficient dimensions provide structural
   bounds. Nonzero problems at a zero threshold remain unavailable rather than
   promising agreement with roundoff-sensitive SVD equality. The zero matrix
   and insufficient row dimension have direct bounds.
3. For structurally nondeficient designs, a local-support witness may certify
   `FullRank`. Columns are grouped by identical sparse row support; each group
   needs exclusive rows and has size one or two. Selected exclusive rows form
   disjoint square blocks. For a one-column block, an outward-rounded magnitude
   lower bound is used. For a two-column block, the lower bound is
   `|det(B)|_lower / ||B||_F,upper`, also evaluated with outward-rounded
   intervals. If `S` contains the selected rows, then
   `Z'Z = S'S + R'R`, so `sigma_min(Z) >= sigma_min(S)`; for the block diagonal
   `S`, the global lower bound is the minimum block bound. Every block bound
   must be strictly greater than the same `RankRequest` threshold upper bound.
   This is a sufficient-only certificate: incomplete coverage, unsupported
   groups, missing exclusive rows, interval uncertainty, or a lower bound at
   or below threshold falls through to the existing SPQR factor inspection and
   reconstruction (or returns unavailable if no factor was supplied); it does
   not imply deficiency.
4. Up to three smallest-pivot triangular directions are checked against the
   original sparse design. Only a direct-action upper bound strictly below the
   threshold lower bound establishes numerical deficiency. Pivot size alone
   never decides rank.
5. Let M(R) have diagonal abs(R_ii) and off-diagonal -abs(R_ij).
   Nonnegative triangular solves with ones bound the infinity and 1-norms of
   R inverse. Their geometric mean bounds its spectral norm. An already
   insufficient comparison bound returns rank-bound-too-wide without an
   expensive futile verification. See
   [Higham's triangular inverse norm bound](https://nhigham.com/2021/03/30/bounds-for-the-norm-of-the-inverse-of-a-triangular-matrix/).
6. For each stored reflector I - tau h h', enclose its exceptional eigenvalue
   1 - tau norm(h)^2. The product of the smaller of its absolute lower bound
   and one bounds sigma-min(Q) below; permutations preserve it.
7. Apply the stored reflectors to each R column using outward-rounded
   intervals, compare with the corresponding Z column, and accumulate a
   Frobenius upper bound delta on E. Only one observation column is retained.

Full rank requires

    q_lower / inverse_norm_upper - delta > threshold_upper

including the native SVD positive-minimum guard. Equality, overlap, overflow,
and budget exhaustion remain unavailable. Every elementary bound expands
outward with nextafter. IEEE binary arithmetic and round-to-nearest are
required; fast-math is rejected. No guessed safety factor or stochastic
estimate is used as a certificate.

Current rank budget defaults are 120 seconds, 100 million charged work entries,
and 256 MiB of workspace. The checked workspace estimate accounts for row
vectors, the sparse support census, and column/group metadata. The workspace
excludes the input design and immutable factor; the process watchdog includes
their storage, construction scratch and allocator overhead. Exported
sparse-array bytes and known shape probes are reported separately; they are
not allocation traces. Benchmark-only budget overrides do not alter these
defaults.

Rank work telemetry names the active certificate stage, including
`local-witness`, and reports the certificate path plus local-witness coverage,
threshold and lower bound. On the SPQR fallback path it reports normalized
design, R-factor, and Householder storage counts. Before reconstruction it
forecasts the exact remaining charge from the stored reflector nonzeros and
observation count: p * (2*nnz(H) + n). Checked size arithmetic makes the
estimate unavailable on overflow without changing the rank decision. The
charged work counter and `rank-work-budget` reason retain their existing
meaning; a failed budget check also reports its active stage.

The prototype allocates no global dense compact, normal matrix, inverse or
singular-vector matrix. It retains sparse QR and its fill-in. Streaming
verification may cost O(p * nnz(H) + n*p), and comparison bounds can be
conservative. This certificate prototype does not claim linear runtime,
complete rank recovery, or large-workflow scalability.

## Current verification and measurements

Permanent sparse tests cover weighted/zero-weight rows, near-degenerate
reference SVD rank, cache invalidation, expired handles, canonical faces,
derivative reuse, rank bounds, equality and overrides, structural deficiency,
resource budgets, factor identity, and dense-oracle parity.

Current fixed-state examples use the same frozen case in separate backend builds:

    python3 tests/integration/joint_benchmark.py \
      --profile fixed --case chain-8 --build-dir build/joint-sparse-eigen \
      --output build/joint-sparse-eigen-fixed.json
    python3 tests/integration/joint_benchmark.py \
      --profile fixed --case chain-8 --build-dir build/joint-sparse-spqr \
      --output build/joint-sparse-spqr-fixed.json

Historical backend results, incomplete resource runs, and unavailable
archives are summarized in the
[canonical historical evidence](joint-component-evidence.md).
