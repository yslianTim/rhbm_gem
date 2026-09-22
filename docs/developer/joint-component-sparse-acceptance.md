# Joint sparse backend acceptance

Baseline: `ceb6155992b9c891ab3fa878e875c6b10aa0590f`.
Candidate production source fingerprint: `9718531067e72cd3ecf180d12ae8033b8e39360b7fe1d49d24ae847de4240fa0`.

## Fixed-state measurements

Times are seconds, medians of three fresh processes for completed cases. Each cell is baseline / SPQR. Peak RSS is the maximum OS-reported process peak among the samples, in MiB.

| Case | Primary | Reference | Derivative preparation | Peak RSS |
| --- | ---: | ---: | ---: | ---: |
| single-128 | 0.121 / 0.033 | 0.820 / 0.465 | 1.067 / 0.676 | 358 / 302 |
| heterogeneous-168 | 0.362 / 0.168 | 1.617 / 0.435 | 2.711 / 1.643 | 318 / 416 |
| single-512 | 3.615 / 1.457 | 68.923 / 54.661 | 83.888 / 78.323 | 704 / 1267 |

6Z6U: one attempt per backend; both stopped during reference at the shared 600-second process limit. Primary reported rank 4314/4334 and `rank-deficient` in both builds. Neither provides an accepted initial state. No endpoint speedup is assigned.

## Complete command

Single 128 analysis → SQLite → JSON/CSV export: 19.20 → 11.78 seconds (median, 1.63×). All six runs completed and passed runtime convergence. All nine cross-backend endpoint pairs retain parameters within 1e-10 and normalized objective within 1e-12.

Single 512: both complete-command attempts reached 600 seconds during analysis;
neither produced a persisted/exported endpoint. Sampled process-tree peaks were
675 MiB (baseline) and 1811 MiB (SPQR). These are censored runs, not equal-time
completion measurements, and no end-to-end speedup is calculated.

## Numerical and packaging validation

Both EIGEN and SPQR Release builds passed 27 CTests covering the normal C++/CLI,
extended and offline groups; the resource-runner CTest was run separately with
permission to inspect processes (nine Python controls). The configuration was
SYSTEM dependencies, OpenMP ON, ROOT/UMAP/Python bindings OFF. This does not
claim coverage of the disabled optional features. The sparse core group contains
62 tests, including pattern/value/policy/scope invalidation, stale handles,
canonical faces, weighted and zero-weight rows, near-degenerate reference SVD,
nonzero residual correction and independent dense controls.

Testing-disabled Eigen installation and SPQR shared/static RHBM-GEM installations
built and ran external consumers. The static RHBM-GEM consumer resolves installed
SuiteSparse dependencies; this is not a fully static operating-system executable.
Disabling SPQR package discovery fails an explicitly requested SPQR configuration;
the default Eigen install build succeeds with that package disabled.

A fresh SPQR testing-disabled source-package build also passed install/consumer
execution from `/private/tmp`, with neither a `.git` directory nor a discoverable
parent repository. Its production source fingerprint matches the tested builds.
Repository guards and the final supplemental core/resource-runner tests passed.

All available fixed-state sample pairs agree in status, rank and active face;
coefficient errors also pass the stricter absolute 1e-10 check. Both backends
retain their independent Guarded trust checks. 6Z6U has matching unavailable
primary states and incomplete reference evidence, so it is not counted as a
successful numerical endpoint. No frozen expectations, rank/convergence gates
or search budgets were changed.

One numerical accommodation is required: nearly cancelled projected width
columns use the existing tiled free-design arithmetic, because QR-ordering
roundoff otherwise exceeded the existing normalized-spectrum comparison. This
is a counted computational branch, not a relaxed acceptance threshold. It is
covered by the frozen active-A cases and a direct cancellation control. None of
the completed fixed-state SPQR measurements triggered it.

## Late diagnostic and delivery decision

The separate 600-second diagnostic stopped at the same limit. Its successful
one-second stack samples show:

| Sample | Active numerical work |
| --- | --- |
| 120 seconds | `SearchProfile → PrepareDerivative → JacobiSVD` |
| 300 seconds | `Profile::Trial → CheckTrust → SolveLinear → JacobiSVD` (reference) |
| 540 seconds | `SearchProfile → PrepareDerivative → JacobiSVD` |

These snapshots identify the remaining free-design and reference SVD work during
search; they do not estimate whole-run percentages or establish assessment as
the bottleneck. Repository lint checks overlapped this diagnostic only; formal
timing had already ended. The [diagnostic receipt](figures/joint-sparse-acceptance/late-profile/receipt.json)
records binary and stack hashes, with the [120](figures/joint-sparse-acceptance/late-profile/stack-120.txt),
[300](figures/joint-sparse-acceptance/late-profile/stack-300.txt) and
[540](figures/joint-sparse-acceptance/late-profile/stack-540.txt) second stacks retained.

The sparse-factorization milestone passes its numerical controls and reduces
measured factorization costs. Single 512 does **not** pass the complete-command
resource gate. Further work should address the measured compact SVD costs while
preserving the rank/trust contract; no operator LM, preconditioner, warm start or
block global assessment was introduced here. SPQR remains explicitly enabled.

## Measurement interpretation and evidence

The host was Apple M1, 16 GiB RAM, macOS arm64, AppleClang 21. SPQR was 4.3.6,
CHOLMOD 5.3.5, SuiteSparse_config 7.14.1 and Apple Accelerate 4.0.0. Numerical
workers were restricted to one. Builds and other benchmarks did not overlap
formal timing. Fixed-state diagnostics and complete-command measurements are
separate processes; the later stack-sampling diagnostic is excluded from both.

Completed fixed states each reused one symbolic analysis for two primary numeric
factorizations and reused the current factor for derivative preparation. The
maximum reported R-fill bounds were 8420, 16200 and 276208 for 128, 168 and 512;
6Z6U's initial factor bound was 915788. These are SPQR's bounds on nonzeros in R,
not exact allocation sizes. Basis/CSC preparation, symbolic, numeric, reference
QR/SVD and phase wall times are retained per sample; nested counters are not
additive. For Single 512, reference SVD alone took about 53 seconds in the first
SPQR sample, while reference sparse QR took about 1.2 seconds.

The fixed-state runner uses shared serialized production initialization values.
Both backends' actual initial arrays match exactly. For 6Z6U, JSON parsing differs
from the in-memory initialization serialization by at most 5.6e-17; actual arrays
are retained in every state record. The complete commands run the production
initializer directly. The original measurement receipt also retains its raw
objective comparison; the final audit separately checks normalized objectives.

- [All measurement samples and configuration fingerprints](figures/joint-sparse-acceptance/measurements.json.gz)
- [All cross-sample numerical comparisons](figures/joint-sparse-acceptance/numerical-comparison.json)
- [Build, test and installation receipt](figures/joint-sparse-acceptance/verification.json)
- [Complete-command JSON/CSV exports and process logs](figures/joint-sparse-acceptance/command-exports.tar.gz)
- [Eigen CTests](figures/joint-sparse-acceptance/eigen-tests.txt) and [SPQR CTests](figures/joint-sparse-acceptance/spqr-final-tests.txt)
- [Benchmark source used for timing](figures/joint-sparse-acceptance/benchmark-source.cpp.txt) and [measurement runner source](figures/joint-sparse-acceptance/runner-source.py.txt)

The current [backend guide](joint-component-sparse-backend.md) describes the build
switch and reproduction commands. Full local build/install logs and generated
inputs remain under `build/joint-sparse-evidence/` and
`build/joint-sparse-measurements/`. Original 6Z6U inputs are identified by hashes
and are not redistributed here.
