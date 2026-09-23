# Joint compact SVD

Derivative preparation checks the column-normalized free design with a
singular-values-only decomposition. Projection coefficients and the full
nonzero-residual correction still use the existing QR and triangular solves.
The independent reference computes its own weighted reduction and retains both
singular-vector factors for its compact free-face least-squares solve.

The private compact helper uses JacobiSVD below 16 singular values and Eigen
BDCSVD otherwise. Both sparse backends use this rule; there is no user-facing
SVD option or additional dependency. The direct dense oracle, block diagnostic,
endpoint spectra, local correction SVD and LM are unchanged.

Rank and solve use the original caller's relative threshold, including the
original observation dimensions and any absolute-threshold override. Compact
dimensions do not replace those dimensions. Eigen's rank/solve boundary and
zero-spectrum semantics are retained.

After an unsuccessful/nonfinite BDCSVD, or when a singular value is within
`64 * epsilon * max(compact.rows, compact.cols) * sigma_max` of the caller's
absolute threshold, the same matrix is recomputed once with JacobiSVD. This
selects arithmetic, not a different rank threshold, and is not an error-bound
certificate. A failed retry remains invalid. The existing near-cancellation
tiled derivative reduction is separate and retains its original trigger.

The implementation uses Eigen 5 template options and its default BDCSVD
switching size. Unsafe floating-point optimization is not used; see the
[Eigen BDCSVD documentation](https://libeigen.gitlab.io/eigen/docs-5.0/classEigen_1_1BDCSVD.html).

## Counters and replay

The existing internal work counters include reference QR, per-face compact
reduction, derivative compact extraction/reduction, reference and free-design
SVD, reference solve, Jacobi retry and cancellation reduction. Counters cover
both EIGEN and SPQR, including early returns and exceptions. An SVD count is a
helper invocation; BDCSVD attempts and Jacobi retries have separate counts.

`derivative_inclusive_seconds` covers the whole preparation. Reference phase
wall time includes reference QR, compact reduction, SVD and solve. Retry time
overlaps its SVD/solve counters; cancellation time overlaps the additional
tiled reduction and SVD. These inclusive and nested times must not be added.
The SVD counters time decomposition, and the reference solve counter times
application of the factors. OS peak RSS is not an allocation-level attribution.

In testing builds, `joint_sparse_benchmark` accepts `--svd-mode legacy`,
`--svd-mode values` or `--svd-mode auto`. Legacy reproduces the previous vector
requests; values removes only the unused derivative U/V; auto uses production
dispatch and its rank-boundary retry. These controls are absent from
testing-disabled builds and never read environment variables or CLI options in
the production executable.

`--audit` records per-call spectra and derivative coefficients, residual
correction and reduced derivative matrices. `--capture DIRECTORY` additionally
saves the actual compact matrices in native-endian binary64 column-major order,
with JSON dimensions, byte order, transformed RHS and threshold inputs. Capture
and audit runs are separate from formal timing. A captured matrix can be replayed:

```sh
build/joint-compact-spqr/bin/joint_sparse_benchmark replay \
  CAPTURE/reference-0.json auto replay.json
```

## Reproduction

Retain EIGEN and SPQR Release builds of `ce58c897` as the baselines. Configure
candidate builds identically with SYSTEM dependencies, OpenMP ON, ROOT/UMAP/
Python bindings OFF, and the joint extended/offline options ON. Run the joint
core, runtime, extended, offline and CLI tests for both backends, the resource
runner tests, repository guards and a testing-disabled build before timing.
Build `tests_all` and `rhbm_gem_cli` in both candidate directories, and explicitly
build `joint_validation` in the SPQR candidate directory for input generation.

```sh
python3 tests/integration/joint_compact_validation.py \
  --baseline-spqr build/joint-sparse-spqr \
  --baseline-eigen build/joint-sparse-eigen \
  --spqr build/joint-compact-spqr --eigen build/joint-compact-eigen \
  --work-dir build/joint-compact-measurements
```

The runner generates shared 128/512 inputs and serialized production widths,
and unpacks the frozen heterogeneous-168 fixture. It compares three fresh
processes per completed fixed state and checks the instrumented legacy mode
against the frozen baseline executable. Full-command runs compare the frozen
SPQR command with production auto dispatch, including SQLite and JSON/CSV
export. There is no production switch for running the intermediate values-only
command; its contribution is isolated in fixed states and exact matrix replay.

Use one numerical worker, with no overlapping builds or other benchmarks.
Each process or analysis/export pipeline is bounded to 600 seconds and a
sampled 4 GiB process-tree limit; the campaign is bounded to 80 minutes. The
watchdog is not an OS hard memory limit. A resource stop ends repetitions for
that case; unfinished and budget-skipped work is explicit. Capture/replay is
diagnostic work and is excluded from timing repetitions. All raw samples,
input and binary fingerprints, build settings and linked libraries are saved.
`--report-only` recomputes gates from the saved receipt and audit outputs without
running estimators. Audit paths are relative to the work directory and their
hashes are checked, so a retained evidence bundle can be relocated.

The compact gate requires numerical controls to pass and Single 512 reference
SVD, free-design SVD and reference-plus-derivative median times each to decrease
by at least 30%. Combined 128/168 fixed-state time may regress at most 10%.
Completed candidate runs must stay within 4 GiB. The 128 full command must
retain converged, numerically equivalent outputs. Single 512 command completion
and runtime convergence are a separate gate; no speedup is assigned to a
timeout or an unavailable endpoint.

See the [acceptance record](joint-component-compact-svd-acceptance.md) for measured
results and the retained evidence. Initial-rank attribution for 6Z6U and
operator/global-assessment changes are outside this work package.
