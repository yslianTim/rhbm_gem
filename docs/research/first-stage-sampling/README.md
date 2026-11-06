# FirstStage sampling and interpolation study

## P0: baseline (2026-11-06)

Initial branch: `develop`; local HEAD and fetched origin/develop:
`ec686c675aee6c953dd9929fc60ae0c241746428`. Working tree was clean.
No branch/worktree was created. The supplied Notion page was read only.

Data flow: deterministic Fibonacci shells (20 radii, 0.05..1.95 A,
10 points/shell) -> neighbor ownership/cone selection -> tricubic map
responses -> selected raw entries -> local alpha training -> identical
log-Gaussian OLS/MDPDE -> FirstStage width/fallback -> qualified target
width calibration -> bounded JCE initialization. Both alpha training and
Gaussian fitting select radii <=1.0 A; positive responses only enter the
log design. No fitting-radius change is planned.

Raw vectors are owned by `LocalPotentialEntry` in model analysis, moved by
`ModelAnalysisEditor::SetAtomLocalRawSamplingEntries`. FirstStage views copy
the selected entries; prepared designs retain source indices only for the
lifetime of the fitting call. JSON stage results serialize the actual raw
coordinates, distances, selection flags and responses. Replacing raw data
clears sample-derived Two-Stage peeling/group state, but does not invalidate
unrelated original-grid Post-Joint results. New sources must use this editor
path, never reuse raw indices as original-grid identities.

QScore uses spiral geometry, radii accumulated from 0.1 through 2.0 A,
8 minimum accepted points, up to 50 attempts with 2 extra candidates each,
all accepted points from the first successful attempt, a 2r neighbor query,
0.9r exclusion, and ignores hydrogen. Failed shells return no observations.
Correlation repeats the center 8 times. QScore uses trilinear interpolation:
half-voxel inclusive margins and coordinate clamping; out-of-range or
nonfinite coordinates throw `out_of_range`. Tricubic uses MapObject's floor
index and clamped stencil nodes; finite out-of-range positions log and use
the existing {-1,-1,-1} sentinel. Its nonfinite-to-integer conversion was
undefined and needs an explicit precondition check in P1.

Baseline data: `baseline.json` (24 isolated synthetic targets, spacing 0.4 A,
known A=1.5+0.02i, tau=0.35+0.005i, truncated normalized Gaussian and
deterministic additive noise 1e-4); `baseline-6z6u.json` (first 120
non-hydrogen targets, complete model neighbors, local 6z6u.cif and
emd_11103_unsharpened.map). Real-map widths have no truth reference.
Synthetic: 24/24 qualified donors, bound 1.8782695995625727 A,
final objective 1.631044078500785e-7. Real: 120/120 donors,
bound 3.1823057120925764 A, final objective 0.012784319782223002.

Validation: Release build passed; all 21 CTest groups passed (37.53 s).
The QScore suite and new fixed tricubic reference passed (22 tests).
Use `OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=4 OPENBLAS_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1 ctest --test-dir build/release --parallel 1
--output-on-failure`. Earlier runs with OMP_THREAD_LIMIT=2 failed the
existing four-thread simulation assertion; OMP_NUM_THREADS=1 failed the
existing Eigen thread restoration assertion. No code was changed for them.

Reproduce baseline with the explicit, noninstalled research target:

```
cmake --build build/release --target first_stage_sampling_study --parallel 4
OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=4 build/release/bin/first_stage_sampling_study output.json 0.4
```

Timing and process high-water RSS are measurements, not deterministic regression
values. Builds use <=4 workers; studies run sequentially with one numerical
worker. No new third-party dependency is used.

## P1: shared interpolation

`MapInterpolation.hpp` now supplies `MapInterpolationMethod`, the independent
`InterpolateTrilinear` and `InterpolateTricubic` kernels, and enum-dispatched
`InterpolateMapValue`. Batch sample construction defaults to Tricubic.
The original trilinear arithmetic was moved without changing operation order.
QScore explicitly selects Trilinear and retains its exception text in a wrapper.
Each finite boundary policy remains unchanged; nonfinite cubic positions now
throw before undefined floating-to-integer conversion.

Release build, sampler and QScore CTest groups passed. All original fields of
synthetic baseline atoms, calibration, complete JCE evidence/state and QScore
are bit-for-bit equal. Tests include constant/affine/quadratic/Gaussian fields,
anisotropic spacing, grid and off-grid points, half-voxel margins, finite
out-of-range sentinel behavior, nonfinite rejection and singleton grids.
The Gaussian cubic reference is an interpolation value, not an assertion that
a finite-grid interpolant exactly reproduces an analytic Gaussian.

## P2: shared QScore-inspired sampling

`QScoreSampling.cpp/.hpp` separates pure spiral geometry from the model-aware
shell policy. Public QScore shell sampling delegates to that single policy.
FirstStage candidate construction uses one center and the same accumulated
shell radii and accepted point ordering. Isolated atoms have 161 candidates,
81 in the fitting domain (including the center). Shell diagnostics retain
minimum requested count, attempts, final candidate and accepted counts, and
success; failed shells contribute no observations even if their final attempt
accepted fewer than the minimum. Existing QScore center weighting is unchanged.

Release build and sampler/QScore CTest passed. Added tests verify one center,
all 20 radii and ordering, 1.0 A fitting boundary, pure geometry reference,
the first successful retry returning 3 rather than 2 points, and all 50 failed
attempts ending at 106 candidates without fabricated observations. Existing
hydrogen, crowded-neighbor, foreign-model and correlation tests passed. F0
atoms, calibration, JCE evidence/state and QScore remain exactly equal to P1.

## P3: internal per-run policies

`JointSeedSamplingPolicy` independently selects Fibonacci/QScore-inspired shells
and Tricubic/Trilinear. `PrepareJointWidthInitialization` and the internal Joint
workflow accept a const per-run policy, default F0. `FitOptions`, public CLI,
Two-Stage sampling and QScore correlation options are unchanged. The common
sampling workflow computes responses in workers, rethrows interpolation errors
on the caller before modifying model analysis, then moves samples through the
same editor. No new sample identity or original-grid index interpretation is
introduced. No estimator or calibration/solver qualification rule changes.

Policies are named in normal Joint sampling logs. Research callers request a
`JointSeedSamplingReport` containing the policy and per-atom shell/sample
diagnostics and save it alongside results. This report is run-owned; there is
no global/thread-local experiment state or public persistence-schema extension.
Persisted raw vectors retain their actual geometry and selection; Post-Joint
uses its separate original-grid identities.

All 21 Release CTest groups passed (39.07 s), including new independent-factor,
legacy-default equality, QScore isolation, serial/four-worker deterministic,
repeated-policy-run and parallel-boundary-failure tests, existing Two-Stage,
joint convergence, persistence and invalidation contracts. All original F0
baseline numerical fields remain exactly equal.
