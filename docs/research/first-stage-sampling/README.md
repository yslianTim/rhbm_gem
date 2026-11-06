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
OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=4 build/release/bin/first_stage_sampling_study output.json coarse:0 F0 production
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

## P4: FirstStage experiment

The explicit research executable and stdlib Python runner are reproducible and
not installed or included in ordinary CTest. The runner saves detailed records
in `build/release/first-stage-study` and compact summaries in `results.json`.
The two P0 baselines pass exact comparison of every originally recorded atom,
calibration, JCE state/evidence and QScore field. No synthetic truth enters
calibration or the estimator.

Run the current study with:

```
cmake --build build/release --target first_stage_sampling_study --parallel 4
python3 tests/research/run_first_stage_sampling_study.py --real-data /Users/yslian/data
# Reaggregate existing detailed outputs without rerunning numerical fits:
python3 tests/research/run_first_stage_sampling_study.py --real-data /Users/yslian/data --summarize-only
```

The executable accepts `output.json scenario[:phase] policy production|bound
[model.cif map.map]`; e.g. `output.json coarse:0 F0 production`. Scenarios are
fine (0.2 A spacing), coarse (0.4 A), noisy (0.4 A, additive sinusoidal noise
amplitude 0.05), crowded (0.4 A, atomic separation 1.4 A), or real. Other
synthetic cases have separation 5.5 A and noise amplitude 1e-4. There are 24
targets with varying known A/tau and three deterministic position/noise phases
per synthetic scenario: 72 observations of atom-level error per cell. These
phases are sensitivity cases, not independent random trials or a confidence
interval. Synthetic map contributions use the same normalized Gaussian and
2.5 A support as JCE. FirstStage remains a single-atom approximation in crowded
maps, which contain signal from neighboring atoms.

MDPDE Width RMSE (A), pooling all 72 atom-level errors per scenario:

| Scenario | F0 | F1 | Q0 | Q1 |
|---|---:|---:|---:|---:|
| Fine | 0.000905 | 0.008654 | 0.000693 | 0.008660 |
| Coarse | 0.011701 | 0.034004 | 0.010327 | 0.034660 |
| Noisy | 0.015587 | 0.034316 | 0.013747 | 0.033833 |
| Crowded | 0.027900 | 0.053872 | 0.072629 | 0.104191 |

MDPDE integrated Amplitude RMSE (same pooling):

| Scenario | F0 | F1 | Q0 | Q1 |
|---|---:|---:|---:|---:|
| Fine | 0.002609 | 0.003461 | 0.001842 | 0.007214 |
| Coarse | 0.032598 | 0.017302 | 0.025432 | 0.033970 |
| Noisy | 0.087842 | 0.043386 | 0.075530 | 0.050453 |
| Crowded | 0.267184 | 0.340042 | 0.931593 | 1.142455 |

The Amplitude ordering differs from Width ordering, particularly under noise;
there is no universally best interpolant for every parameter in this study.

OLS Width RMSE is also retained in `results.json`, along with OLS/MDPDE A/tau
means, bias and RMSE, rank counts, retained counts and distinct radial counts.
All synthetic cells have 24/24 native qualified donors, no refined donors,
unqualified donors or median fallbacks. Designs have rank two. Qualification
was not loosened. The noisy case loses up to three Q0 positive fitting rows;
interpolation therefore changes positive-response selection as well as values.

Representative coarse phase 0, counts summed over 24 targets; time is median
initialization time over the three phases; RSS is maximum process high-water
RSS, including the separate downstream validation:

| Metric | F0 | F1 | Q0 | Q1 |
|---|---:|---:|---:|---:|
| Raw candidates | 4800 | 4800 | 3864 | 3864 |
| Selected candidates | 4800 | 4800 | 3864 | 3864 |
| Positive fitting rows (r <=1 A) | 2400 | 2400 | 1944 | 1944 |
| Distinct fitting radii/atom | 10 | 10 | 11 | 11 |
| Median alpha (phase 0) | 0.30 | 0.25 | 0.35 | 0.20 |
| Qualified donors | 24 | 24 | 24 | 24 |
| Donor coverage | 100% | 100% | 100% | 100% |
| Median fallbacks | 0 | 0 | 0 | 0 |
| Initialization ms | 56.03 | 50.29 | 56.28 | 53.22 |
| Peak RSS MB (decimal) | 46.15 | 46.33 | 46.43 | 46.55 |

For the response **mean Width bias**, the factorial sampling main effect is
((Q0-F0)+(Q1-F1))/2, interpolation effect is ((F1-F0)+(Q1-Q0))/2, and interaction
is (Q1-Q0)-(F1-F0). Effects in A, averaged over phases:

| Scenario | Sampling | Interpolation | Interaction |
|---|---:|---:|---:|
| Fine | -0.000132 | +0.007820 | +0.000203 |
| Coarse | -0.000412 | +0.023307 | +0.002171 |
| Noisy | -0.000525 | +0.021063 | +0.000367 |
| Crowded | +0.037664 | +0.030860 | +0.008378 |

These are complete sampling-policy effects, not Fibonacci-versus-spiral
geometry effects: shells, neighbor rules, counts, retry and center handling
also differ. S0/S1 controls use the F0 shell radii, 10 candidates/shell,
no center, identical Fibonacci ownership/cone exclusion, and the same
estimator/interpolation, replacing only sphere geometry with the shared spiral
generator. Coarse S0 RMSE is 0.011244 A versus F0 0.011701; crowded S0 is
0.025694 versus F0 0.027900. Selected crowded counts differ slightly (1811
versus 1764) because changing directions changes which points pass the common
exclusion. Isolated controls match counts exactly. These small improvements
cannot explain the complete Q0 policy's crowded deterioration.

Q0 crowded sampling has 3948 selected candidates and 1984 fitting rows, versus
F0's 1764 and 1200. Successful shells can return 8..10 points. More retained
rows here accompany more neighbor-contaminated responses and worse truth
recovery; coverage alone is not accuracy. No shell failed in these selected
synthetic/real workloads; deliberately unavailable shells are covered by tests.

Q8 controls repeat the same center eight times with all other Q0 choices fixed.
Coarse Q0/Q8 mean Width biases are 0.009990/0.010581 A and RMSEs
0.010327/0.010967 A; crowded biases are 0.059084/0.059178 A and RMSEs
0.072629/0.071735 A. The median working Width variances change from
7.56e-7 to 4.97e-7 (coarse) and 2.64e-7 to 9.79e-9 (crowded). These are
conditional working covariances from `(X' diag(1/data_covariance) X)^-1`,
transformed to A/tau; they are not independently calibrated statistical
uncertainties. Duplicate centers are perfectly dependent, and map
interpolation also correlates rows. Smaller working variance does not prove
more information. FirstStage retains one center; original QScore retains eight.

P4 conclusion: Tricubic has lower Width bias/RMSE in these synthetic cases.
Q0 has a small isolated/noisy advantage but a substantial crowded regression.
There is insufficient evidence to promote a sampling or interpolation default.
