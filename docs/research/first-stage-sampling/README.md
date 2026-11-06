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
