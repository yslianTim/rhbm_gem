# Fibonacci neighbor selection sensitivity study

## E0: locked baseline (2026-11-06)

Initial and fetched `origin/develop`: `0e9661a81f3480629549361101b1e3609b97686b`.
Branch `develop`, clean working tree; no pre-existing uncommitted work. Fetch
and fast-forward check succeeded (already current). No branch/worktree created.
Notion is not modified. Existing P0–P6 evidence is retained without edits.

Prior pooled MDPDE Width RMSE (72 atom/phase observations per scenario):

| Scenario | F0 | Q0 |
|---|---:|---:|
| Fine | 0.000905 | 0.000693 |
| Coarse | 0.011701 | 0.010327 |
| Noisy | 0.015587 | 0.013747 |
| Crowded | 0.027900 | 0.072629 |

S0 crowded RMSE is 0.025694 A. These agree with the prior README/results.json.
Before edits, fresh coarse:0 F0 and crowded:0 Q0/S0 runs exactly matched all
previous atom, calibration, full JCE evidence/state and QScore fields.
Timing/RSS are excluded from numerical parity. Detailed prior outputs and both
real-map sources are available locally; no baseline data need reconstruction.

Pipeline inspection: `SampleAtomMapValuesFromPoints` queries non-hydrogen
neighbors via `atom.FindNeighborAtoms()` (2.5 A), applies the existing strict
ownership/cone comparisons, then interpolates every raw point, including
unselected points. Fibonacci has 20 shells, 10 ordered points per shell,
radii 0.05..1.95 A. `ModelAnalysisEditor` owns moved raw vectors;
`GetSamplingEntries(First)` copies selected entries. Alpha training and
`PreparedLocalGaussianDesign` retain the existing inclusive 0..1 A fitting
range and positive-response handling. Designs' indices refer to those local
copies, never original-grid observations. Actual raw coordinates, distances,
responses and flags follow the existing persistence contract. Initialization
uses the same prepared-sample estimator, qualification and fallback; calibrated
Width bounds feed unchanged JCE. Raw replacement uses the editor's existing
invalidation rules and leaves unrelated original-grid Post-Joint state intact.

`SampleQScoreShell` queries the owner model at 2r, excludes hydrogen, compares
squared point-neighbor distance strictly against (0.9r)^2, retries at most 50
times (+2 candidates/attempt), and emits all accepted points only on shell
success. Original QScore uses spiral geometry, trilinear interpolation and
eight center observations. Q0 uses tricubic, one center and those shell/retry
rules. It is an external reference, not a selection-only control.

Required invariants: F0 raw positions/order/responses/flags, alpha, OLS/MDPDE,
calibration, QScore, JCE evidence and endpoint must retain exact numerical
parity. Q0/S0 existing numerical evidence must also remain unchanged. FO/FQ
will copy F0 raw samples and modify only `is_selected`, without changing
geometry, interpolation, fitting radius, estimator or Production default.

## E1: shared distance predicate

`sample_filter::IsInsideNeighborExclusion` accepts the already squared exclusion
radius and retains QScore's axis accumulation order and strict `<` comparison.
Geometry validation rejects nonfinite inputs. Discovery, hydrogen exclusion,
retry, shell success, interpolation and center weighting remain in QScore.
48 focused SampleFilter/QScore tests passed, including three new boundary,
off-grid and invalid-geometry tests and existing hydrogen/multiple-neighbor,
empty-neighbor and retry cases. Crowded Q0 atom/calibration/JCE/QScore fields
remain exactly equal after extraction.
