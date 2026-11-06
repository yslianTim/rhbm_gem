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

## E2: fixed-geometry controls

Research-only `FibonacciNeighborSelection.hpp` builds FO (strict ownership with
F0's fixed 2.5 A candidates) and FQ (QScore's per-shell 2r candidates, strict
squared 0.9r exclusion), excluding hydrogen. The explicit study IDs FO/FQ copy
F0 raw observations and change only selection flags through the existing
editor, then call the unchanged prepared-sample FirstStage estimator. No
Production enum, public CLI option, refill or repeated center was added.

Four new GoogleTests cover fixed geometry/responses/order, deterministic and
1/4-worker estimator parity, independent distance checks, dense neighbors,
strict ownership/distance/search boundaries, hydrogen/empty sets and a raw
r>1 A counterexample to universal full-policy inclusion. All 55 focused tests
passed. F0 atom/calibration/JCE/QScore fields still match E0 exactly. Detailed
research output adds per-point masks/responses and per-shell neighbor IDs;
synthetic Gaussian contributions are diagnostics only, outside the estimator.

## E3: synthetic truth and contamination

Completed 36 Production-calibrated runs (four scenarios, three deterministic
phases, F0/FO/FQ) and three crowded phase-0 common-bound runs. Every endpoint
passes the original KKT/Replay/Endpoint checks. All 12 new F0 runs exactly
match the previous detailed numerical records. Fine/coarse/noisy controls
match across all three policies. All designs have rank two; every synthetic
run has 24 qualified donors and no fallback. Four sensitivity tests pass,
including a new truth/contamination regression independent of solver agreement.

Pooled crowded Width RMSE is F0 0.027900, FO 0.050862, FQ 0.065149 A;
Amplitude RMSE is 0.267184, 0.603496, 0.793883. Selected positive fitting rows
per phase rise from 1200 to 2166 and 2250. FO/FQ improve absolute Width error
on only 16.67%/9.72% of matched observations. Mean paired changes are
+0.020625/+0.033929 A; worst regressions +0.058778/+0.084522 A.

Mean analytic neighbor signal on F0 retained fitting rows is 0.029302;
FO/FQ newly accepted rows have 0.096691/0.109163. Their mean contamination
ratios are 0.239340/0.338626 versus F0's 0.078633. FQ's additions beyond FO
have mean ratio 1.480407. Contributions use only the fixture's normalized,
2.5 A supported Gaussian, with no invented charge. Epsilon is 1e-12 times
that atom's largest absolute raw analytic target contribution; the ledger
reports the fraction with abs(target)<=epsilon (zero here), absolute neighbor
signal, per-shell and per-atom summaries. Truth never enters fitting or
selection. More retained rows have not established more independent information.

The independent `neighbor-selection-results.json` retains per-run solver,
qualification, covariance, timing and coverage summaries, per-atom fits,
matched outliers and selection diagnostics. Reproduce with
`python3 tests/research/run_neighbor_selection_study.py` (single numerical
worker), or `--summarize-only` to aggregate saved detailed build outputs.

## E4: real maps, calibration and endpoints

Completed six real Production runs, six real common-bound runs and independent
prediction/objective audits for both datasets. SHA-256 digests match P0–P6's
original sources; the same first 120 non-hydrogen targets and full model
neighbors are used. Both F0 real runs exactly match all prior numerical fields.

| Dataset | Policy | Donors | Invalid Width / fallback | Width p95 A | Bound A | Sweeps | Final objective |
|---|---|---:|---:|---:|---:|---:|---:|
| 6Z6U | F0 | 120 | 0 / 0 | 0.795576 | 3.182306 | 163 | 0.0127843197822230 |
| 6Z6U | FO | 119 | 1 / 1 | 0.836727 | 3.346908 | 165 | 0.0127843197822228 |
| 6Z6U | FQ | 119 | 1 / 1 | 0.871393 | 3.485572 | 163 | 0.0127842669998200 |
| 8RQB | F0 | 119 | 0 / 0 | 0.701763 | 2.807051 | 95 | 0.0169703591348487 |
| 8RQB | FO | 120 | 0 / 0 | 0.698130 | 2.792519 | 95 | 0.0169703591348488 |
| 8RQB | FQ | 120 | 0 / 0 | 0.675597 | 2.702388 | 95 | 0.0169703591348486 |

All KKT, Replay, Endpoint and BCD checks pass without qualification changes.
6Z6U native/refined/unqualified fits are 120/0/0 for all policies, but FO/FQ
have one invalid Width excluded from donors. 8RQB is 119/0/1 for F0 and
120/0/0 for FO/FQ. Initial bound projections are 1/0/0 and 0/0/0 respectively.
Thus native fit success must still be distinguished from valid Width donation.

The unchanged Post-Joint block was extracted into `RunJointPostprocessing`,
so research-prepared endpoints reuse the exact Production peeling/Local/Group
implementation. There is no downstream algorithm change. All six Production
runs retain 120 available Local estimates; Group members are 45 for 6Z6U and
42 for 8RQB. All 33868/141457 original-grid peeling samples are available.
A new 24-donor regression checks all three policies, original-grid downstream
processing and SQLite round trips of actual raw selection flags. Its first
2-atom draft correctly failed the unchanged donor minimum; the fixture was
corrected to 24 donors rather than relaxing calibration. All five sensitivity
tests now pass.

Common bounds are research-only maxima across F0/FO/FQ: 3.485572 A (6Z6U)
and 2.807051 A (8RQB). At 6Z6U atom 221, F0 Production Width is 0.246789 A;
FO Production is 0.103922 A. F0/FO under the common bound are 0.205448/
0.076554 A. These directions have vanishingly small Width-column norms.
FQ instead reaches the active 3.485572 A bound under either mode, A=0,
C=0.001613825, normalized Width-column norm 3.89869e-4. Its objective is
lower by 5.27824e-8 but its maximum normalized prediction difference is
7.19501e-5 and maximum Width difference is 3.238784 A. This is a different
certified endpoint, not evidence of Width truth improvement. At the same
bound, F0 and FO still differ from FQ: initialization effects remain.
Holding F0 initialization while changing its bound also changes atom 221:
calibration effects are separately present. A=0 alone does not prove lack of
Width information when charge depends on Width.

8RQB has no active Width bound and maximum endpoint Width difference <=
1.722e-8 A, normalized prediction difference <=4.923e-10. No real-map Width
truth is available. The opposite donor-coverage changes across datasets do
not establish robust improvement. Reproduce real runs/audits with
`python3 tests/research/run_neighbor_selection_study.py --real-only`.
