# Fibonacci neighbor selection sensitivity study

> Archived research record. The following sections describe the completed studies
> at commit `503a6607674e068e4d3ef3db52aa770abcf3aea4`, including interfaces and
> tools that have since been removed. They are not current execution instructions.
> Production is fixed Fibonacci + Tricubic, Ownership + 30-degree cone, MDPDE
> Width initialization/calibration; QScore independently keeps Spiral + Trilinear.
> No adaptive/hybrid seed policy was promoted. Original JSON/source/runners are
> available through the [evidence archive](evidence-archive.md); current regression
> coverage and verification are in the [cleanup report](cleanup-report.md).
> Historical numerical tables and observations below are unchanged.

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

Full inspection of the 98-run original ledger found one pre-existing prose
discrepancy: the P4 README says all synthetic cells have 24/24 donors, but
`results.json` and detailed outputs agree that crowded Q0 phases 1 and 2
have 23/24 donors (one unqualified fit each; no median fallback). Phase 0
has 24/24. The table below uses the ledger, not that prose generalization.
All 98 prior endpoints have no failed checks. The requested RMSE numbers
remain consistent; original evidence is not edited.

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
matched outliers and selection diagnostics. The historical single-worker
runner and aggregation implementation are preserved in the archived commit.

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
have one invalid Width excluded from donors (atom 172, fitted Width 0,
median fallback). 8RQB is 119/0/1 for F0 and 120/0/0 for FO/FQ; its F0
unqualified atom 23 becomes qualified under both controls. Initial bound projections are 1/0/0 and 0/0/0 respectively.
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
the archived real-only runner invocation.

## E5: decision and acceptance

This is Case D in the tested crowded geometry: both relaxations worsen truth
recovery, and FO is less damaging than FQ. Production remains F0. The 30-degree
cone is beneficial in these tests; this does not establish that 30 degrees is
universally optimal. Ownership-only is not supported as a replacement. The
QScore distance rule increases row coverage but lowers target specificity and
increases Width/Amplitude bias here. FQ reproduces much of Q0's crowded
regression while preserving Fibonacci geometry/responses, supporting selection
as a substantial cause. It does not isolate every Q0 factor or prove selection
is the sole cause: Q0 additionally changes shell radii, geometry, counts,
center and adaptive retry.

No adaptive refill is implemented. The prerequisite truth advantage is absent,
and contamination rather than shortage of observations dominates this case.
No optimizer, fitting radius, qualification, calibration, public policy/CLI,
Two-Stage or downstream numerical rule was changed. Real-map inconsistency and
FQ's bound-sensitive 6Z6U endpoint further argue against promotion.

The runner stores shared mask/neighbor/contamination diagnostics once per
scenario and policy-specific fit/contamination summaries separately. Reaggregation
verified that compaction preserved every numerical field. Detailed observations
and logs remain reproducible build artifacts; only compact evidence is committed.

### Policy definitions and search completeness

For sample p around target x_i, r=|p-x_i|, d_ij=|x_j-x_i| and
 d_j=|p-x_j|, the rejection rules are:

- F0: any fixed-query non-hydrogen neighbor has d_j<r or theta_ij<30 degrees.
- FO: any of the same fixed-query neighbors has d_j<r.
- FQ: any non-hydrogen neighbor returned by the inclusive 2r query has d_j<0.9r.

The implementation preserves stored shell radii for QScore's query/exclusion
and F0's actual coordinate norms for ownership. Their tiny floating-point
representation differences do not change the tested inclusion results.
Equal distance is retained; equal search distance is included. F0's strict
cos(theta)>cos(30 degrees), coincident-position tolerance and neighbor ordering
remain unchanged. Boundary, invalid geometry, hydrogen, retry, deterministic
selection and parallel reproducibility tests are in the ordinary test suite.

For a shared neighbor universe, QScore rejection is a subset of ownership
rejection, which is a subset of F0 rejection (verified on all study samples).
For fitting radii <=1 A, every ownership-rejecting neighbor must have d_ij<2r
by the triangle inequality, and every QScore-rejecting neighbor has d_ij<1.9r.
Thus fixed 2.5 A discovery is complete for the fitting predicates, even though
its actual candidate set differs from each 2r query. All fitting masks match
shared-universe FQ; no full-policy inclusion violation occurred there.
For raw r>1 A, the neighbor sets and masks really differ; this is not a
universal inclusion guarantee. The regression suite includes a counterexample.

### Synthetic numerical tables

Counts are means per 24-target phase; errors pool all 72 target/phase records.
The effective design row count equals positive fitting rows here; it is not an
estimate of independent information. Contamination is mean rho on retained
fitting rows. Q0 uses the immutable prior ledger, with different spiral
geometry/radii/retry/center policy; no Q0 contamination measurement is invented.
Initialization time is the median of three runs, including research mask setup.


**Fine**

| Metric | F0 | FO | FQ | Q0 external reference |
|---|---:|---:|---:|---:|
| Selected fitting rows | 2400 | 2400 | 2400 | 1944 |
| Positive / effective fitting rows | 2400 | 2400 | 2400 | 1944 |
| MDPDE Width bias A | 0.000888 | 0.000888 | 0.000888 | 0.000655 |
| MDPDE Width RMSE A | 0.000905 | 0.000905 | 0.000905 | 0.000693 |
| MDPDE Amplitude RMSE | 0.002609 | 0.002609 | 0.002609 | 0.001842 |
| OLS Width RMSE A | 0.000774 | 0.000774 | 0.000774 | 0.000552 |
| Qualified donors | 24 | 24 | 24 | 24 |
| Median fallback | 0 | 0 | 0 | 0 |
| Initialization ms | 52.811500 | 52.928875 | 53.160292 | 56.035250 |
| Contamination mean rho | 0 | 0 | 0 | unavailable |
| Median local alpha | 0.300000 | 0.300000 | 0.300000 | 0.100000 |
| Median working Width variance | 6.6273e-09 | 6.6273e-09 | 6.6273e-09 | 8.32722e-09 |

**Coarse**

| Metric | F0 | FO | FQ | Q0 external reference |
|---|---:|---:|---:|---:|
| Selected fitting rows | 2400 | 2400 | 2400 | 1944 |
| Positive / effective fitting rows | 2400 | 2400 | 2400 | 1944 |
| MDPDE Width bias A | 0.011487 | 0.011487 | 0.011487 | 0.009990 |
| MDPDE Width RMSE A | 0.011701 | 0.011701 | 0.011701 | 0.010327 |
| MDPDE Amplitude RMSE | 0.032598 | 0.032598 | 0.032598 | 0.025432 |
| OLS Width RMSE A | 0.010902 | 0.010902 | 0.010902 | 0.008328 |
| Qualified donors | 24 | 24 | 24 | 24 |
| Median fallback | 0 | 0 | 0 | 0 |
| Initialization ms | 56.180250 | 56.543625 | 56.615667 | 56.275250 |
| Contamination mean rho | 0 | 0 | 0 | unavailable |
| Median local alpha | 0.200000 | 0.200000 | 0.200000 | 0.100000 |
| Median working Width variance | 7.21996e-07 | 7.21996e-07 | 7.21996e-07 | 8.11615e-07 |

**Noisy**

| Metric | F0 | FO | FQ | Q0 external reference |
|---|---:|---:|---:|---:|
| Selected fitting rows | 2400 | 2400 | 2400 | 1944 |
| Positive / effective fitting rows | 2400 | 2400 | 2400 | 1941.333333 |
| MDPDE Width bias A | 0.012770 | 0.012770 | 0.012770 | 0.012061 |
| MDPDE Width RMSE A | 0.015587 | 0.015587 | 0.015587 | 0.013747 |
| MDPDE Amplitude RMSE | 0.087842 | 0.087842 | 0.087842 | 0.075530 |
| OLS Width RMSE A | 0.009699 | 0.009699 | 0.009699 | 0.007665 |
| Qualified donors | 24 | 24 | 24 | 24 |
| Median fallback | 0 | 0 | 0 | 0 |
| Initialization ms | 84.002542 | 84.444667 | 85.095417 | 83.064958 |
| Contamination mean rho | 0 | 0 | 0 | unavailable |
| Median local alpha | 0.100000 | 0.100000 | 0.100000 | 0.200000 |
| Median working Width variance | 4.02866e-06 | 4.02866e-06 | 4.02866e-06 | 8.53516e-07 |

**Crowded**

| Metric | F0 | FO | FQ | Q0 external reference |
|---|---:|---:|---:|---:|
| Selected fitting rows | 1200 | 2166 | 2250 | 1984 |
| Positive / effective fitting rows | 1200 | 2166 | 2250 | 1984 |
| MDPDE Width bias A | 0.025609 | 0.046234 | 0.059538 | 0.059084 |
| MDPDE Width RMSE A | 0.027900 | 0.050862 | 0.065149 | 0.072629 |
| MDPDE Amplitude RMSE | 0.267184 | 0.603496 | 0.793883 | 0.931593 |
| OLS Width RMSE A | 0.029609 | 0.050819 | 0.064431 | 0.067017 |
| Qualified donors | 24 | 24 | 24 | 23.333333 |
| Median fallback | 0 | 0 | 0 | 0 |
| Initialization ms | 62.552958 | 94.308833 | 101.916750 | 95.182000 |
| Contamination mean rho | 0.078633 | 0.150306 | 0.199963 | unavailable |
| Median local alpha | 0.200000 | 0 | 0 | 0.250000 |
| Median working Width variance | 6.12054e-06 | 2.34092e-05 | 3.53005e-05 | 2.17977e-06 |

All F0/FO/FQ synthetic raw counts are 4800 per phase, with 10 distinct fitting
radii and rank two per atom. Positive retention is 100% in their fitting
domains. Per-atom OLS/MDPDE, alpha, ranks, radii, working covariance, donor
status and selection counts are retained in detailed build JSON; the compact
ledger retains per-atom fits and covariance summaries. The conditional working
covariances neither model interpolation dependence nor establish independent
information. In crowded data, covariance increases under both relaxations.

### Paired Width error changes (candidate minus F0)

| Scenario | Policy | Mean A | Median A | Improved | Unchanged | Worst regression A |
|---|---|---:|---:|---:|---:|---:|
| fine | FO | 0.000000 | 0.000000 | 0.00% | 100.00% | 0.000000 |
| fine | FQ | 0.000000 | 0.000000 | 0.00% | 100.00% | 0.000000 |
| coarse | FO | 0.000000 | 0.000000 | 0.00% | 100.00% | 0.000000 |
| coarse | FQ | 0.000000 | 0.000000 | 0.00% | 100.00% | 0.000000 |
| noisy | FO | 0.000000 | 0.000000 | 0.00% | 100.00% | 0.000000 |
| noisy | FQ | 0.000000 | 0.000000 | 0.00% | 100.00% | 0.000000 |
| crowded | FO | 0.020625 | 0.018253 | 16.67% | 0.00% | 0.058778 |
| crowded | FQ | 0.033929 | 0.030830 | 9.72% | 0.00% | 0.084522 |

Outliers are saved by phase/atom ID. The largest crowded FO regression is
phase 0 atom 19 (+0.058778 A); FQ is phase 2 atom 19 (+0.084522 A). Broad
regression across 60/72 and 65/72 records shows that a single outlier does not
explain the direction. These are matched deterministic sensitivity results,
not random-sample confidence intervals.

### Coverage and contamination detail

Crowded fitting mean absolute neighbor contribution is 0.029302 (F0),
0.059357 (FO), 0.066570 (FQ). Newly accepted observations carry real target
signal too, but have worse target specificity. Mean rho for FQ-only additions
beyond FO is 1.480407, with absolute neighbor signal 0.252583 and target
signal 0.178132. Rejected rows are more contaminated still (FO/FQ mean rho
2.969979/3.804140). All small-denominator fractions are zero. Analytic
contamination versus absolute Width-error Pearson correlations across crowded
atom/phase records are 0.626/0.929/0.883 for F0/FO/FQ; these are descriptive,
not independent-sample tests or a universal causal model. Isolated controls
have zero analytic neighbor contamination, so the correlation is undefined.

| Dataset | Raw F0 selected | FO selected | FQ selected | FQ new vs F0 | FQ rejected vs F0 | Full/shared FQ mask differences |
|---|---:|---:|---:|---:|---:|---:|
| crowded (per phase) | 1764 | 2730 | 3030 | 1266 | 0 | 0 |
| 6z6u | 14763 | 17108 | 17869 | 3434 | 328 | 417 |
| 8rqb | 14677 | 17200 | 18015 | 3660 | 322 | 430 |

The 328/322 newly rejected FQ real-map rows and 417/430 full/shared-mask
differences are outside the fitting domain. Candidate-set differences,
per-shell retained fractions, rejection counts (ownership and cone can
overlap), per-atom changes and FO/FQ differences are explicitly saved. Do not
interpret raw-domain inclusion as guaranteed by fitting-domain completeness.

### Real-map quantiles, runtime and memory

Global calibration donor quantiles (A); the bound is always the unchanged
4*p95 rule. Initial/final endpoint quantiles are separately retained per run.

| Dataset | Policy | p50 | p90 | p95 | p99 | Initial p95 | Final p50 | Final p95 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 6z6u | F0 | 0.632081 | 0.765299 | 0.795576 | 1.294110 | 0.795576 | 0.624846 | 0.785574 |
| 6z6u | FO | 0.645645 | 0.774746 | 0.836727 | 1.097254 | 0.835487 | 0.624846 | 0.785574 |
| 6z6u | FQ | 0.654644 | 0.786847 | 0.871393 | 1.106054 | 0.871349 | 0.625081 | 0.803529 |
| 8rqb | F0 | 0.553407 | 0.623142 | 0.701763 | 0.818569 | 0.701492 | 0.524013 | 0.635541 |
| 8rqb | FO | 0.557619 | 0.620466 | 0.698130 | 0.898623 | 0.698130 | 0.524013 | 0.635541 |
| 8rqb | FQ | 0.563868 | 0.636180 | 0.675597 | 0.941357 | 0.675597 | 0.524013 | 0.635541 |

| Dataset | Policy | Initialization s | Through JCE + input s | JCE s | Initial peak MB | JCE peak MB | Whole-process peak MB |
|---|---|---:|---:|---:|---:|---:|---:|
| 6z6u | F0 | 0.422 | 8.757 | 7.097 | 1496.302 | 1525.105 | 1525.105 |
| 6z6u | FO | 0.510 | 9.187 | 7.476 | 1497.596 | 1523.892 | 1523.892 |
| 6z6u | FQ | 0.540 | 9.116 | 7.283 | 1497.645 | 1523.925 | 1523.925 |
| 8rqb | F0 | 0.424 | 17.906 | 16.284 | 1586.545 | 1586.545 | 1586.545 |
| 8rqb | FO | 0.499 | 18.179 | 16.452 | 1586.725 | 1636.385 | 1636.385 |
| 8rqb | FQ | 0.515 | 18.215 | 16.481 | 1587.905 | 1639.318 | 1639.318 |

These are single-run real timings and process high-water RSS, not precision
speedup estimates. Import dominates RSS; per-policy JCE peak differences are
up to approximately 3.3%. Research FO/FQ initialization includes diagnostic mask
construction and is slower here. Production predicate extraction adds no
per-sample allocation. Whole-process RSS includes separate downstream
validation. F0 reruns its full workflow to check exact agreement; FO/FQ reuse
their captured fit directly through the shared postprocessing block. Therefore
the recorded downstream seconds have different work scopes and are not a
policy speed comparison. Through-JCE time and phase RSS have matching scopes.

### Common-bound numerical controls

| Dataset | Bound A | F0 / FO / FQ sweeps | Objective range |
|---|---:|---|---|
| crowded-r0 | 2.229871 | 124 / 124 / 124 | 2.986211421506671e-08 .. 2.986211434382782e-08 |
| 6z6u | 3.485572 | 163 / 164 / 163 | 0.01278426699981996 .. 0.01278431978222305 |
| 8rqb | 2.807051 | 95 / 95 / 95 | 0.01697035913484863 .. 0.01697035913484876 |

All common-bound runs retain the same observation domain/objective and pass
all original convergence/KKT/Replay/Endpoint checks. The initial objective is
0.397521864730 for 6Z6U and 0.429031571040 for 8RQB for every policy, since
Production initializes A=C=0. This does not assess Width seed quality. The
6Z6U audit above shows why objective agreement alone cannot certify parameter
agreement. JSON audits retain independent objective reconstruction, prediction
differences, changed per-atom Widths, Width-column norms and active bounds.

### Scope and promotion evidence

Evidence is limited to these Gaussian fixtures, three deterministic phases
and two representative 120-target real-map subsets. This is not whole-model
validation or a claim about real-map Width truth. No charge is present in the
synthetic fixture. No interpolation or angular-threshold sweep was performed.
The fixed-geometry comparison supports preserving the cone here, without
establishing an optimal angle. FO/FQ have insufficient promotion evidence on
truth RMSE, contamination, donor robustness or calibration/endpoint stability.
Neither is promoted; no adaptive refill follows from this evidence.

### Final verification and Git provenance

Release and Debug builds passed with at most four workers. Final CTest is
21/21 groups in Release (41.38 s) and 21/21 in Debug (519.19 s). Both
registration guards cover all 1098 enabled GoogleTests (nine new tests,
none disabled). This covers QScore retry/correlation, F0/Two-Stage, Joint
initialization, Width Calibration, BCD/KKT/Replay/Endpoint, persistence and
invalidation. `lint_repo` and `lint_install_smoke` passed. Debug reports three
pre-existing missing-field-initializer warnings in `DataObjectPersistence_test.cpp`;
that unrelated file is unchanged. No numerical rule or timeout was loosened.

CTest uses parallel=1, OMP_NUM_THREADS=2, OMP_THREAD_LIMIT=4,
OPENBLAS_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1; existing Joint-group
OMP_NUM_THREADS=1 settings remain. Study processes are sequential with one
numerical worker. Serial/four-worker tests are reproducibility checks, not
benchmark runs. No third-party dependency was added.

After all implementation changes, seven final baseline reruns (coarse:0 F0,
crowded:0/1/2 Q0 and S0) exactly match every prior policy, scenario, grid,
atom, calibration, initial objective, full JCE evidence/state and QScore field.
In particular Q0's two 23-donor phases reproduce exactly. This is in addition
to the 14 matched F0 study baselines. Final acceptance and SHA-256 checks that
all six original evidence files are unchanged are recorded separately in
`neighbor-selection-validation.json`. Timing/RSS are never parity assertions.

All E0–E5 requirements were executed; there is no unavailable real-map case or
outstanding required test. Cleanup removed generated study bytecode and
transient analysis scripts, deduplicated shared diagnostics, and retained only
the reusable research controls, numerical tests and compact evidence in Git.
Detailed numerical outputs/logs stay under the ignored
`build/release/neighbor-selection-study` directory. No temporary branch or
worktree existed to remove. No Notion document was changed.

Work stayed on the existing `develop`, initially clean at
`0e9661a81f3480629549361101b1e3609b97686b`, matching fetched origin/develop.
Implementation commits, in order:

| Commit | Verification unit |
|---|---|
| 74f258ec | E0 baseline lock |
| 27394cbd | E1 shared QScore predicate and parity |
| 2bfe3690 | E2 fixed Fibonacci FO/FQ controls and mask tests |
| 76724fa6 | E3 synthetic truth/contamination study |
| f43b562f | E4 real maps, common bounds, unchanged downstream reuse/persistence |
| 3c7036a8 | E5 compact evidence and promotion/refill decision |

The final documentation commit and normal `git push origin develop` receipt
are recorded in the delivery message. No force push, detached HEAD, branch,
worktree or destructive reset is used. All source changes are confined to the
shared pure distance predicate, reuse of existing postprocessing, explicit
research tooling and necessary numerical tests. Production remains F0.
