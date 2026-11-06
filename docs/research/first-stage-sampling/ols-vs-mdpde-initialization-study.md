# OLS versus MDPDE Width initialization sensitivity

## S0: baseline and contracts

Work began on the existing `develop` at
`2a1b565e5bf891ec6dcf15767ab45cd372eda3ee`, with a clean working tree.
Fetched `origin/develop` was identical. No branch/worktree is created and
Notion is not modified. Prior evidence is retained unchanged; its SHA-256
digests and real-data source checks are in `ols-vs-mdpde-validation.json`.

F0 remains deterministic Fibonacci shells (20 radii, 0.05..1.95 A, ten
points/shell), Tricubic, Ownership/30-degree cone, positive responses and
inclusive 0..1 A fitting. `PrepareJointWidthInitialization` samples once;
`RunJointFirstStageInitializationFromPreparedSamples` trains one local alpha
and calls one estimator per contributor. `PreparedLocalGaussianDesign` builds
the shared log-response design; `EstimateBetaMDPDE` returns both `beta_ols`
and `beta_mdpde`. Alpha is the MDPDE robustness parameter; OLS has no separate
trained alpha. Failed-only refinement may update MDPDE, preserving OLS.

`record.ols` and `record.mdpde` retain decoded A/Width/offset. The captured
native status, recovery and qualification belong to the MDPDE solver. They
do **not** independently qualify OLS. This study uses that explicitly named
common gate plus source-specific finite positive Width validity. Full parameter
validity (finite A/Width/offset, A>=0, Width>0) is recorded separately.

Production sets `initialization.tau` from MDPDE. `original_tau` retains finite
raw Width (including nonpositive values); `used_tau` is set only for a positive
finite seed. A target with invalid Width can receive the median of target
`reason==valid-width` raw positive finite Widths. This fallback donor pool is
not the calibrated donor pool and may include an unqualified target fit.
Halos never receive this fallback and do not affect `initialization.valid`.
An unqualified Halo has a diagnostic reason even if its Width is finite.
Current Pure Halo Support Exclusion fits targets only and reconstructs Halos
downstream through the unchanged geometry-local target median contract.

`CalibrateGlobalWidth` requires valid target starts, >=20 qualified target
donors and >=95% target coverage. A donor must pass
`IsQualifiedFirstStageEstimate`, have finite positive raw MDPDE Width, and
have neither initialization failure nor median fallback. Native and accepted
refined qualifications are counted separately. Fallback seeds never become
native donors. Type-7 p95 gives the single bound `4*p95`, with overflow/log
checks. `MakeBoundedInitialization` projects only target starts above the bound,
without mutating captured estimates. These Production rules remain unchanged.

The existing synthetic Width RMSEs (72 atom/phase records per scenario) are
confirmed: OLS/MDPDE Fine 0.000774/0.000905, Coarse 0.010902/0.011701,
Noisy 0.009699/0.015587, Crowded 0.029609/0.027900 A. Local truth ordering
alone does not establish a better JCE initialization.

The baseline regression checks all 24 coarse-phase-0 OLS/MDPDE vectors, alpha,
roles, qualification, original/used seeds, calibration, projection and final
Widths/objective against the original committed baseline. Fresh coarse and
both real-map F0 runs additionally compare every prior atom, calibration,
complete JCE evidence/state and QScore field exactly. Time/RSS are excluded.
JSON reference parsing uses precise numbers so decimal round trips do not
introduce a one-ULP comparison artifact.

## S1: derived Width sources

`tests/support/ResearchWidthInitialization.hpp` creates separate, run-owned
MDPDE/OLS selections from a const completed initialization. Each selected
record carries source identity, raw/used Width, Width and full parameter
validity, common qualification gate, fallback reason and donor exclusion.
MDPDE seeds use the exact Production vector/fallback. OLS uses its own finite
positive target Width median; no valid donors means unavailable, never an
invented seed or a switch to MDPDE. Halos retain roles and receive no fallback.
The transient projection bridge explicitly labels its research scope and
retains both original estimator arrays. It is not original FirstStage evidence.

Ten focused baseline/selection/calibration tests pass, covering deterministic
source choice, original-record immutability, nonfinite/zero/negative Widths,
parameter validity versus Width validity, same-source fallback, unavailable
source pools, Target/Halo roles and the existing donor-minimum contract.

## S2: orthogonal execution

The existing explicit research executable accepts the research ID `width`.
It builds one const problem, samples/trains/fits once, then constructs both
source selections and both calibrations. MM uses the actual Production wrapper;
OM/MO/OO call the same Production BCD implementation with an explicit bound.
No public CLI/FitOptions flag is added. The only Production source extraction
is `CalibrateQualifiedWidthDonors`: unchanged minimum/coverage/Type-7/4*p95/log
math, shared after caller-specific eligibility. Production donor selection is
unchanged. OLS calibration is labeled research with the common MDPDE gate.
An intersection-cohort control is produced whenever donor sets differ.

Twelve focused tests pass. Fresh coarse-phase-0 MM matches original atom,
calibration and complete JCE evidence/state exactly. All four cells converge
in two component sweeps; MM/OM share 1.8782695995625727 A and MO/OO share
1.8780533642794233 A bounds. Unprojected seeds match by source. The runner
asserts these factorial invariants and exact repeated common-bound endpoints.
First-sweep diagnostics use the existing observer, without a numerical solver
change. A/C start at zero, so the shared initial objective is uninformative
about Width seed quality. Solver timing starts at `SolverConfigured`, after
calibration/projection, for every cell. MM's wrapper repeats its inexpensive
Production calibration internally to preserve exact Production provenance;
this validation overhead is outside solver timing.

## S3: synthetic truth and JCE

Completed all 48 main runs (four unchanged scenarios, three deterministic
phases, four cells) plus four crowded-phase-0 common-bound controls. Every
endpoint passes the unchanged KKT, Independent Replay and Endpoint Certification
checks. All 12 MM atom/calibration/full JCE/QScore baselines match both prior
studies exactly. All designs retain the same F0 observations; fitting positive
rows are 2400/phase isolated and 1200/phase crowded. Every synthetic source has
24/24 native qualified calibration donors, no invalid Width and no fallback.
No optional contamination scenario was needed to answer the main question.

Width metrics pool 72 matched atom/phase records per scenario; phases are
sensitivity cases, not independent random trials. A is integrated amplitude.

| Scenario | Source | Width bias A | RMSE A | MAE A | A RMSE |
|---|---|---:|---:|---:|---:|
| fine | OLS | 0.000767121 | 0.000774272 | 0.000767121 | 0.001798583 |
| fine | MDPDE | 0.000888351 | 0.000905387 | 0.000888351 | 0.002609013 |
| coarse | OLS | 0.010733798 | 0.010902371 | 0.010733798 | 0.026111037 |
| coarse | MDPDE | 0.011487394 | 0.011700981 | 0.011487394 | 0.032598438 |
| noisy | OLS | 0.008067193 | 0.009699107 | 0.008272966 | 0.043720318 |
| noisy | MDPDE | 0.012769544 | 0.015587107 | 0.013052113 | 0.087842298 |
| crowded | OLS | 0.028517830 | 0.029608840 | 0.028517830 | 0.282168060 |
| crowded | MDPDE | 0.025609353 | 0.027900279 | 0.025609353 | 0.267184002 |

Paired absolute-Width error change is OLS minus MDPDE. Exact-sign fractions
are retained in JSON. The table additionally treats changes <=1e-12 A as ties,
so alpha-zero rounding does not count as a scientific improvement.

| Scenario | Mean delta A | Median A | OLS improved (>1e-12 A) | Tied | Worst OLS regression A | Worst MDPDE regression A |
|---|---:|---:|---:|---:|---:|---:|
| fine | -0.000121230 | -0.000085002 | 69.44% | 29.17% | 0.000031613 | 0.000685217 |
| coarse | -0.000753595 | -0.000289401 | 55.56% | 34.72% | 0.003470604 | 0.005805876 |
| noisy | -0.004779147 | -0.001224145 | 61.11% | 36.11% | 0.001611315 | 0.022895675 |
| crowded | 0.002908477 | 0.000000000 | 13.89% | 37.50% | 0.030169769 | 0.017732279 |

Largest OLS regression is crowded phase 0 atom 5 (+0.030170 A); MDPDE's
largest is noisy phase 1 atom 4 (+0.022896 A). At alpha=0, maximum Width
difference is 2.22e-16 A: the implementation still performs a weighted solve
with unit weights rather than requiring bitwise equality to its OLS solve.
All differences beyond 1e-12 A occur with alpha>0. Per-alpha metrics and
outliers are in JSON. OLS improves isolated/noisy average truth recovery;
MDPDE better handles this fixture's neighbor contamination. This is descriptive,
not a universal robustness theorem. No fitting rows are lost by positive
filtering in these F0 cases, so positive-row differences do not explain this
estimator comparison.

Each table entry is **maximum component sweeps / solver seconds / final
objective**. Runtimes have identical solver scopes and are single observations
for each deterministic case. All entries have KKT/Replay/Endpoint Pass.

| Scenario/phase | MM | OM | MO | OO |
|---|---|---|---|---|
| fine-r0 | 2 / 0.3758 / 1.63334972214e-07 | 2 / 0.3801 / 1.63334972214e-07 | 2 / 0.3757 / 1.63334972214e-07 | 2 / 0.3750 / 1.63334972214e-07 |
| fine-r1 | 2 / 0.3757 / 1.63216263417e-07 | 2 / 0.3766 / 1.63216263417e-07 | 2 / 0.3756 / 1.63216263417e-07 | 2 / 0.3752 / 1.63216263417e-07 |
| fine-r2 | 2 / 0.3778 / 1.63244281161e-07 | 2 / 0.3780 / 1.63244281161e-07 | 2 / 0.3769 / 1.63244281161e-07 | 2 / 0.3778 / 1.63244281161e-07 |
| coarse-r0 | 2 / 0.0469 / 1.6310440785e-07 | 2 / 0.0467 / 1.63104407847e-07 | 2 / 0.0465 / 1.6310440785e-07 | 2 / 0.0466 / 1.63104407847e-07 |
| coarse-r1 | 2 / 0.0469 / 1.63010944828e-07 | 2 / 0.0466 / 1.63010944826e-07 | 2 / 0.0464 / 1.63010944828e-07 | 2 / 0.0464 / 1.63010944826e-07 |
| coarse-r2 | 3 / 0.0485 / 1.63039069686e-07 | 3 / 0.0481 / 1.63039069692e-07 | 3 / 0.0480 / 1.63039069686e-07 | 3 / 0.0476 / 1.63039069692e-07 |
| noisy-r0 | 5 / 0.0709 / 0.0376988288475 | 5 / 0.0706 / 0.0376988288475 | 5 / 0.0704 / 0.0376988288475 | 5 / 0.0703 / 0.0376988288475 |
| noisy-r1 | 5 / 0.0729 / 0.037682939286 | 5 / 0.0721 / 0.037682939286 | 5 / 0.0724 / 0.037682939286 | 5 / 0.0718 / 0.037682939286 |
| noisy-r2 | 5 / 0.0741 / 0.0376909591757 | 5 / 0.0728 / 0.0376909591757 | 5 / 0.0735 / 0.0376909591757 | 5 / 0.0724 / 0.0376909591757 |
| crowded-r0 | 124 / 2.5935 / 2.98621142151e-08 | 124 / 2.6284 / 2.98621142196e-08 | 124 / 2.6203 / 2.98621142151e-08 | 124 / 2.5977 / 2.98621142196e-08 |
| crowded-r1 | 126 / 2.6519 / 2.91632563512e-08 | 126 / 2.6550 / 2.91632563949e-08 | 126 / 2.6563 / 2.91632563512e-08 | 126 / 2.6582 / 2.91632563949e-08 |
| crowded-r2 | 125 / 2.6428 / 2.9522700626e-08 | 125 / 2.6372 / 2.95227006384e-08 | 125 / 2.6351 / 2.9522700626e-08 | 125 / 2.6353 / 2.95227006384e-08 |

Final Width truth RMSE pooled across phases, and largest OM-versus-MM Width
difference. MO equals MM numerically and OO equals OM for all synthetic cases;
all bounds are inactive at their final endpoints. Calibration changes have no
observed synthetic endpoint effect here.

| Scenario | MM/MO final RMSE A | OM/OO final RMSE A | Max OM-MM Width A | Max normalized prediction difference |
|---|---:|---:|---:|---:|
| fine | 2.57107303e-07 | 2.57107277e-07 | 1.86462e-13 | 2.56864e-14 |
| coarse | 5.09087768e-06 | 5.09088787e-06 | 5.82621e-09 | 1.02425e-09 |
| noisy | 0.0025478897 | 0.00254789012 | 1.06571e-08 | 1.70956e-09 |
| crowded | 4.8438526e-06 | 4.84385289e-06 | 2.13893e-11 | 7.87775e-12 |

OM and MM require identical sweeps in all 12 cases: Fine 2/2/2, Coarse
2/2/3, Noisy 5/5/5, Crowded 124/126/125. First-sweep decreases, accepted
updates and search evaluations are retained per run; smaller initial truth
error does not consistently improve these path metrics or final recovery.
Initial A/C=0 objectives are shared and do not measure seed quality. The largest
synthetic endpoint Width change is only 1.07e-8 A. The synthetic truth advantage
therefore does not translate into a demonstrated JCE convergence advantage.

## S4: real maps, calibration and endpoint sensitivity

All eight primary real-map runs and eight common-bound runs completed. Both
120-target subsets use the original full-model neighbor/closure construction,
map values and preprocessing. File byte counts/SHA-256 match both earlier
studies; both MM complete numerical baselines are exactly equal. All 16
endpoints pass unchanged BCD/KKT/Replay/Endpoint rules. No real Width truth is
available, so neither a smaller bound nor a smaller Width proves accuracy.

The two source donor cohorts are identical: 120 common donors for 6Z6U;
119 common donors and one common unqualified fit (atom 23) for 8RQB. No
source-only donor occurs, so quantile changes here are estimator effects,
not cohort effects. Both sources have 120 valid Widths per map, zero fallback,
zero refined donors. In 8RQB valid Width and initialization availability do
not turn the unqualified atom into a donor. Matched-cohort controls are covered
by tests and are automatically emitted by the runner if cohorts differ.

Raw positive finite Width distribution, including unqualified estimates:

| Map | Source | p50 A | p90 A | p95 A | p99 A |
|---|---|---:|---:|---:|---:|
| 6z6u | OLS | 0.623572175 | 0.711822687 | 0.765333377 | 0.915289949 |
| 6z6u | MDPDE | 0.632081042 | 0.765298567 | 0.795576428 | 1.294110099 |
| 8rqb | OLS | 0.538074256 | 0.601859143 | 0.647571343 | 0.814430454 |
| 8rqb | MDPDE | 0.553382740 | 0.621134288 | 0.701491808 | 0.818278114 |

**6Z6U**. The table p95 is the qualified calibration quantile.

| Metric | MM | OM | MO | OO |
|---|---:|---:|---:|---:|
| Qualified donors | 120 | 120 | 120 | 120 |
| Coverage | 100.0000% | 100.0000% | 100.0000% | 100.0000% |
| Median fallback | 0 | 0 | 0 | 0 |
| Width p95 A | 0.795576428 | 0.795576428 | 0.765333377 | 0.765333377 |
| Width upper bound A | 3.182305712 | 3.182305712 | 3.061333508 | 3.061333508 |
| Bounded projections | 1 | 0 | 1 | 0 |
| Sweeps | 163 | 165 | 164 | 164 |
| Final objective | 0.012784319782223002 | 0.012784319782222695 | 0.012784319782222858 | 0.012784319782223001 |
| Global AC / Width KKT | 1.4486e-08 / 9.6824e-10 | 1.3184e-08 / 9.0081e-10 | 1.3736e-08 / 9.405e-10 | 1.4423e-08 / 9.8516e-10 |
| Replay / Endpoint | Pass / Pass | Pass / Pass | Pass / Pass | Pass / Pass |
| JCE seconds | 7.084179 | 7.170618 | 7.123300 | 7.128432 |
| Process peak RSS MB | 1469.563 | 1469.563 | 1469.563 | 1469.563 |

OLS bound minus MDPDE bound is -0.120972204 A (-3.8014%); the p95 difference is -0.030243051 A.

**8RQB**. The table p95 is the qualified calibration quantile.

| Metric | MM | OM | MO | OO |
|---|---:|---:|---:|---:|
| Qualified donors | 119 | 119 | 119 | 119 |
| Coverage | 99.1667% | 99.1667% | 99.1667% | 99.1667% |
| Median fallback | 0 | 0 | 0 | 0 |
| Width p95 A | 0.701762643 | 0.701762643 | 0.647987068 | 0.647987068 |
| Width upper bound A | 2.807050572 | 2.807050572 | 2.591948271 | 2.591948271 |
| Bounded projections | 0 | 0 | 0 | 0 |
| Sweeps | 95 | 95 | 95 | 95 |
| Final objective | 0.016970359134848675 | 0.016970359134849095 | 0.016970359134848675 | 0.016970359134849095 |
| Global AC / Width KKT | 6.2032e-08 / 9.9243e-10 | 6.2323e-08 / 9.9407e-10 | 6.2032e-08 / 9.9243e-10 | 6.2323e-08 / 9.9407e-10 |
| Replay / Endpoint | Pass / Pass | Pass / Pass | Pass / Pass | Pass / Pass |
| JCE seconds | 16.192022 | 16.640994 | 16.213977 | 16.228951 |
| Process peak RSS MB | 1579.418 | 1579.418 | 1579.418 | 1579.418 |

OLS bound minus MDPDE bound is -0.215102301 A (-7.6629%); the p95 difference is -0.053775575 A.

6Z6U's largest source discrepancy is atom 172: OLS 0.710901 A versus
MDPDE 13.830992 A. Both are finite positive estimates under the shared gate;
this is not an invalid-Width case. MDPDE seeds are projected by 10.648687 A
under the MDPDE bound and 10.769659 A under the OLS bound. OLS seeds need no
projection. 8RQB's largest source difference is atom 6 (OLS 0.557500,
MDPDE 0.794736 A); no seed projects under either bound. Per-atom outliers,
parameter validity, qualification and source/fallback evidence are retained.

### Common-bound control and independent audit

The common bound is the maximum available calibrated bound for that exact
dataset. Each constituent calibration independently meets >=20 donors and
>=95% coverage; the maximum provides a common feasible domain without
shrinking any calibrated domain. This is explicitly research-only and does
not acquire Production provenance. Identical-source repeated runs under the
common bound must reproduce the full endpoint and numerical evidence exactly.

| Dataset | Common bound A | MM/OM/MO/OO sweeps | Objective range |
|---|---:|---|---|
| crowded-r0 | 124 / 2.5935 / 2.98621142151e-08 | 124 / 2.6284 / 2.98621142196e-08 | 124 / 2.6203 / 2.98621142151e-08 | 124 / 2.5977 / 2.98621142196e-08 |
| 6z6u | 3.182305712 | 163/165/163/165 | 0.012784319782222695 .. 0.012784319782223002 |
| 8rqb | 2.807050572 | 95/95/95/95 | 0.016970359134848675 .. 0.016970359134849095 |

All common endpoints pass the original numerical criteria. Independent
support-based reconstruction audits final A/C/Width, prediction, objective,
active Width bound and normalized d(prediction)/d(log Width) columns. The
largest real objective reconstruction error is below 1e-15. No final Width
bound is active in any main/common run. Absence of final bound activity does
not prove the optimization path was independent of the bound.

6Z6U atom 221 starts with OLS/MDPDE Width 0.596216742/0.632147427 A.
Its nearest retained observation is 1.961676 A away; A is zero in every run.
The charge coefficient remains approximately 0.000149081, and its normalized
signal norm is approximately 8.587e-5. The signal is not zero, but its Width
derivative is tiny. The table's objective contribution is F(full prediction)
minus F(prediction with this atom removed), retaining other parameters.

| Run | Width A | A | C | Width-column norm | Objective contribution |
|---|---:|---:|---:|---:|---:|
| MM | 0.246788791 | 0.0 | 0.000149080861419 | 3.1144e-18 | -3.6867e-09 |
| OM | 0.250733151 | 0.0 | 0.000149081612166 | 8.2206e-18 | -3.6868e-09 |
| MO | 0.129805491 | 0.0 | 0.000149081412693 | 7.919e-54 | -3.6868e-09 |
| OO | 0.095005598 | 0.0 | 0.000149081108279 | 1.1192e-96 | -3.6867e-09 |
| MM-common | 0.246788791 | 0.0 | 0.000149080861419 | 3.1144e-18 | -3.6867e-09 |
| OM-common | 0.250733151 | 0.0 | 0.000149081612166 | 8.2206e-18 | -3.6868e-09 |
| MO-common | 0.246788791 | 0.0 | 0.000149080861419 | 3.1144e-18 | -3.6867e-09 |
| OO-common | 0.250733151 | 0.0 | 0.000149081612166 | 8.2206e-18 | -3.6868e-09 |

At identical MDPDE bounds, MM versus OM changes atom 221 by 0.003944 A,
with maximum normalized prediction change 7.045e-10. This is a direct seed
sensitivity, and both endpoints are certified. Holding the OLS seed fixed,
OM versus OO changes it by 0.155728 A, with normalized prediction change
6.323e-10. Their bounded OLS initial vectors are identical and unprojected:
this isolates a bound/path effect beyond initial clipping. With MDPDE seeds,
MM versus MO changes it by 0.116983 A; there the bound also changes the initial
projection of atom 172. Thus bound sensitivity is substantially larger than
seed sensitivity for this particular weak direction, without demonstrating
better scientific accuracy. A=0 alone does not prove unidentifiability, because
charge depends on Width; the directly computed Width-column norms provide
the relevant local sensitivity evidence.

Common-bound MM/MO and OM/OO reproduce exactly by seed source. 8RQB has maximum
seed endpoint Width difference 1.313e-8 A and normalized prediction difference
3.198e-10; changing the bound has no observed endpoint effect. Objective
agreement must not be described as universal parameter endpoint identity.

### Downstream and persistence

MM uses the unchanged original-grid peeling/Post-Joint Local/Group workflow:
6Z6U has all 33868 peeling rows, 120 Local results and 45 Group members;
8RQB has all 141457 rows, 120 Local and 42 Group members. Research-source and
common-bound results remain in independent JSON. The current persistence
validator requires Production calibration evidence and a recognized
initialization data scope. The projection-only bridge intentionally has a
research scope; raw solver qualification is insufficient to publish it as
Production evidence. No provenance is forged and no validator is bypassed.
A regression verifies rejection of raw research solver persistence, plus
MM downstream processing and exact codec round-trip. Existing persistence,
invalidation and downstream tests remain enabled for final CTest.

## S5: timing and acceptance scope

An optional internal `FirstStageTimingReport` measures the existing alpha and
estimator calls separately. It is absent for Production callers, adds no policy
or public CLI flag, and changes no estimator input/order/decision. Reported call
durations are sums across workers; this study uses one worker, so they are
sequential durations. Bookkeeping is the residual of the enclosing FirstStage
wall time. Only completed calls contribute; all study calls complete. Regression
checks both timed and default paths against the exact same original baseline.
All primary/common cases are revalidated after this measurement-only addition.

Sampling, alpha training, the complete OLS+MDPDE estimator call, initialization
adapter, calibration, bounded projection, solver and downstream have separate
scopes. Solver time is `SolverConfigured` through return, including search,
numerical qualification and assembly. It excludes initial calibration and
projection for every cell, including MM. The serial, fixed cell order may
introduce cache/order effects. Real timings are single observations, not
statistically reliable speedups. RSS is process high-water memory, cumulative
across all cells (including earlier MM downstream); it cannot rank cell-specific
peak allocations. Import dominates real-map memory. The adapter is linear in
the small contributor count and does not duplicate map storage.

This is a seed-source switch benchmark after a complete FirstStage fit. OLS
is read from the same call that computes MDPDE and uses the trained alpha;
no MDPDE or alpha-training cost has been avoided. A future OLS-only workflow
could have a different cost, but its savings are unimplemented and are **not**
included in these measurements. It would also need an independent qualification
contract, since the study currently depends on the MDPDE solve gate.

### Performance observations

Synthetic entries are medians across three deterministic phases; real entries
are single observations. Adapter/calibration times construct both sources;
FirstStage fitting is the complete estimator call, including OLS, MDPDE and
any existing failed-only refinement, without separating an invented OLS-only cost.

| Dataset | Sampling s | Alpha training s | FirstStage fitting s | Bookkeeping s | Adapter us | Both calibrations us |
|---|---:|---:|---:|---:|---:|---:|
| fine | 0.004304 | 0.047931 | 0.000545 | 0.000074 | 7.334 | 2.500 |
| coarse | 0.003971 | 0.051613 | 0.000586 | 0.000075 | 6.250 | 2.417 |
| noisy | 0.003886 | 0.080205 | 0.000690 | 0.000079 | 7.375 | 2.792 |
| crowded | 0.003905 | 0.058438 | 0.000539 | 0.000066 | 6.500 | 2.750 |
| 6z6u | 0.025500 | 0.392232 | 0.004059 | 0.000414 | 24.750 | 6.583 |
| 8rqb | 0.018454 | 0.401691 | 0.004535 | 0.000413 | 27.417 | 7.875 |

| Dataset | MM solver s | OM solver s | MO solver s | OO solver s | MM downstream s | Process peak MB |
|---|---:|---:|---:|---:|---:|---:|
| fine | 0.375778 | 0.378023 | 0.375696 | 0.375151 | 0.313180 | 146.162 |
| coarse | 0.046899 | 0.046662 | 0.046512 | 0.046617 | 0.082739 | 46.105 |
| noisy | 0.072946 | 0.072060 | 0.072440 | 0.071805 | 0.087185 | 46.055 |
| crowded | 2.642767 | 2.637150 | 2.635076 | 2.635316 | 0.081071 | 46.170 |
| 6z6u | 7.084179 | 7.170618 | 7.123300 | 7.128432 | 0.258470 | 1469.563 |
| 8rqb | 16.192022 | 16.640994 | 16.213977 | 16.228951 | 0.501686 | 1579.418 |

Bounded projection is measured separately per run (1.417..21.709 us for main runs). These tiny durations are observations, not a meaningful performance ranking.

The major FirstStage cost in these measurements is alpha cross-validation,
not the final estimator call. Both remain executed for every source-switch
comparison. OLS seed does not reduce sweeps, and the small timing differences
across deterministic phases do not establish a reliable speedup. Single real
OM observations are slower than MM on both maps; no statistical claim follows.

Alpha summary (one shared trained MDPDE alpha per atom; OLS does not have a
separately trained robustness parameter). All values span 0..1. Both decoded
parameter vectors are valid in every observed main-study atom record.

| Dataset | Atom/phase records | Alpha median | Alpha=0 records | Both Widths valid |
|---|---:|---:|---:|---:|
| fine | 72 | 0.30 | 21 | 72 |
| coarse | 72 | 0.20 | 25 | 72 |
| noisy | 72 | 0.10 | 26 | 72 |
| crowded | 72 | 0.20 | 27 | 72 |
| 6z6u | 120 | 0.30 | 23 | 120 |
| 8rqb | 120 | 0.30 | 22 | 120 |

### Early progress and solver work

Each entry is **first-sweep objective / total accepted local updates / total
search evaluations** (phase 0 for synthetic). First-sweep objectives sum each
component after its first complete sweep plus parameter-independent rows,
using the same parent observation scale. They do not alter solver decisions.

| Dataset | MM | OM | MO | OO |
|---|---|---|---|---|
| fine | 1.63335067989e-07 / 48 / 96 | 1.63334999358e-07 / 48 / 96 | 1.63335067989e-07 / 48 / 96 | 1.63334999358e-07 / 48 / 96 |
| coarse | 1.66508035448e-07 / 48 / 96 | 1.65387082733e-07 / 48 / 96 | 1.66508035448e-07 / 48 / 96 | 1.65387082733e-07 / 48 / 96 |
| noisy | 0.0376990360367 / 91 / 182 | 0.0376989301857 / 91 / 182 | 0.0376990360367 / 91 / 182 | 0.0376989301857 / 91 / 182 |
| crowded | 0.0224335564184 / 248 / 496 | 0.0226159841262 / 248 / 496 | 0.0224335564184 / 248 / 496 | 0.0226159841262 / 248 / 496 |
| 6z6u | 0.0273555717449 / 1610 / 4321 | 0.027394802789 / 1618 / 4374 | 0.0273558844713 / 1609 / 4366 | 0.0273948016104 / 1629 / 4366 |
| 8rqb | 0.0306223025663 / 929 / 2090 | 0.0307603703582 / 940 / 2090 | 0.0306223025663 / 929 / 2090 | 0.0307603703582 / 940 / 2090 |

OLS has a slightly lower first-sweep objective on isolated synthetic phase 0,
but a higher one on crowded phase 0 and both real maps. These small early
advantages do not reduce synthetic sweeps/evaluations. Real accepted-update
counts can differ even when the final objective agrees; there is no consistent
work reduction from the OLS seed.

### Production decision

**Retain F0 + MDPDE seed + MDPDE calibration.** The evidence combines a
scenario-dependent truth result (Outcome E) with insufficient promotion evidence
(retain Production, Outcome A). No adaptive/hybrid Production policy is added.

| Question | Evidence-based answer |
|---|---|
| Is OLS FirstStage Width more accurate? | On these isolated/fine/coarse/noisy fixtures, average bias/RMSE/MAE improve. Crowded truth favors MDPDE. There is no real-map Width truth. |
| Is OLS a better JCE seed (MM vs OM)? | No demonstrated advantage: identical synthetic sweeps, 165 vs 163 on 6Z6U, and 95 vs 95 on 8RQB. Final truth recovery is essentially unchanged. |
| Does OLS destabilize numerical qualification? | No observed failure in this matrix. Every available endpoint passes the original criteria. This is not a guarantee for untested maps. |
| Is OLS calibration better (MM vs MO, OM vs OO)? | It yields smaller bounds with identical donor cohorts, but no truth/convergence superiority. 6Z6U's weak endpoint changes substantially. |
| Are endpoints the same? | Synthetic and 8RQB agree closely, but are not forced bitwise identical by estimator. 6Z6U has different certified Widths along a weak direction despite almost identical predictions/objectives. |
| What comes from seed? | MM vs OM at the identical bound isolates the 6Z6U 0.003944 A weak-Width difference and two additional sweeps. |
| What comes from bound? | OM vs OO has identical bounded OLS seeds but a 0.155728 A weak-Width change. MM vs MO also includes altered initial clipping. This direction is more bound-sensitive than seed-sensitive here. |
| Is full OLS justified (MM vs OO)? | No. A scenario-dependent local truth benefit, no consistent JCE benefit and no realized fitting-cost reduction do not justify promotion. |
| Is adaptive/hybrid follow-up warranted? | A limited, preregistered follow-up is reasonable because noise/isolation and crowding favor different local estimates. It must evaluate contamination diagnostics, independent OLS qualification, additional maps/charge cases and repeated timings; these data do not establish a winning rule. |
| What should remain? | The source-selection adapter, shared calibration math, exact-baseline/factorial/fallback/projection/persistence regressions, optional call timer, sequential runner, compact source/endpoint evidence and this report. |

### Reproduction and limitations

```
cmake --build build/release --target first_stage_sampling_study --parallel 4
python3 -B tests/research/run_width_initialization_study.py
# Optional scopes: --scope synthetic or --scope real
# Reaggregate already saved detailed build outputs:
python3 -B tests/research/run_width_initialization_study.py --summarize-only
```

The runner uses committed FirstStage/JCE summaries even in a fresh checkout;
when old detailed build outputs are present, it additionally compares complete
raw sampling, JCE evidence/state and calibration records. The original committed
coarse and 6Z6U full baseline snapshots are always compared for those datasets.
Detailed observations and logs are reproducible ignored build artifacts; compact
JSON retains per-atom source evidence, final parameter vectors, qualification,
independent audits, truth errors and timing. Prior evidence is never overwritten.

Limitations are the unchanged normalized-Gaussian, no-charge synthetic truth,
three deterministic phases, two first-120-target real subsets, common MDPDE
qualification rather than independent OLS qualification, no optional impulse
contamination trial, single ordered real timings and cumulative RSS. Independent
prediction comparisons cover retained JCE fitting observations. Small local
Width-column norms diagnose weak local sensitivity, not a population-wide
statistical identifiability theorem. Production/source policies, fitting radius,
interpolation, neighbor rules, BCD/LM/KKT/Replay/Endpoint thresholds and downstream
numerical rules are unchanged.

### Final verification and Git provenance

Release and Debug full builds passed with at most four workers. CTest passed
21/21 groups in Release (43.12 s) and Debug (529.18 s), including FirstStage,
Width source selection, donor eligibility, median fallback, calibration,
projection, BCD convergence, KKT, Independent Replay, Endpoint Certification,
persistence/invalidation and Post-Joint downstream. Registration guards cover
all 1105 enabled GoogleTests, including seven new Width source tests, with
none disabled or missing. Focused source/calibration tests pass 13/13.
`lint_repo`, `lint_install_smoke`, runner syntax and `git diff --check` pass.
Three pre-existing missing-field-initializer warnings remain in
`DataObjectPersistence_test.cpp`; that unrelated file is unchanged.

All tests use the existing thresholds/timeouts. Builds use parallel=4 and CTest
parallel=1, with OMP_NUM_THREADS=2, OMP_THREAD_LIMIT=4, OPENBLAS_NUM_THREADS=1,
VECLIB_MAXIMUM_THREADS=1 and unchanged per-group overrides. Numerical study
processes run sequentially with one worker. After the optional timer addition,
all 68 endpoints were rerun: every non-time/RSS field of all 14 complete dataset
ledgers matches the pre-timer run exactly. All 14 MM runs also retain full
original F0 parity, and all nine prior research artifacts retain their SHA-256.
Log digests, source digests and acceptance counts are saved separately in
`ols-vs-mdpde-validation.json`.

Cleanup removes the temporary pre-timer comparison copies after recording
parity; final detailed observations/logs remain in the ignored build directory.
There are no experimental public flags, new dependencies, unused helpers,
temporary solver instrumentation or new Production numerical branches. The
optional FirstStage call timer is retained for the requested reproducible
performance decomposition, with the default path unchanged. Research endpoints
are not relabeled as Production persistence evidence. No unrelated work was
overwritten; no branch/worktree, detached HEAD, force push, reset or Notion edit
was used.

All changes are directly on the existing `develop`. Implementation commits:

| Commit | Verified unit |
|---|---|
| 78717218 | S0 baseline and contracts |
| 330bb501 | S1 immutable source selection |
| cee14b98 | S2 orthogonal seed/calibration controls |
| ba86ad80 | S3 synthetic truth and convergence |
| 557381c0 | S4 real maps, common bounds and persistence guard |
| 74c99985 | S5 separate alpha/fitting measurements and timer parity |
| 95a00747 | Fresh-checkout reproduction and tooling cleanup |

The final documentation commit and normal push receipt are reported in the
delivery message. Fetched origin/develop still equals the initial
`2a1b565e5bf891ec6dcf15767ab45cd372eda3ee` before that final commit/push;
the completed history is a safe fast-forward from it.
