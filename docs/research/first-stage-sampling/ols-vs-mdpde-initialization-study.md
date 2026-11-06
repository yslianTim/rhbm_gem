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
| fine-r0 | 2 / 0.3768 / 1.63334972214e-07 | 2 / 0.3754 / 1.63334972214e-07 | 2 / 0.3763 / 1.63334972214e-07 | 2 / 0.3746 / 1.63334972214e-07 |
| fine-r1 | 2 / 0.3753 / 1.63216263417e-07 | 2 / 0.3762 / 1.63216263417e-07 | 2 / 0.3761 / 1.63216263417e-07 | 2 / 0.3756 / 1.63216263417e-07 |
| fine-r2 | 2 / 0.3784 / 1.63244281161e-07 | 2 / 0.3786 / 1.63244281161e-07 | 2 / 0.3782 / 1.63244281161e-07 | 2 / 0.3779 / 1.63244281161e-07 |
| coarse-r0 | 2 / 0.0473 / 1.6310440785e-07 | 2 / 0.0466 / 1.63104407847e-07 | 2 / 0.0466 / 1.6310440785e-07 | 2 / 0.0464 / 1.63104407847e-07 |
| coarse-r1 | 2 / 0.0468 / 1.63010944828e-07 | 2 / 0.0467 / 1.63010944826e-07 | 2 / 0.0463 / 1.63010944828e-07 | 2 / 0.0463 / 1.63010944826e-07 |
| coarse-r2 | 3 / 0.0484 / 1.63039069686e-07 | 3 / 0.0476 / 1.63039069692e-07 | 3 / 0.0480 / 1.63039069686e-07 | 3 / 0.0474 / 1.63039069692e-07 |
| noisy-r0 | 5 / 0.0710 / 0.0376988288475 | 5 / 0.0705 / 0.0376988288475 | 5 / 0.0703 / 0.0376988288475 | 5 / 0.0703 / 0.0376988288475 |
| noisy-r1 | 5 / 0.0730 / 0.037682939286 | 5 / 0.0720 / 0.037682939286 | 5 / 0.0742 / 0.037682939286 | 5 / 0.0741 / 0.037682939286 |
| noisy-r2 | 5 / 0.0738 / 0.0376909591757 | 5 / 0.0727 / 0.0376909591757 | 5 / 0.0735 / 0.0376909591757 | 5 / 0.0726 / 0.0376909591757 |
| crowded-r0 | 124 / 2.5959 / 2.98621142151e-08 | 124 / 2.5945 / 2.98621142196e-08 | 124 / 2.5935 / 2.98621142151e-08 | 124 / 2.6043 / 2.98621142196e-08 |
| crowded-r1 | 126 / 2.6570 / 2.91632563512e-08 | 126 / 2.6576 / 2.91632563949e-08 | 126 / 2.6587 / 2.91632563512e-08 | 126 / 2.6558 / 2.91632563949e-08 |
| crowded-r2 | 125 / 2.6408 / 2.9522700626e-08 | 125 / 2.6363 / 2.95227006384e-08 | 125 / 2.6365 / 2.9522700626e-08 | 125 / 2.6357 / 2.95227006384e-08 |

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
