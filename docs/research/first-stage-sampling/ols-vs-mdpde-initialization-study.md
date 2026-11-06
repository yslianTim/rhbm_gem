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
