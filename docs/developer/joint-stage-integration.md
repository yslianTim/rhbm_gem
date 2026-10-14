# Joint stage integration

This guide describes the current Joint-to-group workflow and its data ownership.
Joint numerical fitting and qualification are defined in the
[Runtime guide](joint-component-runtime.md); SQLite payload structure is
documented in the [data I/O architecture](architecture/dataobject-io-architecture.md).

## Result ownership

The estimator-neutral `FittingStage::Second` record remains the Joint endpoint:
A, B, C, target/halo role, component identity, run identity, convergence,
qualification and source provenance. It is used to separate overlapping
signals and to retain Joint diagnostics. Downstream group inference does not
consume its A/B values as measurements.

Each selected target stores a separate `PostJointLocalResult`, keyed by the
same atom and exact `EstimateSource`. It contains local alpha, OLS and MDPDE
results, local solver diagnostics, the fixed Joint charge profile, sample
counts and failures, status, reason, and compact `RHBMInformation` (precision
and information moment). The local result has its own lifecycle and never
writes into the Joint endpoint or objective. Halo atoms contribute to peeling
but do not receive local or group results.

## Current data flow

1. `RunJointComponentWorkflow` fits the retained selected targets and eligible
   halos with the existing FullABC Joint estimator and publishes their Second
   stage endpoints.
2. `BuildPostFitPeelingSamples` uses raw sample entries and their saved
   geometry. For each target sample it subtracts the tricubic interpolation
   of every other retained contributor on the same Joint observation rows:
   `y_peel = y_raw - I_k[sum(j != i) fhat_j]`. Each `PostFitPeelingResult`
   keeps one entry per raw sample, including an unavailable reason when the
   stencil leaves the Joint domain or a contributor state is missing.
   `ApplyPostFitPeelingsBatch` validates target, run and sample identities
   before publishing the peeling results.
3. `BuildPostJointLocalSamples` makes a per-target `PostJointLocalSamples`
   view. It keeps raw and peeled values paired by raw sample index and
   geometry. It computes the target's own fixed nuisance profile
   `q_i^J(r) = C_i^J K(r; B_i^J)` on Joint grid support and uses the same
   tricubic interpolation convention as peeling. The Gaussian response is
   `z_ik = y_peel_ik - I_k[q_i^J]`. The profile's Joint width remains fixed
   even when the fitted local width differs.
4. `FitPostJointLocalTarget` passes those charge-corrected samples to the
   existing `TrainAlphaR`, `PreparedLocalGaussianDesign`, and
   `EstimateLocalGaussian`. The design is the same log-Gaussian design used by
   first-stage MDPDE. Nonpositive responses are retained in the paired sample
   record and excluded by the existing log-response dataset builder.
   Local rank, positive sample count, positive A/B, native solver status,
   variance/weights, source identity and numerical diagnostics all contribute
   to local eligibility. Failed or unqualified fits remain unavailable with
   explicit reasons.
5. `RunJointGroupPotentialFitting` trains group alpha from eligible
   post-Joint local MDPDE points. It uses those same points with the compact
   information stored by each `PostJointLocalResult`; it does not repeat local
   MDPDE fitting. `EstimateMuMDPDE` and `EstimateWEBFromInformation` provide
   group means, prior/posterior parameters, uncertainty and outlier flags.
   Only selected targets enter group statistics.

The group information for local design `X_i`, log-response `y_i`, and
working diagonal covariance `Sigma_y_i` is
`P_i = X_i^T Sigma_y_i^-1 X_i` and
`h_i = X_i^T Sigma_y_i^-1 y_i`. It comes from local MDPDE fitting data and
weights, not a Joint-endpoint conditional least-squares covariance. The
persisted information is six doubles per eligible target; per-sample
weights/covariances are transient and are not serialized.

## Invalidation and provenance

All derived records carry the Joint source and run. Editors reject mixed Joint
runs and mismatched atom, component or role identities.

- A Joint endpoint change clears its peeling, post-Joint local result, related
  group posterior and Joint snapshot for that run.
- A raw sample change clears peeling, local result and group inference while
  retaining a still-valid Joint numerical snapshot.
- A peeling change clears the local result and group inference.
- A local MDPDE result change clears affected group inference only.
- A group alpha change clears that group's posterior.
- Halo results remain peeling contributors and are excluded from local fitting,
  group alpha training and group membership statistics.

These changes are batched by affected group so each group is invalidated once.
Missing samples never become zero-valued observations.

## Resource and numerical limits

Peeling and charge correction use bounded per-target preparation and avoid a
dense observations-by-contributors matrix. `PostJointLocalSamples` is temporary;
only compact fit results, diagnostics, fixed charge profile and 2-by-2
information are retained. Group WEB works from those per-target summaries.

The fixed charge correction is an approximation to the FullABC Joint model:
it holds `C_i^J` and `B_i^J` fixed while the local Gaussian fit estimates its
own A/B. Tricubic interpolation, finite Joint support and the log transform's
positive-response requirement remain part of the data contract. The local
measurement covariance is a working model derived from the local MDPDE
observations; it does not propagate uncertainty in other Joint contributors or
cross-atom errors.

## Persistence and coverage

SQLite schema v20 stores neutral stage-result JSON format 3. It round-trips
the Joint snapshot, paired peeling, post-Joint local result, fixed charge
profile, compact information, status/reason, and group posterior. Version-2
Joint group results based on the previous measurement source are invalidated
during migration and require the post-Joint local workflow to run again.
Legacy stage formats remain readable through their existing migration path.

Permanent coverage exercises exact grid peeling, fixed-charge correction,
local MDPDE regression, group information equivalence, identity and
invalidation checks, persistence, and two-stage compatibility. The complete
test lane is described in [tests/README.md](../../tests/README.md).
