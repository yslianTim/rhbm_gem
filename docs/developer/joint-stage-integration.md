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

1. The production `RunJointComponentWorkflow` builds
   `BuildPureHaloSupportExclusionProblem`. Supported non-target atoms are used
   only to exclude their support rows from the fitting domain; they are not
   Joint parameters and no FirstStage halo `(A, B, 0)` snapshot is required.
   Joint fits the selected targets' FullABC `(A, B, C)` endpoints and publishes
   their Second-stage diagnostics. The explicit
   `EstimateJointComponentsFixedHalo` / `BuildFixedHaloResearchProblem` APIs
   remain historical Fixed-Halo research/compatibility paths.
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

`RunJointComponentWorkflow` creates one immutable `JointPostprocessingContext`
after publishing the trusted Joint snapshot and target endpoints. It builds
the observation-grid-index to row-index lookup and atom-ID to contributor-index
lookup once for that Joint run. Construction checks the snapshot contributor
IDs and complete observation row IDs against the same `JointProblem` once;
`BuildPostJointLocalSamples` does not repeat that observation-sized identity
comparison. It validates the selected target, endpoint and peeling provenance,
sample pairing and geometry, and that target's support stencils. The context
reads saved component representatives and endpoint provenance without
rebuilding Joint stage estimates. Its borrowed map, model endpoints, problem
input and snapshot must outlive the context and remain unchanged through local
preparation. It is scoped to one Joint run and is not reused across runs;
publishing peeling/local derived state does not change those inputs. Each target
creates only its own support charge lookup and sample workspace. The context
does not hold a dense observations-by-contributors matrix or a per-target copy
of contributor state.

The group information for local design `X_i`, log-response `y_i`, and
working diagonal covariance `Sigma_y_i` is
`P_i = X_i^T Sigma_y_i^-1 X_i` and
`h_i = X_i^T Sigma_y_i^-1 y_i`. It comes from local MDPDE fitting data and
weights, not a Joint-endpoint conditional least-squares covariance. The
persisted information is six doubles per eligible target; per-sample
weights/covariances are transient and are not serialized.

## Evidence and posterior representation

`AssessPostJointLocalResult` is the shared contract used when fitting local
evidence, publishing a batch, selecting group observations, and validating
stage persistence. It checks target/source/run/component identity, Joint
qualification, paired peeling and raw sample counts, fixed charge equality,
sample failure indices, finite local A/B and alpha, native MDPDE diagnostics,
effective sample count, and finite positive-definite compact information.
Unavailable evidence remains excluded from group inference. Publication
validates the whole batch before changing any atom or group result.

Group posterior Gaussian parameters and the fixed Joint charge profile are
stored separately. Each `GroupGaussianMemberResult` contains posterior A/B
with Gaussian offset zero and that member's original `(C_i^J, B_i^J)` profile.
`EvaluateGroupMemberResponseAtDistance` evaluates
`A_post G(r; B_post) + C_i^J K(r; B_i^J)`. Group means, priors, and MDPDE
summaries contain Gaussian parameters only. A group charge median is a
descriptive statistic and does not define a physical group-level FullABC
response. Plots label peeled samples as including the target-own Joint charge
and group priors/means as Gaussian-only with `C = 0`. Their overlay is a view of
different signal representations, not a measurement-versus-fit comparison.
Raw-map and neighbor-signal curves are labeled separately. Plots do not
reconstruct charge-corrected samples from stored peeling data; a full member
response uses that member's persisted fixed charge profile.

## Invalidation and provenance

All derived records carry the Joint source and run. Editors reject mixed Joint
runs and mismatched atom, component or role identities.

- A Joint endpoint change clears the run's peeling, post-Joint local results,
  group inference and Joint snapshot. Any contributor can affect every target
  whose peeling subtracts its signal.
- A raw sample change clears peeling, local result and group inference while
  retaining a still-valid Joint numerical snapshot.
- A peeling change clears the local result and group inference.
- A local MDPDE result change clears affected group inference only.
- A group alpha change clears that group's posterior.
- Halo results remain peeling contributors and are excluded from local fitting,
  group alpha training and group membership statistics.

These changes are batched by affected group so each group is invalidated once.
Missing samples never become zero-valued observations.

## Qualification and downstream availability

Initialization validity, Joint search completion, endpoint verification,
independent replay, and target qualification are evaluated before any
qualified-dependent downstream stage. An initialization failure or an
unqualified Joint endpoint is saved with its diagnostics and stop reason; Halo
median reconstruction, peeling, Post-Joint Local, and Group Fitting are not
started. An unqualified candidate never becomes a formal target estimate.

For a qualified Pure Exclusion result, Halo median reconstruction is a
separate original-domain step used only to build peeling responses. If no
qualified target donor can produce a valid median, the Joint endpoint remains
qualified but the reconstruction evidence is `Unavailable`; downstream
reconstruction-dependent stages are skipped with that explicit reason. This
state is distinct from Joint numerical non-convergence. A qualified result
with no Halo contributors does not require a median snapshot.

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

SQLite schema v20 stores neutral stage-result JSON format 4. It round-trips
the Joint snapshot, paired peeling, post-Joint local result, fixed charge
profile, compact information, status/reason, and group posterior. Version-3
stage documents are decoded with their per-member fixed charge reconstructed
from available local evidence, then emitted as format 4. Version-2
Joint group results based on the previous measurement source are invalidated
during migration and require the post-Joint local workflow to run again.
Legacy stage formats remain readable through their existing migration path.

Permanent coverage exercises shared mapping counts, exact grid peeling,
fixed-charge correction, local MDPDE regression, group information
equivalence, available/unavailable evidence, atomic publication, identity and
cross-target invalidation, persistence, and two-stage compatibility. The complete
test lane is described in [tests/README.md](../../tests/README.md).
