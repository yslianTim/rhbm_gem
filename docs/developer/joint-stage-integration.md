# Joint Stage Integration

Joint Components is the default potential estimator. The
[command contract](commands/potential_analysis.md) also exposes `two-stage`.
Potential fitting excludes hydrogen from targets and eligible contributors.
The Joint workflow builds the Original and Fitting Domains, initializes
Target widths, runs the [Joint estimator](joint-component-runtime.md), and
captures a snapshot. A qualified snapshot feeds Halo reconstruction,
Original-grid peeling, Post-Joint Local fitting and Group Fitting.

## FirstStage initialization

The Joint workset resolves selected Target identities from the Fitting
Domain. FirstStage uses 20 deterministic shells at 0.05, 0.15, …, 1.95 Å,
with ten points per shell. C++ callers can set `FitOptions::first_stage_sampling_method`
to `SphereSamplingMethod::SpiralDeterministic`. The default is
`SphereSamplingMethod::FibonacciDeterministic`. Both methods apply the
strict Ownership and 30° cone exclusions for non-hydrogen neighbors.
Tricubic samples feed alpha training and OLS/MDPDE on selected positive
responses within 1 Å.
FirstStage can sample outside the Joint voxel domain. Raw FirstStage records
retain their sample coordinates and selection independently of grid results.

Joint uses finite positive FirstStage MDPDE widths as seeds and the Target
median to fill missing or invalid seeds when donors exist. FirstStage uses
the native local fit and records its qualification. The
[width calibration](joint-component-runtime.md#initialization-and-global-width-bound)
uses a stricter donor gate than seed fallback. See
[Production Fitting](production-fitting.md) for the shared fitting and
recovery contract.

## Joint result ownership

The snapshot retains the target-only state, Original and Fitting mappings,
available rows, component and Target evidence, conditioning, and candidate
diagnostics. `IsQualifiedJointAnalysisResult` requires valid initialization,
completed search and Passed runtime and Target convergence. Unqualified
snapshots can be saved with diagnostics and reasons. The workflow skips
publication and downstream fitting for them.

The stage adapter exposes qualified Target A/tau/C and convergence in
estimator-neutral Second-stage records. Source identity includes method,
atom, component, run and role. Missing or unqualified points remain
unavailable. Halo records and reconstructed values have separate provenance.
Local and Group results live alongside Joint state and never rewrite it.

## Target Point Quality

`target-point-validation-v2` records component convergence, parameter
validity, Local Conditioning, conditioning block IDs, quality and reason.
Parameter validity checks dimensions, finite A ≥ 0, signed finite C,
finite tau > 0 and agreement of physical and log widths. Blocks must have
nonempty unique IDs and valid unique Target membership. Missing membership
in a nonempty diagnostic set is invalid. An empty set is NotAssessed.

The first matching rule determines quality.

| Evidence | Point quality | Reason |
| --- | --- | --- |
| Component convergence did not pass | Unavailable | `component-not-converged` |
| Invalid or unavailable parameters | Unavailable | `invalid-or-unavailable-parameters` |
| Invalid diagnostic block mapping | Warning | `invalid-diagnostic-block-mapping` |
| Singular, IllConditioned, ActiveFace or Unavailable conditioning | Warning | `conditioning-warning` |
| NotAssessed conditioning | NotAssessed | `quality-diagnostics-not-assessed` |
| Passed numerical checks and WellConditioned evidence | Available | `numerically-converged-point` |

ActiveFace flags possible active-set switching. The
[conditioning calculation](joint-component-numerics.md#compact-spectra-and-local-conditioning)
is diagnostic and does not establish statistical identifiability or
precision. `TargetRuntimeConvergence` combines numerical convergence,
parameter and selection checks, and at assembly also requires complete
fitting-row availability and mapping consistency. A Point Quality warning
can coexist with a qualified endpoint.

## Halo reconstruction

During fitting, Halo support supplies only the exclusion mask. After
qualification, Halo receives a geometry-local median from up to four
nearest eligible Target points by Euclidean distance in Å, with ties
ordered by Atom ID. Donors must match the snapshot and selection, have
Passed component convergence and parameter validity, have `Adequate`
Local Conditioning, and provide finite valid A/tau/C and positions.
Chemical labels play no role. The algorithm takes the median of each
parameter separately and preserves signed charge.

The `geometry-local-target-median-v1` record carries Halo identity, donor
IDs, count, run and donor policy. The global donor median is recorded as
evidence. Reconstruction remains outside the Joint objective. If retained
Halo atoms have no valid donors, the workflow records
`halo-median-no-valid-targets`, publishes the qualified Target endpoints and
skips peeling and subsequent fits. No Halo is a valid case with no required
reconstruction record.

## Original-grid peeling

For each Target, collect Original support rows in increasing Original-row
order at `r² ≤ 4 Å²`, including rows excluded from fitting. Subtract the
direct contributions of every other supported non-hydrogen atom from the
raw voxel response, using qualified Target or reconstructed Halo values
and their actual 2.5 Å support. Each sparse contribution is evaluated once
and shared across Targets. Direct non-self summation preserves the own
signal. Missing contributor state makes that response unavailable.

For Local fitting subtract the Target's own charge,
`y_Gaussian = y_peeled - C_joint * K(r; tau_joint)`. All response and
prediction evaluations here use the Original grid without interpolation.
The fixed-charge profile preserves the Joint C and width independently of
the Gaussian estimates. Each observation records Original row, global
grid index, squared distance, voxel position, raw response and either
a peeled response or an unavailable reason.

## Post-Joint Local and Group Fitting

Post-Joint Local uses selected, converged Joint Targets and current peeling
from the same source. From the 2 Å candidates it retains finite corrected
responses, then strictly positive responses within the inclusive 1 Å
fitting radius. Finite nonpositive corrections are valid observations
excluded from the fit. Counts distinguish candidates, corrected responses,
positive responses and effective solver rows.

One selected log-response Gaussian design serves alpha training, OLS/MDPDE,
working covariance and information. Availability requires more than two
effective rows, rank two, finite nonnegative alpha, qualified successful
MDPDE, valid positive Gaussian A and width, and finite positive working
variance, weights and covariance. The compact information is
`Xᵀ V⁻¹ X` and `Xᵀ V⁻¹ y`, with positive-definite precision for working
covariance `V`. Unavailable or ineligible records carry reasons and
supply no measurement information.

Group Fitting uses only eligible current Local results and their compact
information within a Joint run. Alpha training uses main-chain groups
with at least ten eligible members. Group inference requires at least
two eligible members and a qualified group covariance and posterior
solve. The Group posterior is Gaussian and preserves each member's
separate fixed Joint charge. Local information and Group covariance
are conditional working quantities.

## Invalidation

Changing a Joint stage endpoint invalidates its run's peeling, Local and
Group results and snapshot. Replacing or clearing the snapshot invalidates
dependent grid evidence. Peeling changes invalidate Local and Group
results, and Local changes invalidate affected Groups. Independent FirstStage
raw sample edits leave grid evidence intact. Atom position edits clear
analysis and spatial caches.

Published live results bind to the source map's observation revision.
Replacing, clearing or normalizing that map invalidates the bound analysis
on next access, including model copies. Map copies have independent
revisions. Map geometry is immutable. Destroying the map leaves recorded
evidence available. Loading a saved model restores evidence without a
live map binding. A new workflow collects fresh observations. Runtime
postprocessing contexts borrow their immutable inputs.

## Persistence

[DataObject I/O](architecture/dataobject-io-architecture.md) owns version
acceptance and the database lifecycle. SQLite v21 stores Joint Result
JSON schema 24 and Stage Result JSON version 5. Joint JSON requires
`target-only-pure-halo-support-exclusion-v1`, the current point validation
contract, and `strict-v2` or `practical-v2` convergence criteria. Its state
keys are `ac`, `tau`, `log_tau`, `eta_gradient` and `objective`.
Initialization uses `tau`, `original_tau` and `used_tau`. Candidate
diagnostics use `candidate_eta_gradient_inf_norm` and
`global_eta_gradient_threshold` (the projected threshold). Units use
`joint-kernel-map-units-v2` with tau in Å. Joint CSV orders parameters
`A,tau,C`.

Every saved Joint result, including pre-search failures, requires the
configured solver identity. The reader rejects missing, unknown or mixed
solver fields.

```json
{
  "search_method": "gauss-seidel-bcd",
  "contract_version": "joint-solver-provenance-v6",
  "sparse_backend": "EIGEN",
  "bcd_core_atoms": 12,
  "bcd_policy_version": "gauss-seidel-bcd-production-v1",
  "bcd_maximum_sweeps": 300,
  "bcd_order": "forward",
  "bcd_local_search": "active-set-lm"
}
```

Stage JSON retains FirstStage samples and estimates, Second-stage records,
grid peeling, `post_joint_local` and Group results. Grid peeling requires
its snapshot. Available Halo-subtracted responses require reconstruction
provenance. Validation checks ordered unique `original_row` and `grid`
identities, snapshot row mapping and voxel geometry, response/reason pairs,
source/run identity, and Local failure-row identities and counts. It
replays the correction for available Local records. `fixed_charge` must
match the Joint endpoint and `fitting_radius` must be exactly 1.0 Å for
Available, Unavailable and Ineligible records. Persisted Group members
must match their Local source and fixed charge. Nullable diagnostics,
missing states and unavailable reconstruction are representable outcomes.

See [PotentialFittingWorkflow.cpp](../../src/core/detail/PotentialFittingWorkflow.cpp),
[FirstStageInitialization.cpp](../../src/core/detail/FirstStageInitialization.cpp),
[TargetEvidence.cpp](../../src/core/detail/joint_component/TargetEvidence.cpp),
[HaloMedianReconstruction.cpp](../../src/core/detail/joint_component/HaloMedianReconstruction.cpp),
[PostFitPeeling.cpp](../../src/core/detail/joint_component/PostFitPeeling.cpp),
[GroupPotentialFitting.cpp](../../src/core/detail/GroupPotentialFitting.cpp)
and [PostJointLocalValidation.cpp](../../src/data/detail/PostJointLocalValidation.cpp).
