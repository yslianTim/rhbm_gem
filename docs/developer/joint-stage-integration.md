# Joint Stage Integration

The [Potential Analysis command](commands/potential_analysis.md) uses the Joint
Component Estimator.
Potential fitting excludes hydrogen from targets and eligible contributors.
The workflow builds the Original and Fitting Domains before running three
execution stages.

1. **InitialStage** samples contributors and fits native OLS/MDPDE estimates,
   including local alpha training and width seeds.
2. **JointStage** runs the [Joint Component Estimator](joint-component-runtime.md)
   and certifies the endpoint. A qualified snapshot feeds contributor reconstruction,
   Original-grid peeling, charge correction and Post-Joint Local MDPDE.
3. **GroupingStage** uses eligible Post-Joint Local evidence for group alpha
   training, group MDPDE, WEB posterior inference and publication. It runs only
   after JointStage completes downstream processing.

The analysis API exposes separate atom-local records for InitialStage and JointStage.

## InitialStage

The initialization workset resolves contributor identities and Output or Buffer roles
in Fitting Domain order. InitialStage uses 20 deterministic shells at 0.05, 0.15, …, 1.95 Å,
with ten points per shell. The internal `FitOptions::initial_stage_sampling_method`
can select `SphereSamplingMethod::SpiralDeterministic`. The default is
`SphereSamplingMethod::FibonacciDeterministic`. Both methods apply the
strict Ownership and 30° cone exclusions for non-hydrogen neighbors.
Tricubic interpolation of stored voxel Map Values produces Sampling Values.
`GetSelectedSamplingEntries()` supplies only selected interpolated samples to
alpha training and OLS/MDPDE on positive responses within 1 Å.
`GetSamplingEntriesView()` exposes all stored samples, including excluded points,
to EstimatorTester and Joint contributor initialization. The latter copies them
through `SetAtomLocalSamplingEntries()` without changing order or selection flags.
The internal `RunInitialStageLocalAlphaTraining` and
`RunInitialStageLocalFitting` support the current selection or an explicit atom list.
InitialStage can sample outside the Joint voxel domain. Raw InitialStage records
retain their sample coordinates and selection independently of grid results.

Use `GetInitialStageEstimate()` for the atom's published estimate record and
`HasInitialStageModel()` before reading its point with `GetInitialStageModel()`.
`GetInitialStageOLS()`, `GetInitialStageMDPDE()` and
`GetInitialStageGaussianResult()` expose native fit evidence. Native MDPDE
can retain a seed or an unqualified model while the published point is absent.

Joint uses finite positive InitialStage MDPDE widths as seeds. Missing or
invalid seeds for either fitting role can use the median from Output donors.
InitialStage uses
the native local fit and records its qualification. The
[width calibration](joint-component-runtime.md#initialization-and-global-width-bound)
uses a stricter donor gate than seed fallback. See
[Production Fitting](production-fitting.md) for the production fitting workflow.

## Joint result ownership

The snapshot retains the complete fitted state for Output and Buffer Targets, Original and Fitting mappings,
available rows, component and Target evidence, conditioning, and candidate
diagnostics. `IsQualifiedJointAnalysisResult` requires valid initialization,
completed search and Passed runtime and Target convergence. Unqualified
snapshots can be saved with diagnostics and reasons. The workflow skips
publication and downstream fitting for them.

The analysis editor publishes qualified A/tau/C and convergence for both fitting roles with
`ApplyJointResult`. Read the estimate and its evidence with
`GetJointStageEstimate()`. Use `HasJointStageModel()` before reading the point
with `GetJointStageModel()`. Source identity includes method, atom, component,
run and role. Missing or unqualified points remain unavailable. Contributor reconstruction
uses separate provenance and preserves both fitting roles. Local and Group results are stored
alongside the JointStage estimate and consume evidence from the same run.

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

## Contributor reconstruction

A Fitting Target can supply its own neighbor profile and donate its Joint point
when the Joint, Target and component convergence and parameter-validity checks
pass, A/tau/C are finite with A ≥ 0 and tau > 0, Local Conditioning is
`Adequate`, and the accepted log width is interior. The [Width KKT boundary
assessment](joint-component-numerics.md#projected-width-boundaries) classifies
each atom independently of its block summary. Boundary and conditioning
reasons can coexist. Point Quality and the accepted endpoint retain their
existing meaning.

A qualified fitting endpoint that fails either reconstruction quality check,
or an External Halo, uses the componentwise median of up to four eligible
original Joint points. Selection uses the contributor's own coordinates,
Euclidean distance in Å, and Atom ID as the deterministic distance tie-breaker.
The search has no radius or chemical-label filter. An immutable KD-tree serves
the entire run, including QScore requests outside the Original Domain.
A replaced fitting point retains its Output or Buffer role.

The `quality-qualified-contributor-reconstruction-v1` snapshot records the
qualified donor count and each replacement's identity, role, reasons, run,
parameters and ordered donor IDs. The `rejection_reasons` bit mask assigns
1 to an External Halo, 2 to the lower width bound, 4 to the upper width
bound, and 8 to conditioning below Adequate. Direct contributors reuse
their Joint points. If a replacement is needed and the donor pool is empty, its record
reports unavailable with no parameters. The workflow preserves the qualified
Joint endpoint, records `no-qualified-reconstruction-donors` and skips
Peeling, Post-Joint Local and Grouping. A run whose contributors are all direct
has no replacement records.

## Original-grid peeling

For each Output Target, collect Original support rows in increasing Original-row
order at `r² ≤ 4 Å²`, including rows excluded from fitting. Subtract the
direct contributions of every other supported non-hydrogen atom from the
raw voxel response using each neighbor’s reconstruction profile
and its recorded 2.5 Å support. Each sparse contribution is evaluated once
and shared across Targets. Direct non-self summation preserves the own
signal. Missing contributor state makes that response unavailable.

For Local fitting subtract the Target's own charge,
`y_Gaussian = y_peeled - C_joint * K(r; tau_joint)`. All response and
prediction evaluations here use the Original grid without interpolation.
The fixed-charge profile preserves the Joint C and width independently of
the Gaussian estimates. Each observation records Original row, global
grid index, squared distance, voxel position, `PostJointGridObservation::map_value`
(the unpeeled stored voxel value) and either
a peeled response or an unavailable reason.

`JointProblemInput::observations` is a neutral solver field populated directly
from `map.GetMapValue(index)`, without interpolation. InitialStage's raw
interpolated samples and these raw voxel observations have different sources.
Peeling, charge correction and Gaussian responses are derived quantities;
they are not Sampling Values. Shared response containers and correlation
helpers retain neutral names when they accept more than one source.

`GetPeeledObservationEntries()` projects only available Original-grid peeled
responses into the existing `LocalPotentialSampleList` for plotting and response
statistics. It preserves observation positions, `sqrt(squared_distance)` and
ordering without interpolation. The projection's `is_selected=true` is only a
container convention. Use `GetPostFitPeeling()` for complete observations,
Original row/global voxel identity, unavailable reasons and provenance.

JointStage uses `AssessLocalFittingPeelingRatio(PeelingRegion::Tail)` for
coverage and ratio availability diagnostics after publishing Peeling.
`umap_embedding` uses the same assessment for Signal `[0, 1]` Å and
Tail `(1, 2]` Å. Classification compares squared distances from the stored
observations. A complete region contributes every binary64 raw/peeled pair
to `sum(|raw - peeled|) / (sum(|raw - peeled|) + sum(|peeled|))`, with
independently scaled component sums. This derived diagnostic is available
when Joint evidence qualifies, coverage is complete, and the sum of component
magnitudes is positive. It uses the pairs before own-charge correction,
independently of Local solver row selection and fitting status. Ratio
calculation leaves the Joint, Local and Group qualification rules in place.

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

For exactly the eligible members `E_g` used by Group MDPDE/WEB,
`C_g = median({post_joint.fixed_charge.coefficient : i in E_g})`.
Each coefficient and its fixed width must match the same-source Joint
endpoint. Charges may be positive, negative or zero. The existing
`descriptive_charge_median` stores this descriptive Joint statistic.
`descriptive_mean` retains the arithmetic means of Local A and tau;
Group `mean`, `mdpde` and WEB `prior` retain their decoded Gaussian A/tau
estimates. All four published models use C_g:

```text
f_g(r) = A_g G(r; tau_g) + C_g K(r; tau_g)
```

`GaussianModel3D::offset` is a charge-kernel coefficient, not a constant
background. C_g does not enter the two-dimensional Group optimization or
covariance. Prior A/tau standard deviations are unchanged and its C standard
deviation is NaN (not inferred). Member posteriors remain Gaussian-only with
C=0, `charge_inferred=false`, a two-dimensional parameter covariance and
NaN C standard deviation. Each member separately retains its own Joint
`fixed_charge = (C_i, tau_i)`; Group tau never overwrites that profile.
One eligible member has a descriptive mean and median but no inference;
zero eligible members have neither. Failed group inference retains its
descriptive evidence without publishing a prior.

Group curves evaluate `C_g K(r; tau_g)`, while Original-grid peeling retains
the individual target's `C_i K(r; tau_i)`. Local fitting subtracts that
individual profile once. Parameter heterogeneity, model approximation and
sample selection can therefore leave differences between a representative
group curve and peeled observations; exact agreement is not required.
Component-average plots retain equal-weight arithmetic parameter averages
across available groups, so their C is the mean of group charge medians,
not the pooled eligible-atom median.

## Post-Joint QScore

Post-Joint QScore uses the same reconstruction records for Original-domain
neighbors and queries the same donor population once per Extra-domain Halo.
The QScore sample strategy and standard reference are retained. Each Output
Target's own Joint reference uses its accepted A/tau/C. Missing neighbors
make that comparison unavailable. Stage JSON records the ordered donor IDs of
Extra-domain Halos under the `joint-peeled-qscore-v2` contract.

## Invalidation

Changing a Joint stage endpoint invalidates its run's Peeling, QScore, Local and
Group results and snapshot. Replacing or clearing the snapshot invalidates
dependent grid evidence. Peeling changes invalidate Local and Group
results, and Local changes invalidate affected Groups. Independent InitialStage
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
acceptance and the database lifecycle. SQLite v27 stores Joint Result
JSON schema 27 and Stage Result JSON version 11. Joint JSON requires
`fixed-output-domain-buffer-expansion-v1`, the current point validation
contract, and `strict-v2` or `practical-v2` convergence criteria. Its state
keys are `ac`, `tau`, `log_tau`, `eta_gradient` and `objective`.
Initialization uses `tau`, `original_tau` and `used_tau`. Its atom records
store `initial_stage_fit_qualified`, `initial_stage_qualification` and
`initial_stage_recovery`. Initialization roles derive from the canonical
`output_indices` and `buffer_indices` through the Fitting-to-Original mapping.
Stage JSON stores the derived roles with the stage records and checks JointStage
points against the fitted snapshot. Candidate
diagnostics use `candidate_eta_gradient_inf_norm` and
`global_eta_gradient_threshold` (the projected threshold). Units use
`joint-kernel-map-units-v2` with tau in Å. Joint CSV keeps accepted fitting parameters in `A,tau,C`
and exports replacement parameters in dedicated reconstruction columns.

Every saved Joint result, including pre-search failures, requires the
configured solver identity. The reader rejects missing, unknown or mixed
solver fields. Results that reach the solver require passing
`global-width-calibration` evidence with policy `global-initial-stage-q50-k4`,
the qualified Output donor count in the reason, the donor median as the
value and `4 × q50` as the threshold. The codec validates the policy and
threshold and checks each saved endpoint and candidate against that bound.

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

Stage JSON stores InitialStage estimates in `initial`, native OLS/MDPDE
diagnostics in `initial_native`, and samples and coordinates in `raw` and
`geometry`. The `joint` record stores JointStage points, convergence and
source identity. The optional source flag is `fixed_from_initial_stage`.
The readers require the current versions and field names.
`peeling` stores the paired observations and neighbor count under `original-grid-v2`.
`post_joint_local` stores charge-corrected fit evidence, and `posterior` stores
Group member inference. Grid peeling requires
its snapshot. Available responses require usable contributor reconstruction.
JSON validates each replacement against accepted endpoint evidence, then SQLite
save and load verify nearest donor order against the model's atom coordinates. Validation checks ordered unique `original_row` and `grid`
identities, snapshot row mapping and voxel geometry, response/reason pairs,
source/run identity, and Local failure-row identities and counts. It
replays the correction for available Local records. `fixed_charge` must
match the Joint endpoint and `fitting_radius` must be exactly 1.0 Å for
Available, Unavailable and Ineligible records. Persisted Group members
must match their Local source and fixed charge. Group validation independently
rebuilds eligible membership from the selected GroupKey members and current
Local evidence, checks counts and run identities, and recomputes C_g. The
stored charge median and all four published model coefficients must equal
that value exactly. Stage JSON uses precise double parsing and SQLite stores
the JSON text, preserving endpoint coefficients and their median on
round-trip. Empty groups cannot acquire a zero fallback or successful
inference; member posterior C remains zero and its charge uncertainty,
like the prior's, remains NaN. Unsupported database/JSON versions are rejected
without migration. Rerun `potential_analysis` into a new database. Nullable diagnostics,
missing states and unavailable reconstruction are representable outcomes.

See [PotentialFittingWorkflow.cpp](../../src/core/estimator/PotentialFittingWorkflow.cpp),
[InitialStage.cpp](../../src/core/estimator/InitialStage.cpp),
[TargetEvidence.cpp](../../src/core/estimator/joint_component/TargetEvidence.cpp),
[ContributorReconstruction.cpp](../../src/core/estimator/joint_component/ContributorReconstruction.cpp),
[PostFitPeeling.cpp](../../src/core/estimator/joint_component/PostFitPeeling.cpp),
[GroupingStage.cpp](../../src/core/estimator/GroupingStage.cpp)
and [PostJointLocalValidation.cpp](../../src/data/detail/PostJointLocalValidation.cpp).
