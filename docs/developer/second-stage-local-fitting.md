# Second-Stage Local Fitting

`rhbm_gem::core::RunSecondStageLocalFitting` refines selected atom-local
Gaussian estimates after the first-stage per-atom fit. It is a
workflow-internal lifecycle seam declared in
[`src/core/detail/GaussianEstimatorStages.hpp`](/src/core/detail/GaussianEstimatorStages.hpp)
and implemented in
[`src/core/GaussianEstimator.cpp`](/src/core/GaussianEstimator.cpp); it is not
part of the installed `GaussianEstimator.hpp` API.

The stage is an outer fixed-point iteration over selected atoms. Each iteration
solves joint offsets independently for every active sample-contributor cluster,
refits active atoms with selected-neighbor signals removed, assembles a damped
or Anderson-accelerated candidate, and accepts cluster-local progress through
an objective backtracking gate.

## Inputs and State

At entry, the function reads `model_object.GetSelectedAtoms()` and builds a
`SecondStageLocalFittingContext` in the same selected-atom order. The context
stores:

- the selected atom pointer;
- original local sampling entries from `GetSamplingEntries(false)`;
- selected-neighbor indexes within `kNeighborAtomSearchRange`;
- per-sample selected-neighbor contributions within
  `kNeighborContributionDistanceMax`;
- the per-atom `alpha_r` used by local refitting; and
- a prior width for objective penalties. The prior-width preference is a finite
  group prior, then the selected atoms' finite group-median width, then the
  atom's current finite width, and finally `1.0`.

The fixed-point state is a `LocalFittingState`, an ordered vector of per-atom
`LocalGaussianResult` values. Its MDPDE model is the single parameter source:

```text
state[i].mdpde.GetModel() = per-atom [amplitude, width, offset]
```

Snapshots and transformed Anderson coordinates are derived from that model only
at their use boundaries, so result and parameter storage cannot diverge. The
stage does not keep a parallel editor list or mutate atom-local results while it
iterates. `ApplyLocalFittingState` obtains an editor for each context atom only
when an exit path commits the selected state.

`previous_state` is replaced only after an accepted non-terminal iteration.
Rejected iterations keep the previous state.

Before constructing `previous_state`, the stage validates every local MDPDE
model with the same rule used by refit and objective evaluation: amplitude and
width must be finite and positive, offset must be finite, and the transformed
coordinates must be representable. An invalid local MDPDE shape is repaired
only in the internal second-stage state. The shape
source order is atom-specific group posterior, group prior, local OLS,
same-group valid-model parameter median, then global valid-model parameter
median. A finite original local offset is retained; otherwise the selected
source offset is used. Direct sources retain their uncertainty, while median
sources use zero uncertainty. Repair changes only `mdpde` in the internal seed;
the other `LocalGaussianResult` fields remain those read at second-stage entry.

The group and global median pools use each atom's valid local MDPDE model when
available; otherwise they use that atom's first valid direct source in the same
posterior, group-prior, then local-OLS order. They include only models satisfying
the validity rule and are separate from the later third-stage group-median
builder. If any selected atom cannot obtain a valid repaired seed, the second
stage emits one warning in non-quiet mode and returns without applying an
internal state.
Non-quiet mode reports repair-source counts; debug logging additionally reports
each repaired atom's original and replacement parameters.

## Outer Iteration

Each loop performs the following high-level steps:

```text
previous state
    -> collect active atoms from the freeze tracker
    -> build active clusters and reconcile cluster-local state
    -> solve each cluster's offsets and refit the active fixed-point map G(x)
    -> mark suspicious offset rollbacks and temporary ridge multipliers
    -> try localized Anderson damping per cluster
    -> try damped fixed-point fallback for remaining clusters
    -> jointly polish A/B/C for each healthy non-suspicious candidate
    -> accept the polished or original candidate through the local objective gate
    -> update cluster quality state, trust radius, Anderson history, ridge, freeze/thaw
    -> exit, retry, or continue
```

The maximum outer iteration count is `kLocalFittingMaximumIterations` (`200`).
`RunSecondStageLocalFitting` owns this outer lifecycle. Anderson-first attempts,
fixed-point fallback, damping, objective evaluation, and accepted-state assembly
are contained in `SelectLocalFittingClusterCandidates`. The same selection
boundary suppresses rejected Anderson histories and releases suppression after
accepted fixed-point progress.

## Joint Offset Step

`RunLocalFittingIteration` first builds a selected-atom Gaussian snapshot from
`previous_state`. Frozen selected atoms remain in the snapshot as fixed signal
contributors, but they do not become active columns in the joint offset system.
Unselected atoms are outside this second-stage snapshot.

The active sample-contributor components built by the outer loop also partition
the joint-offset work. Every cluster solve reads the same immutable previous
snapshot. Its solved offsets and local coupling edges are mapped back to the
complete active-atom order, then all offsets are applied to one refit snapshot
at once. Cluster iteration order therefore cannot expose one cluster to another
cluster's newly solved offset.

For each original sample on each active atom, the joint system:

1. subtracts the target atom's current zero-offset signal;
2. adds the target atom's offset basis as a row entry;
3. subtracts full fitted responses from frozen selected neighbors; and
4. subtracts active selected-neighbor zero-offset signals and adds their offset
   bases as row entries.

The sparse system is solved with weighted ridge regression. The ridge is
relative to the previous offsets, so weakly constrained columns stay close to
their prior offsets. Its effective multiplier is the maximum of three sources:

- a global `ridge_ratio`;
- suspicious-offset retry multipliers on affected atoms; and
- cluster-local objective-backtracking multipliers.

The joint offset builder also applies a per-solve collinearity guard. If two
active offset-basis columns have normalized overlap at least
`kJointOffsetCollinearityOverlapThreshold`, both receive a local ridge
multiplier for that solve. A second full-matrix guard column-normalizes the
sparse design, factorizes its Gram matrix, and applies the same proactive ridge
to the complete cluster when factorization fails, a column is empty, or the
normalized LDLT pivot ratio is at most `1.0e-8`. This catches multi-column
dependencies that do not exceed the pairwise threshold.

Robust-loss weights are updated by IRLS. The source-local
`kSecondStageRobustLossKind` policy is currently Cauchy and is shared with the
cluster objective gate. The IRLS loop stops when the weighted-ridge surrogate
objective deteriorates, normalized offset movement is small, or the robust-loss
iteration limit is reached. `JointOffsetSolveStatus` distinguishes converged,
system-build failure, empty-system, initial-solve failure, IRLS-solve failure,
objective deterioration, and IRLS-iteration-limit exits. Build, empty-system,
initial-solve, and IRLS-solve failures use the previous offsets. Objective
deterioration and the IRLS iteration limit retain the last valid finite offsets.
The rest of the local fitting iteration still runs so its
amplitude/width refits can pass through the normal cluster objective gate.
Each cluster retains its own status and effective ridge multipliers; one failed
solve does not make a disconnected, converged cluster unhealthy.

Joint-offset status has separate progress and stationarity meanings.
`Converged` is eligible for both. `IrlsObjectiveDeteriorated` and
`IrlsMaximumIterationsReached` are finite deterministic outputs eligible for
solver-ridge retry and Anderson progress, but not freeze or convergence.
System-build, empty-system, initial-solve, and IRLS-solve failures are hard
failures eligible for neither.

## Refit and Rollback

After joint offsets are attached to the snapshot, the iteration checks every
active atom for suspicious offsets before refit. The check first preserves the
hard finite-sample guard: if the previous model can build finite zero-offset
samples but the current joint-offset model cannot, the offset is suspicious.
A baseline-independent magnitude guard also rejects a candidate whose center
offset response grows beyond the previous signal, offset, and sampled-profile
scale. This protects the Cauchy solve from finite but physically implausible
offset compensation even when a multi-radius baseline cannot be constructed.
When both profiles are finite, the remaining checks compare the fit-range
zero-offset profile against the previous accepted profile and flag strong center
sign flips, radial rebound or multi-peak shape, sudden width growth, and large
amplitude/offset compensation. These shape checks are conservative and only run
when the previous profile is usable as a baseline.

Suspicious atoms seed bounded rollback propagation through the active coupling
graph produced by the joint offset system. Propagation follows finite edges with
overlap at least `kSuspiciousOffsetClusterMinimumOverlap` and stops after
`kSuspiciousOffsetClusterMaxDepth` graph steps. All seeds in a pass are expanded
with one multi-source traversal, producing the same union of bounded
neighborhoods. This bounded propagation applies only before refit. Every newly
reached active atom is rolled back once to its previous model and skips refit.

After pre-refit rollback, the snapshot is fixed for the complete refit pass.
Each remaining active atom reads that same immutable snapshot and is refit by:

1. subtracting fitted selected-neighbor responses from its original samples;
2. using the joint-offset snapshot model as the fixed-offset model;
3. calling `EstimateLocalGaussian` for amplitude and width; and
4. accepting the result only if its Gaussian model is valid, zero-offset sample
   construction stays finite, and the profile/parameter plausibility gate does
   not mark it suspicious.

If refit fails or returns an invalid Gaussian, the previous validated amplitude
and width are reused with the joint offset when that fallback remains finite and
passes the same suspicious-offset gate. This is an explicit unhealthy refit
fallback. If that fallback is itself suspicious, the atom seeds post-refit
rollback. Post-refit rollback uses the complete active sample-contributor
cluster rather than the bounded joint-offset graph: every atom in the affected
cluster restores its previous validated model and result. This guarantees that
no retained refit was built from a neighbor snapshot that no longer exists.
Other active clusters keep their provisional refits.

Refit health separates progress from stationarity. A finite, plausible result
with `SUCCESS` or `MAX_ITERATIONS_REACHED` is progress-eligible, but only
`SUCCESS` is stationarity-eligible. A result with
`NUMERICAL_FALLBACK`, `INSUFFICIENT_DATA`, or `SINGLE_MEMBER`, a missing fit
result, an exception fallback, or reuse of the previous result is eligible for
neither. Progress-eligible results may pass the objective gate and retain their
cluster history, while stationarity-ineligible results cannot freeze or provide
convergence evidence.

Post-refit affected clusters clear and suppress Anderson history before
candidate construction, so their raw fallback cannot be replaced by an
extrapolation from earlier iterations. They remain subject to the normal
objective gate and persistent suspicious policy; five accepted no-progress
iterations are still required before terminal fallback.

## Clusters and Acceleration

Active clusters are rebuilt every outer iteration from sample-level
contributors. For each selected sample, the active target atom and active
selected neighbors that affect the sample residual are connected. Connected
components define both Anderson history scope and objective sample ownership.
Samples without active contributors do not participate in the cluster objective
gate for that iteration. A canonical cluster work map stores each component's
objective sample references and an optional accepted source (`Anderson` or
`FixedPoint`; unset while pending); the same keys reconcile the acceleration,
objective-quality, trust-radius, and ridge managers. Cluster construction
records one representative active contributor per sample while joining all of
that sample's contributors, then assigns objective samples after the connected
components are complete.

The raw refit result is treated as the fixed-point output `G(x)`. For each
cluster with compatible history, localized Anderson acceleration proposes:

```text
candidate = sum(gamma_i * G(x_i)), where sum(gamma_i) = 1
```

Compatibility includes a cluster-local regime signature, not only the cluster
key. The signature records the joint-offset solve status, the global joint-offset
ridge ratio and, in canonical cluster-key order, each atom's effective ridge
multiplier. The effective multiplier is the maximum of the cluster objective-
ridge multiplier, the suspicious-retry multiplier, and the proactive
collinearity multiplier actually used to build the joint-offset system.
Signatures use exact comparison, so a status or control change starts a new
fixed-point regime.

Each joint-offset system build that produces its cluster's complete
effective-multiplier list can form a signature, including a build whose later
solve does not converge. The signature map is therefore partial when another
cluster fails during system construction. With the default health policy
enabled, a hard-failure or local-refit-fallback cluster clears and suppresses
only its own Anderson history before candidate construction and cannot commit
its raw map or signature. Progress-eligible joint-offset results may commit
history only while their complete signature remains unchanged. A mismatch
clears only the affected cluster. Any
pre-refit or post-refit suspicious rollback clears its containing cluster before
candidate construction. Compatible healthy remote clusters retain their
histories.

Anderson inputs and fixed-point outputs are encoded as log peak height, log
width, and offset-to-peak ratio. The coefficient solve includes L2
regularization and rejects candidates whose
coefficient L1 norm or maximum absolute coefficient exceeds the configured
limits. A rejected coefficient solve is not fatal; that cluster simply has no
Anderson candidate for the attempt and can continue through damped fixed-point
fallback. Candidate construction is structural; transformed candidates must be
finite and decode to finite models with positive amplitude and width.

Each outer iteration tries damping values `1.0`, `0.5`, `0.25`, `0.125`, and
`0.0625`. Anderson attempts run first for clusters that have a localized
candidate. Pending clusters then try the same damping sequence as fixed-point
fallback. Before candidate construction, the requested damping is capped by the
cluster's adaptive trust radius:

```text
step norm = max over atoms and parameters of
            abs(candidate - previous) / [0.50, 0.35, 1.00]
effective damping = min(requested damping, trust radius / step norm)
damped = previous + effective damping * (candidate - previous)
```

The scales correspond to log peak height, log width, and offset-to-peak ratio.
Identical effective damping values caused by the cap are evaluated only once.
An Anderson cluster clears and suppresses its history only when every damping
attempt fails. If a larger damping fails but a smaller damping is accepted, the
history remains active and receives the normal committed fixed-point sample.

Both candidate kinds use the same `BuildLocalFittingCandidateState` path. The
builder performs the displayed interpolation in transformed coordinates,
decodes the result to the canonical raw `[amplitude, width, offset]` state, and
returns an `std::optional<LocalFittingState>` containing the complete
candidate only when every changed model is valid. Exact transformed no-op
candidates retain the previous raw model without a round-trip conversion.
Fixed-point candidates retain raw refit MDPDE uncertainty; Anderson candidates
retain MDPDE uncertainty from the previous accepted state. Candidate building
starts from `previous_state` and replaces only `mdpde`, so other
`LocalGaussianResult` fields are retained from the previous accepted state. In
particular, the raw refit's `fit_result` is used to classify progress and
stationarity health but is not committed into the candidate. Accepted clusters
copy only their active atoms from that candidate into the assembled state.
Rejected or structurally invalid candidates produce no state, so their clusters
keep the previous atom parameters for this iteration.

### Joint A/B/C Polishing

Every structurally valid Anderson or fixed-point candidate in a cluster with a
converged joint-offset solve, stationarity-eligible refits, and no suspicious
atom is polished before it can contribute freeze or convergence evidence. The
polisher uses the same dimensionless transformed coordinates as acceleration:

```text
q = [log peak height, log width, offset-to-peak ratio]
```

It collects the cluster's existing objective samples and evaluates one joint
residual/Jacobian system with three columns per active atom. Target and active-
neighbor derivatives occupy their cluster columns; frozen and terminal selected
atoms remain fixed contributors in the candidate snapshot. The analytic
Jacobian follows the same center-distance branch as `GaussianModel3D`, including
the limiting offset-basis response at distances below `1.0e-5`.

One Cauchy weight vector is calculated from the candidate residual and its MAD
scale. A single weighted-ridge Gauss-Newton direction is then solved around zero
movement. The ridge uses the current global ratio, the atom's effective
cluster/suspicious multiplier, column squared norms, and the full-matrix
conditioning guard. This is a bounded polishing step, not an inner IRLS loop.
Directions below `kLocalFittingTransformedChangeTolerance` are stationary and
the original candidate proceeds directly to the objective gate.

For a non-stationary direction, the normal damping sequence is applied from the
original candidate toward the joint solution. The polish substep is clipped so
the total displacement from the outer previous state remains inside the same
cluster trust radius. A candidate already on the boundary cannot polish farther
outward, but an inward polish remains eligible. Polished variants are tried first
through the existing cluster objective gate. Rejected polished variants do not
commit scale or quality state and do not independently trigger solver-ridge
backtracking. If none is accepted, the original candidate receives its normal
single objective attempt. An accepted original candidate is retained as a
polishing fallback; an accepted polished or stationary candidate supplies valid
stationarity evidence. Thus every cluster still commits at most one quality
candidate per outer iteration, and polishing does not become a third candidate
kind or change Anderson-before-fixed-point ordering. Polished candidates retain
the original candidate's MDPDE uncertainty and other local-result fields.

## Scientific Objective, Trust Radius, and Ridge Retry

Objective scoring is cluster-local. The residual term uses the same robust-loss
policy as joint-offset IRLS, then adds conservative parameter plausibility
penalties for width drift from the group prior or local group median, offset
dominance over the local peak. Candidate, previous, best, committed, and global
audit states all use this same scientific objective; movement is not an
objective term. Each cluster owns its objective sample refs,
residual-scale tracker, a tracked previous candidate that pairs its objective
stats and samples, best local objective stats, solver-ridge multiplier, and
trust radius.
During scale warm-up,
candidate, previous, and best objective values are scored with the same
provisional scale before backtracking is evaluated. Stored objective samples
contain residuals, an atom/model snapshot for the active cluster, and the
scalar residual/response-derived scale sample. Raw response values exist only
while the samples are collected and are discarded after that scale is calculated. This
snapshot is sufficient for provisional rescoring
without retaining a response list or reconstructing models from a separate
fitting state.

Unavailable scale references and unavailable objective values are represented
as optional values rather than separate presence flags plus infinity sentinels.
The damping loop owns candidate retry. The quality manager returns whether a
candidate was accepted and immediately commits that cluster's scale, previous,
best, and objective-sample state on acceptance; rejected candidates leave the
manager state unchanged.

A cluster candidate is rejected when a local reference exists and the candidate
has no finite objective, or when its objective deteriorates beyond
`kLocalFittingConvergenceObjectiveRelativeTolerance` relative to the cluster's
previous or best tracked state. Objective ties are broken by the maximum
scale-consistent transformed percentile movement.

Rejected Anderson clusters clear and suppress only their own histories. A
suppressed cluster can use Anderson again after accepted fixed-point progress.
Healthy accepted clusters commit the raw fixed-point map and the signature used
to generate it together. Suspicious and terminal clusters do not commit either.
Ridge changes made after candidate selection do not rewrite that signature; the
next iteration detects the new regime, clears the old samples, and requires a
fresh fixed-point commit before Anderson can resume.

The trust radius starts at `1.0`, is bounded to `[0.0625, 4.0]`, shrinks by
`0.5`, and grows by `2.0`. Exact cluster keys retain their radius across outer
iterations; a split or merge starts at the initial radius. A cluster that
remains rejected shrinks its radius. A candidate that improves the scientific
objective beyond the tie tolerance and uses at least 80% of the radius grows it.
Unavailable or tied objectives leave it unchanged.

Trust-radius and solver-ridge retry are staged:

- rejected clusters first shrink their trust radius;
- only clusters already at the minimum radius raise their local solver-ridge
  multiplier, then reset their radius because the generated direction changed;
- partial accepted iterations lower solver ridge for progress-eligible,
  accepted, non-suspicious clusters; and
- the global `ridge_ratio` increases only when all rejected cluster-local
  solver-ridge multipliers are saturated, then resets the affected radii.

Accepted iterations without objective rejection shrink the global `ridge_ratio`
toward `kJointOffsetRidgeRatioMin` when every active cluster is progress-
eligible. Soft incomplete IRLS results can therefore continue the numerical
trajectory without being treated as stationarity evidence. Hard-failure or
local-refit-fallback clusters do not lower their local solver ridge;
eligible remote clusters retain their local ridge decrease.

## Freeze, Thaw, and Convergence

All movement and stability decisions use the same dimensionless transformed
coordinates:

```text
log peak height = log(A * (2 pi B^2)^(-3/2))
log width       = log(B)
offset ratio    = C * OffsetBasis(0) / peak height
```

An atom's accepted change is the component-wise absolute difference between its
assembled and previous transformed coordinates. Its raw fixed-point residual is
the same difference between the undamped refit output `G(previous)` and the
previous state. Scaling amplitude and offset by the same map-intensity factor
leaves both quantities unchanged. Using peak height rather than amplitude also
prevents amplitude-width compensation from appearing as center-signal movement.

Non-positive amplitude or width, and coordinates that cannot be represented as
finite values, produce infinite change or residual. Invalid accepted states
cannot freeze, converge, win an objective movement tie, or count as persistent
no-progress; dependency thaw treats them conservatively as movement. An invalid
raw state independently prevents freeze and parameter convergence.

Freeze tracking is updated only for atoms in clusters that accepted progress, so
rejected clusters are not frozen because their assembled movement is zero and
now reset their accumulated freeze stability. For eligible atoms, freeze
evidence is the component-wise maximum of accepted change and raw fixed-point
residual, so damping, Anderson, or an objective gate cannot make a large raw map
look stationary.
An unhealthy cluster resets only its own active freeze-stability counters.
An accepted polishing fallback likewise resets only its cluster's stability and
is omitted from that iteration's freeze update. This preserves its objective-
accepted parameters without allowing an unpolished candidate to freeze.
Accepted parameter movement from either healthy or unhealthy clusters may still
thaw frozen neighbors.

Atoms freeze after the maximum component of their combined accepted/raw
freeze evidence remains below `1.0e-4` for
`kLocalFittingFreezeStableIterations`. Frozen atoms can thaw when an active
selected neighbor's maximum accepted transformed change reaches `1.0e-3` times
the dependency-thaw hysteresis multiplier. Dependency thaw can repeat without a
per-atom count limit or permanent lock; threshold growth, its maximum multiplier,
and frozen-state decay still suppress rapid oscillation. Suspicious rollback
atoms are force-thawed after freeze tracking so they can retry with their
temporary ridge multiplier.

A single persistent terminal-failure tracker handles accepted no-progress
clusters. A suspicious rollback has priority and uses the complete expanded
suspicious atom set as its failure reason. An accepted cluster without a
suspicious rollback instead uses its hard-failure `JointOffsetSolveStatus` as
the reason. The tracker advances only while the cluster key and exact failure
reason remain unchanged and its assembled transformed movement stays below the
convergence tolerance. Rejection, effective movement, solver recovery, a reason
or status change, or a cluster-key change resets the count.

After five consecutive accepted no-progress iterations with the same reason,
the complete active cluster restores the current iteration's `previous_state`
and becomes terminal for the remainder of the second stage. Its atoms no longer
participate in offset solving, acceleration, ridge updates, freeze/thaw, or
convergence statistics, but remain in the fitted snapshot as fixed contributors
while disconnected healthy clusters continue fitting. Terminal diagnostics
retain separate suspicious-rollback and joint-offset-failure summaries,
including the hard-failure status breakdown.

Parameter convergence requires:

- no suspicious offset rollback in the accepted iteration;
- every active cluster to have a converged joint-offset solve and a
  stationarity-eligible local refit;
- no accepted cluster to have used the unpolished fallback;
- no cluster objective rejection in the accepted iteration;
- all active cluster objective scale references to be locked when present; and
- each active-set accepted transformed change and raw fixed-point residual to
  have p99 below `1.0e-4` and maximum below `1.0e-3` for every transformed
  component.

Raw residual is used only for freeze and parameter convergence. Dependency
thaw, objective tie-breaking, persistent suspicious or solver no-progress, and
ridge adjustment retain their accepted-change semantics.

Hard-failure and local-refit-fallback clusters clear and suppress their own
Anderson history before candidate selection. Soft incomplete status changes
also start a new history regime. In non-quiet mode the existing progress or
retry line reports every joint status as
`joint-offset statuses clusters/atoms = status:C/A`, reports strict health as
`health-unhealthy clusters/atoms = C/A`, counts non-stationary joint reasons in
enum order such as `joint-offset = system-build-failed:1`, and reports refit fallbacks as
`local-refit-fallback clusters/atoms = C/A`. Stationarity-ineligible refits are
reported as `local-refit-nonstationary clusters/atoms = C/A`. Candidate summaries
report `joint-ABC polish clusters accepted/stationary/fallback = C/C/C`. No separate
per-iteration warning or fallback summary is emitted.

This end-to-end health lifecycle is always applied by
`RunSecondStageLocalFitting`; it is not separately configurable through the
internal stage interface, `FitOptions`, commands, CMake, or environment
variables.

## Maximum-Iteration Audit State

Maximum-iteration fallback uses a global audit that is independent of
cluster-local objective history. Its sample refs contain every original
second-stage sample in selected-atom order, and its model set contains every
selected atom regardless of active, frozen, or terminal status. The first state
with a finite complete audit locks one robust scale for the remainder of the
stage; the initial state is considered first.

The fixed audit objective uses the normal Cauchy residual, width prior, and
offset plausibility penalty. It excludes trust-radius, ridge, and freeze state.
Its residual is a per-sample mean, while its width
and offset penalties are per-atom means. This keeps the audit ranking invariant
when the same scientific case is replicated to a different atom count.

Cluster-local objectives use the same cardinality-independent aggregation:
width-prior and offset-plausibility penalties are each averaged over the active
atoms in that cluster before their weights are applied. Candidate, committed
previous, best, and audit states therefore have identical objective semantics.

Every complete assembled state with at least one objective-gated accepted
cluster can compete, including partial or unhealthy iterations. The tracker
keeps both the global best, which may be the initial state, and the best state
drawn only from accepted iterations. A candidate replaces either tracked state
only when its finite objective improves beyond the objective-tie tolerance;
ties retain the earlier state.

When atoms become terminal, their validated fallback models and results are
overlaid onto the stored audit best and the fixed objective is recomputed. This
keeps historical improvements in remote atoms without allowing an
iteration-limit exit to restore superseded terminal parameters. If the
reconciled state cannot be audited, it is discarded and the current assembled
state can establish a new best under the already fixed scale.

## Exit Paths

The loop exits through one of five terminal cases:

- **All atoms frozen:** apply the previous state when the loop starts with no
  active atoms, or apply the current accepted assembled state when the last
  accepted update freezes all remaining atoms.
- **Parameter convergence:** apply the current accepted assembled state after
  the convergence conditions pass.
- **Terminal fallback:** when persistent suspicious or joint-offset-failure
  clusters have been removed and all remaining active clusters freeze or
  converge, apply the assembled state and emit one warning with reason-specific
  terminal cluster/atom counts instead of reporting whole-stage convergence.
- **Objective backtracking failure:** when all acceleration and fixed-point
  attempts are rejected at the minimum trust radius and maximum global solver
  ridge, apply the best validated
  accepted-iteration audit state when available, otherwise apply the previous
  state. If the final rejection instead reaches the outer iteration limit,
  apply the global best validated audit state when available.
- **Maximum iterations reached:** after an accepted final iteration, apply the
  best validated audit state when available; otherwise preserve the current
  accepted assembled-state fallback.

Progress logs report iteration, accepted Anderson/fixed-point cluster counts,
joint-A/B/C accepted/stationary/fallback cluster counts,
active/frozen atom counts, and raw/accepted 99th-percentile
offset-to-peak-ratio change. Objective retry lines report the raw offset
statistics because no candidate was accepted. Accepted and retry lines also
report objective-gate accepted/rejected cluster and atom counts before health,
suspicious, or terminal filtering, plus every joint-offset status as cluster/
atom counts. Active counts exclude all terminal atoms;
progress reports
`terminal-suspicious atoms` and `terminal-joint-offset-failure atoms`
separately. Final warnings report reason-specific cluster/atom counts and the
terminal joint-offset status breakdown. Convergence logs retain the percentile
log-peak-height, log-width, and offset-to-peak-ratio changes. Maximum-iteration
warnings instead report whether the best validated audit state or the applicable
previous/current accepted-state fallback was applied, the audit source and
objective breakdown when available, and the offset distribution of the state
actually applied.

At debug verbosity (`-v4`), every cluster that remains rejected after all
candidate attempts emits an objective diagnostic after the normal progress or
retry line. The cluster is identified by atom count and the first/last atom
index in its canonical key. Attempts are reported in execution order, including
Anderson candidates when available and distinct effective fixed-point damping
values. Requested damping `1.0` is labeled `raw`.

Each scored attempt reports the provisional objective scale and candidate,
previous, and best breakdowns in
`residual/width/offset/total` order. It also reports requested and effective
damping, trust radius, normalized step, and polish damping when applicable.
`rejected-by` identifies whether the existing tolerance comparison failed
against previous, best, or both. Invalid damped models and unavailable
objectives are labeled without fabricated component values. An
invalid-model entry also reports the first failing atom, whether the failure is
parameter-size, non-finite-parameter, invalid-transformed-coordinates, or
non-positive-width, and the candidate `u/v/q` values when their structure is
valid. These
diagnostics observe the existing gate only; they do not change objective
weights, scale tracking, trust radius, damping, or acceptance.

## Workflow Context

The installed estimator API exposes `RunPotentialFittingWorkflow`, while the
local-alpha, fixed-offset first/third pass, and second-stage entry points remain
private lifecycle seams in `src/core/detail/GaussianEstimatorStages.hpp`.
Keeping these stage hooks private prevents callers from bypassing their required
workflow ordering.

In `RunPotentialFittingWorkflow`, the second stage runs after initial local-alpha
training, seed-model initialization, first-stage local fitting, and an initial
group-alpha/group-potential fit. After second-stage local fitting, the workflow:

1. retrains group alpha;
2. rebuilds updated samples from group-median Gaussians and refits group
   potentials;
3. rebuilds updated samples from the fitted group Gaussians;
4. retrains local alpha on those updated samples and runs third-stage local
   fitting; and
5. retrains group alpha and group potentials once more.

Second-stage local fitting updates selected atom-local Gaussian results only.
Group Gaussian results are handled by `RunGroupPotentialFitting`. Info-level
stage-boundary summaries report selected-atom offsets after second stage, group
prior offsets before third stage, and selected-atom offsets after third stage.
Third-stage fitting intentionally keeps the existing group-prior median offset
as its fixed offset; the summaries make that propagation visible without
changing it.
