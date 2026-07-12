# Second-Stage Local Fitting

`rhbm_gem::core::RunSecondStageLocalFitting` refines selected atom-local
Gaussian estimates after the first-stage per-atom fit. It is a
workflow-internal lifecycle seam declared in
[`src/core/detail/GaussianEstimatorStages.hpp`](/src/core/detail/GaussianEstimatorStages.hpp)
and implemented in
[`src/core/GaussianEstimator.cpp`](/src/core/GaussianEstimator.cpp); it is not
part of the installed `GaussianEstimator.hpp` API.

The stage is an outer fixed-point iteration over selected atoms. Each iteration
solves active joint offsets, refits active atoms with selected-neighbor signals
removed, assembles a damped or Anderson-accelerated candidate, and accepts
cluster-local progress through an objective backtracking gate.

## Inputs and State

At entry, the function reads `model_object.GetSelectedAtoms()`, creates one
`AtomLocalPotentialEditor` per selected atom, and builds a
`SecondStageLocalFittingContext` in the same selected-atom order. The context
stores:

- original local sampling entries from `GetSamplingEntries(false)`;
- selected-neighbor indexes within `kNeighborAtomSearchRange`;
- per-sample selected-neighbor contributions within
  `kNeighborContributionDistanceMax`;
- the per-atom `alpha_r` used by local refitting; and
- a prior width for objective penalties. The prior-width preference is a finite
  group prior, then the selected atoms' finite group-median width, then the
  atom's current finite width, and finally `1.0`.

The fixed-point state is a `GaussianFittingState` with aligned result and model
vectors:

```text
result_list      = per-atom LocalGaussianResult
estimation_list  = per-atom MDPDE model vector [amplitude, width, offset]
```

`previous_state` is replaced only after an accepted non-terminal iteration.
Rejected iterations keep the previous state.

## Outer Iteration

Each loop performs the following high-level steps:

```text
previous state
    -> collect active atoms from the freeze tracker
    -> build active clusters and reconcile cluster-local state
    -> solve/refit the active fixed-point map output G(x)
    -> mark suspicious offset rollbacks and temporary ridge multipliers
    -> try localized Anderson damping per cluster
    -> try damped fixed-point fallback for remaining clusters
    -> accept any cluster that passes its local objective gate
    -> update cluster quality state, Anderson history, ridge, freeze/thaw
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
multiplier for that solve.

Robust-loss weights are updated by IRLS. The source-local
`kSecondStageRobustLossKind` policy is currently Cauchy and is shared with the
cluster objective gate. The IRLS loop stops when the weighted-ridge surrogate
objective deteriorates, normalized offset movement is small, or the robust-loss
iteration limit is reached. `JointOffsetSolveStatus` distinguishes converged,
system-build failure, empty-system, initial-solve failure, IRLS-solve failure,
objective deterioration, and IRLS-iteration-limit exits. Build, empty-system,
and initial-solve failures use the previous offsets; later IRLS failures use the
last valid offsets. The rest of the local fitting iteration still runs so its
amplitude/width refits can pass through the normal cluster objective gate.

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
4. accepting the result only if zero-offset sample construction stays finite
   and the profile/parameter plausibility gate does not mark it suspicious.

If refit fails, the previous atom result is reused with the joint offset when
that fallback remains finite and passes the same suspicious-offset gate. If that
fallback is itself suspicious, the atom seeds post-refit rollback. Post-refit
rollback uses the complete active sample-contributor cluster rather than the
bounded joint-offset graph: every atom in the affected cluster restores its
previous validated model and result. This guarantees that no retained refit was
built from a neighbor snapshot that no longer exists. Other active clusters keep
their provisional refits.

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
objective-quality, and ridge managers. Cluster construction
records one representative active contributor per sample while joining all of
that sample's contributors, then assigns objective samples after the connected
components are complete.

The raw refit result is treated as the fixed-point output `G(x)`. For each
cluster with compatible history, localized Anderson acceleration proposes:

```text
candidate = sum(gamma_i * G(x_i)), where sum(gamma_i) = 1
```

Compatibility includes a cluster-local ridge regime signature, not only the
cluster key. The signature records the global joint-offset ridge ratio and,
in canonical cluster-key order, each atom's effective ridge multiplier. The
effective multiplier is the maximum of the cluster objective-ridge multiplier,
the suspicious-retry multiplier, and the proactive collinearity multiplier
actually used to build the joint-offset system. Signatures use exact comparison,
so any real change to one of these controls starts a new fixed-point regime.

Any joint-offset system build that produces the complete effective-multiplier
list can form a signature, including builds whose later solve does not converge.
With the default health policy enabled, only a converged solve may use or commit
that signature: before Anderson candidate construction, an unhealthy solve
clears and suppresses every current cluster. A signature mismatch otherwise
clears and suppresses only the affected cluster. Any pre-refit or post-refit
suspicious rollback likewise clears its containing cluster before candidate
construction. Compatible remote clusters retain their histories.

Residuals are scaled per parameter using the normalized-change scale floor.
The coefficient solve includes L2 regularization and rejects candidates whose
coefficient L1 norm or maximum absolute coefficient exceeds the configured
limits. A rejected coefficient solve is not fatal; that cluster simply has no
Anderson candidate for the attempt and can continue through damped fixed-point
fallback. Candidate construction is structural; the damped candidate that would
actually be applied must still have finite active parameters and positive active
widths.

Each outer iteration tries damping values `1.0`, `0.5`, `0.25`, `0.125`, and
`0.0625`. Anderson attempts run first for clusters that have a localized
candidate. Pending clusters then try the same damping sequence as fixed-point
fallback:

```text
damped = previous + damping * (candidate - previous)
```

Both candidate kinds use the same `BuildLocalFittingCandidateState` path. The
builder copies the previous fitting state, applies the damped values, and
returns an `std::optional<GaussianFittingState>` containing the complete
candidate only when every changed model is finite and has positive width.
Fixed-point candidates retain raw refit uncertainty;
Anderson candidates retain uncertainty from the previous accepted state.
Accepted clusters copy only their active atoms from that candidate into the
assembled state. Rejected or structurally invalid candidates produce no state,
so their clusters keep the previous atom parameters for this iteration.

## Objective Gate and Ridge Retry

Objective scoring is cluster-local. The residual term uses the same robust-loss
policy as joint-offset IRLS, then adds conservative parameter plausibility
penalties for width drift from the group prior or local group median, offset
dominance over the local peak, and single-step movement from the atom/model
snapshot stored with the previous objective candidate. The movement-reference
snapshot is validated against the current active atom indexes before use. The
movement penalty is used only for the current acceptance gate; committed
previous and best objective references keep the residual, width-prior, and
offset-plausibility terms so a one-step movement surcharge does not permanently
pollute future comparisons. Each cluster owns its objective sample refs,
residual-scale tracker, a tracked previous candidate that pairs its objective
stats and samples, best local objective stats, and objective ridge multiplier.
During scale warm-up,
candidate, previous, and best objective values are scored with the same
provisional scale before backtracking is evaluated. Stored objective samples
contain residuals, an atom/model snapshot for the active cluster, and the
scalar residual/response-derived scale sample. Raw response values exist only
while the samples are collected and are discarded after that scale is calculated. This
snapshot is sufficient for provisional rescoring and movement comparison
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

Ridge retry is staged:

- partial accepted iterations lower objective ridge for accepted clusters and
  raise it for rejected clusters;
- if every cluster is rejected, rejected clusters first raise their local
  objective ridge multipliers; and
- the global `ridge_ratio` increases only when all rejected cluster-local
  objective ridge multipliers are saturated.

Accepted iterations without objective rejection shrink the global `ridge_ratio`
toward `kJointOffsetRidgeRatioMin`. An iteration whose joint-offset status is
not converged does not lower its accepted cluster ridge or global ridge. Rejected
clusters retain the existing ridge-increase behavior.

## Freeze, Thaw, and Convergence

All movement and stability decisions use the same dimensionless transformed
coordinates:

```text
log peak height = log(A * (2 pi B^2)^(-3/2))
log width       = log(B)
offset ratio    = C * OffsetBasis(0) / peak height
```

An atom's change is the component-wise absolute difference between its current
and previous transformed coordinates. Scaling both amplitude and offset by the
same map-intensity factor leaves this change unchanged. Using peak height rather
than amplitude also prevents amplitude-width compensation from appearing as
center-signal movement.

Non-positive amplitude or width, and coordinates that cannot be represented as
finite values, produce infinite change. Such states cannot freeze, converge, win
an objective movement tie, or count as persistent no-progress; dependency thaw
treats them conservatively as movement.

Freeze tracking is updated only for atoms in clusters that accepted progress, so
rejected clusters are not frozen because their assembled movement is zero.
If the joint-offset solve is not healthy, all active freeze-stability counters
are reset instead. Accepted parameter movement may still thaw frozen neighbors.

Atoms freeze after their maximum transformed change remains below `1.0e-4` for
`kLocalFittingFreezeStableIterations`. Frozen atoms can thaw when an active
selected neighbor's maximum transformed change reaches `1.0e-3` times the
dependency-thaw hysteresis multiplier. Dependency thaw retains its capped thaw
count and hysteresis decay. Suspicious rollback atoms are force-thawed after
freeze tracking so they can retry with their temporary ridge multiplier.

An accepted cluster enters persistent suspicious rollback tracking only when it
contains the same expanded suspicious atom set as the previous iteration and
its assembled transformed movement is below the convergence tolerance.
Rejected clusters, clusters with effective movement, changed suspicious sets,
and changed cluster keys reset that tracking state. After five consecutive
accepted no-progress iterations, the complete active cluster keeps its previous
validated state and becomes terminal for the remainder of this second stage.
Its atoms no longer participate in offset solving, acceleration, ridge updates,
freeze/thaw, or convergence statistics, but remain in the fitted snapshot as
fixed contributors for other active clusters.

Parameter convergence requires:

- no suspicious offset rollback in the accepted iteration;
- a converged joint-offset solve status;
- no cluster objective rejection in the accepted iteration;
- all active cluster objective scale references to be locked when present; and
- each active-set transformed change percentile below
  `kLocalFittingTransformedChangeTolerance`.

Unhealthy joint-offset iterations clear and suppress current cluster Anderson
history before candidate selection and do not commit their raw fixed-point map
or a ridge regime signature. In non-quiet mode their
existing progress line includes the status, for example
`joint-offset = system-build-failed`; no separate per-iteration warning or
fallback summary is emitted.

Developers can disable only this health policy through the non-installed stage
interface:

```cpp
RunSecondStageLocalFitting(
    model,
    fit_options,
    SecondStageLocalFittingInternalOptions{ false });
```

The default is enabled, including every production workflow call that omits the
internal options. When disabled, solver status is still computed and the current
safe offset fallback remains in force, but status does not reset freeze
stability, block freeze/ridge/history updates or convergence, or appear in
progress output. Ridge-regime compatibility still applies when the system build
produced a signature. A system-build failure has no signature and follows the
legacy history path without health-driven clearing. This switch is intended
only for controlled A/B diagnostics; it is not exposed through `FitOptions`,
commands, CMake, or environment variables.

## Exit Paths

The loop exits through one of five terminal cases:

- **All atoms frozen:** apply the previous state when the loop starts with no
  active atoms, or apply the current accepted assembled state when the last
  accepted update freezes all remaining atoms.
- **Parameter convergence:** apply the current accepted assembled state after
  the convergence conditions pass.
- **Terminal suspicious fallback:** when persistent clusters have been removed
  and all remaining active clusters freeze or converge, apply the assembled
  state and emit one warning with terminal cluster/atom counts instead of
  reporting whole-stage convergence.
- **Objective backtracking failure:** when all acceleration and fixed-point
  attempts are rejected and no further local or global ridge retry is available,
  apply the previous state.
- **Maximum iterations reached:** apply the current accepted assembled state at
  the iteration limit.

Progress logs report iteration, acceleration kind, damping, and active/frozen
atom counts. When terminal suspicious atoms exist, active counts exclude them
and the progress line reports their count separately. Terminal logs report the
stop reason and, for convergence or maximum-iteration exits, the percentile
log-peak-height, log-width, and offset-to-peak-ratio changes.

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
Group Gaussian results are handled by `RunGroupPotentialFitting`.
