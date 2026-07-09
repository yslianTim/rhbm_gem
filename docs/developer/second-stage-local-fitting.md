# Second-Stage Local Fitting

`rhbm_gem::core::RunSecondStageLocalFitting` refines selected atom-local
Gaussian estimates after the first-stage per-atom fit. The implementation is in
[`src/core/GaussianEstimator.cpp`](/src/core/GaussianEstimator.cpp).

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
- selected-neighbor indexes within `kNeighborAtomSearchRange`; and
- per-sample selected-neighbor contributions within
  `kNeighborContributionDistanceMax`.

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

Robust-loss weights are updated by IRLS. The current internal policy defaults
to Huber, with Cauchy available through the same source-local loss setting. The
IRLS loop stops when the weighted-ridge surrogate objective deteriorates,
normalized offset movement is small, or the robust-loss iteration limit is
reached. If the system cannot be built, is empty, or cannot be solved, the
offset step uses the previous offsets and the rest of the local fitting
iteration still runs.

## Refit and Rollback

After joint offsets are attached to the snapshot, the iteration checks every
active atom for suspicious offsets before refit. The check first preserves the
hard finite-sample guard: if the previous model can build finite zero-offset
samples but the current joint-offset model cannot, the offset is suspicious.
When both profiles are finite, the same check also rejects clearly non-physical
local profiles. It compares the fit-range zero-offset profile against the
previous accepted profile and flags strong center sign flips, radial rebound or
multi-peak shape, sudden width growth, and large amplitude/offset compensation.
These profile checks are conservative and only run when the previous profile is
usable as a baseline.

Suspicious atoms seed bounded rollback propagation through the active coupling
graph produced by the joint offset system. Propagation follows finite edges with
overlap at least `kSuspiciousOffsetClusterMinimumOverlap` and stops after
`kSuspiciousOffsetClusterMaxDepth` graph steps. Every reached active atom is
rolled back to its previous model for this iteration and skips refit.

Each remaining active atom is refit by:

1. subtracting fitted selected-neighbor responses from its original samples;
2. using the joint-offset snapshot model as the fixed-offset model;
3. calling `EstimateLocalGaussian` for amplitude and width; and
4. accepting the result only if zero-offset sample construction stays finite
   and the profile/parameter plausibility gate does not mark it suspicious.

If refit fails, the previous atom result is reused with the joint offset when
that fallback remains finite and passes the same suspicious-offset gate. If that
fallback is itself suspicious, the atom is kept at its previous state and can
seed the same bounded rollback pass after the refit loop.

## Clusters and Acceleration

Active clusters are rebuilt every outer iteration from sample-level
contributors. For each selected sample, the active target atom and active
selected neighbors that affect the sample residual are connected. Connected
components define both Anderson history scope and objective sample ownership.
Samples without active contributors do not participate in the cluster objective
gate for that iteration.

The raw refit result is treated as the fixed-point output `G(x)`. For each
cluster with compatible history, localized Anderson acceleration proposes:

```text
candidate = sum(gamma_i * G(x_i)), where sum(gamma_i) = 1
```

Residuals are scaled per parameter using the normalized-change scale floor.
Candidate construction is structural; the damped candidate that would actually
be applied must still have finite active parameters and positive active widths.

Each outer iteration tries damping values `1.0`, `0.5`, and `0.25`.
Anderson attempts run first for clusters that have a localized candidate.
Pending clusters then try the same damping sequence as fixed-point fallback:

```text
damped = previous + damping * (candidate - previous)
```

Accepted clusters copy only their active atoms into the assembled state.
Rejected clusters keep their previous atom parameters for this iteration.

## Objective Gate and Ridge Retry

Objective scoring is cluster-local and uses the same robust-loss policy as
joint-offset IRLS. Each cluster owns its objective sample refs, residual-scale
tracker, previous objective samples, best local objective stats, and objective
ridge multiplier. During scale warm-up, candidate, previous, and best objective
values are scored with the same provisional scale before backtracking is
evaluated.

A cluster candidate is rejected when a local reference exists and the candidate
has no finite objective, or when its objective deteriorates beyond
`kLocalFittingConvergenceObjectiveRelativeTolerance` relative to the cluster's
previous or best tracked state. Objective ties are broken by the maximum
normalized percentile parameter movement.

Rejected Anderson clusters clear and suppress only their own histories. A
suppressed cluster can use Anderson again after accepted fixed-point progress.
Accepted clusters commit quality state and Anderson history when not suppressed.

Ridge retry is staged:

- partial accepted iterations lower objective ridge for accepted clusters and
  raise it for rejected clusters;
- if every cluster is rejected, rejected clusters first raise their local
  objective ridge multipliers; and
- the global `ridge_ratio` increases only when all rejected cluster-local
  objective ridge multipliers are saturated.

Accepted iterations without objective rejection shrink the global `ridge_ratio`
toward `kJointOffsetRidgeRatioMin`.

## Freeze, Thaw, and Convergence

Accepted assembled candidates compute absolute parameter changes for freeze/thaw
and normalized parameter changes for objective tie-breaking and convergence.
Freeze tracking is updated only for atoms in clusters that accepted progress, so
rejected clusters are not frozen because their assembled movement is zero.

Atoms freeze after their maximum absolute parameter movement remains below the
freeze threshold for `kLocalFittingFreezeStableIterations`. Frozen atoms can
thaw when an active selected neighbor changes enough to pass the dependency
thaw threshold. Dependency thaw uses per-atom hysteresis and a capped thaw count
to limit repeated freeze/thaw cycling. Suspicious rollback atoms are force-thawed
after freeze tracking so they can retry with their temporary ridge multiplier.

Parameter convergence requires:

- no suspicious offset rollback in the accepted iteration;
- no cluster objective rejection in the accepted iteration;
- all active cluster objective scale references to be locked when present; and
- active-set normalized percentile changes below
  `kLocalFittingNormalizedChangeTolerance`.

## Exit Paths

The loop exits through one of four terminal cases:

- **All atoms frozen:** apply the previous state when the loop starts with no
  active atoms, or apply the current accepted assembled state when the last
  accepted update freezes all remaining atoms.
- **Parameter convergence:** apply the current accepted assembled state after
  the convergence conditions pass.
- **Objective backtracking failure:** when all acceleration and fixed-point
  attempts are rejected and no further local or global ridge retry is available,
  apply the previous state.
- **Maximum iterations reached:** apply the current accepted assembled state at
  the iteration limit.

Progress logs report iteration, acceleration kind, damping, and active/frozen
atom counts. Terminal logs report the stop reason and, for convergence or
maximum-iteration exits, normalized terminal movement.

## Workflow Context

In `RunPotentialFittingWorkflow`, this stage runs after first-stage local
fitting initializes selected atom-local Gaussian results. Later workflow stages
train group alpha values, fit group potentials, rebuild updated sampling
entries, retrain local alpha on updated samples, and run third-stage local
fitting.

Second-stage local fitting updates selected atom-local Gaussian results only.
Group Gaussian results are handled by `RunGroupPotentialFitting`.
