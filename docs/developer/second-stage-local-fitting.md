# Second-Stage Local Fitting

`rhbm_gem::core::RunSecondStageLocalFitting` refines selected atom-local
Gaussian estimates after the first-stage per-atom fit. The implementation is in
[`src/core/GaussianEstimator.cpp`](/src/core/GaussianEstimator.cpp).

This stage is a fixed-point iteration across the selected atoms. Each iteration
estimates a joint offset update for the active selected atoms, refits those atoms
with selected-neighbor contributions removed, then uses percentile normalized
parameter movement to decide whether the active set has converged.

## Inputs and State

The function receives:

- a `ModelObject`, used to read the selected atom list and edit each atom's
  local analysis entry; and
- `FitOptions`, which supplies fit distance range, thread
  count for the lower-level Gaussian estimator, and logging mode.

At entry, it reads `model_object.GetSelectedAtoms()`, builds one
`AtomLocalPotentialEditor` per selected atom, and builds a
`SecondStageLocalFittingContext` from the same selected-atom order. The context
stores original local sampling entries via `GetSamplingEntries(false)`, so this
stage does not apply the sample selection filter and does not use updated
sampling entries produced later in the workflow. The previous iteration state is
stored as a `GaussianFittingState` with two aligned vectors:

```text
previous_result_list      = current per-atom LocalGaussianResult
previous_estimation_list  = current MDPDE model vector [amplitude, width, offset]
```

These vectors are the fixed-point state. They are replaced at the end of each
non-terminal iteration.

## Per-Iteration Flow

Each loop first asks the freeze tracker for the active atom indexes. If no atom
is active, the previous state is applied and the stage exits. Otherwise the
active indexes are passed to `RunLocalFittingIteration`, and the updated model
vectors are checked for freezing and convergence.

```text
previous state
    -> build active atom set
    -> exit if active set is empty
    -> build selected-atom fitted Gaussian snapshot
    -> estimate joint offsets for active atoms
    -> refit active atoms with selected-neighbor contributions removed
    -> roll back any suspicious offset clusters found before or during refit
    -> build an Anderson-accelerated or damped fixed-point candidate
    -> compute active-atom p95 absolute and normalized parameter changes
    -> backtrack against the objective; retry with lower damping when needed
    -> update candidate ranking and Anderson history
    -> freeze stable atoms and thaw changed selected neighbors
    -> exit, fallback, or continue
```

The maximum iteration count is `kLocalFittingMaximumIterations` (`200`).

## Joint Offset Step

`RunLocalFittingIteration` first converts the full previous selected-atom
vectors into a snapshot aligned with the second-stage context atom index.
`EstimateJointOffsets` then solves one sparse linear system for active atom
offsets. Frozen selected atoms remain in the snapshot, but they do not become
columns in the linear system. Unselected atoms are not present in this snapshot,
so second-stage neighbor subtraction does not use them as fixed contributors.

For each original sampling entry on each active atom, with the sample selection
filter disabled:

1. subtract the target atom's current zero-offset signal;
2. add the target atom's offset basis as the row entry for that atom;
3. for selected neighbors present in the snapshot and within
   `kNeighborContributionDistanceMax`, subtract the full fitted response when
   the neighbor is not in the active solve; and
4. for active selected neighbors, subtract the current zero-offset signal and add
   the neighbor's offset basis as another row entry.

The resulting system is solved with weighted ridge regression. The ridge term is
relative to the previous offsets, so weakly constrained columns stay close to
their prior values. Its global ratio starts at `kJointOffsetRidgeRatio`
(`1.0e-3`) and is adjusted across outer fixed-point iterations:
objective-backtracking rejections that reject every active cluster increase the
ratio only after the rejected clusters have exhausted their cluster-local
objective ridge multipliers. Accepted iterations without cluster rejections
gradually decrease the global ratio. Cluster-local objective rejections use
per-atom ridge multipliers first, so a rejected cluster can be retried without
raising the baseline ridge for unrelated clusters.

The joint offset builder also applies a proactive local ridge guard before the
solve. While assembling the sparse design matrix, it accumulates active-column
cross products and converts each pair to a normalized overlap. If two active
offset-basis columns overlap by at least `0.98`, both atoms receive a local
ridge multiplier for that solve only. This catches near-collinear neighboring
atoms before they can produce an extreme offset update; it is independent from
the global objective-backtracking ridge controller and is recomputed each time
the active set changes.

Atoms that previously participated in a suspicious joint offset rollback can
also receive a temporary per-atom ridge multiplier, which keeps their next joint
offset solve closer to the previous offset without changing public fitting
options. Suspicious-offset and objective-backtracking multipliers are tracked
separately and combined with `max()` when the joint system is built, so accepting
one source does not accidentally clear the other. Huber weights are then updated
from residual median absolute deviation. The IRLS loop stops when the weighted-ridge surrogate objective would
deteriorate, when the maximum normalized offset movement drops below
`kJointOffsetIrlsNormalizedChangeTolerance`, or when the Huber iteration limit is
reached.

If the joint system cannot be built, is empty, or cannot be solved during the
initial or robust solve, the step falls back to the previous offsets. The rest of
the local fitting iteration still runs, and the fallback is counted for the final
diagnostic summary.

## Per-Atom Refit

After joint offsets are attached to the snapshot, each active atom is refit.
Frozen atoms are left in the iteration state copied from the previous state. The
second-stage active-atom loop itself is sequential in `RunLocalFittingIteration`;
`FitOptions::thread_size` is passed through to the lower-level fixed-offset
Gaussian estimator.

Before the per-atom refit loop runs, each active atom's joint-offset snapshot
model is checked against the atom's raw sampling entries. If the previous model
can build finite zero-offset samples but the joint-offset model cannot, the atom
seeds suspicious offset rollback propagation. The propagation graph comes from
active columns that co-occur in the sparse joint-offset system, with each edge
weighted by normalized joint-offset column overlap. Rollback expands only across
finite edges whose overlap is at least `0.05`, and only within two topological
steps from the original suspicious atom. This is narrower than all spatial
neighbors and avoids rolling back a large connected component through distant
weak links. Frozen selected neighbors remain fixed snapshot contributors;
unselected neighbors are outside the second-stage snapshot.
Every reached atom has its snapshot entry rolled back to the previous model and
skips refit for that iteration, so strongly coupled nearby active atoms see a
synchronous rollback rather than a one-sided update.

For each atom:

1. build a sample list by subtracting fitted selected-neighbor responses from
   the atom's unfiltered local sampling entries;
2. use the joint-offset snapshot model as the fixed offset model;
3. call `EstimateLocalGaussian` to fit amplitude and width with that fixed
   offset model and `FitOptions` distance limits; and
4. accept the candidate only if `CanBuildFiniteZeroOffsetSamples` confirms that
   subtracting the candidate offset remains finite and inside the `float`
   response range.

If fitting throws or the finite-sample check fails, the previous atom result is
kept, but its OLS and MDPDE models are rewritten with the joint offset when that
fallback model still builds finite zero-offset samples. If the forced-sync
fallback itself would become invalid while the previous model was valid, the atom
is marked as suspicious and the previous result is kept unchanged. After the
refit loop, any suspicious atom found this way is expanded through the same
bounded weighted active coupling graph, and every reached atom is rolled back to
the previous iteration state. This keeps post-refit suspicious-offset handling
synchronous with the pre-refit joint-offset check without letting atoms reached
only by propagation become new expansion seeds.

## Acceleration, Ranking, and Convergence

The raw iteration result is treated as the fixed-point map output `G(x)`.
Second-stage fitting keeps short internal Anderson Acceleration histories for
active-atom clusters. Clusters are rebuilt each outer iteration from selected
atom samples: the active target atom and active selected neighbors that
contribute to the same sample residual are connected, and connected components
keep independent histories. When a cluster has at least one compatible prior
pair `(x, G(x))`, with residuals `G(x) - x`, the stage solves a constrained least
squares problem over scaled residuals for that cluster and proposes:

```text
candidate = sum(gamma_i * G(x_i)), where sum(gamma_i) = 1
```

Residuals are scaled per parameter using the same normalized-change scale floor
used by convergence checks, so amplitude does not dominate width and offset.
Active-set changes no longer clear unrelated histories globally: clusters whose
active members are unchanged keep their history, while merged, split, missing,
or newly created clusters start fresh. Candidates with non-finite active
parameters, non-positive active widths, invalid coefficients, or excessive
extrapolation are discarded for that cluster.

Each outer iteration tries the Anderson candidate with damping values `1.0`,
`0.5`, and `0.25`. A damping value applies as:

```text
damped = previous + damping * (candidate - previous)
```

The first accepted Anderson attempt is logged as `acceleration = aa`; lower
damping is logged as `acceleration = damped-aa`. A full iteration candidate
starts from the previous state. Anderson damping is tried per cluster first; each
cluster that passes its local objective gate writes only its active atoms into
the assembled candidate. Clusters with missing, invalid, or rejected Anderson
candidates remain pending and then try the damped fixed-point sequence logged as
`acceleration = damped-fixed-point`. If no cluster can produce an Anderson
candidate, the stage starts directly with fixed-point damping. Anderson history
candidate construction only checks structural compatibility; the damped
candidate that would actually be applied must still be finite and have positive
active widths before it can reach objective scoring.

Objective samples are owned by clusters. A selected sample belongs to the
connected component containing the active target atom and active selected
neighbors that contribute to that sample residual. Samples without active
contributors are ignored by the cluster gate for that iteration, and a sample is
never counted by more than one cluster. Each cluster keeps its own residual-scale
tracker, previous objective samples, best local objective stats, and objective
ridge multiplier. During residual-scale warm-up, that cluster's previous,
current, and best objective values are re-scored with the same provisional scale
before the backtracking decision is made.

If a cluster candidate lacks a finite objective while its local reference exists,
or is worse than that cluster's previous or best tracked state by more than
`kLocalFittingConvergenceObjectiveRelativeTolerance`, only that cluster is
rejected. Other accepted clusters still update the assembled candidate and become
the next previous state. Rejected Anderson clusters clear and suppress only their
own history; a suppressed cluster must receive accepted fixed-point progress
before Anderson is enabled again. Rejected clusters increase only their
cluster-local objective ridge multiplier for the next outer iteration. The
global ridge ratio increases only when no active cluster is accepted and the
rejected clusters' local objective ridge multipliers are already saturated.
Suspicious-offset rollback clears and suppresses only clusters containing
rolled-back atoms; unrelated accepted clusters may continue to commit history.

Each accepted assembled candidate computes absolute and normalized parameter
movement for amplitude, width, and offset for every selected atom. Active atoms
are summarized by the 95th percentile, and the accepted candidate's normalized
movement drives convergence checks. Freeze tracking is updated only for active
atoms in clusters that accepted progress; rejected clusters keep their previous
state and cannot be frozen merely because their accepted movement is zero.
Parameter convergence requires all three active-set normalized percentile
changes to be below `kLocalFittingNormalizedChangeTolerance`, no suspicious
offset rollback, no cluster objective rejection in that iteration, and every
active cluster objective scale with a reference to be locked.

Atoms are frozen when their maximum absolute parameter movement stays below
`sqrt(kLocalFittingParameterChangeTolerance) * 0.1` for three consecutive active
iterations. A frozen atom can be thawed again when a currently active selected
neighbor changes by at least `sqrt(kLocalFittingParameterChangeTolerance)`.
Dependency thawing applies per-atom hysteresis: each dependency thaw raises that
atom's next dependency-thaw threshold, up to a capped multiplier, and the
multiplier decays back toward the base threshold while the atom remains frozen.
To prevent flat-region freeze/thaw thrashing from consuming the full iteration
budget, dependency thawing is also capped per atom within one second-stage run:
after five successful neighbor-triggered thaws, later dependency-thaw requests
for that atom are denied and the atom stays frozen at its current local state.
Frozen atoms do not participate in the joint offset solve or per-atom refit while
they remain frozen, but their fitted Gaussian remains in the snapshot so active
neighbors can subtract them as fixed signal contributions.
Suspicious offset cluster members are thawed after freeze tracking so the next
iteration can retry them with the temporary per-atom ridge multiplier; those
forced retry thaws are not counted against the dependency thaw cap.

Objective scoring is cluster-local. At second-stage entry, each cluster seeds
its residual normalization scale from the cluster-owned objective samples when
they are finite. Each residual scale sample is floored by a small fraction of
the robust response scale from the same objective samples, then by the absolute
Huber scale minimum, so a near-perfect entry fit cannot create an overly
sensitive denominator. A local candidate is better when its normalized robust
objective improves beyond the tie tolerance; objective ties are broken by the
maximum of the three normalized percentile parameter changes. Global progress
and terminal logs report parameter movement, acceleration, active/frozen/thawed
atom counts, and terminal reasons; they do not compute a full-state objective.

## Exit Paths

The loop has four terminal cases:

- **All atoms frozen:** apply the previous state if the loop starts with no
  active atoms. If the last accepted update froze the remaining active atoms,
  apply the current accepted accelerated result, then log an info message when logging is
  enabled.
- **Parameter convergence:** apply the current accepted accelerated iteration result when the
  iteration had no suspicious offset rollback, no cluster objective rejection,
  all three active-set normalized percentile changes are below
  `kLocalFittingNormalizedChangeTolerance`, and all active cluster objective
  scales with references have finished warm-up, then log an info message when
  logging is enabled.
- **Objective backtracking failure:** if all clusters are still rejected after the
  acceleration damping sequence, rejected cluster-local ridge multipliers are
  saturated, and global-ridge retry is no longer possible, apply the current
  previous state and log a warning when logging is enabled.
- **Maximum iterations reached:** apply the current accepted assembled candidate
  and log a warning when logging is enabled.

On non-terminal accepted iterations, the accepted estimation and result vectors
become the previous state for the next loop iteration. A rejected iteration does
not update the fixed-point state or the freeze tracker.

When logging is enabled, the function also emits a warning summary after the
loop if any iteration used the joint-offset fallback, any atom refit fallback, or
any suspicious offset rollback. The summary reports joint-offset fallback
iterations, refit fallback atom-events, distinct atoms that used the refit
fallback, suspicious offset atom-events, and distinct atoms that were marked for
suspicious offset rollback.

## Related Notes

- In `RunPotentialFittingWorkflow`, `RunLocalAlphaTraining` runs before local
  fitting, then `InitializeLocalFittingSeedModels` and `RunFirstStageLocalFitting`
  seed the per-atom local Gaussian results before this stage runs. After this
  stage, the workflow runs group alpha training, builds updated sampling entries
  from group-median local Gaussians, runs group potential fitting, rebuilds
  updated sampling entries from fitted group Gaussians, retrains local alpha on
  updated samples, then runs `RunThirdStageLocalFitting`. Group alpha training
  and group potential fitting run again after third-stage local fitting.
- [`estimate-local-gaussian-with-offset.md`](/docs/developer/estimate-local-gaussian-with-offset.md)
  documents the single-atom fixed-offset estimator used by related local fitting
  paths.
- Second-stage local fitting does not update group Gaussian results; group
  fitting is handled by `RunGroupPotentialFitting`.
