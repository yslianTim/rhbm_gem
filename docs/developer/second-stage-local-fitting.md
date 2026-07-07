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
    -> apply under-relaxation
    -> compute active-atom p95 absolute and normalized parameter changes
    -> backtrack against the objective; shrink beta and retry when needed
    -> update candidate ranking and adaptive relaxation
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
objective-backtracking rejections increase the ratio for the next recomputed
joint solve, while accepted iterations without backtracking gradually decrease
it.

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
options. Huber weights are then updated from residual median absolute deviation.
The IRLS loop stops when the weighted-ridge surrogate objective would
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

## Relaxation, Ranking, and Convergence

The raw iteration result is under-relaxed before convergence is checked:

```text
relaxed = beta * current + (1 - beta) * previous
```

`beta` starts from the internal `kAdaptiveRelaxationInitialBeta` value (`0.5`)
and is clamped to the local adaptive relaxation range `[0.05, 1.0]`. After each
iteration, the controller looks at the maximum of the three 95th-percentile
normalized parameter changes. A change more than 1% larger than the previous
iteration for three consecutive accepted iterations halves `beta`. Two
consecutive changes more than 1% smaller than the previous iteration allow
`beta` to grow by `1.2x`, up to `1.0`. This is a trend-based adaptive relaxation
rule, not Anderson acceleration. The relaxed vector replaces the candidate MDPDE
model while preserving its standard-deviation model.

When an objective reference is available, the relaxed candidate must produce a
finite objective and pass the objective quality gate before it can update
ranking, freezing, thawing, or the next iteration's previous state. If no
objective reference can be built at all, this quality gate is skipped and
candidate ranking falls back to normalized parameter movement. During
residual-scale warm-up the objective reference is provisional and can include the
candidate's own scale sample; after warm-up it is locked. If the candidate lacks
an objective while a reference exists, or is worse than the previous state or the
best tracked state by more than
`kLocalFittingConvergenceObjectiveRelativeTolerance`, the attempt is rejected.
When another attempt is available and `beta` is above its minimum, the stage
shrinks `beta` and rebuilds the relaxed candidate from the same raw iteration
result. Each outer iteration tries at most three relaxation candidates. If all
attempts fail and `beta` is still above the local minimum, the raw iteration is
rejected and the next outer iteration retries from the unchanged previous state
with the smaller `beta`. Any rejected relaxation attempt also increases the
dynamic joint offset ridge ratio, but that ridge change does not trigger an
immediate refit; it applies when the next outer iteration rebuilds the
joint-offset system.

Each relaxation attempt computes absolute and normalized parameter movement for
amplitude, width, and offset for every selected atom. Active atoms are summarized
by the 95th percentile, and the accepted attempt's normalized movement drives
relaxation updates and convergence checks. Parameter convergence requires all
three active-set normalized percentile changes to be below
`kLocalFittingNormalizedChangeTolerance`, and no suspicious offset rollback may
have occurred in that iteration. Because only candidates accepted by objective
backtracking reach this point, convergence never applies a rejected candidate.

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

The best fixed-point candidate is tracked separately. At second-stage entry,
the initial residuals seed the residual normalization scale when they are finite.
During warm-up, each accepted objective-scored candidate contributes its residual
median absolute deviation to a moving-average scale; after five such accepted
candidates, that average is locked for the rest of the fitting run. Each
residual scale sample is floored by a small fraction of the robust response scale
from the same objective samples, then by the absolute Huber scale minimum, so a
near-perfect entry fit cannot create an overly sensitive denominator. During
warm-up, every previous, current, and best candidate involved in a quality
comparison is re-scored with the same provisional scale before backtracking or
ranking uses the objective values. A candidate is better when its normalized
robust objective improves beyond the tie tolerance; objective ties are broken by
the maximum of the three normalized percentile parameter changes.
Objective acceptance is enforced before convergence is evaluated, when objective
values are available. Parameter convergence is not allowed to terminate the stage
while an available objective scale is still in warm-up.

## Exit Paths

The loop has four terminal cases:

- **All atoms frozen:** apply the previous state if the loop starts with no
  active atoms. If the last accepted update froze the remaining active atoms,
  apply the current relaxed result, then log an info message when logging is
  enabled.
- **Parameter convergence:** apply the current relaxed iteration result when the
  iteration had no suspicious offset rollback, all three active-set normalized
  percentile changes are below `kLocalFittingNormalizedChangeTolerance`, and any
  available objective scale has finished warm-up, then log an info message when
  logging is enabled.
- **Objective backtracking failure:** if a candidate is still rejected after the
  retry limit and `beta` is already at the local minimum, or the maximum
  iteration limit prevents retrying, apply the best tracked candidate when
  available, otherwise apply the unchanged previous state, then log a warning
  when logging is enabled.
- **Maximum iterations reached:** apply the best fixed-point candidate found so
  far and log a warning when logging is enabled.

On non-terminal accepted iterations, the relaxed estimation and result vectors
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
