# Second-Stage Local Fitting

`rhbm_gem::core::RunSecondStageLocalFitting` refines selected atom-local
Gaussian estimates after the first-stage per-atom fit. The implementation is in
[`src/core/GaussianEstimator.cpp`](/src/core/GaussianEstimator.cpp).

This stage is a fixed-point iteration across atoms. Each iteration estimates a
joint offset update for the active atoms, refits those atoms with neighbor
contributions removed, then uses percentile parameter movement to decide whether
the active set has converged.

## Inputs and State

The function receives:

- a `ModelObject`, used to edit each atom's local analysis entry;
- the selected `AtomObject` list to refit; and
- `FitOptions`, which supplies fit distance range, relaxation factor, thread
  count, and logging mode.

At entry, it builds one `AtomLocalPotentialEditor` per selected atom and reads
the current `LocalGaussianResult` from each atom. The previous iteration state is
stored as a `LocalFittingState` with two aligned vectors:

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
    -> build fitted Gaussian snapshot
    -> estimate joint offsets for active atoms
    -> refit active atoms with neighbor contributions removed
    -> apply under-relaxation
    -> compute active-atom p95 parameter changes
    -> freeze stable atoms
    -> exit, fallback, or continue
```

The maximum iteration count is `kLocalFittingMaximumIterations` (`100`).

## Joint Offset Step

`RunLocalFittingIteration` first converts the full previous per-atom vectors
into a snapshot keyed by atom pointer. `EstimateJointOffsets` then solves one
sparse linear system for active atom offsets. Frozen atoms remain in the
snapshot, but they do not become columns in the linear system.

For each unfiltered sampling entry on each active atom:

1. subtract the target atom's current zero-offset signal;
2. add the target atom's offset basis as the row entry for that atom;
3. subtract selected neighbors' current zero-offset signals when those neighbors
   are present in the snapshot and their distance is within
   `kNeighborContributionDistanceMax`; and
4. add each active selected neighbor's offset basis as another row entry.

The resulting system is solved with weighted ridge regression. The ridge term is
relative to the previous offsets, so weakly constrained columns stay close to
their prior values. Huber weights are then updated from residual median absolute
deviation until the maximum offset movement is below `kHuberSlopeTolerance` or
the Huber iteration limit is reached.

If the joint system cannot be built, is empty, or cannot be solved, the step
falls back to the previous offsets. The rest of the local fitting iteration still
runs.

## Per-Atom Refit

After joint offsets are attached to the snapshot, each active atom is refit.
Frozen atoms are left in the iteration state copied from the previous state. The
active loop may run under OpenMP using `FitOptions::thread_size`.

For each atom:

1. build a sample list by subtracting fitted neighbor responses from the atom's
   unfiltered local sampling entries;
2. use the joint-offset snapshot model as the fixed offset model;
3. call `EstimateLocalGaussianWithOffsetModel` to fit amplitude and width with
   that fixed offset and `FitOptions` distance limits; and
4. accept the candidate only if `CanBuildFiniteZeroOffsetSamples` confirms that
   subtracting the candidate offset remains finite and inside the `float`
   response range.

If fitting throws or the finite-sample check fails, the previous atom result is
kept, but its OLS and MDPDE models are rewritten with the joint offset. This
keeps the iteration state aligned with the joint offset update even when a local
shape refit cannot be accepted.

## Relaxation, Ranking, and Convergence

The raw iteration result is under-relaxed before convergence is checked:

```text
relaxed = beta * current + (1 - beta) * previous
```

`beta` starts from `FitOptions::relaxation_factor` and is clamped to the local
adaptive relaxation range `[0.05, 1.0]`. After each iteration, the controller
looks at the maximum of the three 95th-percentile parameter changes. A change
more than 1% larger than the previous iteration immediately halves `beta`.
Two consecutive changes more than 1% smaller than the previous iteration allow
`beta` to grow by `1.2x`, up to `1.0`. This is a trend-based adaptive relaxation
rule, not Anderson acceleration. The relaxed vector replaces the candidate MDPDE
model while preserving its standard-deviation model.

The stage then computes absolute parameter movement for amplitude, width, and
offset for every selected atom. Active atoms are summarized by the 95th
percentile. The candidate is considered converged only when all three active-set
squared percentile changes are below `kLocalFittingParameterChangeTolerance`.

Atoms are frozen when their maximum absolute parameter movement stays below
`sqrt(kLocalFittingParameterChangeTolerance) * 0.1` for three consecutive active
iterations. Frozen atoms are not thawed. They no longer participate in the joint
offset solve or per-atom refit, but their fitted Gaussian remains in the
snapshot so active neighbors can subtract them as fixed signal contributions.

The best fixed-point candidate is tracked separately. A candidate is better when
the maximum of its three percentile changes is lower than the previous best.

## Exit Paths

The loop has three terminal cases:

- **All atoms frozen:** apply the previous state if the loop starts with no
  active atoms, or the current relaxed iteration result if the last update froze
  the remaining active atoms, then log an info message when logging is enabled.
- **Parameter convergence:** apply the current relaxed iteration result when all
  three active-set squared percentile changes are below
  `kLocalFittingParameterChangeTolerance`, then log an info message when logging
  is enabled.
- **Maximum iterations reached:** apply the best fixed-point candidate found so
  far and log a warning when logging is enabled.

On non-terminal iterations, the relaxed estimation and result vectors become the
previous state for the next loop iteration.

## Related Notes

- `RunFirstStageLocalFitting` seeds the per-atom local Gaussian results before
  this stage runs.
- [`estimate-local-gaussian-with-offset.md`](/docs/developer/estimate-local-gaussian-with-offset.md)
  documents the single-atom fixed-offset and residual-offset estimators used by
  related local fitting paths.
- Second-stage local fitting does not update group Gaussian results; group
  fitting is handled by `RunGroupPotentialFitting`.
