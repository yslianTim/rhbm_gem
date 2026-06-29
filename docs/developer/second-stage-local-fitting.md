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

At entry, it builds one `AtomLocalPotentialEditor` per atom in the input list and
reads the current `LocalGaussianResult` from each atom. The previous iteration
state is stored as a `GaussianFittingState` with two aligned vectors:

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
    -> reject, shrink beta, and retry if the objective deteriorates
    -> update candidate ranking and adaptive relaxation
    -> freeze stable atoms and thaw changed selected neighbors
    -> exit, fallback, or continue
```

The maximum iteration count is `kLocalFittingMaximumIterations` (`1000`).

## Joint Offset Step

`RunLocalFittingIteration` first converts the full previous per-atom vectors
into a snapshot keyed by atom pointer. `EstimateJointOffsets` then solves one
sparse linear system for active atom offsets. Frozen atoms remain in the
snapshot, but they do not become columns in the linear system.

For each unfiltered sampling entry on each active atom:

1. subtract the target atom's current zero-offset signal;
2. add the target atom's offset basis as the row entry for that atom;
3. for selected neighbors present in the snapshot and within
   `kNeighborContributionDistanceMax`, subtract the full fitted response when
   the neighbor is not in the active solve; and
4. for active selected neighbors, subtract the current zero-offset signal and add
   the neighbor's offset basis as another row entry.

The resulting system is solved with weighted ridge regression. The ridge term is
relative to the previous offsets, so weakly constrained columns stay close to
their prior values. Huber weights are then updated from residual median absolute
deviation. The IRLS loop stops when the weighted-ridge surrogate objective would
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

When a fixed objective reference is available, the relaxed candidate must pass
the objective quality gate before it can update ranking, freezing, thawing, or
the next iteration's previous state. If the candidate is worse than the previous
state or the best tracked state by more than
`kLocalFittingConvergenceObjectiveRelativeTolerance`, the stage shrinks `beta`
and rebuilds the relaxed candidate from the same raw iteration result. Each
outer iteration tries at most three relaxation candidates. If all attempts fail
and `beta` is still above the local minimum, the raw iteration is rejected and
the next outer iteration retries from the unchanged previous state with the
smaller `beta`.

The stage then computes absolute parameter movement for amplitude, width, and
offset for every input atom. Active atoms are summarized by the 95th percentile.
Parameter convergence requires all three active-set squared percentile changes
to be below `kLocalFittingParameterChangeTolerance`. Because only candidates
accepted by objective backtracking reach this point, convergence never applies a
rejected candidate.

Atoms are frozen when their maximum absolute parameter movement stays below
`sqrt(kLocalFittingParameterChangeTolerance) * 0.1` for three consecutive active
iterations. A frozen atom can be thawed again when a currently active selected
neighbor changes by at least `sqrt(kLocalFittingParameterChangeTolerance)`.
Frozen atoms do not participate in the joint offset solve or per-atom refit while
they remain frozen, but their fitted Gaussian remains in the snapshot so active
neighbors can subtract them as fixed signal contributions.

The best fixed-point candidate is tracked separately. At second-stage entry,
the initial residuals define one fixed residual scale for this fitting run.
Every previous, current, and best candidate is scored by first normalizing its
residuals by that entry-time scale, then applying the same dimensionless Huber
cutoff. This keeps objective values comparable across iterations and across map
value scales. A candidate is better when its normalized robust objective
improves beyond the tie tolerance; objective ties are broken by the maximum of
the three percentile parameter changes.
Convergence is accepted only if the current objective is not worse than both the
previous candidate and the best candidate by more than
`kLocalFittingConvergenceObjectiveRelativeTolerance`, when those objective values
are available.

## Exit Paths

The loop has four terminal cases:

- **All atoms frozen:** apply the previous state if the loop starts with no
  active atoms. If the last accepted update froze the remaining active atoms,
  apply the current relaxed result, then log an info message when logging is
  enabled.
- **Parameter convergence:** apply the current relaxed iteration result when all
  three active-set squared percentile changes are below
  `kLocalFittingParameterChangeTolerance`, then log an info message when logging
  is enabled.
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
loop if any iteration used the joint-offset fallback or any atom refit fallback.
The summary reports joint-offset fallback iterations, refit fallback atom-events,
and distinct atoms that used the refit fallback.

## Related Notes

- `RunFirstStageLocalFitting` seeds the per-atom local Gaussian results before
  this stage runs.
- [`estimate-local-gaussian-with-offset.md`](/docs/developer/estimate-local-gaussian-with-offset.md)
  documents the single-atom fixed-offset and residual-offset estimators used by
  related local fitting paths.
- Second-stage local fitting does not update group Gaussian results; group
  fitting is handled by `RunGroupPotentialFitting`.
