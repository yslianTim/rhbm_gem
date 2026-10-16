# Joint component estimator runtime

Include `<rhbm_gem/core/JointComponentEstimator.hpp>` and link `RHBM_GEM::rhbm_gem`.
The installed API uses the production Joint estimator directly.

```cpp
model.SelectAllAtoms();
auto result = rhbm_gem::core::EstimateJointComponents(map, model);
if (result.problem && result.initialization.valid) {
    auto repeated = rhbm_gem::core::FitJointComponentsInitialized(
        *result.problem, result.initialization);
}
```

`BuildJointProblem` creates the raw observation domain and contributor mapping;
it does not create a qualified fixed-halo snapshot. A width vector by itself
cannot initialize the supported halos. `FitJointComponents` accepts one
initial width per selected target and requires a problem that already owns the
complete qualified FirstStage `(A, B, 0)` snapshot. `FitJointComponentsInitialized`
checks the supplied FirstStage records and attaches that snapshot before fitting.
FirstStage initialization has its own qualification result, independent of
whether the subsequent Joint target fit converges. `EstimateJointComponents`
is the convenience API that performs problem construction, FirstStage
initialization, and target-only fitting in sequence.

`JointProblem::WithFixedHalos` accepts exactly the supported non-target
contributors. Each value must carry qualified FirstStage provenance and C=0;
the returned snapshot stays constant during optimization. Structural
components, FixedNeighbor blocks, KKT, width gradients, conditioning, and
stability cover selected targets only. Full prediction, independent replay,
and PostFitPeeling include the fixed halo signal.

## Objective and production search

For each observation, Joint uses the FullABC contribution
`A_i G(r; B_i) + C_i K(r; B_i)`, with `A_i >= 0`, `B_i > 0`, and signed `C_i`.
`C_i` multiplies the charge kernel; it is not a constant background. The fixed
observation-domain objective is

```text
0.5 * ||prediction - observations||^2 / max(1, ||observations||)^2
```

There is no global regularizer, prior, pseudo-observation, or added width bound.
The objective and FixedNeighbor update strategy are unchanged by the diagnostic
architecture.

For `first-stage-fixed-halo-v1` selected-target problems, only selected targets
have free Joint parameters. Every supported halo contributes its qualified
FirstStage `(A, B, 0)` snapshot as a fixed response; the solver fits targets to
the raw observations after subtracting that response, with the scale still
computed from the raw observations. Numerical convergence therefore certifies
the target optimization conditional on those fixed halo values. It does not
establish statistical accuracy of the FirstStage halo parameters. Halo
parameter validity and target point quality are reported separately; halo model
error requires truth-controlled simulation evidence.

Production uses disjoint mutable cores of 12 atoms, forward serial
Gauss-Seidel sweeps, at most 200 complete sweeps, and at most one trusted
Profile LM update per block visit. Neighbor parameters stay fixed during each
block solve. Every accepted update passes the existing trust checks and a full
global-objective monotonic replay. The local constrained A/C solve uses
SparseQR and keeps the nonnegative-amplitude active set.

When a block visit produces no trusted local profile, FixedNeighbor preserves
the current accepted component state, records the failure, and completes the
remaining blocks in that sweep. It retries that block once on the next sweep,
rebuilding the effective response from the then-current neighbor state. A
second failed visit stops with `block-recovery-exhausted`; reaching the sweep
budget with a pending retry stops with `block-recovery-sweep-budget`. Any sweep
containing a failed block is excluded from stationarity confirmation.

## Global numerical convergence

A component becomes a trusted scientific state only after a confirmed
stationary sweep and endpoint verification. Production still defaults to the
practical global criteria (practical-v1): projected A/C KKT `<= 1e-7`,
log-width gradient infinity norm `<= 1e-9`, and confirmed log-width change
`<= 1e-6`. Strict research reference criteria (strict-v1) remain available at
`1e-10`, `1e-12`, and `1e-10`, respectively. Result candidate diagnostics save
the criteria identifier, all three thresholds, and the final complete-sweep
log-width change alongside the achieved KKT and width-gradient values. Both
policies require finite, nonnegative measures and inclusive thresholds;
confirmation requires a complete prior sweep.
Changing the global criteria does not change the local SparseQR/active-set
qualification, profile trust checks, or replay-consistency enclosures. At the
accepted endpoint, a primary streaming evaluation and an independent streaming
replay check feasibility, prediction, residual, objective, KKT and
width-gradient consistency under the selected global stationarity criteria.

For fixed-halo problems, the final endpoint check also rebuilds the complete
contributor prediction from the original observations and immutable halo
snapshots. It verifies the full/effective residual and objective identities,
original row and contributor mappings, and target-only projected KKT and
log-width gradients. Fixed halo gradients are not convergence conditions.

Endpoint replay alone does not establish convergence. A valid replay with a
recoverable numerical failure can be retried at a bounded interval while
stationarity holds. An invalid primary state or an endpoint-verification
resource limit stops the search. Reaching the sweep limit is budget exhaustion;
it cannot mark the search converged. If the last sweep has not been verified,
one final bounded verification can run, but it passes only when that sweep also
meets the complete stationarity criteria.

`JointComponentResult::RuntimeConvergence()` depends on numerical convergence,
endpoint verification, and a valid state. An assembled result additionally
checks component and observation mappings, constant rows, prediction and
objective consistency, and convergence of every required component. Missing
component states are not filled with zeros. Local conditioning and parameter
stability do not reject a numerically converged state.

## Local conditioning

After convergence, each FixedNeighbor block is assessed at the final accepted
A/C/log-width values, with its external neighbors held fixed. The assessment
reports separate evidence for:

- Free A/C design conditioning on the active face.
- Conditional projected-width conditioning.
- The full profiled width derivative, including the nonzero-residual
  coefficient correction.

The reported ranks, singular values, conditions, active-face notes, and
unavailable reasons describe a conditional local block. They do not certify or
rule out global parameter identifiability. Zero columns, exact or near rank
deficiency, active amplitudes, nonfinite values, and unavailable decompositions
remain explicit diagnostic outcomes. The accepted state and global objective
are read-only during this assessment.

For fixed-halo problems, these blocks contain targets only; fixed halo columns
are not included in conditioning or parameter-stability probes.

Historical Free-Halo FullABC result JSON remains readable under its saved
schema and parameterization contract. Loading that data does not reinterpret
the free halo estimates as a FirstStage fixed-halo snapshot; a new fixed-halo
fit must pass its own FirstStage qualification.

Each design and Jacobian is column-normalized before its spectrum is summarized,
so the reported condition numbers are dimensionless within their own matrix.
The provisional warning classifies a full-rank matrix as ill-conditioned when
`kappa * epsilon_double >= 1e-8` (roughly `kappa >= 4.5e7`). This is a binary64
roundoff-amplification indicator, not a bound on parameter error or statistical
uncertainty. Exact rank loss remains `Singular`; an active amplitude boundary
remains `ActiveFace`; threshold-near cases retain their rank margin and
threshold-sensitive flag. These classifications do not enter convergence.

## Parameter stability and point quality

Cross-sweep telemetry retains only the most recent three complete-sweep
changes. Fewer than three sweeps is reported as insufficient history; the
trajectory summary is evidence, not a convergence threshold.

Blocks are ranked after local conditioning by class (singular, ill-conditioned,
active-face, then well-conditioned), descending maximum available
`log10(condition number)`, and stable block identity. Exact free A/C rank loss
is skipped when the local constrained refit cannot support a width probe. At
most eight blocks are probed per fit, with at most one well-conditioned
reference block per component.

Within a selected block, unique directions are tried in this order: the weakest
right singular direction of the column-normalized projected-width Jacobian,
the coordinate with the largest component in that vector, and the normalized
common mode. Equivalent directions are deduplicated; a missing weak-width
direction stays unavailable. A/C singular vectors are never converted into
width directions. Each direction gets positive and negative `1e-3` log-width
seeds, with outside neighbors fixed and temporary local state. Each seed is
limited to 12 profile evaluations including its final endpoint evaluation, six
accepted updates, and the remaining portion of a two-second CPU budget per
block. The total caps are 24 directions, 48 seeds, 576 profile evaluations,
and 16 CPU seconds across eight selected blocks. Memory preflight remains
block-local.

Parameter and prediction changes use separate normalized summaries and a
provisional `1e-8` warning tolerance; the objective difference is recorded
separately, while objective consistency continues to mean the independent
replay check. `Stable` requires every available direction to finish and both
parameter summaries and prediction change to stay within their tolerances.
Parameter change with similar predictions is `Sensitive`; changed predictions
are `Unstable`. Failed or unselected probes remain `Unavailable` or
`NotAssessed`, never `Stable`. These probes measure conditional local behavior,
not global identifiability.

Target point validation builds compact atom-to-block indexes for conditioning
and stability once, then maps points and selected targets through those
indexes. Its diagnostic mapping work is `O(P + B + M)`, where `P` is point
count, `B` is block count, and `M` is the total number of block-core
memberships. Duplicate identities, duplicate memberships, missing mappings,
and out-of-range atom indices are reported as mapping evidence; they do not
become convergence gates.

The `target-point-validation-v1` evidence reports component convergence,
parameter validity, local conditioning, parameter stability, point quality,
and a reason. A conditioning warning or unstable parameter probe is a point
quality warning and does not change the fitted estimate. Prediction stability
and parameter stability are separate: similar predictions do not imply similar
parameters. This evidence makes no global identifiability claim.

## Memory planning

`ResourceWork` preserves preflights for problem construction, structural
construction, FixedNeighbor preparation, global streaming evaluation, endpoint
verification, reference solve, local conditioning, and parameter stability.
The planner uses checked byte arithmetic, actual support and affected-row
counts, the largest block, resident bytes, and a 3,200,000,000-byte estimated
peak budget below the 4,000,000,000-byte hard limit. These are allocation
estimates, not observed process RSS. A memory rejection is reported with its
stage and reason.

Production endpoint verification streams over support and does not allocate a
global dense design or global rank factor. Local conditioning allocates only
bounded block workspace. Parameter stability has its own local workspace
estimate and bounded solve counts. Target evidence stores two compact
atom-to-block indexes proportional to points and core memberships; it does not
build a global rank matrix. Construction and FixedNeighbor preparation remain
separate possible large-input bottlenecks.

## Persistence and downstream use

Current point evidence uses `target-point-validation-v1`. Historical
`target-estimability-v1` through `target-estimability-v4` keeps its saved
meaning and is never converted into point-quality evidence. Historical rank
records remain decoding-only provenance; they are not reinterpreted as local
conditioning or stability. Joint results store the FullABC state for selected
targets and fixed FirstStage halo parameters separately under
`target-only-fixed-halo-v1` (Joint JSON schema 14). Earlier schemas retain their
original parameterization and are not reinterpreted. Stage records mark fixed
halos as FirstStage-sourced and not Joint-optimized. Post-fit peeling
reconstructs the complete target-plus-halo prediction; Post-Joint local MDPDE
and group fitting remain limited to selected targets and keep their own
qualification contracts.

The runtime regression harness is `tests/integration/joint_component_runtime.py`.
The ordinary Joint CTest groups cover component contracts, workflow,
FixedNeighbor numerics, local conditioning, stability, persistence, and
historical decoding. Optional offline tests cover high-precision local
numerical controls; retired global-rank and weak-direction certification tools
are not part of the production or test entry points.
