# Fixed-state factor and normal actions

This increment keeps `SearchMethod::LegacyCompact` as the production default.
It changes the internal operator path, not the model, active-set A/C solver,
independent reference, trust-region policy, or uncertainty calculation.

See the [bounded acceptance report](joint-fixed-actions-validation.md) for measured results.

## Rank contract

`EvaluateRank` accepts a `RankRequest`: original `RankPolicy`, decision-column
count, absolute override, and boundary rule. Decision dimensions can exceed a
component compact's dimensions. Operator and tiled derivative preparation and
endpoint X/W/normalized-W/J evidence share this entry point. Normalization and
parent residual scaling remain at the existing callers.

`SvdNative` preserves Eigen's cutoff (including its positive-minimum guard);
`StrictGreater` preserves the explicit `sigma > threshold` spectrum counts.
These intentionally differ at equality. The BDC dispatch, Jacobi retry near the
boundary, solution and right-vector requests are unchanged. Valid deficient
rank is evidence, while decomposition failure is unavailable evidence. Neither
can qualify full rank. No iterative spectrum approximation is introduced.

## Fixed factors and normal action

Only operator-owned SPQR factors use the public one-shot R/E/H/HPinv/HTau API.
The exported sparse R and Householder data belong to that immutable factor;
compact extraction restores the original column order from R. Primary A/C
workspaces still use expert symbolic/numeric factorization and generation checks.
Independent reference factors are not shared with either path.

For normalized free design Z, projection P=ZZ+, raw width derivative D and
owner-form residual contraction T, the operator computes

`J'J v = (D'(I-P)D v + T'(Z'Z)^-1 T v) / s^2`.

`NormalSolve` uses two triangular solves with the factor's column permutation.
No normal matrix or inverse is formed. The two derivative terms lie in
orthogonal subspaces in exact arithmetic; nonzero residual correction is kept.
Finite-precision parity, symmetry, curvature and weak directions are tested.
PCG now uses this action. Its tolerance, refresh policy, damping, metric,
predicted reduction from `Apply(step)`, and failure behavior are unchanged.

## Cost interpretation

Preparation includes design preparation, factor acquisition and rank work.
Rank work includes compact extraction and SVD. SVD retry time includes its
fallback work. Action timers and factor-operation timers are nested; do not
sum them. One-shot SPQR factorization has one combined timer, not guessed
symbolic/numeric components. Eigen's opaque least-squares solve is timed as
`eigen_least_squares_seconds`; its Q and triangular call counts are recorded,
but its combined time is not assigned to the Q-only timer.

R nonzeros are not full factor memory. `exported_factor_bytes` describes only
the exported R/H/permutation arrays, excluding the retained sparse design,
allocator overhead and factorization scratch; unavailable values are null.
Process-tree RSS and OS peaks remain separate measurements. Dense probes are
known allocations, not an allocator trace. A p-by-p compact still exists.

## Reproduction

Configure both EIGEN and SPQR Release builds with BUILD_TESTING=ON,
RHBM_GEM_DEP_PROVIDER=SYSTEM, RHBM_GEM_ENABLE_UMAP=OFF,
RHBM_GEM_ROOT_MODE=OFF and RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON.
Build `rhbm_tests`, `joint_sparse_benchmark`, `joint_component_runtime`,
`joint_validation` and `rhbm_gem_cli` for candidate regression and measurement.

Export pristine `ec24f3056427408cb1978cff004ed477908fb0d0` to a separate source
tree and apply `figures/joint-fixed-actions/baseline-instrumentation.patch`.
This patch only adds counters/timers; it retains the original numerical actions.
Copy the current benchmark driver and `tests/support/JointFixedDiagnostic.hpp`
into that tree and add `PR4_BASELINE_DRIVER` to the benchmark target definitions.
Build `joint_sparse_benchmark` for both backends and `joint_validation` for Eigen.
The runner verifies the frozen instrumented source hash and identical wrappers.

Run `tests/integration/joint_fixed_validation.py` with `--work-dir`,
`--input-dir`, `--baseline-eigen`, `--baseline-spqr`, `--eigen`, and `--spqr`.
Input files may reuse frozen PR0/PR1 inputs; their hashes are saved. A new
campaign directory is mandatory. `--compare --work-dir ...` only regenerates
comparison output, never executions. The original PR0-PR3 receipts are untouched.

The dedicated benchmark modes are `--fixed freeze|composed|normal` and
`--state FILE`, with the existing initial/fixture/synthetic-fixed input forms.
Freeze uses the baseline Eigen profile and replay. Every subsequent process
reconstructs that exact beta/eta and context with `EvaluateState`, avoiding
backend-dependent re-profiling. No reference solve, derivative reduction,
assessment, search, or uncertainty runs in the measured diagnostic path.

A uses baseline factors and composed J'/J; B uses new factors and composed
J'/J; C uses new factors and ApplyNormal. Each run measures actions and one
fixed-damping solve (mu=1e-3) for identity, diagonal and Schwarz. Per-step totals
include operator preparation, separately timed gradient setup, metric and applicable partition/model/factor
construction, plus solve/prediction. Input/basis setup and audit output are
separate. A fixed step is not a trust-region acceptance or convergence claim.

The single campaign has a 1,800-second deadline including state preparation,
300 seconds per process, 4 GiB sampled process-tree RSS and one numerical thread.
Sampling can overshoot. Runs are serial, ordered chain-8, single-128,
heterogeneous-168, then optional single-512. Three independent processes are
required per backend/mode; A/B/C order rotates each repetition. A resource stop
ends that group's remaining repetitions. There is no restart or budget increase.

Comparisons check same-state residual/objective/gradient, rank/spectrum,
J/J'/normal actions, all three steps and true residuals. Missing/invalid results
never qualify. Only three completed numerical matches, a lower median and
three faster paired times establish an observed phase improvement. A faster
phase with a slower total is not an overall speedup. This diagnostic never
promotes the production search backend or claims large-component scalability.

The archived campaign predates the gradient-timer correction. Its total_seconds
fields are explicitly treated as subtotals by postprocessing; see the acceptance
report for the retained limitation and conservative tracked-phase bounds.
