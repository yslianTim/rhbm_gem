# Joint local stability practicalization

## Baseline contract

Local stability now reports four separate pieces of evidence for the accepted
global slice: numerical validity, conditional local optimality, independent
replay, and diagnostic eligibility. The original accepted global slice remains
the baseline; fixed-width A/C reprofiling remains an audit comparison and does
not replace the Joint estimates.

An accepted slice is eligible for local perturbation diagnostics when it is
finite and feasible, its local A/C projected KKT is at most the global
practical `1e-7`, its bound-projected Width KKT is at most the global practical
`1e-9`, and independent replay passes with those same practical criteria.
The result is labeled `PracticalOnly` when it meets those criteria but not the
strict local A/C (`1e-10`) and Width (`1e-12`) criteria. This is diagnostic
eligibility only; it does not mark the strict A/C solver certificate as
passed, alter global qualification, or alter runtime convergence.

Conditional reprofiling and perturbation endpoints still require the existing
strict constrained A/C certificate. Perturbation endpoints also keep the
existing strict independent replay. The explicit local projected Width check
uses `1e-9`; its scaling is the same normalized gradient used by global
qualification. Bound feasibility is required in both the primary state and
independent replay. The `1e-8` parameter and prediction warning thresholds are
unchanged.

## Upper-bound directions

For each signed perturbation, the diagnostic stores requested and effective
log-width displacement separately. If upper-bound projection removes the
entire perturbation, the seed is `ZeroDisplacement` and no profile solve is
started. If projection would change a requested outward direction, it is
reported as `BoundInfeasibleDirection`; that altered projected direction is
not substituted for the requested one. The opposite feasible inward direction
can still be solved and is labeled `OneSidedAssessed`. A block cannot become
Stable when every direction has zero effective displacement or no perturbation
solve was performed.

## Regression evidence

Focused regression coverage includes a practical-qualified baseline whose
local strict A/C KKT fails, bounded and unbounded perturbations, active-bound
inward/outward and zero-displacement handling, scale-invariant projected
Width KKT, baseline replay rejection, local A/C solve failure, persistence
round-trip, unchanged estimates, and unchanged `RuntimeConvergence`.

The P0 6Z6U unsharpened block remains the numerical reference: local A/C KKT
`3.5208046721135755e-10` and projected Width KKT
`5.4610760358286825e-12` are within practical criteria but outside the former
post-fit strict checks; practical independent replay passes. Its conditional
reprofile has strict A/C KKT `1.3877787807814457e-17`, objective difference
`4.1257784491995526e-19`, and relative prediction difference
`1.1157094118308093e-8`. The map-level rerun after practicalization is recorded
with P2's probe-coverage comparison.

The P1 6Z6U simulation rerun retained the previous fit: 1,418 targets,
180,478 observations, objective `1.2330854370811315e-13`, global A/C KKT
`9.434574186999167e-8`, projected Width KKT `8.776510069988752e-10`, 83
sweeps, 22,908 profile evaluations, and `RuntimeConvergence=Passed`. The selected
block remained `1/block-117` and `Stable`; its three directions used six
replayed perturbation solves and 24 profile evaluations. The audit reprofile
objective difference was `3.037728304414033e-22`. Measured Peak RSS was
`618,905,600 B` and wall time was `303.731 s`. This is within the previously
observed run-to-run resource variation and shows no numerical regression.

The P1 6Z6U unsharpened rerun preserved the endpoint: 1,418 targets,
189,377 observations, objective `0.018417539368903946`, global A/C KKT
`8.174157334426418e-8`, projected Width KKT `9.402726930396454e-10`, two
bound-active targets, 124 sweeps, 34,848 profile evaluations,
`RuntimeConvergence=Passed`, and independent replay passed. Its selected block
`1/block-5` is now numerically valid and baseline-eligible with
`PracticalOnly` conditional optimality: A/C KKT `3.5208046721135755e-10`,
projected Width KKT `5.4610760358286825e-12`, and practical independent replay
passed. Conditional reprofiling remains strictly A/C-qualified with objective
difference `4.1257784491995526e-19` and relative prediction change
`1.1157094118308093e-8`.

The block attempted three perturbation solves and 31 profile evaluations
(`0.080683` process CPU seconds), then reported `Unavailable` because the
existing strict independent replay rejected each first-side perturbation.
The opposite sides were left `NotAssessed`; they were not interpreted as
stable. This is the intended separation: practical baseline eligibility
revealed the block to diagnostics, while the unchanged strict perturbation
replay still prevented an unsupported stability result. Peak RSS was
`1,434,451,968 B`; wall time was `435.112 s`.

The P1 build used `cmake --build build --target rhbm_tests
joint_width_repro_no_instrumentation --parallel 4`, followed by rebuilding
ABI-dependent candidate executables with `cmake --build build --target
joint_component_runtime joint_production_bounded_synthetic --parallel 4`.
The focused stability/persistence filter passed 12/12 tests. The focused CTest
set passed 8/8, including numerical, contract, workflow, bounded synthetic,
component regression, and persistence contracts.
