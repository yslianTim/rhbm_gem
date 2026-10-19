# Bounded Joint production workflow validation — 2026-10-20

## Decision

The bounded internal Joint entry point now has an end-to-end CLI regression
covering calibration, Active-set LM, endpoint qualification, persistence,
reload, PostFitPeelingSamples, Post-Joint Local MDPDE, and Group Fitting. A
separate one-target run confirms that calibration failure is persisted as
unqualified and downstream estimates are not published. The public
`EstimateJointComponents` default remains unchanged.

## Test seam and path

The integration test reuses the existing `JointWorkflowFitHook` test seam. It
replaces the in-memory Joint result with the result of
`FitJointComponentsInitializedBoundedWithProgress`, using the actual
FirstStage-derived initialization and problem created by the normal workflow.
The same `RunJointComponentWorkflow` then captures, persists, reloads, peels,
fits local MDPDE results, and performs group fitting. The candidate result
comes from the bounded internal entry point; the ordinary result produced
before the test hook is discarded and is not persisted. No public CLI option
or second downstream implementation was added.

The success fixture contains 25 selected CA targets. The test invokes the real
CLI `potential_analysis --estimator joint-components`, reloads the saved model
from `DataRepository`, and asserts:

- Joint result qualification and `TargetRuntimeConvergence` survive reload.
- The contributor policy is Pure Halo-Support Exclusion and all 25 target
  mappings remain present.
- Solver metadata identifies the FixedNeighbor outer solver. One
  `active-set-lm-subproblems` evidence item passes.
- `global-width-calibration` passes, its evidence records P95 and Bmax with
  `Bmax = 4 × P95`, and qualified donor coverage meets its persisted threshold.
- Every target has a Second-stage Joint estimate and a persisted peeling and
  Post-Joint local result. At least one local result is Available and at least
  one Group Fitting member result is persisted.

The one-target CLI fixture takes the same bounded candidate route. It fails the
minimum donor-count calibration preflight, reloads with an Unavailable
`global-width-calibration` check, and remains unqualified. The test verifies
that no Second-stage Joint point, peeling, local result, or Group Fitting
member is published. There is no silent unbounded fallback.

## Failure lifecycle evidence

The focused numerical suite exercises the bounded failure surfaces without
changing the qualification gates:

- Calibration: too few qualified donors, poor target coverage, invalid target
  starts, malformed initialization, and unavailable/overflowed donor data.
- Active-set subproblems: upper-face multiplier sign, releases, multiple
  bounds, rank deficiency, infeasible starts, fixed coordinates, and nonfinite
  coefficients.
- Endpoint path: bound feasibility, projected Width KKT, independent replay
  disagreement, subproblem failures, unavailable components, and incomplete
  sweep budget.
- Persistence/workflow: qualified bounded results reload as qualified; an
  unavailable bounded calibration reloads as unqualified and skips downstream
  publication.

The full focused CTest command passed 7/7, including `rhbm_tests_integration_command`,
the Joint numerical, contract and workflow groups, bounded synthetic,
persistence contract, and component runtime runner. The dedicated bounded
CLI integration case also passed individually. Compilation used at most four
jobs and numeric tests set OpenMP and BLAS internal threads to one.

## Scope and remaining gates

This validates the internal bounded candidate through the current Production
CLI and downstream contracts while preserving the current default. It does not
promote the default and it does not validate 37,406-target resource use. Those
remain separate gates. No Notion documents or public CLI options were changed.
