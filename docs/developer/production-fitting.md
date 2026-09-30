# Production fitting

## Inputs and statistical contract

Production fitting consumes sampled map responses, selected atomic geometry,
explicit fitting options, and previously estimated state. It never loads the
simulation manifest, generated charges, truth widths, or scoring output. Atomic
number is not an amplitude truth initializer. `potential_analysis --simulation`
continues to skip normalization; `-r supplies model metadata` only.

The shape solve uses the production MDPDE fixed-point equations, and selected
offsets use joint-offset IRLS. Parameter accuracy alone does not establish
convergence. Inner solver failures retain their native status and remain
reproducible.

The [Second-stage local fitting specification](second-stage-local-fitting.md)
owns outer iteration, candidate acceptance, recovery, convergence, and
finalization.

## Solver endpoints and qualification

Each nominal shape endpoint stores its MDPDE status, iteration count, terminal
squared beta change, relative variance change, and variance. Each offset
endpoint stores its joint solve status, iteration count, normalized change,
and robust scale. Quarantine may require a separate unrestricted shape solve;
its status cannot be borrowed from the constrained proposal. Every nominal
endpoint must be qualified, including coordinates excluded from accepted
movement statistics.

## Failed-only refinement

Second-stage shape solves run the native OLS/MDPDE fixed-point solve first.
`Native SUCCESS` results are returned unchanged. Only a non-success result
from a prepared second-stage shape fit enters the refinement eligibility
checks. Invalid, underdetermined, or rank-deficient data remain unqualified and
cannot be promoted. First-stage fitting, alpha training, group solvers, and
generic local/MDPDE calls remain native.

The production switch `--second-stage-failed-only-refinement false` selects
native-only behavior. Refinement uses the same Powell hybrid method in
`(beta0, log(beta1), log(variance))`, a central-difference Jacobian,
`xtol=1e-12`, and a 128-equation budget.

Acceptance requires finite positive variance, a valid Gaussian, full weighted
rank, a positive denominator, and fresh floor-inclusive scaled equations whose
maximum norm is at most `1e-8`. The candidate must agree with a continuation
from the native endpoint: that reference reaches residual `1e-10` within
10,000 total fixed-point updates, including native iterations; transformed
coordinates and weights agree within the existing `1e-6` thresholds; and
floor masks match.

An accepted result receives fresh weights and the existing covariance formula.
Rejection preserves the native parameters, weights, covariance, iteration
diagnostics, and status. Budget exhaustion, invalid numerical state, rank
deficiency, reference failure, and wrong-root or branch mismatch remain
explicitly unqualified.

Permanent owners are `MDPDERegression_test.cpp` and
`ProductionFitting_test.cpp`. They cover native-success bypass and unchanged
values, failed-only invocation, ineligible rank-deficient endpoints, wrong-root
branch rejection, fresh weights, covariance, budget exhaustion, invalid
variance/state, rank deficiency, and preservation of native status.
`MDPDETestSupport.*`, `EndpointRefinementTestSupport.*`, and
`SolverFailureCapture.*` provide numerical references, endpoint probes, and
captured solver-failure replay.

## Solver failure replay

Match the compiler, feature flags, worker count, verbosity, and quiet setting
when comparing runs. Debug verbosity affects scheduling, so a `-v 4`
diagnosis is not the same benchmark run as `-v 3`. Compare work, commits,
terminal state, parameters, and peeling. A reproducible failure case
demonstrates repeatability; it does not establish convergence.

Testing builds can save the first failure of each solver/status pair:

    RHBM_TEST_SOLVER_CAPTURE_DIR=build/production-fitting/failures \
      build/production-fitting/on/bin/RHBM-GEM potential_analysis [normal arguments]
    RHBM_TEST_REPLAY_DIR=build/production-fitting/failures \
      build/production-fitting/on/bin/RHBM-GEM-TEST \
      --gtest_filter=ProductionFittingTest.ReplaysCapturedSolverFailureWithoutSimulationInputs

These environment variables exist only in `BUILD_TESTING` binaries. Capture is
independent of the observer and never changes a result. The validation harness
must verify that requested artifacts were produced. No map or manifest is
needed to replay a captured subproblem. Retain executable/library hashes and
build/source provenance beside the fixtures; replay uses the same compiled
solver constants.

Fixture version 1 is whitespace-delimited with round-trip double precision:

- `shape 1`: alpha, threads, iteration limit, tolerance, and data-weight floor;
  dense X and y; expected status/variance; expected OLS and MDPDE beta vectors.
- `offset 1`: sparse X dimensions/nonzero count and row/column/value entries
  (including explicit zeros); y, anchor, and ridge vectors; expected status and
  offset vector. The compiled IRLS constants remain unchanged.
- Dense matrices and vectors contain row and column counts followed by
  row-major values. Replay requires exact parameters and matching terminal
  status.

The MDPDE fixtures remain checked in because permanent tests use them.
