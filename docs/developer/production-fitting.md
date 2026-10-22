# Production fitting

## Inputs and statistical contract

Production fitting consumes sampled map responses, selected atomic geometry,
explicit fitting options, and previously estimated state. It never loads the
simulation manifest, generated charges, truth widths, or scoring output. Atomic
number is not an amplitude truth initializer. `potential_analysis --simulation`
continues to skip normalization; `-r supplies model metadata` only.

The shape solve uses the production MDPDE fixed-point equations, and selected
offsets use joint-offset IRLS. Parameter accuracy alone does not establish
convergence. Inner solver failures retain their native status and are handled
by production qualification and recovery rules.

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
cannot be promoted. This failed-only refinement is fixed on for Second Stage.
First-stage fitting, Joint initialization, alpha training, group solvers, and
generic local/MDPDE calls retain their existing native behavior.

Potential-fitting workflows exclude hydrogen from targets, sample-neighbor
masks and contributors. `FitOptions` exposes only `estimator`, `thread_size`
and `quiet_mode`; sampling and Second Stage fixed policies are owned by their
implementations. Refinement uses the same Powell hybrid method in
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

`MDPDERegression_test.cpp` retains deterministic checks for MDPDE equations
and numerical boundaries. Its failed-only cases use generated inputs to verify
that native success stays unchanged, rank-deficient endpoints remain ineligible,
and an eligible generated failure can be refined. `ProductionFitting_test.cpp`
checks that every nominal endpoint is qualified and that nominal evaluation is
read-only and uses the unrestricted solve status. The failed-only cases generate
their input data in the test itself; MDPDE equation helpers are declared in
`MDPDETestSupport.hpp`.

## Joint result persistence

Joint component estimation uses the separate production contract in the
[Joint runtime guide](joint-component-runtime.md). Its saved data supports
SQLite schema v20, Stage Result JSON version 4, and Joint Result JSON schema
v19 only. These formats do not migrate or decode older formats.
