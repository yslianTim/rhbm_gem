# Production failed-only refinement

Second-stage shape solves run the native OLS/MDPDE fixed-point solve first.
Production then applies the current failed-only refinement policy. Native
`SUCCESS` results are returned unchanged. A non-success result enters the
refinement eligibility checks; invalid, underdetermined or rank-deficient data
remain explicitly unqualified and cannot be promoted. This policy applies only
to prepared second-stage shape fits. First-stage fitting, alpha training, group
solvers and generic local/MDPDE calls remain native.

The existing production switch remains available through
`--second-stage-failed-only-refinement false` for native-only behavior. This PR
does not change that policy or any solver, root, convergence, objective, weight,
covariance or tolerance rule.

## Numerical contract

The refinement uses the same Powell hybrid method in
`(beta0, log(beta1), log(variance))`, central-difference Jacobian and `xtol=1e-12`.
The candidate has a 128-equation budget. Acceptance requires finite positive
variance, a valid Gaussian, full weighted rank, positive denominator, and fresh
floor-inclusive scaled equations with maximum norm at most `1e-8`.

Before acceptance, the candidate is compared with a continuation from the native
endpoint. That reference must reach residual `1e-10` within 10,000 total fixed
point updates, including native iterations. Transformed coordinates and weights
must agree within the existing `1e-6` thresholds, and floor masks must match.
The accepted result receives fresh weights and the existing covariance formula.
Rejected candidates preserve the native parameters, weights, covariance,
iteration diagnostics and status. Budget exhaustion, invalid numerical state,
rank deficiency, reference failure and wrong-root/branch mismatch remain
explicitly unqualified.

## Permanent owners

`MDPDERegression_test.cpp` and `ProductionFitting_test.cpp` cover native-success
bypass and unchanged values, failed-only invocation, ineligible rank-deficient
endpoints, wrong-root branch rejection, fresh weights, covariance, budget
exhaustion, invalid variance/state, rank deficiency and preserved native status.
Captured shape fixtures remain checked in because permanent tests use them.
`MDPDETestSupport.*`, `EndpointRefinementTestSupport.*`, and
`SolverFailureCapture.*` remain for numerical references, endpoint policy
probes, and captured failure replay. The MDPDE forward-comparison campaign and
the Fold-168 external regression have been retired; permanent tests do not need
external model, map, or manifest inputs.

The historical endpoint and failed-only campaign runners have been retired.
Their exact work totals and fold-wide outputs are not permanent correctness
contracts. The current refinement policy is owned by this document,
`production-fitting.md`, and the permanent C++ tests above. Historical artifact
provenance remains recorded in the retirement manifest and Git history.
