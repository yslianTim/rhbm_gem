# Joint second-stage integration acceptance

Baseline: `adf16940493bf46de4a467ade719d15ad7cb9b75`.
Each numbered stage is validated and committed before the next stage starts.

## 1. Estimator-neutral stage contracts

The new stage view distinguishes an absent estimate from a zero amplitude, keeps
method provenance separate from fitting stage, and does not expose Joint values
as OLS/MDPDE diagnostics. Initialization seeds are not published as final points.
Existing method-specific computation is unchanged.

Baseline: all eight related CTest groups passed (core estimator, sampler,
commands, joint component, data runtime, data schema, HRL and command integration).

After the change: the same eight groups, core contracts and the frozen
`joint_component_regression` passed (10/10). `tests_all` built successfully;
`git diff --check` passed. The availability test exercises seed/unavailable,
real zero amplitude, method-specific rejection and transient-state clearing.
