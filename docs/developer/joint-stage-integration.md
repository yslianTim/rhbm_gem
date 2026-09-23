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

## 2. Shared sampling and First initialization

The map-aware fitting workflow prepares one immutable Joint problem, initializes
contributors using explicit identities, and calls `FitJointComponents` directly.
The standalone wrapper delegates to the same First implementation on its private
copy and preserves its historical writeback contract. Two-stage keeps its existing
sampling and numerical path. The model-only workflow rejects Joint requests.

Validation: seven focused CTest groups passed, including frozen Joint regression,
command/CLI smoke, sampler and estimator regression. A final targeted rerun after
adding First target/halo provenance passed both Joint and CLI smoke. The new
instrumented test verifies one raw sampling and one formal First per contributor,
unchanged selections, and exact equality with a direct fit of the same problem/B0.

## 3. Joint endpoints and Second summary

A data-only adapter maps component-local A/C/B through contributor identities to
neutral Second estimates. It retains target/halo roles, component convergence and
one source ID per workflow; missing states replace any previous point with an
explicit unavailable reason. Summary reads the common point interface before
group fitting, excludes halo and labels C and between-atom dispersion correctly.

Validation: all six related tests/groups passed (Joint, estimator, commands,
command integration, frozen regression and CLI smoke). The new adapter test
permutes identities/mappings and verifies target-only summary and stale-state
removal. No solver settings or numeric algorithms changed.
