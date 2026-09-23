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

## 4. Grid-consistent post-fit peeling

Sampling and model subtraction now share the original nested tricubic arithmetic
and clamped stencil. The pure conversion subtracts fitted neighbors on the fixed
Joint rows, including halo, and preserves each raw sample. Per-sample responses
are optional, with explicit outside-domain or missing-contributor reasons.
Ratios require complete paired coverage in their interval and never clamp signed
results. Neighbor counts include contributors affecting interpolation nodes.

Validation: sampler, estimator, data runtime, commands, command integration,
frozen Joint regression and CLI smoke passed. The final Joint group rerun also
passed after correcting a new fixture to use binary-exact spacing: tiny nonzero
weights at decimal-grid nodes correctly retain conservative coverage rejection.
Tests cover signed C, zero A, negative responses, zero denominator, distant halo,
clamping, partial coverage, missing state, and unchanged raw/Second parameters.
The existing independent observation-stencil oracle passes after extraction.

## Stage 5 — parameter evidence and group inference

- Full-component raw `(A,C,log B)` Jacobian is column-scaled, reduced by tiled QR, and factored by SVD. Marginal blocks retain charge and neighbor coupling; residual variance uses `RSS/(N-3m)` and original dimensions set the rank threshold.
- Uncertainty requires runtime convergence, interior positive amplitudes, full rank and positive finite residual variance/degrees of freedom. Zero amplitude, missing state, nonconvergence, rank and variance failures retain Second points and explicit reasons.
- Group inference consumes only eligible parameter evidence and its covariance. Information-form WEB shares the original sample-domain core; no local MDPDE runs in this route. The correlation approximation is `block-diagonal-by-atom`. Posterior and unchanged Second are distinct; C remains descriptive with no inferred C uncertainty.
- Alpha training uses eligible members. Descriptive statistics retain all target points. Single-member and singular group covariance yield no substitute posterior. Results are written by atom identity, including exclusions in the middle of a group.
- Validation: seven affected CTest groups passed (Joint, estimator, data runtime, HRL, core commands, command integration and frozen Joint regression). Final Joint group passed after adding degeneracy tests. Added dense full-component covariance reference (including nuisance C), zero variance/rank/df/boundary/convergence gates, sample-vs-information WEB equivalence, and posterior independence from raw samples with sensitivity to changed evidence. Empty samples are sufficient for the new group route. `git diff --check` passed.
