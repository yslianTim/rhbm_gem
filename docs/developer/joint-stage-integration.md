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

## Stage 6 — persistence and downstream analysis

- SQLite v18 stores neutral stage/source/role, uncertainty, evidence, posterior,
  paired peeling coverage and actual sample geometry separately from legacy
  method columns. Joint Second never occupies MDPDE/OLS columns. Snapshot and
  common Second identities/values are checked within the save transaction.
- v17 reads do not mutate the database. First write upgrades and saves in one
  transaction; failed validation rolls back schema and records. Legacy Joint
  snapshots expose recorded Second points without inventing peeling, uncertainty
  or posterior. Legacy samples explicitly lack geometry.
- Display, painters, Gaussian/position/outlier exports and feature construction
  consume common results and respect fitted target/halo roles. Curves identify
  Joint, charge coefficient and available uncertainty correctly. Dataset-specific
  Demo figures are explicitly skipped when their required named inputs are absent.
- UMAP retains its three features and standardization. Missing required features
  exclude rows with reasons; fewer than three valid rows fail. Ancillary missing
  fields remain empty, and a metadata JSON records estimator, features, peeling
  mode, normalization and exclusions.
- New tests cover precise round-trip of points/covariance/sample coordinates and
  unavailable states, snapshot/Second mismatch rollback, byte-identical v17 reads,
  successful and failed migration, legacy snapshot adaptation, and persisted group
  posterior identity. Command tests delete original map/model files before display
  and export, exercise all painter choices and missing components, and verify Joint
  UMAP target/coverage exclusions with six valid embedded rows.

Final acceptance (2026-09-23):

| Configuration/check | Result |
| --- | --- |
| EIGEN, ROOT/UMAP disabled, complete general CTest suite | 23/23 passed |
| ROOT and UMAP enabled, complete general CTest suite | 23/23 passed |
| SPQR: Joint component, frozen regression, estimator and HRL groups | 4/4 passed |
| Frozen Joint regression | Passed in both general suites and SPQR |
| `lint_repo` | Passed |
| Installed consumer smoke (`lint_install_smoke`) | Passed with ROOT/UMAP build |
| `git diff --check` | Passed |

The general suite excludes separately labelled offline/extended research runs;
those are not claimed as final acceptance evidence. Solver objective, support,
search budget and convergence thresholds were not changed.

Unavailable data remain intentional and inspectable: fixed-domain stencil gaps
or missing contributor states block individual peeling samples; old sample blobs
lack geometry; nonconvergence, boundary amplitudes, rank/df/variance failures block
uncertainty; insufficient eligible members or singular group covariance block
posterior. Historical unrecorded derivatives stay NotRun. No replacement local
fit, fabricated zero, inferred charge error bar or substitute posterior is used.
