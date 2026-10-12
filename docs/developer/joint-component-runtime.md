# Joint component estimator and regression workflow

Include `<rhbm_gem/core/JointComponentEstimator.hpp>` and link
`RHBM_GEM::rhbm_gem`. The installed API needs no testing or offline audit target.

```cpp
model.SelectAllAtoms();
auto result = rhbm_gem::core::EstimateJointComponents(map, model);
auto problem = rhbm_gem::core::BuildJointProblem(map, model);
auto repeated = rhbm_gem::core::FitJointComponents(problem, initial_b);
```

## Runtime contract

`JointProblem` owns an immutable snapshot of observations, identities and
squared-distance memberships. V1 uses equal-weight least squares over its fixed
observation domain, with A nonnegative, C signed, and log-B widths. The public
objective is `0.5 * ||prediction - observations||² / ObservationScale()²`.

Production fitting uses the FixedNeighbor route described below. Each accepted
local Profile LM update is replayed against the global state and passes the
roundoff-safe monotone acceptance rule before it is kept. At the end of each
complete sweep, global stationarity checks identify a candidate endpoint. A
confirmed candidate is independently certified before search can complete;
reference coefficients never replace the accepted search coefficients. A valid
endpoint that fails trust keeps its assessment and may continue through the
same block-coordinate updates. While all stationarity requirements remain
confirmed, another certification is attempted at most once every 10 sweeps. A
failed KKT, width-gradient or eta-change confirmation resets retry eligibility;
after the full stationarity conditions are confirmed again, certification may
be retried immediately. An invalid primary endpoint stops immediately. At
budget exhaustion, if the final sweep has no certification, one bounded
assessment certifies the final accepted state. It does not change accepted
parameters, add a sweep or mark the search complete. If the final sweep already
has a certification, that assessment is not repeated. A passing final trust
check may preserve a usable state while the search remains incomplete; a failed
check preserves its evidence and leaves state unavailable. Whenever a final
assessment is completed, its provenance is `endpoint_certification_sweep ==
search_sweeps`; zero-sweep searches and earlier unrecoverable failures do not
invent certification provenance. Search completion, certification,
state availability and runtime convergence remain separate outcomes.

`AssessProfileSearchForDiagnostics()` is a test-support diagnostic for Profile
search history and trusted historical states. Production components use
`CertifyFixedNeighborEndpoint()` for independent endpoint certification. Both
paths use the shared evaluation, replay and assessment primitives.

### Production FixedNeighbor route

The sole production route for `JOINT_COMPONENTS` is `FixedNeighbor`. The
estimator, workflow API, CLI and Python API do not expose a Joint route
selector; `--estimator joint-components` is sufficient.

The production FixedNeighbor policy is intentionally small and fixed: disjoint
mutable cores of 12 atoms, at most 200 forward Gauss-Seidel sweeps, and one
trusted accepted local Profile LM update per block visit. Outside-core
neighbors remain fixed while a block is solved. Each accepted block is replayed
globally, checked for monotone roundoff-safe acceptance, and contributes to the
existing global A/C KKT (`1e-10`), width-gradient (`1e-12`) and eta-confirmation
(`1e-10`) criteria. Meeting these checks makes the endpoint a candidate; the
existing independent trust assessment must also pass for successful stopping.
Valid trust failures remain eligible for bounded retry; invalid primary
assessments stop, and an untrusted endpoint at the sweep limit reports budget
exhaustion without a usable state. The at-most-one accepted-update contract is
intrinsic; diagnostic controls such as reverse order and trajectory capture are
not production configuration.

The local numerical primitive is named **Profile LM** (Levenberg–Marquardt
optimization of the profiled width objective). `SearchProfile()` optimizes
`η = log(B)` with `InstrumentedLM<Profile>`. At every profile evaluation,
`EvaluateProfile()` conditionally solves the A/C coefficients with the existing
constrained `SolveLinear()` call, which uses EIGEN SparseQR:

```text
β*(η) = arg min 1/2 ||X(η)β - y_eff||²
        subject to A_i >= 0 and C_i unrestricted
```

The full profile derivative drives the LM width update. `AssessProfile()` keeps
its assessment responsibility. Each block visit accepts at most one trusted
update and replays the candidate against the global state.

The production search returns an internal `FixedNeighborComponentResult` with
the search state and independent endpoint certification. Shared orchestration
turns that into the certified public `JointComponentResult`, records component
cost/status and progress, then releases the larger search result. The batch
retains those public component results, aggregate costs and only the optional
single-component assessment that can be reused by exact-match FullABC assembly;
it does not hold a second vector of solver results. The route returns through
`JointFitResult` and the standard post-processing and persistence path. It uses
the defined observation domain, contributor/halo closure, parameterization, rank
semantics, parent scale, objective and `RuntimeConvergence()` semantics. Its
progress event includes sweep, block, accepted-local-update, KKT, width-gradient,
eta-change and objective fields.

The production FixedNeighbor implementation uses `MakeComponentProblemView()`
to prepare an immutable indexed view over the shared root snapshot. The view
keeps explicit parent/local atom and row mappings, identities and rank shape;
support rows and squared distances remain in the root snapshot. Structural
partitioning then creates one prepared block view per block, including its core
and affected-row mappings, indexed `Domain` and context template. A sweep still
rebuilds the state-dependent `eta`, `beta`, effective response and numerical
design; it does not rebuild those structural views. No numeric factor, active
set or coefficient solution is carried across sweeps.

`JointProblem::ParameterLayout()` is the parameterization authority.
`OrchestrateFixedNeighborComponents()` projects that layout onto each component
and owns component traversal, width-seed validation, search, independent
certification and component progress. The presence of nuisance groups selects
observable assembly; callers do not pass a separate mode flag. FullABC results
keep `component.layout` absent, while observable results retain each projected
layout for reconstruction. The orchestration has no assembly callback. FullABC
assembly consumes the certified public component states and performs global
reconstruction and consistency assessment.
Observable assembly runs after component traversal: it reconstructs each
component's full-atom contributions, analytically profiles each singleton-halo
nuisance amplitude as `lambda = observation - current prediction`, updates the
available-row mask, and checks full prediction/objective consistency. Both paths
emit the established assembly progress phases once assembly begins. Target
evidence is added once after assembly finishes.

`FitFixedNeighborComponentsImpl()` owns the Eigen thread scope, solver
provenance and FullABC/observable dispatch; `FitFixedNeighborComponents()` adds
target evidence once after implementation and assembly finish. Public
`FitJointComponents()` and `EstimateJointComponents()` keep their existing
contracts and do not expose assembly details to callers.

Production search retains aggregate counters, final sweep state and the latest
endpoint assessment needed by the component route. Persisted component counters
record actual sweeps, first and confirmed stationarity, and the last certification
sweep. A completed final assessment at budget exhaustion records the same sweep
as `search_sweeps`; no certification sweep is recorded when no assessment ran.
Detailed block records, local trial trajectories, reverse order and local
certification remain opt-in experiment/test diagnostics. The preparation and
mapping counters are reported as implementation attribution; they are not
convergence criteria.

The whole-problem `RunFixedNeighborDiagnostic()` facade is compiled into test
support and called only by component tests and the fixed-neighbor experiment.
The library retains lower-level telemetry and observer hooks used by that
support; they are implementation instrumentation. Production receives only the
qualified core-size, sweep-count and block-order controls through
`FixedNeighborSearchPolicy`, plus the progress observer.

Within one immutable problem, a full single component can share its assessment
with assembly only when observations, structural support, identities, state,
scale and linear/rank policy match exactly. This avoids repeating the reference
evaluation and assessment for the same returned state. Offline audit directions
do not affect runtime assessment. Constant rows or different contexts require
separate assessments. This is scoped result reuse, not a persistent cache. The
test-support Profile-search diagnostic assesses a selected historical trial at
its saved coefficients and widths. Assessment work and timing are reported
separately from search work. Sparse factorization details can change
floating-point trajectories within the fixed endpoint contract.

### Selection-domain contract

The Map/Model builder uses selected non-hydrogen atoms as targets S. Their
2.5 Angstrom spheres fix the unique observation rows V. Selected targets remain
in the target mapping even if clipping leaves one with no observed support; no
rows are invented for it. The builder intersects every eligible non-hydrogen
atom's support with V, retaining every intersecting unselected atom as a halo
contributor. Contributors follow model order. Halo never adds rows or
recursively adds its neighbors. Atom centers outside the map are eligible if
their support intersects V. Empty non-hydrogen selection is rejected. The
existing `sphere-fma-v1` arithmetic and negative/zero observations are retained.
Structural partition runs after closure; a shared or bridge halo has one
parameter vector and can connect separate target regions.

`JointProblemInput.selection_domain` is optional for hand-built frozen inputs
(null means selection not recorded). The Map/Model builder always records
`fixed-selected-voxel-closure-v1`, ordered unique `target_indices`, both radii
(2.5 Angstrom), and `all-non-hydrogen` contributor eligibility. Halo indices are
the complement within the problem's atom identities. This establishes completeness
only for the supplied catalogue and truncated support model. The saved outcome
retains this domain and atom ordering; stage and CSV metadata report target,
halo, or not-recorded roles. No role is inferred when selection was not recorded.

`EstimateJointComponents` and the production workflow share
`PrepareJointWidthInitialization()`: it samples only FullABC contributors with
the deterministic Fibonacci sampler, then trains local alpha and computes
first-stage MDPDE width seeds for all required contributors. The workflow first
initializes the selected model analysis and applies the prepared results there.
`EstimateJointComponents` runs the same preparation on a full model copy, then
copies raw samples and first-stage results back only for successfully initialized
targets. Only B initializes joint fitting. Its initialization data scope is
`contributor-local-sampling-may-read-outside-target-domain`: interpolation and
sampling can read beyond V, while the joint objective cannot. Supplied widths
use `caller-provided-widths`. Global map normalization is also separate from V.
The standalone copy/writeback contract preserves halo history, original
selection, second-stage and group results. Initialization records retain
atom-specific reasons and nonfinite diagnostics as unavailable, not replacement
coefficients.

An incorrectly sized initial-width vector remains a whole-input failure.
Otherwise invalid widths/initialization exceptions block only their structural
component, which retains mappings, stop reason and unavailable evidence; other
components continue. `initialization.valid` describes all widths, not whether
any component can run. Missing component states still prevent global assembly.
Peeling results, truth, prior fits and certificates do not initialize the
estimator.
`costs.construction_seconds` measures the Map/Model builder; initialization timing
includes the model copy and setup. Direct frozen-input fitting has zero builder
and initializer costs.

A/C order is `[A0, C0, A1, C1, ...]`; widths follow atom identities. Public
component, assembled-state and result objectives are
`0.5 * ||prediction - observations||² / ObservationScale()²`. Components use the
same parent scale `max(1, ||parent observations||₂)`; constant rows contribute
once to the global objective.
Width gradients are with respect to log-B. Internal raw objectives retain
half-RSS semantics. This public objective correction changes returned
scale, not the optimizer or its acceptance thresholds.

Component `state` is the actual last trusted state used by assembly. Search
completion, usable state and evidence status are separate. A missing component
prevents complete prediction/objective; the row mask retains available rows,
including constant rows. Assembly does not reprofile or zero-fill missing states.
Two-step derivative, Richardson, multiprecision, boundary and regular certificates
are `NotRun` in runtime results. Offline tools produce separate evidence without
mutating the runtime result. Operational runtime checks are not regular
certification.

`component.RuntimeConvergence()` and `result.RuntimeConvergence()` derive their
status from the returned evidence; they do not run numerical work or store another
success flag. Component convergence requires an actual trusted state and all five
checks: inner/reference agreement, A/C KKT (including feasibility), log-B gradient,
local correction and numerical identifiability. Existing thresholds remain
`1e-10` for coefficient agreement, KKT and local correction, and `1e-12` for both
primary/reference width gradients. Full profile derivatives, including the
nonzero-residual correction, remain in the numerical core.

Global convergence additionally requires every component to converge, complete
assembly, global evidence and independently reprofiled assembly consistency.
Missing required states or assembly give `Unavailable`. Otherwise, required
checks combine in this order: any `Failed`, any missing/`Unavailable`, any
`NotRun`, then `Passed`. Only evidence from the required scope participates.
Search completion is independent: a budget/untrusted-trial stop can leave a
trusted endpoint whose runtime checks pass; its stop reason is still reported.

Runtime convergence does not imply an offline regular certificate.

## Immutable storage and component views

The input snapshot owns observations, memberships and identities once. Numerical
observations are read-only Eigen maps; contexts and component support views retain
shared ownership, including after the original problem handle is destroyed.
Partition reverse mappings are shared once across components. Local support is
mapped on access without copying squared distances. Contiguous observation ranges
are borrowed; noncontiguous ranges use one scoped gather buffer for search and
assessment. Census and offline bundles materialize legacy mappings only at their
serialization boundary.

Component, profiled-component and prepared-block views each keep their local
atom selection and precomputed row-to-local mapping over the parent view. Support
iteration applies those mappings directly in constant work per view layer, and
FixedNeighbor reuses prepared block mappings across sweeps.

Uncertainty preparation uses the same `MakeComponentProblemView()` for ordinary
component covariance and observable target-quotient covariance. It selects the
needed observation values for numerical work but reads support membership and
identity through that indexed view; it does not construct local support ownership
or synthetic `PartitionMappings`. The target-quotient and full-component
covariance calculations remain separate because their rank, identifiability and
degrees-of-freedom contracts differ.

The Map builder collects and sorts relevant voxel indices instead of allocating
an index array for the entire map. Voxel ordering and support arithmetic remain
unchanged. Current Joint search, solve, workflow, postprocessing and command
measurements are run through the single workflow in the
[benchmark guide](joint-benchmark.md).

## Internal numerical path

The current numerical path is:

```text
FixedNeighbor
  -> local Profile LM over η = log(B)
  -> conditional A/C profile solve / EIGEN SparseQR
  -> full profile derivative and LM width update
  -> global replay, endpoint certification, RuntimeConvergence
```

These stages preserve the original parent rank context and full residual norm.
They do not cut component edges or use a size-based route. The
[sparse-backend guide](joint-component-sparse-backend.md) documents the fixed
EIGEN implementation.

Permanent tests check the EIGEN path and numerical contract: trusted-state
availability and required runtime checks must agree, missing states remain
explicit, and the production KKT, width-gradient, eta-confirmation, replay,
endpoint-certification and `RuntimeConvergence` thresholds remain unchanged.

## Routine regression

```sh
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target tests_all -j4
ctest --test-dir build -L joint:runtime --output-on-failure
```

The runtime test lane protects the estimator behavior and public API. The
integration regression runner and its fixture catalog are documented in
[tests/README.md](../../tests/README.md). Offline certification runs are a
separate lane and do not change the runtime result.

Runtime search/assembly reports use schema version 2: `runtime_convergence`,
`runtime_checks` and `runtime_failure` replace the mixed `joint_qualified` and
`qualification_*` fields. Search-endpoint identity is retained in diagnostics.
Missing evidence never means passing. Legacy `joint_qualified` is produced only
by offline assessment and retains its derivative-dependent meaning.

`joint_component_runtime` provides commands for public API runs, physical-input
smokes, and fixture cases. The Python runner exposes run, summary, comparison,
regression and physical-smoke operations. Current entry points are maintained in
`tests/integration/joint_component_runtime.py`. It is a validation harness,
not a benchmark or experiment entry point. Current experiments and validation
tools are listed in [joint-experiments.md](joint-experiments.md).

## Extended and offline checks

Both options default to `OFF` and require `BUILD_TESTING=ON`:

```sh
cmake -S . -B build -DBUILD_TESTING=ON \
  -DRHBM_GEM_ENABLE_JOINT_EXTENDED_TESTS=ON \
  -DRHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON
cmake --build build --target tests_all -j4
ctest --test-dir build -L 'joint:extended|joint:offline' --output-on-failure
```

Extended and offline lanes are opt-in and stay separate from the default
runtime regression.

Local-audit preparation and precision code are linked only into
`joint_component_audit`, `joint_offline_diagnostic` and `joint_offline_tests`,
never the library or ordinary test executable. Kernel or Jacobian changes
require the applicable derivative
audits; constraint changes also require boundary controls. Sparse-factorization
changes require rank and precision controls. Recorded campaign results and their limits
belong in the [canonical evidence index](joint-component-evidence.md).

```sh
python3 tests/integration/joint_component_audit.py \
  --executable build/bin/joint_component_audit --work-dir build/joint-audits
build/bin/joint_component_audit local-bundle INPUT_JSON OUTPUT_DIRECTORY
```

A standalone bundle carries immutable parent observations/context, component
memberships and initial widths. Local directions come from the actual component
endpoint; sibling/global endpoints are not inputs. Modified parent context or
invalid widths are rejected. Parent normalization, actual-state preservation,
local weak-direction evidence and independently reprofiled consistency remain
required; a local certificate cannot promote a missing or failed component.

Current sparse-factorization behavior and benchmark commands are documented by
the sparse-backend and benchmark guides.

## Saved production outcomes

`potential_analysis --estimator joint-components` stores a `JointAnalysisResult`
on the model and saves it through `DataRepository`. The default estimator remains
two-stage. See the [command contract](commands/potential_analysis.md#joint-component-opt-in).

Direct `EstimateJointComponents` callers keep the return-only joint contract.
To explicitly retain a result, use `CaptureJointAnalysisResult(fit, metadata)` and
`model.EditAnalysis().SetJointResult(...)`, then `DataRepository::SaveModel`.
Read it through `model.GetAnalysisView().GetJointResult()`. Include
`<rhbm_gem/data/io/JointAnalysisFileIO.hpp>` to export a saved result directly with
`WriteJointAnalysisResult(result, json_path, csv_path)`.

The snapshot, prediction vector and numerical workspaces are not captured. Only
result values and identity/mapping information survive. Global and component
runtime-convergence statuses are captured from the existing runtime methods;
storage validates the document structure but neither recomputes convergence nor
runs an audit. These saved documents are not standalone numerical audit bundles.
Model copying retains the outcome; analysis `Clear()` and `ClearJointResult()`
remove it, while `ClearTransientFitStates()` does not.

### Saved provenance and units

Production outcome JSON writes schema 5 and retains schema 3 and 4 reading
(distinct from fixture and offline report schemas). The [observable-halo contract](joint-observable-halo.md)
describes parameter layouts, nuisance contributions and seed provenance. See the
[metadata contract](commands/potential_analysis.md#provenance-and-map-units-joint-json-schemas-3-4-and-5).
`JointAnalysisMetadata` contains optional `JointMapNormalization`, input SHA-256
values, `JointSoftwareProvenance` and optional `JointSolverProvenance`.
`JointSolverProvenance` records the current FixedNeighbor estimator and
preserves solver and local-work fields found in saved records.
Current writers record the fixed EIGEN sparse backend and the FixedNeighbor policy
fields; current FixedNeighbor records use `fixed-neighbor-production-v2` with
intrinsic OneAccepted behavior and the `joint-solver-provenance-v3` contract.
Saved `joint-solver-provenance-v2` records with LegacyCompact and OperatorPcg
fields remain readable.
`CaptureJointAnalysisResult`
records the current library's version/source/configuration/build identity. The
decoder and exporter preserve saved identity without recalculation. Direct
in-memory callers leave unknown input hashes and normalization null; capture
does not open paths.

`MapObject::MapValueArrayNormalization()` returns `std::optional<double>`: the
actual divisor if applied, or nullopt when zero SD skips normalization. Existing
callers may ignore the return. A caller retaining provenance records its requested
flag, whether the returned divisor exists, and `divisor.value_or(1)`; a command
that skips normalization records `applied=false, divisor=1`.

B and geometry use Angstrom, A uses fitted-map-unit*Angstrom^3, and C uses
fitted-map-unit*Angstrom. Multiply A/C by the recorded divisor to recover input-map
scale. These are kernel coefficient units, not an assertion of physical map
calibration; C is not a constant background offset. The objective's parent
observation scale is separate. Unknown normalization cannot support conversion.

The command's default estimator remains two-stage. Joint LS is opt-in and its
runtime checks do not establish an offline regular certificate.

Permanent [partial-selection tests](../../tests/core/joint_component/PartialSelection_test.cpp)
cover the fixed observation domain, nonrecursive halo closure, bridge
parameterization, unavailable target states, and the numerical policy.

## Joint validation tiers

Build with at most four jobs, then run the bounded gates from the build tree:

```sh
cmake --build build/debug-tests -j4
ctest --test-dir build/debug-tests -L 'joint:contract'
ctest --test-dir build/debug-tests -L 'joint:workflow'
ctest --test-dir build/debug-tests -L 'joint:numerical'
```

Contract covers deterministic component, linear algebra, and small assessment
checks, including all 14 `JointComponentChecksTest` cases. Workflow covers
runtime, persistence, CLI, and smoke paths. Numerical covers bounded
FixedNeighbor, stationarity, order, prepared-block, KKT, and derivative/reference
qualification. Each gate has a CTest timeout and sets `OMP_NUM_THREADS=1`.

Ordinary Joint cleanup or refactoring accepts passing all three gates. A
numerical solver change also runs each relevant scalability case. Release or
performance qualification runs Tier 4 for its campaign scope:

```sh
ctest --test-dir build/debug-tests -L 'joint:scalability'
ctest --test-dir build/debug-tests -L 'joint:offline'
```

Tier 4 gives `joint_component_gram_scalability`, the two 256-atom FixedB cases,
and optional offline audits their own CTest entries with a one-hour timeout.
These entries are outside ordinary acceptance. Offline audit entries are
registered when configured with
`-DRHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON`.

Report test states as `passed`, `failed`, `timeout`, `not-run`, or
`not-completed`. A manually stopped run without a test result is
`not-completed`; it is neither a pass nor a failure.
