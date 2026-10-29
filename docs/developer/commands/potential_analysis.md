# Potential Analysis Command

## Registration and Implementation Files

The command membership entry lives in
the internal command catalog in
[`src/core/command/detail/CommandCatalog.hpp`](/src/core/command/detail/CommandCatalog.hpp).

Public request and entrypoint:

- [`include/rhbm_gem/core/CommandSystem.hpp`](/include/rhbm_gem/core/CommandSystem.hpp)

Internal schema and wiring:

- [`src/core/command/detail/CommandCatalog.hpp`](/src/core/command/detail/CommandCatalog.hpp)
- [`src/core/command/CommandSystem.cpp`](/src/core/command/CommandSystem.cpp)
- [`src/python/CommandSystemBindings.cpp`](/src/python/CommandSystemBindings.cpp)

Concrete implementation:

- [`src/core/command/PotentialAnalysisCommand.cpp`](/src/core/command/PotentialAnalysisCommand.cpp)
- [`src/core/MapSampler.cpp`](/src/core/MapSampler.cpp)

## Request Surface

`PotentialAnalysisRequest` is a plain DTO that inherits shared fields from `CommandRequestBase`.

Shared fields:

- `job_count`
- `verbosity`
- `output_dir`

Command-specific fields:

- `estimator` (default `PotentialEstimator::JOINT_COMPONENTS`; `TWO_STAGE`
  remains explicitly selectable)
- `database_path`
- `model_file_path`
- `map_file_path`
- `selected_chain_ids` (empty means all chains)
- `simulation_flag`
- `simulated_map_resolution`
- `saved_key_tag`

The request catalog generates both CLI options and Python bindings. C++ and
Python requests share the same estimator default. Potential analysis always
excludes hydrogen from its target, sampling mask and contributor domains. The
sampler owns the deterministic Fibonacci sampling default; generic sampling
APIs retain explicit `SphereSamplingMethod` selection. Second Stage refinement
and dependency-polish settings are fixed production policies, not request
fields. Alpha training is enabled, with initial local/group alpha values set
internally to `0.0`.

`FitOptions` contains only the estimator and execution settings:

```cpp
struct FitOptions
{
    PotentialEstimator estimator{ PotentialEstimator::JOINT_COMPONENTS };
    int thread_size{ 1 };
    bool quiet_mode{ false };
};
```

## Execution Contract

`rhbm_gem::core::RunCommand(PotentialAnalysisRequest{...})` returns `CommandResult`.

Expected result contract:

- `result.succeeded == true` when execution completes and saves its result,
  including a Joint diagnostic result whose fitting did not qualify
- `result.succeeded == false` when normalization, validation, preflight, or execution stops the command
- `result.outcome` is `Completed` for ordinary command completion,
  `JointInitializationFailed` when Joint setup failed, `JointNotConverged` when
  fitting did not qualify, `JointDownstreamUnavailable` when Joint qualified but
  Halo reconstruction is unavailable, or `JointQualified` after all runtime
  gates passed
- `result.issues` reports public validation diagnostics as option/message pairs

The CLI uses the following exit codes:

| Exit | Meaning |
| ---: | --- |
| 0 | Qualified Joint result or ordinary command success |
| 1 | Execution or persistence failure |
| 2 | Joint initialization failure was saved |
| 3 | Joint numerical non-convergence was saved |
| 4 | Joint qualified, but downstream Halo reconstruction was unavailable |

Thus an automation can distinguish saved diagnostics and unavailable downstream
outputs from a qualified Joint estimate without treating every saved file as a
successful fit.

## Command Behavior

The anonymous-namespace `NormalizeAndValidateRequest(...)` phase handles field validation:

- validates required model and map paths
- trims and deduplicates selected chain IDs and rejects empty IDs
- rejects invalid scalar inputs that this command cannot safely recover from
- rejects an empty `saved_key_tag`

`ValidatePreparedRequest(...)` performs semantic checks after normalization:

- `--simulation` requires positive simulated resolution
- loads the model once, verifies requested chain IDs and requires at least one
  non-hydrogen target before map preprocessing or fitting

`ExecutePreparedRequest(...)`:

- reuses the validated model and loads the map through `ReadMap(...)`
- optionally switches the model object into simulation mode
- skips map normalization for simulation requests and always attempts it for
  other requests
- computes Q-scores and fits the target selection without removing atoms from
  the model
- delegates sampling to `MapSampler`, whose default uses deterministic
  Fibonacci points
- excludes hydrogen from all potential-fitting targets, neighbor masks and
  contributor sets
- applies the fixed Second Stage refinement and final-polish policies
- persists the prepared model through `DataRepository::SaveModel(...)`
- writes to the repository using `request.saved_key_tag` as the persisted key
- clears sampled local-potential distance/value buffers after persistence to keep runtime state lean
- returns a numerical `CommandOutcome` after persistence so the API and CLI
  distinguish saved diagnostics from a qualified Joint result

`CommandRunner` creates `output_dir` during filesystem preflight when needed.

## Target chains and fixed policies

With no chain selection, all non-hydrogen atoms in the complete model are
targets, including side-chain atoms and chains that share an entity with other
chains. `--select-chain` narrows the target set while preserving the full
`ModelObject`. Chain IDs are compared case-sensitively using
`AtomObject::GetChainID()`, which for mmCIF input prefers `label_asym_id` and
falls back to `auth_asym_id`. Unknown IDs, empty CSV tokens and selections with
no non-hydrogen targets fail before Q-score preprocessing and fitting. Repeated
IDs are deduplicated; token whitespace is trimmed.

Direct C++ calls to `RunPotentialFittingWorkflow` remove hydrogen only from the
currently selected targets, preserving the caller's non-hydrogen subset. They
do not select previously unselected atoms, and fail clearly if no non-hydrogen
target remains. Hydrogen neighbors do not affect sample masks, Second Stage
neighbor counts, frozen background or peeling.

C++ callers can assign `request.selected_chain_ids = {"A", "B", "C"}`;
Python callers can assign `request.selected_chain_ids = ["A", "B", "C"]`.
An empty vector or list selects every non-hydrogen atom in the model.

For Joint Components, selected chains define targets and their voxel
observations. Other non-hydrogen atoms remain eligible as halo contributors
when their support overlaps those observations. Chain selection does not crop
the contributor model or change the target/halo distinction.
Production Joint fitting uses bounded Active-set LM.

Simulation maps skip normalization. Every non-simulation map is normalized by
its standard deviation when that operation succeeds; no mean is subtracted. A
zero-SD map records `requested: true`, `applied: false`, and `divisor: 1`.
Map Simulation is a separate workflow; `MapSimulationRequest::exclude_hydrogen`
continues to control whether hydrogen contributes to a simulated map.

```sh
# Default: Joint Components for all non-hydrogen atoms in all chains
rhbm_gem potential_analysis \
    --model model.cif \
    --map map.mrc

# One target chain
rhbm_gem potential_analysis \
    --estimator joint-components \
    --model model.cif \
    --map map.mrc \
    --select-chain A

# Several target chains
rhbm_gem potential_analysis \
    --estimator joint-components \
    --model model.cif \
    --map map.mrc \
    --select-chain A,B,C

# Explicitly select the Two-Stage estimator
rhbm_gem potential_analysis \
    --estimator two-stage \
    --model model.cif \
    --map map.mrc
```

## Failed-only endpoint refinement

Second Stage always enables failed-only refinement for its shape refits. Native
successes remain unchanged; only eligible failed endpoints can be refined.
First Stage, Joint initialization, alpha training, group estimation and generic
`EstimateLocalGaussian` / `EstimateBetaMDPDE` entrypoints retain their existing
native behavior. Final dependency polish is enabled with the production
10-round limit and runs only after the existing eligible convergence stop.

A native successful solve is used unchanged, without a fresh-equation evaluation.
A failed solve may be refined from its own endpoint. Acceptance requires fresh
scaled equations at `1e-8`, a valid positive-variance Gaussian, weighted rank and
denominator checks, and agreement with a fixed-point continuation from the
original endpoint. The reference must reach `1e-10`; transformed coordinates and
weights must agree within `1e-6`, with identical weight-floor activation.

The endpoint-refinement candidate budget is 128 equation evaluations including
verification. Reference continuation is separate work, capped at 10,000 total
fixed-point updates including the native iterations. These constants are not
public tuning options. Exact-fit variance boundaries are not promoted to
successful solves. Rejection preserves the original numerical result; the
reference is never used as a fallback.

`RHBMBetaEstimateResult.status` and `diagnostics` retain native solver history.
Use `Qualification()` to distinguish `NativeSuccess`, `RefinedSuccess` and
`Unqualified`. The optional `refinement` diagnostics record acceptance, rejection
reason, residuals, branch differences and candidate/reference work. Audit shape
records retain effective `status` for existing consumers and additionally
expose `native_status`, `qualification` and optional refinement/recovery
evidence.

Recovery continues to require qualified inner endpoints. Refinement does not
relax outer convergence, offset IRLS or the best-objective bound, and does not
itself establish convergence. See the [production fitting contract](../production-fitting.md).

## Internal fitting ranges

Fitting ranges are internal constants in `src/core/detail/gaussian_fit/FittingRanges.hpp`.
Signal fitting and alpha training use `[0, 1.0]` Å; the second-stage tail
objective uses `[1.2, 2.0]` Å. Both boundaries are inclusive. The two ranges
are independent and may overlap or leave a gap when the constants are changed.
An overlapping sample contributes to both objective terms, each with its own
scale and sample-count normalization. Samples in neither region contribute to
neither residual objective term.

Offset fitting still uses all raw samples, and joint polish keeps its existing
sample sources. The gap is therefore not excluded from every estimation step;
tail is an objective constraint region, not a strictly held-out validation set.
Neither PotentialAnalysis nor RHBMTest accepts `--fit-min` or `--fit-max`.

## Tests to Update When Behavior Changes

- [`tests/core/command/CommandScenarios_test.cpp`](/tests/core/command/CommandScenarios_test.cpp)
- [`tests/core/contract/CommandExecutionContract_test.cpp`](/tests/core/contract/CommandExecutionContract_test.cpp)
- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)
- [`tests/integration/python_bindings_runtime_smoke.py`](/tests/integration/python_bindings_runtime_smoke.py)

## Joint component production

The map-aware `RunPotentialFittingWorkflow(map, model, options)` dispatches
`JOINT_COMPONENTS` to the production target-only Pure Halo-Support Exclusion
workflow. Selected non-hydrogen atoms are Targets. Halo support is used only to
exclude fitting rows within 2.5 Å; Halos are not free Joint parameters.

Qualified FirstStage Target widths initialize the fit and determine the global
width bound `Bmax = 4 × P95`. Calibration needs at least 20 qualified Target
donors and 95% Target coverage. If calibration is unavailable, the saved
result retains that status and does not retry with an unbounded search.

Production uses FixedNeighbor core size 12, forward serial Gauss-Seidel, at
most 200 sweeps, and no more than one trusted accepted update per block visit.
Bounded Active-set LM updates widths. Constrained A/C subproblems use Eigen
SparseQR; A is nonnegative and C is signed. Global convergence, complete-sweep
stationarity, endpoint verification, Independent Replay, Local Conditioning,
and Parameter Stability remain explicit evidence. See the
[Joint runtime guide](../joint-component-runtime.md).

Qualified target A/B/C points provide geometry-local kNN-4 Halo Median
Reconstruction. Original-domain PostFitPeeling uses available Target and
reconstructed Halo points, followed by Post-Joint Local MDPDE and Group
Fitting for qualified Targets. These downstream results retain their own
availability and qualification states.

Initialization failure, unavailable width calibration, search failure,
unqualified endpoint, missing component state, absent Halos, unavailable Halo
reconstruction, and incomplete downstream results are valid saved outcomes.
`CommandResult.succeeded` means command execution and persistence completed; it
does not mean the Joint estimate qualified. Invalid input or a persistence
failure remains a command failure.

```cpp
rhbm_gem::core::PotentialAnalysisRequest request;
request.estimator = rhbm_gem::core::PotentialEstimator::JOINT_COMPONENTS;
request.model_file_path = "model.cif";
request.map_file_path = "map.mrc";
request.selected_chain_ids = {"A"};
request.database_path = "joint.sqlite"; // SQLite schema v20
request.saved_key_tag = "example";
auto completed = rhbm_gem::core::RunCommand(request);
```

`result_dump --printer joint` writes one Joint JSON and one atom CSV. The CSV
columns are `AtomID,ComponentID,A,B,C,StateAvailable,SearchCompleted,StopReason,RuntimeConvergence,SelectionRole,Parameterization,ContributionGroupID,TargetRuntimeConvergence,PointQuality`.
`ContributionGroupID` remains an empty compatibility column in the current
CSV layout. No singleton contribution file is written. Available unqualified
Target states remain distinguishable from unavailable point estimates through
convergence and point-quality columns.

## Joint persistence formats

| Data | Accepted format |
| --- | --- |
| SQLite database | Schema v20 only |
| Stage Result JSON | Version 4 only |
| Joint Result JSON | Schema v21 only |

An empty database is initialized as v20. An existing database must have the
exact v20 schema. Stage JSON accepts only version 4; Joint JSON always writes
schema 21 and reads only schema 21. No migrations, upgrades, historical
fallbacks, or automatic repairs are performed. Unsupported inputs are rejected
without modifying the database. See the
[data I/O architecture](../architecture/dataobject-io-architecture.md).

## Tests to update when behavior changes

- [`tests/data/DataObjectPersistence_test.cpp`](/tests/data/DataObjectPersistence_test.cpp)
- [`tests/data/DataObjectSchemaLifecycle_test.cpp`](/tests/data/DataObjectSchemaLifecycle_test.cpp)
- [`tests/data/DataObjectSchemaValidation_test.cpp`](/tests/data/DataObjectSchemaValidation_test.cpp)
- [`tests/core/joint_component`](/tests/core/joint_component)
- [`tests/integration/joint_workflow_cli_smoke.py`](/tests/integration/joint_workflow_cli_smoke.py)
