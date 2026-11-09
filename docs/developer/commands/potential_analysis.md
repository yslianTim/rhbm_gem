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

| Request field | CLI flags |
| --- | --- |
| `job_count` | `-j, --jobs` |
| `verbosity` | `-v, --verbose` |
| `output_dir` | `-o, --folder` |
| `estimator` | `--estimator` |
| `database_path` | `-d, --database` |
| `model_file_path` | `-a, --model` |
| `map_file_path` | `-m, --map` |
| `selected_chain_ids` | `--select-chain` |
| `simulation_flag` | `--simulation` |
| `simulated_map_resolution` | `-r, --sim-resolution` |
| `saved_key_tag` | `-k, --save-key` |

The request catalog supplies the CLI and Python bindings. All surfaces default
to Joint Components (`PotentialEstimator::JOINT_COMPONENTS`), with
`two-stage` (`PotentialEstimator::TWO_STAGE`) available explicitly. Joint uses
the job count for FirstStage, capped at four workers. Two-Stage passes it
through as its fitting thread count.

## Execution Contract

`rhbm_gem::core::RunCommand(PotentialAnalysisRequest{...})` returns `CommandResult`.

Expected result contract:

- `result.succeeded == true` when execution completes and saves its result,
  including a Joint diagnostic result whose fitting did not qualify
- `result.succeeded == false` when normalization, validation, preflight, or execution stops the command
- `result.outcome` is `Completed` for ordinary command completion,
  `JointInitializationFailed` when the saved initialization is invalid,
  `JointNotConverged` when fitting did not qualify (including unavailable
  calibration after valid initialization), `JointDownstreamUnavailable` when
  Joint qualified but Halo reconstruction is unavailable, or `JointQualified`
  after the runtime gates passed
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

## Fixed fitting policies

[Production Fitting](../production-fitting.md) owns the shared per-atom
MDPDE and failed-only refinement contract.
[Two-Stage local fitting](../second-stage-local-fitting.md) owns its outer
iteration, signal and tail ranges, and final dependency polish.

## Tests to Update When Behavior Changes

- [`tests/core/command/CommandScenarios_test.cpp`](/tests/core/command/CommandScenarios_test.cpp)
- [`tests/core/contract/CommandExecutionContract_test.cpp`](/tests/core/contract/CommandExecutionContract_test.cpp)
- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)
- [`tests/integration/python_bindings_runtime_smoke.py`](/tests/integration/python_bindings_runtime_smoke.py)

## Joint component workflow

See [Runtime](../joint-component-runtime.md) for the Target-only fitting
problem and endpoint qualification, [Numerics](../joint-component-numerics.md)
for the block solver, and [Stage Integration](../joint-stage-integration.md)
for downstream data and result availability.

```cpp
rhbm_gem::core::PotentialAnalysisRequest request;
request.estimator = rhbm_gem::core::PotentialEstimator::JOINT_COMPONENTS;
request.model_file_path = "model.cif";
request.map_file_path = "map.mrc";
request.selected_chain_ids = {"A"};
request.database_path = "joint.sqlite"; // SQLite schema v21
request.saved_key_tag = "example";
auto completed = rhbm_gem::core::RunCommand(request);
```

`result_dump --printer joint` writes one Joint JSON and one atom CSV. The CSV
columns are `AtomID,ComponentID,A,tau,C,StateAvailable,SearchCompleted,StopReason,RuntimeConvergence,SelectionRole,Parameterization,ContributionGroupID,TargetRuntimeConvergence,PointQuality`.
`ContributionGroupID` remains an empty compatibility column in the current
CSV layout. No singleton contribution file is written. Available unqualified
Target states remain distinguishable from unavailable point estimates through
convergence and point-quality columns.

The [data I/O architecture](../architecture/dataobject-io-architecture.md)
defines which database and JSON versions the persistence entrypoints accept.
Joint workflow test groups are indexed in the [test guide](../../../tests/README.md#joint-component-regression).
