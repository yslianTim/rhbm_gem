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

## Request Surface

`PotentialAnalysisRequest` is a plain DTO that inherits shared fields from `CommandRequestBase`.

Shared fields:

- `job_count`
- `verbosity`
- `output_dir`

Command-specific fields:

- `database_path`
- `model_file_path`
- `map_file_path`
- `simulation_flag`
- `map_normalization_flag`
- `exclude_hydrogen`
- `simulated_map_resolution`
- `saved_key_tag`
- `training_alpha_flag`
- `asymmetry_flag`
- `sampling_method`
- `fit_range_min`
- `fit_range_max`
- `alpha_r`
- `alpha_g`

`sampling_method` uses the shared `SphereSamplingMethod` enum and is exposed
through the `--sampling-method` CLI flag.

## Execution Contract

`rhbm_gem::core::RunCommand(PotentialAnalysisRequest{...})` returns `CommandResult`.

Expected result contract:

- `result.succeeded == true` when the command completes
- `result.succeeded == false` when normalization, validation, preflight, or execution stops the command
- `result.issues` reports public validation diagnostics as option/message pairs

## Command Behavior

`PotentialAnalysisCommand::NormalizeAndValidateRequest()` handles request normalization and field validation:

- validates required model and map paths
- leaves map normalization as an execution-time choice: simulation requests always skip it,
  and non-simulation requests use `map_normalization_flag`
- coerces invalid scalar inputs back to command defaults when the command is designed to recover
- resets empty `saved_key_tag` values back to `"model"` and records a validation issue

`PotentialAnalysisCommand::ValidatePreparedRequest()` performs semantic checks after normalization:

- `--simulation` requires positive simulated resolution
- fit ranges must be ordered correctly

`PotentialAnalysisCommand::ExecuteImpl()`:

- builds command-owned data objects through `BuildDataObject()`
- optionally switches the model object into simulation mode
- optionally runs map normalization, then runs model preprocessing
- performs atom sampling, classification, fitting, and optional alpha training
- saves the prepared model through `SavePreparedModel()`

`PotentialAnalysisCommand::BuildDataObject()`:

- loads the model and map through `ReadModel(...)` and `ReadMap(...)`
- assigns command-local key tags before storing the loaded objects
- stores loaded objects directly in typed command-owned members
- wraps load failures with command-specific error context

`PotentialAnalysisCommand::SavePreparedModel()`:

- persists the prepared model through `DataRepository::SaveModel(...)`
- writes to the repository using `request.saved_key_tag` as the persisted key
- clears sampled local-potential distance/value buffers after persistence to keep runtime state lean

`CommandBase` also creates `output_dir` during filesystem preflight when needed.

## Tests to Update When Behavior Changes

- [`tests/core/command/CommandScenarios_test.cpp`](/tests/core/command/CommandScenarios_test.cpp)
- [`tests/core/contract/CommandExecutionContract_test.cpp`](/tests/core/contract/CommandExecutionContract_test.cpp)
- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)
- [`tests/integration/python_bindings_runtime_smoke.py`](/tests/integration/python_bindings_runtime_smoke.py)
