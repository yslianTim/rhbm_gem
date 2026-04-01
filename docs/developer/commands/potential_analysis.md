# Potential Analysis Command

## Registration and Implementation Files

The command membership entry lives in
[`src/core/command/CommandManifest.def`](/src/core/command/CommandManifest.def).

Public request and entrypoint:

- [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp)

Internal schema and wiring:

- [`src/core/command/detail/CommandRequestSchema.hpp`](/src/core/command/detail/CommandRequestSchema.hpp)
- [`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp)
- [`src/core/command/CommandCli.cpp`](/src/core/command/CommandCli.cpp)
- [`src/python/CommandApiBindings.cpp`](/src/python/CommandApiBindings.cpp)

Concrete implementation:

- [`src/core/command/PotentialAnalysisCommand.hpp`](/src/core/command/PotentialAnalysisCommand.hpp)
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
- `simulated_map_resolution`
- `saved_key_tag`
- `training_report_dir`
- `training_alpha_flag`
- `asymmetry_flag`
- `sampling_size`
- `sampling_range_min`
- `sampling_range_max`
- `sampling_height`
- `fit_range_min`
- `fit_range_max`
- `alpha_r`
- `alpha_g`

## Execution Contract

`RunPotentialAnalysis(...)` returns `CommandResult`.

Expected result contract:

- `result.succeeded == true` when the command completes
- `result.succeeded == false` when normalization, validation, preflight, or execution stops the command
- `result.issues` reports public validation diagnostics as option/message pairs

## Command Behavior

`PotentialAnalysisCommand::NormalizeRequest()` handles parse-phase normalization and validation:

- validates required model and map paths
- coerces invalid scalar inputs back to command defaults when the command is designed to recover
- resets empty `saved_key_tag` values back to `"model"` and records a parse-phase issue

`PotentialAnalysisCommand::ValidateOptions()` performs prepare-phase semantic checks:

- `--simulation` requires positive simulated resolution
- sampling and fit ranges must be ordered correctly

`PotentialAnalysisCommand::ExecuteImpl()`:

- builds command-owned data objects through `BuildDataObject()`
- optionally switches the model object into simulation mode
- runs map preprocessing and model preprocessing
- performs atom sampling, classification, fitting, and optional alpha training
- runs the experimental bond workflow when enabled
- saves the prepared model through `SavePreparedModel()`

`PotentialAnalysisCommand::BuildDataObject()`:

- attaches `DataRepository` to `request.database_path`
- loads the model and map through `LoadInputFile<T>(...)`, which delegates to `FileIO`
- stores loaded objects in the command-local runtime cache owned by `CommandBase`
- wraps load failures with command-specific error context

`PotentialAnalysisCommand::SavePreparedModel()`:

- persists the prepared model through `SaveStoredObject(...)`
- writes to the repository using `request.saved_key_tag` as the persisted key
- clears sampled local-potential distance/value buffers after persistence to keep runtime state lean

`CommandBase` also creates `output_dir` during filesystem preflight when needed.

## Tests to Update When Behavior Changes

- [`tests/core/command/CommandValidationScenarios_test.cpp`](/tests/core/command/CommandValidationScenarios_test.cpp)
- [`tests/core/command/CommandWorkflowScenarios_test.cpp`](/tests/core/command/CommandWorkflowScenarios_test.cpp)
- [`tests/core/contract/CommandExecutionContract_test.cpp`](/tests/core/contract/CommandExecutionContract_test.cpp)
- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)
- [`tests/integration/python_bindings_runtime_smoke.py`](/tests/integration/python_bindings_runtime_smoke.py)
