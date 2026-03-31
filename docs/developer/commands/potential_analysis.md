# Potential Analysis Command

This page documents the current `potential_analysis` command implementation.

## Registration and Implementation Files

The command membership entry lives in
[`src/core/command/CommandManifest.def`](/src/core/command/CommandManifest.def).

Public request and entrypoint:

- [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp)

Internal schema and wiring:

- [`src/core/internal/command/CommandRequestSchema.hpp`](/src/core/internal/command/CommandRequestSchema.hpp)
- [`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp)
- [`src/core/command/CommandCli.cpp`](/src/core/command/CommandCli.cpp)
- [`src/python/CommandApiBindings.cpp`](/src/python/CommandApiBindings.cpp)

Concrete implementation:

- [`src/core/internal/command/PotentialAnalysisCommand.hpp`](/src/core/internal/command/PotentialAnalysisCommand.hpp)
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

Expected outcomes:

- `Succeeded` when the command completes
- `ValidationFailed` when normalization, validation, or output-directory preflight stops execution
- `ExecutionFailed` when preparation succeeds but execution still returns false

Validation details are reported through `result.issues`.

## Command Behavior

`PotentialAnalysisCommand::NormalizeRequest()` handles parse-time normalization:

- validates required model and map paths
- normalizes numeric ranges and sampling options
- resets invalid `saved_key_tag` values back to `"model"` with a parse-phase issue

`PotentialAnalysisCommand::ValidateOptions()` performs prepare-phase semantic checks:

- `--simulation` requires positive simulated resolution
- sampling and fit ranges must be ordered correctly

Generic preparation in `CommandBase` also creates `output_dir` when needed.

## Tests to Update When Behavior Changes

- [`tests/core/command/CommandValidationScenarios_test.cpp`](/tests/core/command/CommandValidationScenarios_test.cpp)
- [`tests/core/command/CommandWorkflowScenarios_test.cpp`](/tests/core/command/CommandWorkflowScenarios_test.cpp)
- [`tests/core/contract/CommandExecutionContract_test.cpp`](/tests/core/contract/CommandExecutionContract_test.cpp)
- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)
- [`tests/integration/python_bindings_runtime_smoke.py`](/tests/integration/python_bindings_runtime_smoke.py)
