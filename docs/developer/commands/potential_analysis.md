# Potential Analysis Command

This page documents the existing `potential_analysis` command implementation.

## Registration and public surface

The command is declared in
[`include/rhbm_gem/core/command/CommandList.def`](/include/rhbm_gem/core/command/CommandList.def):

```cpp
RHBM_GEM_COMMAND(
    PotentialAnalysis,
    "potential_analysis",
    "Run potential analysis")
```

That manifest entry is expanded into:

- [`include/rhbm_gem/core/command/CommandContract.hpp`](/include/rhbm_gem/core/command/CommandContract.hpp) for shared metadata and `CommandId`
- [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp) for `PotentialAnalysisRequest` and `RunPotentialAnalysis(...)`
- [`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp) for the `RunPotentialAnalysis(...)` definition
- [`src/core/command/CommandOptionSupport.cpp`](/src/core/command/CommandOptionSupport.cpp) for CLI registration
- [`src/python/CommandApiBindings.cpp`](/src/python/CommandApiBindings.cpp) for Python request and `Run*` bindings

The concrete implementation lives in:

- [`src/core/internal/command/PotentialAnalysisCommand.hpp`](/src/core/internal/command/PotentialAnalysisCommand.hpp)
- [`src/core/command/PotentialAnalysisCommand.cpp`](/src/core/command/PotentialAnalysisCommand.cpp)

`potential_analysis` exposes these common and command-specific CLI options:

- `-j,--jobs`
- `-v,--verbose`
- `-d,--database`
- `-o,--folder`

## Request schema

`PotentialAnalysisRequest` is declared in
[`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp).
Its `VisitFields(...)` definition is the schema used by both CLI registration and Python bindings.

The request fields are:

- Common options: `common.thread_size`, `common.verbose_level`, `common.folder_path`
- Database option: `database_path`
- Input paths: `model_file_path`, `map_file_path`
- Runtime flags: `simulation_flag`, `simulated_map_resolution`, `saved_key_tag`, `training_report_dir`, `training_alpha_flag`, `asymmetry_flag`
- Sampling controls: `sampling_size`, `sampling_range_min`, `sampling_range_max`, `sampling_height`
- Fitting controls: `fit_range_min`, `fit_range_max`
- Alpha controls: `alpha_r`, `alpha_g`

`RunPotentialAnalysis(...)` constructs `PotentialAnalysisCommand`, calls `ApplyRequest(...)`,
calls `PrepareForExecution()`, and then executes the command only if preparation succeeds.

## Request normalization and validation

`PotentialAnalysisCommand::NormalizeRequest()` applies per-field normalization. The command uses
`CommandBase` helpers for most request fields:

- `model_file_path` and `map_file_path` use `SetRequiredExistingPathOption(...)`
- `simulated_map_resolution` uses `SetFiniteNonNegativeScalarOption(...)`
- `sampling_size` uses `SetNormalizedScalarOption(...)`
- `sampling_range_min` and `sampling_range_max` use `SetFiniteNonNegativeScalarOption(...)`
- `sampling_height` uses `SetFinitePositiveScalarOption(...)`
- `fit_range_min` and `fit_range_max` use `SetFiniteNonNegativeScalarOption(...)`
- `alpha_r` and `alpha_g` use `SetFinitePositiveScalarOption(...)`

Command-specific normalization details:

- `saved_key_tag` treats an empty string as a parse-phase validation error, resets the value to
  `"model"`, and records the issue on `--save-key`
- `sampling_size` falls back to `1500` with a parse-phase warning when the value is not positive

`PotentialAnalysisCommand::ValidateOptions()` performs prepare-phase checks:

- `--simulation true` requires `simulated_map_resolution > 0.0`
- `--sampling-min <= --sampling-max`
- `--fit-min <= --fit-max`

`training_report_dir` is optional. The command does not validate or create it during generic
preparation. Its parent directory is created only if training-report emission is attempted.

## Execution flow

`PotentialAnalysisCommand::ExecuteImpl()` runs this sequence:

1. `BuildDataObject()`
2. `RunMapObjectPreprocessing()`
3. `RunModelObjectPreprocessing()`
4. `RunAtomMapValueSampling()`
5. `RunAtomGroupClassification()`
6. `RunAtomAlphaTraining()` or `RunLocalAtomFitting(...)`
7. `RunAtomPotentialFitting()`
8. `RunExperimentalBondWorkflowIfEnabled()`
9. `SaveDataObject()`

`BuildDataObject()`:

- calls `m_data_manager.SetDatabaseManager(request.database_path)`
- loads the model through `command_data_loader::ProcessModelFile(...)`
- loads the map through `command_data_loader::ProcessMapFile(...)`
- updates the model resolution metadata when `simulation_flag` is enabled

`RunMapObjectPreprocessing()` calls `NormalizeMapObject(...)`.

`RunModelObjectPreprocessing()` calls `PrepareModelObject(...)` with:

- atom selection enabled
- bond selection enabled
- symmetry filtering enabled for atoms and bonds unless `asymmetry_flag` bypasses it
- local atom and bond entry initialization enabled

If symmetry filtering removes every atom, the command logs a warning suggesting
`--asymmetry true`.

## Data persistence and reports

The command uses two internal transient keys:

- `m_model_key_tag = "model"`
- `m_map_key_tag = "map"`

`SaveDataObject()` saves the processed model object from the internal `"model"` key to the
user-facing `saved_key_tag`.

After saving, the command clears sampled distance and map-value lists for selected atoms that
still have a local potential entry.

Alpha-training PDF output is optional. Report emission is routed through
`detail::EmitTrainingReportIfRequested(...)`, which:

- returns early when `training_report_dir` is empty
- creates the report directory only for the requested report path
- logs warnings if directory creation or painting fails

## Experimental bond workflow

`RunExperimentalBondWorkflowIfEnabled()` is compiled only when
`RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` is enabled.

When that feature flag is enabled, the method runs only if both the model and map objects were
built successfully, then delegates to
[`src/core/experimental/PotentialAnalysisBondWorkflow.cpp`](/src/core/experimental/PotentialAnalysisBondWorkflow.cpp).

When the feature flag is disabled, the method compiles to a no-op.

## Tests to update when behavior changes

- [`tests/core/command/PotentialAnalysisCommand_test.cpp`](/tests/core/command/PotentialAnalysisCommand_test.cpp)
- [`tests/core/command/CommandApi_test.cpp`](/tests/core/command/CommandApi_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)
- [`tests/integration/python_bindings_runtime_smoke.py`](/tests/integration/python_bindings_runtime_smoke.py)
- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)

## Related references

- [`docs/developer/adding-a-command.md`](/docs/developer/adding-a-command.md)
- [`docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md)
- [`resources/README.md`](/resources/README.md)
