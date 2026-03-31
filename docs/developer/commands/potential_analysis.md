# Potential Analysis Command

This page documents the current `potential_analysis` command implementation.

## Registration and Implementation Files

The command is declared in
[`include/rhbm_gem/core/command/CommandList.def`](/include/rhbm_gem/core/command/CommandList.def):

```cpp
RHBM_GEM_COMMAND(
    PotentialAnalysis,
    "potential_analysis",
    "Run potential analysis")
```

That manifest entry expands into:

- `CommandId::PotentialAnalysis` and catalog metadata in
  [`include/rhbm_gem/core/command/CommandContract.hpp`](/include/rhbm_gem/core/command/CommandContract.hpp)
- `PotentialAnalysisRequest` and `RunPotentialAnalysis(...)` in
  [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp)
- the `RunPotentialAnalysis(...)` definition and CLI registration in
  [`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp)
- Python request and `Run*` bindings in
  [`src/python/CommandApiBindings.cpp`](/src/python/CommandApiBindings.cpp)

The concrete implementation lives in:

- [`src/core/internal/command/PotentialAnalysisCommand.hpp`](/src/core/internal/command/PotentialAnalysisCommand.hpp)
- [`src/core/command/PotentialAnalysisCommand.cpp`](/src/core/command/PotentialAnalysisCommand.cpp)

## Request Surface

`PotentialAnalysisRequest` is declared in
[`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp).
Its `VisitFields(...)` schema drives both CLI registration and Python bindings.

Shared fields under `request.common`:

- `thread_size`
- `verbose_level`
- `folder_path`

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

Primary CLI flags exposed by this request:

- `-j,--jobs`
- `-v,--verbose`
- `-o,--folder`
- `-d,--database`
- `-a,--model`
- `-m,--map`
- `--simulation`
- `-r,--sim-resolution`
- `-k,--save-key`
- `--training-report-dir`
- `--training-alpha`
- `--asymmetry`
- `-s,--sampling`
- `--sampling-min`
- `--sampling-max`
- `--sampling-height`
- `--fit-min`
- `--fit-max`
- `--alpha-r`
- `--alpha-g`

## Request Application and Lifecycle

`RunPotentialAnalysis(...)` follows the shared command entrypoint flow:

1. construct `PotentialAnalysisCommand`
2. call `ApplyRequest(request)`
3. call `Run()`
4. return `ExecutionReport`

`ApplyRequest(...)` first coerces common options, then runs `PotentialAnalysisCommand::NormalizeRequest()`.
`Run()` performs validation, output-folder preflight, and `ExecuteImpl()`.

## Normalization

`PotentialAnalysisCommand::NormalizeRequest()` uses the current `CommandBase` helper API.

Path normalization and parse-time path checks:

- `model_file_path` uses `ValidateRequiredPath(...)`
- `map_file_path` uses `ValidateRequiredPath(...)`

Scalar normalization with `CoerceScalar(...)`:

- `simulated_map_resolution` must be finite and non-negative
- `sampling_size` must be positive
- `sampling_range_min` and `sampling_range_max` must be finite and non-negative
- `sampling_height` must be finite and positive
- `fit_range_min` and `fit_range_max` must be finite and non-negative
- `alpha_r` and `alpha_g` must be finite and positive

Command-specific normalization behavior:

- `saved_key_tag` is handled manually instead of through `CoerceScalar(...)`
- when `saved_key_tag` is empty, the command resets it to `"model"` and records a parse-phase
  error on `--save-key`
- when `sampling_size <= 0`, the command resets it to `1500` and records an auto-corrected
  parse-phase warning on `--sampling`

## Prepare-Phase Validation

`PotentialAnalysisCommand::ValidateOptions()` performs prepare-phase checks with
`RequireCondition(...)`.

Current command-specific prepare rules:

- `--simulation true` requires `simulated_map_resolution > 0.0`
- `--sampling-min <= --sampling-max`
- `--fit-min <= --fit-max`

Generic preparation from `CommandBase` also creates the output folder from `request.common.folder_path`
when needed.

`training_report_dir` is not validated or created during generic preparation. The command only
tries to create directories for report output when report emission is actually requested.

## Execution Flow

`PotentialAnalysisCommand::ExecuteImpl()` runs this sequence:

1. `BuildDataObject()`
2. `RunMapObjectPreprocessing()`
3. `RunModelObjectPreprocessing()`
4. `RunAtomMapValueSampling()`
5. `RunAtomGroupClassification()`
6. `RunAtomAlphaTraining()` or `RunLocalAtomFitting(request.alpha_r)`
7. `RunAtomPotentialFitting()`
8. `RunExperimentalBondWorkflowIfEnabled()`
9. `SaveDataObject()`

Current behavior of the major steps:

- `BuildDataObject()` sets the database manager, loads the model and map files, and updates the
  model metadata when `simulation_flag` is enabled
- `RunMapObjectPreprocessing()` normalizes the map value array
- `RunModelObjectPreprocessing()` selects atoms and bonds, applies symmetry filtering unless
  `asymmetry_flag` bypasses it, initializes local potential entries, and warns if symmetry
  filtering leaves no selected atoms
- `training_alpha_flag` selects the alpha-training branch; otherwise the command uses
  `RunLocalAtomFitting(request.alpha_r)`
- `RunAtomPotentialFitting()` uses trained per-group `alpha_g` values only when
  `training_alpha_flag` is enabled; otherwise it uses the request value

## Data Persistence and Reports

The command uses two internal transient keys:

- `"model"`
- `"map"`

`SaveDataObject()` writes the processed model object from the internal `"model"` key to the
user-facing `saved_key_tag`.

After saving, the command clears sampled distance and map-value vectors for selected atoms that
still have a local potential entry.

Training-report generation is optional and goes through the internal helper
`detail::EmitTrainingReportIfRequested(...)`.

Current report behavior:

- returns immediately when `training_report_dir` is empty
- creates the report directory only for the requested output path
- logs warnings when directory creation fails
- logs warnings when painting fails or when no file is produced

## Experimental Bond Workflow

`RunExperimentalBondWorkflowIfEnabled()` is compiled only when
`RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` is enabled.

When the feature flag is enabled, the workflow runs only if both the model and map objects were
built successfully. It delegates to
[`src/core/experimental/PotentialAnalysisBondWorkflow.cpp`](/src/core/experimental/PotentialAnalysisBondWorkflow.cpp).

When the feature flag is disabled, the method compiles to a no-op.

## Tests to Check When Behavior Changes

- [`tests/core/command/CommandValidationScenarios_test.cpp`](/tests/core/command/CommandValidationScenarios_test.cpp)
- [`tests/core/command/CommandWorkflowScenarios_test.cpp`](/tests/core/command/CommandWorkflowScenarios_test.cpp)
- [`tests/core/contract/CommandExecutionContract_test.cpp`](/tests/core/contract/CommandExecutionContract_test.cpp)
- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)
- [`tests/integration/python_bindings_runtime_smoke.py`](/tests/integration/python_bindings_runtime_smoke.py)

## Related References

- [`docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md)
- [`docs/developer/adding-a-command.md`](/docs/developer/adding-a-command.md)
- [`resources/README.md`](/resources/README.md)
