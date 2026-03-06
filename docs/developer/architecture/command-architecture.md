# Command Architecture

This document describes the current command system in this repository.
Use it as the command-layer developer manual, not as change history.

Related guides:

- [`../development-guidelines.md`](../development-guidelines.md)
- [`../adding-a-command.md`](../adding-a-command.md)

## 1. Scope and ownership

The command layer is the orchestration boundary between input surfaces (CLI/Python) and domain/data logic.

Command responsibilities:

- accept CLI11 or Python setter input
- store normalized options
- validate options in parse/prepare phases
- prepare runtime prerequisites (log level, output folders, database parent folder)
- orchestrate workflow execution in `ExecuteImpl()`
- report user-facing issues through `Logger` and `ValidationIssue`

Non-responsibilities:

- file-format parsing internals
- database persistence internals
- numerical kernel internals
- painter/printer implementation internals

Those responsibilities stay in `DataObjectManager`, data/file modules, and domain helpers.

## 2. Built-in command source of truth

Built-in command definitions are centralized in:

- `src/core/BuiltInCommandList.def`

This list is consumed by:

- `src/core/BuiltInCommandCatalog.cpp` to build `BuiltInCommandCatalog()` for CLI registration and metadata checks
- `src/core/BuiltInCommandBindingInternal.hpp` to map `CommandType -> python binding name`

`CommandDescriptor` fields (defined in `src/core/BuiltInCommandCatalogInternal.hpp`):

- `id` (`CommandId`)
- `name` (CLI subcommand name)
- `description`
- `common_options` (`CommonOptionMask`)
- `python_binding_name`
- `factory`

Rules enforced by implementation:

- `id` and `common_options` are derived from concrete command type (`kCommandId`, `kCommonOptions`)
- built-in descriptors must provide a non-empty `python_binding_name`
- built-in order in the catalog defines CLI subcommand help order
- there is no plugin/self-registration path for commands

### Built-in command manifest

<!-- BEGIN GENERATED: built-in-command-manifest -->
1. `potential_analysis`
2. `potential_display`
3. `result_dump`
4. `map_simulation`
5. `map_visualization`
6. `position_estimation`
7. `model_test`
<!-- END GENERATED: built-in-command-manifest -->

## 3. CLI startup and registration flow

Startup path:

1. `src/main.cpp` creates `CLI::App`.
2. `Application` (`src/core/Application.cpp`) requires exactly one subcommand.
3. `Application::RegisterAllCommands()` iterates `BuiltInCommandCatalog()`.
4. For each descriptor:
   - create concrete command via `descriptor.factory()`
   - register subcommand (`name`, `description`)
   - call `CommandBase::RegisterCLIOptions()`
   - register callback that calls `Execute()`
5. Callback throws `CLI::RuntimeError(1)` when `Execute()` returns `false`.

There is no second built-in CLI registration path in the project.

## 4. Concrete command contract

Standard command structure:

1. define `Options` struct derived from `CommandOptions`
2. derive from `CommandWithProfileOptions<OptionsT, CommandId::..., Profile>`
3. implement command-specific setters
4. register command-local CLI options in `RegisterCLIOptionsExtend(...)`
5. keep cross-field checks in `ValidateOptions()`
6. clear transient runtime state in `ResetRuntimeState()`
7. keep `ExecuteImpl()` orchestration-focused

Base options from `CommandOptions`:

- `thread_size`
- `verbose_level`
- `database_path`
- `folder_path`

Template aliases in `include/core/CommandBase.hpp`:

- `CommandWithOptions<...>` for explicit mask
- `CommandWithProfileOptions<...>` for profile-driven mask

`CommonOptionProfile` mapping (in `CommandMetadata.hpp`):

- `FileWorkflow` -> `Threading | Verbose | OutputFolder`
- `DatabaseWorkflow` -> `Threading | Verbose | Database | OutputFolder`

## 5. Shared option surface

Shared options are registered in `CommandBase::RegisterCLIOptionsBasic(...)` from `GetCommonOptionsMask()`:

- `-j,--jobs` (`CommonOption::Threading`)
- `-v,--verbose` (`CommonOption::Verbose`)
- `-d,--database` (`CommonOption::Database`)
- `-o,--folder` (`CommonOption::OutputFolder`)

Concrete commands should only register command-local options in `RegisterCLIOptionsExtend(...)`.

### Shared policy matrix

<!-- BEGIN GENERATED: command-surface-matrix -->
| Command | Uses database at runtime | Uses output folder |
| --- | --- | --- |
| `potential_analysis` | yes | yes |
| `potential_display` | yes | yes |
| `result_dump` | yes | yes |
| `map_simulation` | no | yes |
| `map_visualization` | no | yes |
| `position_estimation` | no | yes |
| `model_test` | no | yes |
<!-- END GENERATED: command-surface-matrix -->

## 6. Execution lifecycle and validation model

`Execute()` is the single execution entry for CLI callbacks and Python bindings.

Behavior:

- if command is not prepared, `Execute()` calls `PrepareForExecution()`
- if already prepared, `Execute()` runs `ExecuteImpl()` directly
- after every execution attempt, prepared state is cleared

`PrepareForExecution()` steps:

1. `BeginPreparationPass()`
   - apply log level
   - call `ResetRuntimeState()`
   - clear `m_data_manager` cache (`ClearDataObjects()`)
   - invalidate prepared state
2. `RunValidationPass()`
   - call `ValidateOptions()`
   - abort when any error-level issue exists
3. `RunFilesystemPreflight()`
   - create parent folder for `database_path` when database option is enabled
   - create output folder for `folder_path` when output-folder option is enabled
   - emit issues and finalize prepared state

Validation phases:

| Phase | Typical source | Typical purpose |
| --- | --- | --- |
| `Parse` | setter callbacks | single-field checks, enum/path validation, normalization |
| `Prepare` | `ValidateOptions()` + preflight | cross-field checks, mode-dependent checks, directory creation failures |
| `Runtime` | reserved API phase | runtime failures currently surface through `Execute()` failure + logs |

Prepared-state invalidation rule:

- all setter mutations must go through `MutateOptions(...)` or helpers built on top of it
- this clears prepared state and stale `Prepare` issues

## 7. Setter and option-binding APIs

Core setter/validation APIs in `CommandBase`:

- `MutateOptions(...)`
- `AddValidationError(...)`
- `AddNormalizationWarning(...)`
- `ResetParseIssues(...)`
- `ResetPrepareIssues(...)`

Convenience setter helpers:

- `SetRequiredExistingPathOption(...)`
- `SetOptionalExistingPathOption(...)`
- `SetNormalizedScalarOption(...)`
- `SetFinitePositiveScalarOption(...)`
- `SetFiniteNonNegativeScalarOption(...)`
- `SetPositiveScalarOption(...)`
- `SetValidatedEnumOption(...)`

Data/output helpers on command boundary:

- `RequireDatabaseManager()`
- `BuildOutputPath(...)`

CLI registration helpers in `include/core/CommandOptionBinding.hpp`:

- `command_cli::AddScalarOption(...)`
- `command_cli::AddStringOption(...)`
- `command_cli::AddPathOption(...)`
- `command_cli::AddEnumOption(...)`

## 8. Data access boundary for commands

Use command-facing boundaries; do not bypass into persistence internals.

Primary command-side APIs:

- `m_data_manager.ProcessFile(...)`
- `m_data_manager.LoadDataObject(...)`
- `m_data_manager.SaveDataObject(...)`
- `m_data_manager.GetTypedDataObject(...)`
- `m_data_manager.ForEachDataObject(...)`
- `command_data_loader::*` helpers in `src/core/CommandDataLoaderInternal.hpp`

Workflow helpers used by current commands:

- map/model preprocessing via `DataObjectWorkflowOps` (for example `NormalizeMapObject`, `PrepareModelObject`)
- map sampling via `SampleMapValues(...)`
- command-local workflows (for example `PotentialAnalysisCommandWorkflow.cpp`, `ResultDumpCommandWorkflow.cpp`)

`DataObjectDispatch` (`include/data/DataObjectDispatch.hpp`) is optional when a command iterates generic `DataObjectBase` and needs explicit runtime type dispatch.

## 9. Python binding contract

Built-in Python command classes are manually bound in `bindings/CoreBindings.cpp`.

Current binding model:

- built-in membership/name mapping comes from `BuiltInCommandList.def` + `BuiltInCommandBindingInternal.hpp`
- method exposure is explicit in `CoreBindings.cpp` (not auto-generated)
- Python execution path calls the same `Execute()` / `PrepareForExecution()` contract as C++

### Python built-in surface

<!-- BEGIN GENERATED: built-in-python-command-surface -->
### Built-in Python command classes
- `PotentialAnalysisCommand`
- `PotentialDisplayCommand`
- `ResultDumpCommand`
- `MapSimulationCommand`
- `MapVisualizationCommand`
- `PositionEstimationCommand`
- `HRLModelTestCommand`

### Shared diagnostics types
- `LogLevel`
- `ValidationPhase`
- `ValidationIssue`

### Shared diagnostics methods on built-in Python commands
- `PrepareForExecution()`
- `HasValidationErrors()`
- `GetValidationIssues()`
<!-- END GENERATED: built-in-python-command-surface -->

## 10. Built-in command change checklist

When adding or significantly modifying a built-in command, update the same change set:

1. command header under `include/core/`
2. command implementation under `src/core/`
3. built-in list entry in `src/core/BuiltInCommandList.def`
4. metadata/registration behavior (if needed) in `BuiltInCommandCatalog*`
5. Python binding exposure in `bindings/CoreBindings.cpp`
6. tests (contract + command-specific)
7. this architecture document and related developer docs

## 11. Anti-patterns to avoid

- bypassing built-in registration flow and hard-coding command wiring elsewhere
- reintroducing static self-registration/plugin-style registration
- writing directly to `m_options` outside `MutateOptions(...)`
- pushing obvious single-field validation into `ExecuteImpl()`
- creating directories in setters
- bypassing `DataObjectManager` to access persistence internals
- adding a new top-level command when an enum mode extension is sufficient

## 12. Key reference files

Core command architecture:

- `src/main.cpp`
- `include/core/Application.hpp`
- `src/core/Application.cpp`
- `include/core/CommandBase.hpp`
- `src/core/CommandBase.cpp`
- `include/core/CommandMetadata.hpp`
- `src/core/BuiltInCommandList.def`
- `src/core/BuiltInCommandCatalogInternal.hpp`
- `src/core/BuiltInCommandCatalog.cpp`
- `src/core/BuiltInCommandBindingInternal.hpp`
- `include/core/CommandOptionBinding.hpp`
- `src/core/CommandDataLoaderInternal.hpp`
- `include/core/DataObjectManager.hpp`
- `src/core/DataObjectManager.cpp`
- `bindings/CoreBindings.cpp`

Representative concrete commands:

- `src/core/PotentialAnalysisCommand.cpp`
- `src/core/PotentialDisplayCommand.cpp`
- `src/core/ResultDumpCommand.cpp`
- `src/core/MapSimulationCommand.cpp`
- `src/core/MapVisualizationCommand.cpp`
- `src/core/PositionEstimationCommand.cpp`
- `src/core/HRLModelTestCommand.cpp`
