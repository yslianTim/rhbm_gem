# Command Architecture

This document describes the current command system in this repository.
Use it as the command-layer developer manual, not as change history.

Related guides:

- [`../development-guidelines.md`](../development-guidelines.md)
- [`../adding-a-command.md`](../adding-a-command.md)

## 1. Scope and ownership

The command layer is the orchestration boundary between input surfaces
(CLI/Python/GUI execution requests) and domain/data logic.

Command responsibilities:

- accept CLI11 callbacks, Python setter calls, and GUI request-adapter setter calls
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

- `src/core/internal/BuiltInCommandList.def`

This list is consumed by:

- `src/core/command/BuiltInCommandCatalog.cpp` to build `BuiltInCommandCatalog()` for CLI registration and metadata checks
- `src/core/internal/BuiltInCommandBindingInternal.hpp` to map `CommandType -> python binding name`

`CommandDescriptor` fields (defined in `src/core/internal/BuiltInCommandCatalogInternal.hpp`):

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

1. `src/main.cpp` creates `CLI::App` and one `DataIoServices` instance.
2. `Application` (`src/core/command/Application.cpp`) receives `DataIoServices` and requires exactly one subcommand.
3. `Application::RegisterAllCommands()` iterates `BuiltInCommandCatalog()`.
4. For each descriptor:
   - create concrete command via `descriptor.factory(data_io_services)`
   - register subcommand (`name`, `description`)
   - call `CommandBase::RegisterCLIOptions()`
   - register callback that calls `Execute()`
5. Callback throws `CLI::RuntimeError(1)` when `Execute()` returns `false`.

There is no second built-in CLI registration path in the project.

GUI integration path:

- `GuiCommandExecutor` (`src/gui/GuiCommandExecutor.cpp`) does not register subcommands
- it creates concrete command instances directly, applies request DTOs through setters, then runs `PrepareForExecution()`/`Execute()`
- this path reuses the same command lifecycle and validation APIs as CLI/Python

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

Template aliases in `include/rhbm_gem/core/command/CommandBase.hpp`:

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

`Execute()` is the single execution entry for CLI callbacks, Python bindings, and GUI adapters.

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

CLI registration helpers in `include/rhbm_gem/core/command/CommandOptionBinding.hpp`:

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
- `command_data_loader::*` helpers in `src/core/internal/CommandDataLoaderInternal.hpp`

Workflow helpers used by current commands:

- map/model preprocessing via `DataObjectWorkflowOps` (for example `NormalizeMapObject`, `PrepareModelObject`)
- map sampling via `SampleMapValues(...)`
- command-local workflows (for example `PotentialAnalysisTrainingWorkflow.cpp`, `PotentialAnalysisAnalysisWorkflow.cpp`, `PotentialAnalysisReportWorkflow.cpp`, `ResultDumpCommandWorkflow.cpp`)

`DataObjectDispatch` (`include/rhbm_gem/data/dispatch/DataObjectDispatch.hpp`) is optional when a command iterates generic `DataObjectBase` and needs explicit runtime type dispatch.

## 9. Python binding contract

Built-in Python command classes are bound through per-command binding units under `bindings/`,
with module assembly in `bindings/CoreBindings.cpp`.

Current binding model:

- built-in membership/name mapping comes from `BuiltInCommandList.def` + `BuiltInCommandBindingInternal.hpp`
- method exposure is explicit in `bindings/*Bindings.cpp` (not auto-generated)
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

1. command header under `include/rhbm_gem/core/command/`
2. command implementation under `src/core/command/` (and `src/core/workflow/` when introducing/extracting workflows)
3. built-in list entry in `src/core/internal/BuiltInCommandList.def`
4. metadata/registration behavior (if needed) in `BuiltInCommandCatalog*`
5. Python binding exposure in `bindings/*Bindings.cpp` (and wire module init in `bindings/CoreBindings.cpp` if adding a new binding unit)
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
- `include/rhbm_gem/core/command/Application.hpp`
- `src/core/command/Application.cpp`
- `include/rhbm_gem/core/command/CommandBase.hpp`
- `src/core/command/CommandBase.cpp`
- `include/rhbm_gem/core/command/CommandMetadata.hpp`
- `src/core/internal/BuiltInCommandList.def`
- `src/core/internal/BuiltInCommandCatalogInternal.hpp`
- `src/core/command/BuiltInCommandCatalog.cpp`
- `src/core/internal/BuiltInCommandBindingInternal.hpp`
- `include/rhbm_gem/core/command/CommandOptionBinding.hpp`
- `src/core/internal/CommandDataLoaderInternal.hpp`
- `include/rhbm_gem/data/io/DataObjectManager.hpp`
- `src/data/io/DataObjectManager.cpp`
- `include/rhbm_gem/gui/GuiCommandExecutor.hpp`
- `src/gui/GuiCommandExecutor.cpp`
- `bindings/BindingHelpers.hpp`
- `bindings/CommonBindings.cpp`
- `bindings/CoreBindings.cpp`

Representative concrete commands:

- `src/core/command/PotentialAnalysisCommand.cpp`
- `src/core/command/PotentialDisplayCommand.cpp`
- `src/core/command/ResultDumpCommand.cpp`
- `src/core/command/MapSimulationCommand.cpp`
- `src/core/command/MapVisualizationCommand.cpp`
- `src/core/command/PositionEstimationCommand.cpp`
- `src/core/command/HRLModelTestCommand.cpp`
