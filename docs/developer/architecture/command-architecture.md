# Command Architecture

Current command runtime and extension contracts.

Related guides:

- [`../adding-a-command.md`](../adding-a-command.md)
- [`../development-guidelines.md`](../development-guidelines.md)

## 1. Runtime topology

```mermaid
flowchart LR
    A["CLI (main.cpp + CLI11)"] --> B["Application::RegisterAllCommands()"]
    B --> C["CommandCatalog() from CommandList.def"]
    C --> D["Concrete CommandBase instances"]

    E["Python bindings (pybind11)"] --> D
    F["GUI requests (GuiCommandExecutor)"] --> D

    D --> G["PrepareForExecution()"]
    G --> G1["ResetRuntimeState + ClearDataObjects"]
    G --> G2["ValidateOptions()"]
    G --> G3["Filesystem preflight<br/>database parent/output folder"]

    D --> H["Execute() / ExecuteImpl()"]
    H --> I["DataObjectManager + workflows + I/O"]
```

## 2. Source of truth

Top-level command membership is defined in:

- `src/core/internal/CommandList.def`

The manifest drives:

- `CommandCatalog()` construction (`src/core/command/CommandCatalog.cpp`)
- CLI subcommand registration order (`src/core/command/Application.cpp`)
- Python module command registration (`bindings/CoreBindings.cpp`)
- generated catalog/CMake/docs artifacts (`scripts/generate_command_artifacts.py`)

`CommandDescriptor` fields (`src/core/internal/CommandCatalogInternal.hpp`):

- `id`
- `name`
- `description`
- `common_options`
- `python_binding_name`
- `factory`

### Command manifest

<!-- BEGIN GENERATED: command-manifest -->
1. `potential_analysis`
2. `potential_display`
3. `result_dump`
4. `map_simulation`
5. `map_visualization`
6. `position_estimation`
7. `model_test`
<!-- END GENERATED: command-manifest -->

## 3. Registration paths

CLI path:

1. `src/main.cpp` creates `CLI::App`.
2. `Application` requires exactly one subcommand.
3. `RegisterAllCommands()` iterates `CommandCatalog()`.
4. For each descriptor, it creates a command instance, registers CLI options, and binds callback to `Execute()`.
5. Callback throws `CLI::RuntimeError(1)` when execution fails.

Python path:

- `bindings/CoreBindings.cpp` expands `CommandList.def` and calls `BindCommand<...>()` per command.

GUI path:

- `GuiCommandExecutor` creates concrete command objects directly, applies request DTOs via setters, then calls `PrepareForExecution()` and `Execute()`.

## 4. Concrete command contract

Command shape:

1. Define `Options` derived from `CommandOptions`.
2. Derive command class from `CommandWithProfileOptions<...>` or `CommandWithOptions<...>`.
3. Implement command-local setters.
4. Register command-local CLI options in `RegisterCLIOptionsExtend(...)`.
5. Keep cross-field checks in `ValidateOptions()`.
6. Reset transient runtime fields in `ResetRuntimeState()`.
7. Keep `ExecuteImpl()` focused on orchestration.

Base `CommandOptions` fields:

- `thread_size`
- `verbose_level`
- `database_path`
- `folder_path`

`CommonOptionProfile` (`CommandMetadata.hpp`) maps to:

- `FileWorkflow` -> `Threading | Verbose | OutputFolder`
- `DatabaseWorkflow` -> `Threading | Verbose | Database | OutputFolder`

## 5. Shared options

`CommandBase::RegisterCLIOptionsBasic(...)` exposes shared flags by `common_options` mask:

- `-j,--jobs`
- `-v,--verbose`
- `-d,--database`
- `-o,--folder`

Concrete commands should only add command-specific flags in `RegisterCLIOptionsExtend(...)`.

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

## 6. Lifecycle and validation

`Execute()` is the single execution entry point across CLI/Python/GUI.

Behavior:

- If not prepared, `Execute()` calls `PrepareForExecution()`.
- If already prepared, `Execute()` runs `ExecuteImpl()` directly.
- Prepared state is cleared after each execution attempt.

`PrepareForExecution()` sequence:

1. `BeginPreparationPass()`
2. `RunValidationPass()`
3. `RunFilesystemPreflight()`

Validation phases:

| Phase | Typical source | Typical purpose |
| --- | --- | --- |
| `Parse` | setters | single-field validation/normalization |
| `Prepare` | `ValidateOptions()` + preflight | cross-field checks and runtime preflight failures |

Prepared-state invalidation rule:

- All option mutations must go through `MutateOptions(...)` or wrappers built on it.

## 7. Command helper APIs

Core APIs (`CommandBase`):

- `MutateOptions(...)`
- `AddValidationError(...)`
- `AddNormalizationWarning(...)`
- `ResetParseIssues(...)`
- `ResetPrepareIssues(...)`

Convenience setters:

- `SetRequiredExistingPathOption(...)`
- `SetOptionalExistingPathOption(...)`
- `SetNormalizedScalarOption(...)`
- `SetFinitePositiveScalarOption(...)`
- `SetFiniteNonNegativeScalarOption(...)`
- `SetPositiveScalarOption(...)`
- `SetValidatedEnumOption(...)`

Execution boundary helpers:

- `RequireDatabaseManager()`
- `BuildOutputPath(...)`

CLI binding helpers (`CommandOptionBinding.hpp`):

- `command_cli::AddScalarOption(...)`
- `command_cli::AddStringOption(...)`
- `command_cli::AddPathOption(...)`
- `command_cli::AddEnumOption(...)`

## 8. Data boundary

Commands should use command-facing APIs and helpers:

- `DataObjectManager` (`ProcessFile`, `LoadDataObject`, `SaveDataObject`, `ForEachDataObject`, typed getters)
- `command_data_loader::*` (`src/core/internal/CommandDataLoaderInternal.hpp`)
- typed workflow helpers in `src/core/workflow/DataObjectWorkflowOps.*`
- `SampleMapValues(...)` for map sampling

Use `DataObjectDispatch` only when dispatching generic `DataObjectBase` at runtime.

## 9. Python binding contract

Per-command binding units are explicit (`bindings/*Bindings.cpp`), module wiring is manifest-driven (`bindings/CoreBindings.cpp`).

Binding rules:

- command class name comes from `CommandPythonBindingName(...)`
- expose command-local setters + `Execute`
- include shared bindings via `BindCommonCommandSetters(...)` and `BindCommandDiagnostics(...)`
- execution contract is the same `PrepareForExecution()` / `Execute()` path as C++

### Python command surface

<!-- BEGIN GENERATED: command-python-surface -->
### Python command classes
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

### Shared diagnostics methods on Python commands
- `PrepareForExecution()`
- `HasValidationErrors()`
- `GetValidationIssues()`
<!-- END GENERATED: command-python-surface -->

## 10. Update checklist

When adding or changing a command:

1. Update command implementation (`include/.../command/*.hpp`, `src/core/command/*.cpp`).
2. Update membership in `src/core/internal/CommandList.def`.
3. Refresh generated artifacts (`python3 scripts/generate_command_artifacts.py`).
4. Update command bindings (`bindings/*Bindings.cpp`).
5. Update command and contract tests.
6. Keep this document and `docs/developer/adding-a-command.md` in sync.

## 11. Key files

- `src/main.cpp`
- `include/rhbm_gem/core/command/Application.hpp`
- `src/core/command/Application.cpp`
- `include/rhbm_gem/core/command/CommandBase.hpp`
- `src/core/command/CommandBase.cpp`
- `include/rhbm_gem/core/command/CommandMetadata.hpp`
- `src/core/internal/CommandList.def`
- `src/core/internal/CommandCatalogInternal.hpp`
- `src/core/command/CommandCatalog.cpp`
- `include/rhbm_gem/core/command/CommandOptionBinding.hpp`
- `src/core/internal/CommandDataLoaderInternal.hpp`
- `include/rhbm_gem/data/io/DataObjectManager.hpp`
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
