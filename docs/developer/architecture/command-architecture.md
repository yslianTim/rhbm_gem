# Command Architecture

Current command registration and execution flow.

Related guides:

- [`../adding-a-command.md`](../adding-a-command.md)
- [`../development-guidelines.md`](../development-guidelines.md)

## 1. Source of truth

Top-level command membership is defined in:

- `include/rhbm_gem/core/command/CommandList.def`

Each entry uses:

- `RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)`

This manifest is expanded directly with X-macros by:

- `include/rhbm_gem/core/command/CommandMetadata.hpp`
- `include/rhbm_gem/core/command/CommandApi.hpp`
- `src/core/command/CommandApi.cpp`
- `src/core/command/CommandCatalog.cpp`
- `bindings/CommandApiBindings.cpp`

The manifest does not generate request structs or command-specific CLI field bindings.

Current command list:

1. `potential_analysis`
2. `potential_display`
3. `result_dump`
4. `map_simulation`
5. `map_visualization`
6. `position_estimation`
7. `model_test`

## 2. Execution surfaces

All entrypoints converge on the public `Run*` functions in `CommandApi`.

```mermaid
flowchart LR
    A["CLI (src/main.cpp)"] --> B["ConfigureCommandCli(...)"]
    B --> C["RegisterCommandSubcommands(...)"]
    C --> D["Run* in CommandApi"]

    E["Python (pybind11)"] --> D
    F["GUI (subset in MainWindow)"] --> D

    D --> G["Construct concrete command with profile"]
    G --> H["ApplyRequest(...)"]
    H --> I["PrepareForExecution()"]
    I --> J["Execute()"]
```

## 3. Registry shape

`CommandCatalog()` returns metadata only:

- `id`
- `name`
- `description`
- `profile`

`RegisterCommandSubcommands(CLI::App &)` is the only CLI registration entrypoint.
It creates the request object, binds common options from the profile, binds command-specific
options, and wires the callback to the matching `Run*` function.

There is no exported runtime-binder function object layer anymore.

## 4. Public contract and request surface

Shared command contract types live in `include/rhbm_gem/core/command/CommandContract.hpp`.

This header owns:

- default command data/database paths
- `ValidationPhase`
- `ValidationIssue`
- `ExecutionReport`

Public requests and `Run*` entrypoints live in `include/rhbm_gem/core/command/CommandApi.hpp`.

Shared request fields:

- `thread_size`
- `verbose_level`
- `database_path`
- `folder_path`

`CommonOptionProfile` in `CommandMetadata.hpp` controls shared CLI/common-option behavior:

- `FileWorkflow` -> `Threading | Verbose | OutputFolder`
- `DatabaseWorkflow` -> `Threading | Verbose | Database | OutputFolder`

## 5. Concrete command contract

Concrete command classes are internal types under `src/core/command/`.

Current pattern:

1. Define `Options` deriving from `CommandOptions`.
2. Derive the command from `CommandWithOptions<Options>`.
3. Construct the command with a profile supplied by the caller (`Run*` from the manifest).
4. Implement `ApplyRequest(const XxxRequest &)`.
5. Call `ApplyCommonRequest(request.common)` inside `ApplyRequest(...)`.
6. Keep cross-field validation in `ValidateOptions()`.
7. Reset transient execution state in `ResetRuntimeState()`.
8. Keep `ExecuteImpl()` focused on orchestration.

Useful `CommandBase` helpers:

- `MutateOptions(...)`
- `AddValidationError(...)`
- `AddNormalizationWarning(...)`
- `ResetParseIssues(...)`
- `ResetPrepareIssues(...)`
- `SetRequiredExistingPathOption(...)`
- `SetOptionalExistingPathOption(...)`
- `SetNormalizedScalarOption(...)`
- `SetFinitePositiveScalarOption(...)`
- `SetFiniteNonNegativeScalarOption(...)`
- `SetPositiveScalarOption(...)`
- `SetValidatedEnumOption(...)`
- `BuildOutputPath(...)`

## 6. Lifecycle and validation

`Run*` functions follow this sequence:

1. Construct the concrete command with the manifest profile.
2. Call `ApplyRequest(...)`.
3. Call `PrepareForExecution()`.
4. Return early with validation issues if preparation fails.
5. Call `Execute()`.
6. Return an `ExecutionReport`.

`PrepareForExecution()` runs:

1. `BeginPreparationPass()`
2. `RunValidationPass()`
3. `RunFilesystemPreflight()`

Preflight now stays command-local:

- resets transient runtime state
- clears loaded `DataObjectManager` state
- runs `ValidateOptions()`
- creates the output folder when needed

Database parent directory creation happens when the database layer is actually opened, not during
generic command preparation.

## 7. Command support helpers

Cross-command model/map helper logic lives in:

- `src/core/command/CommandDataSupport.hpp`
- `src/core/command/CommandDataSupport.cpp`

This module groups:

- typed file/database loaders (`command_data_loader::*`)
- map normalization
- model preparation/selection
- atom collection and simulation helpers
- atom/bond context construction

## 8. Python and GUI integration

Python bindings are split across:

- `bindings/CoreBindings.cpp`
- `bindings/CommonBindings.cpp`
- `bindings/CommandApiBindings.cpp`

`bindings/CommonBindings.cpp` exposes shared enums and diagnostics.
`bindings/CommandApiBindings.cpp` exposes request structs, `ExecutionReport`, and all `Run*`
functions via the manifest X-macro.

The GUI is not manifest-driven for forms. `src/gui/MainWindow.cpp` still builds request objects
manually for its supported subset, but reads names/profile metadata from `CommandCatalog()`.
