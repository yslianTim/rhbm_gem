# Command Architecture

This document describes how command classes are structured in this repository, how they are wired into the CLI, and what implementation rules new commands are expected to follow.

Use this document together with [`../development-guidelines.md`](../development-guidelines.md), especially Section 7. Editable diagram sources for this area live under [`./diagrams/`](./diagrams/).

## 1. Purpose

The command layer is the orchestration boundary of the application:

- It translates CLI or binding inputs into validated command options.
- It builds prerequisites such as parsed files, database-backed objects, and runtime helpers.
- It invokes domain logic from `core`, `data`, and `utils`.
- It emits user-facing logs, output files, and persistence side effects.

The command layer is not the place for file format parsing, DAO details, or low-level algorithms.

## 2. Main Runtime Topology

```mermaid
flowchart TD
    A["main.cpp<br/>CLI::App app"] --> B["Application"]
    B --> C["CommandRegistry::Instance()"]
    C --> D["Registered CommandInfo entries"]

    D --> E["Application::RegisterCommand(...)"]
    E --> F["Factory creates concrete CommandBase"]
    F --> G["CommandBase::RegisterCLIOptions(...)"]
    G --> H["CommandBase::RegisterCLIOptionsBasic(...)"]
    G --> I["Derived::RegisterCLIOptionsExtend(...)"]

    J["CLI11 parse"] --> K["Option setter callbacks"]
    K --> L["Command options updated / normalized"]

    E --> M["Subcommand callback"]
    M --> N["CommandBase::PrepareForExecution()"]
    N --> O["Logger::SetLogLevel(...)"]
    O --> P["ResetRuntimeState() / ClearDataObjects()"]
    P --> Q["ValidateOptions() + filesystem preflight"]
    Q --> R{"Has validation errors?"}
    R -- true --> S["Report issues and abort"]
    R -- false --> T["Derived::Execute()"]

    T --> U["BuildDataObject*() phase"]
    U --> V["DataObjectManager"]
    V --> W["FileProcessFactoryRegistry / DatabaseManager"]
    T --> X["Workflow phase"]
    X --> Y["Painters / writers / algorithms / persistence"]
```

## 3. Registration Model

### 3.1 Static self-registration

Every CLI command is registered locally in its own `.cpp` file through a namespace-scope `CommandRegistrar<T>` object.

Current commands:

- `potential_analysis`
- `potential_display`
- `result_dump`
- `map_simulation`
- `map_visualization`
- `position_estimation`
- `model_test`

This keeps command discovery localized: adding a new command does not require editing a central switch statement.

### 3.2 Registry responsibilities

`CommandRegistry` stores:

- the command name
- the user-facing description
- a factory returning `std::unique_ptr<CommandBase>`

`Application` reads the registry and creates one CLI subcommand per registered entry.

## 4. Base Class Contract

All CLI commands derive from `CommandBase`.

### 4.1 Required virtual contract

Each command must provide:

- `bool Execute()`
- `void RegisterCLIOptionsExtend(CLI::App * command)`
- `const CommandOptions & GetOptions() const`
- `CommandOptions & GetOptions()`

### 4.2 Shared base features

`CommandBase` already owns:

- `DataObjectManager m_data_manager`
- a `ValidationIssue` collection for preflight diagnostics
- base CLI options in `CommandOptions`

Shared base options:

- `thread_size`
- `verbose_level`
- `database_path`
- `folder_path`

Base setters are not passive assignment helpers. They also normalize and validate:

- `SetThreadSize()` clamps invalid values to `1`
- `SetVerboseLevel()` restores an allowed log level when needed
- `SetDatabasePath()` normalizes the path but does not create directories during parse
- `SetFolderPath()` normalizes the path but does not create directories during parse

### 4.3 Common option capability mask

`CommandBase` also exposes `GetCommonOptionMask()` so each concrete command can
choose which shared base options it actually needs:

- `Threading`
- `Verbose`
- `Database`
- `OutputFolder`

File-only commands can still register a hidden deprecated `--database` alias
for CLI compatibility without exposing it in help output or using it at
runtime.

### 4.4 Options pattern

Each concrete command follows the same pattern:

1. Define `struct Options : public CommandOptions`.
2. Store one `Options m_options`.
3. Return `m_options` from both `GetOptions()` overloads.
4. Keep command-specific validation inside setter methods.

This is the preferred extension point for new command parameters.

## 5. Standard Execution Lifecycle

Most commands in this repository follow the same three-step shape:

```mermaid
sequenceDiagram
    participant User as User / CLI
    participant CLI as CLI11
    participant App as Application
    participant Cmd as Concrete Command
    participant DM as DataObjectManager

    User->>CLI: invoke subcommand with options
    CLI->>Cmd: setter callbacks populate / normalize Options
    CLI->>App: run subcommand callback
    App->>Cmd: PrepareForExecution()
    Cmd->>Cmd: ResetRuntimeState()
    Cmd->>Cmd: ValidateOptions()
    alt validation failure
        Cmd-->>User: structured validation issues
    else valid options
        App->>Cmd: Execute()
        Cmd->>Cmd: BuildDataObject*()
        Cmd->>DM: ProcessFile / LoadDataObject / SaveDataObject
        Cmd->>Cmd: Run workflow steps
        Cmd-->>App: true / false
    end
```

Recommended command shape:

1. `RegisterCLIOptionsExtend()` binds CLI options to setter functions.
2. Setters normalize user input and register option-local validation issues when needed.
3. `PrepareForExecution()` resets runtime state, clears command-local data objects, validates cross-field constraints, and performs filesystem preflight.
4. `Execute()` starts by ensuring preparation has happened, then builds prerequisites.
5. `Execute()` then runs the main workflow in clear phases.
6. Any persistence or file output happens after the main computation is ready.

## 6. Data Boundary

`DataObjectManager` is the command layer's gateway for moving data between files, memory, visitors, and the database.

Core operations used by commands:

- `ProcessFile(path, key_tag)` for parsing model or map files into in-memory data objects
- `LoadDataObject(key_tag)` for loading previously saved objects from SQLite
- `SaveDataObject(key_tag, renamed_key_tag)` for persistence
- `GetTypedDataObject<T>(key_tag)` for typed retrieval
- `Accept(visitor, key_tag_list)` for visitor-based traversal

Design intent:

- commands decide *which* objects are needed and *when*
- the data layer decides *how* files and persistence are implemented

When a command needs database-backed objects, call `SetDatabaseManager()` before `LoadDataObject()` or `SaveDataObject()`.

## 7. Current Command Families

### 7.1 File-driven analysis commands

These commands ingest files directly through `ProcessFile(...)`:

| Command | Main inputs | Main phases | Main outputs |
| --- | --- | --- | --- |
| `potential_analysis` | model file, map file, database path | build objects, preprocess, sample, classify, fit, save | updated `ModelObject` persisted to database |
| `map_simulation` | model file | build atom list, simulate maps | generated map files |
| `map_visualization` | model file, map file | preprocess, sample around one atom, visualize | rendered plot/PDF output |
| `position_estimation` | map file | threshold voxels, KD-tree, iterative convergence, deduplication | ChimeraX point file |

### 7.2 Database-driven presentation / export commands

These commands primarily load previously saved `ModelObject` instances from the database:

| Command | Main inputs | Main phases | Main outputs |
| --- | --- | --- | --- |
| `potential_display` | model key list, optional reference groups | load models, apply selection, dispatch painter | painter-specific output files |
| `result_dump` | model key list, optional map file | load models, collect selected atoms, dispatch dump mode | CSV, CIF, CMM and related dump files |

### 7.3 Standalone test harness command

| Command | Main inputs | Main phases | Main outputs |
| --- | --- | --- | --- |
| `model_test` | tester mode and fitting parameters | choose simulation scenario, run `HRLModelTester` workflows | logs and optional ROOT-based plots |

`model_test` is still a `CommandBase` subclass for CLI consistency, but it does not currently rely on `DataObjectManager`.

### 7.4 Command surface matrix

| Command | Uses database at runtime | Uses output folder | Exposed to Python | Hidden deprecated `--database` alias |
| --- | --- | --- | --- | --- |
| `potential_analysis` | yes | yes | yes | no |
| `potential_display` | yes | yes | yes | no |
| `result_dump` | yes | yes | yes | no |
| `map_simulation` | no | yes | yes | yes |
| `map_visualization` | no | yes | no | yes |
| `position_estimation` | no | yes | no | yes |
| `model_test` | no | yes | no | yes |

## 8. Concrete Command Notes

### 8.1 `potential_analysis`

This is the most representative end-to-end command in the repository.

Architecture pattern:

1. Parse model and map files.
2. Optionally adjust model metadata for simulated maps.
3. Normalize map values.
4. Prepare selected atoms and bonds.
5. Sample map values around atoms.
6. Run classification and fitting logic.
7. Save the resulting `ModelObject` back to the database under a caller-controlled key.

If you are adding another analysis-style command, this class is the closest reference implementation.

### 8.2 `potential_display`

This command demonstrates the "database object -> selection -> strategy output" pattern:

1. Load one or more `ModelObject` instances from SQLite.
2. Apply `AtomSelector` rules.
3. Choose a `PainterBase` implementation from `PainterType`.
4. Feed the selected data objects into the chosen painter.

This is the preferred pattern for future visualization or reporting commands that operate on saved analysis results.

### 8.3 `result_dump`

This command demonstrates the "database object -> mode switch -> writer output" pattern:

1. Load one or more saved `ModelObject` instances.
2. Optionally parse a map file if map-based dumping is requested.
3. Build a selected-atom cache.
4. Dispatch to a dump mode selected by `PrinterType`.

If a new export mode is needed, extending this command is usually a better fit than creating a brand new top-level command.

### 8.4 `map_simulation`

This command is intentionally file-driven and does not need a database-backed workflow for its core behavior.

It shows that `CommandBase` can still be reused even when the command only needs:

- common CLI options
- logging level control
- output folder management

### 8.5 `map_visualization`

This command is structurally similar to `potential_analysis`, but it narrows the workflow to:

- load model and map
- normalize map
- pick one atom context
- generate a focused visualization artifact

It is a good reference for single-purpose exploratory commands.

### 8.6 `position_estimation`

This command shows a map-only workflow:

- parse one map
- normalize values
- build a filtered voxel list
- create a KD-tree
- iteratively update candidate points
- deduplicate and export the final positions

It is a good reference for algorithm-heavy commands that still fit the shared CLI lifecycle.

## 9. CLI vs Python Binding Surface

The CLI surface is driven by:

- `main.cpp`
- `Application`
- `CommandRegistry`
- concrete `CommandBase` subclasses

The Python surface is separate and currently exposes only a subset of commands through `bindings/CoreBindings.cpp`:

- `PotentialAnalysisCommand`
- `PotentialDisplayCommand`
- `ResultDumpCommand`
- `MapSimulationCommand`

Implication for future work:

- adding a new CLI command does not automatically expose it to Python
- if a command is intended to be public from Python, bindings must be updated explicitly
- bound command instances call `Execute()` directly, so `Execute()` must remain safe for repeated calls on the same instance
- CLI-only commands should remain documented as such

## 10. Implementation Rules For New Commands

When adding a new command, follow this checklist:

1. Add the public interface under `include/core/`.
2. Add the implementation under `src/core/`.
3. Derive from `CommandBase`.
4. Define `Options : CommandOptions`.
5. Implement validation in setter methods, not in scattered workflow code.
6. Implement `RegisterCLIOptionsExtend()` using setter callbacks.
7. Keep `Execute()` phase-oriented: rely on `PrepareForExecution()`, then build prerequisites, then run the workflow.
8. Use `DataObjectManager` as the boundary for file parsing and persistence.
9. Add a namespace-scope `CommandRegistrar<YourCommand>` in the `.cpp` file.
10. Update bindings, examples, tests, and user-facing docs if the command is part of a supported public workflow.

## 11. What Future Contributors Should Avoid

Avoid these anti-patterns:

- bypassing `CommandRegistry` and hard-coding new CLI subcommands in `Application`
- putting format-specific parsing logic directly inside a command
- delaying basic validation until deep inside `Execute()`
- mutating the filesystem during CLI parse instead of during preflight
- mixing CLI parsing concerns with algorithm implementation details
- letting commands manipulate database internals directly instead of using `DataObjectManager`
- introducing multiple unrelated responsibilities into one command when an existing mode switch or strategy object is a better extension point

## 12. Recommended Reference Files

For future command work, inspect these files first:

- `include/core/CommandBase.hpp`
- `src/core/CommandBase.cpp`
- `include/core/CommandRegistry.hpp`
- `src/core/CommandRegistry.cpp`
- `include/core/Application.hpp`
- `src/core/Application.cpp`
- `include/core/PotentialAnalysisCommand.hpp`
- `src/core/PotentialAnalysisCommand.cpp`
- `include/core/PotentialDisplayCommand.hpp`
- `src/core/PotentialDisplayCommand.cpp`
- `include/core/ResultDumpCommand.hpp`
- `src/core/ResultDumpCommand.cpp`
- `include/core/DataObjectManager.hpp`
- `src/core/DataObjectManager.cpp`
