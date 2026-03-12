# Adding a Command

Current workflow for adding a new top-level command.

Related references:

- [`architecture/command-architecture.md`](architecture/command-architecture.md)
- [`development-guidelines.md`](development-guidelines.md)

## 1. Decide command scope

Create a new top-level command when:

- the workflow has a distinct entry contract
- options or validation rules differ materially from existing commands
- shared option profile (`FileWorkflow` vs `DatabaseWorkflow`) is different

Prefer extending an existing command when behavior is only a mode/option variant.

## 2. Required file changes

Manual updates for a new command:

- `include/rhbm_gem/core/command/<YourCommand>.hpp`
- `src/core/command/<YourCommand>.cpp`
- `src/core/internal/CommandList.def` (register command)
- `bindings/<YourCommandStem>Bindings.cpp`
- `tests/core/command/<YourCommand>_test.cpp`
- docs impacted by command surface changes

Optional/manual depending on design:

- `src/core/workflow/*` additions
- `src/gui/GuiCommandExecutor.*` if GUI needs this command
- contract tests under `tests/core/contract/`

Generated artifacts from `CommandList.def`:

- `include/rhbm_gem/core/command/internal/CommandIdEntries.generated.inc`
- `src/core/internal/CommandCatalogIncludes.generated.inc`
- `src/core/internal/CommandCatalogEntries.generated.inc`
- `bindings/internal/CommandBindingNames.generated.inc`
- `src/CommandSources.generated.cmake`
- `bindings/CommandBindingUnits.generated.cmake`
- `tests/cmake/CoreCommandTests.generated.cmake`
- generated sections in `docs/developer/architecture/command-architecture.md`

Generate/update:

```bash
python3 scripts/generate_command_artifacts.py
```

## 3. Scaffold (optional)

Create skeleton files:

```bash
python3 scripts/command_scaffold.py --name Example --profile FileWorkflow
```

Create skeletons + wire manifest + regenerate artifacts:

```bash
python3 scripts/command_scaffold.py --name Example --profile FileWorkflow --wire
```

Strict mode (`--wire` required):

```bash
python3 scripts/command_scaffold.py --name Example --profile FileWorkflow --wire --strict
```

Note: scaffolded binding files do not include command-specific setters or `BindCommonCommandSetters(...)`; add them manually.

## 4. Implement command class

Standard structure:

1. `Options` derives from `CommandOptions`.
2. Command derives from `CommandWithProfileOptions<...>` (or `CommandWithOptions<...>`).
3. Setters use `MutateOptions(...)` or `CommandBase` helper setters.
4. CLI-only option wiring goes in `RegisterCLIOptionsExtend(...)`.
5. Cross-field checks go in `ValidateOptions()`.
6. Runtime cache/state cleanup goes in `ResetRuntimeState()`.
7. `ExecuteImpl()` coordinates workflow; avoid embedding parse-time validation.

Implementation rules:

- use `command_cli::Add*Option(...)` helpers for CLI wiring
- do not write `m_options` directly outside mutation helpers
- do not create filesystem directories in setters

## 5. Register command

Add manifest entry in `src/core/internal/CommandList.def`:

```cpp
RHBM_GEM_COMMAND(
    Example,
    ExampleCommand,
    "example_command",
    "Run example command",
    FileWorkflow,
    "ExampleCommand")
```

Then run:

```bash
python3 scripts/generate_command_artifacts.py
```

## 6. Add Python bindings

Bindings are explicit per command in `bindings/*Bindings.cpp`.

Required binding elements:

1. `template <> void BindCommand<YourCommand>(py::module_ & module)`
2. `BindCommandClass<YourCommand>(module)`
3. `.def(py::init<>())`
4. `.def("Execute", &YourCommand::Execute)`
5. `.def(...)` for command-local setters
6. `BindCommonCommandSetters(...)`
7. `BindCommandDiagnostics(...)`

`bindings/CoreBindings.cpp` already loads all commands via `CommandList.def`; no manual module entrypoint edits are needed.

## 7. Tests and docs

Command tests:

- add `tests/core/command/<YourCommand>_test.cpp`
- regenerate artifacts so test CMake list is updated

Contract tests that usually need manual updates for a new command:

- `tests/core/contract/CommandCatalog_test.cpp` (explicit expected metadata/order)
- `tests/core/contract/PublicHeaderSurface_test.cpp` (approved public header list)

`tests/core/contract/DocsSync_test.cpp` usually does not need manual edits; it validates generated blocks in:

- `docs/developer/architecture/command-architecture.md`

## 8. Validation checklist

Before merge:

1. command added to `src/core/internal/CommandList.def` with complete metadata (`COMMAND_ID`, `COMMAND_TYPE`, `CLI_NAME`, `DESCRIPTION`, `PROFILE`, `PYTHON_BINDING_NAME`).
2. generated artifacts refreshed (`python3 scripts/generate_command_artifacts.py`).
3. command bindings complete (`bindings/<YourCommandStem>Bindings.cpp`).
4. command tests and required contract tests updated.
5. docs are synced with final command surface.

Recommended checks:

1. `python3 scripts/check_command_sync.py`
2. `ctest --output-on-failure -R "Command|DocsSync"`
