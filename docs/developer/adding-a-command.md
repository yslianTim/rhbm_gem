# Adding a Command

Checklist for adding a new top-level command to the current command system.

Related references:

- [`architecture/command-architecture.md`](architecture/command-architecture.md)
- [`development-guidelines.md`](development-guidelines.md)

## 1. Decide whether it should be a new command

Create a new top-level command when the workflow needs its own request type, CLI surface, or
validation/execution path.

Prefer extending an existing command when the change is only a new mode, option, or minor
behavior branch.

## 2. Files that usually need manual changes

Concrete command implementation:

- `src/core/command/<YourCommand>.hpp`
- `src/core/command/<YourCommand>.cpp`

Public command API:

- `include/rhbm_gem/core/command/CommandContract.hpp` only if shared diagnostics/defaults change
- `include/rhbm_gem/core/command/CommandApi.hpp`
- `src/core/command/CommandApi.cpp`

Manifest and CLI registration:

- `include/rhbm_gem/core/command/CommandList.def`
- `src/core/command/CommandCatalog.cpp`

Python bindings:

- `bindings/CommandApiBindings.cpp`

Tests and docs:

- `tests/core/command/<YourCommand>_test.cpp`
- `src/CMakeLists.txt` to add `src/core/command/<YourCommand>.cpp`
- `tests/CMakeLists.txt` to add `tests/core/command/<YourCommand>_test.cpp`
- any integration/contract tests impacted by the public surface
- developer/user docs touched by the command surface

Optional:

- `src/gui/MainWindow.cpp` and `src/gui/MainWindow.hpp` if the GUI should expose the command
- `docs/developer/commands/<your-command>.md` for command-specific notes

## 3. Manifest shape

Add an entry to `include/rhbm_gem/core/command/CommandList.def`:

```cpp
RHBM_GEM_COMMAND(
    Example,
    "example",
    "Run example command",
    FileWorkflow)
```

Parameters:

1. `COMMAND_ID`: `CommandId` enum token.
2. `CLI_NAME`: subcommand token.
3. `DESCRIPTION`: CLI description.
4. `PROFILE`: `FileWorkflow` or `DatabaseWorkflow`.

## 4. Scaffold (optional)

Create skeleton files:

```bash
python3 scripts/developer/command_scaffold.py --name Example --profile FileWorkflow
```

Create skeleton files, append the manifest entry, and update the command/test CMake source lists:

```bash
python3 scripts/developer/command_scaffold.py --name Example --profile FileWorkflow --wire
```

## 5. Implement the command class

Current command classes are internal implementation types under `src/core/command/`.

Use this shape:

1. Define `Options` deriving from `CommandOptions`.
2. Derive the command from `CommandWithOptions<...>`.
3. Construct the command with a `CommonOptionProfile`.
4. Implement `ApplyRequest(const XxxRequest &)`.
5. Call `ApplyCommonRequest(request.common)` inside `ApplyRequest(...)`.
6. Keep per-field normalization in setters or `ApplyRequest(...)`.
7. Keep cross-field validation in `ValidateOptions()`.
8. Reset transient execution state in `ResetRuntimeState()`.
9. Keep `ExecuteImpl()` focused on orchestration.

Useful base helpers from `CommandBase`:

- `MutateOptions(...)`
- `SetRequiredExistingPathOption(...)`
- `SetOptionalExistingPathOption(...)`
- `SetFinitePositiveScalarOption(...)`
- `SetFiniteNonNegativeScalarOption(...)`
- `SetPositiveScalarOption(...)`
- `SetValidatedEnumOption(...)`
- `BuildOutputPath(...)`

## 6. Add the public request and run entrypoint

Add the request struct to `include/rhbm_gem/core/command/CommandApi.hpp`.

Shared diagnostics and `ExecutionReport` live in
`include/rhbm_gem/core/command/CommandContract.hpp`, so adding a normal command usually does not
require contract changes.

`Run*` declarations/definitions are expanded from `CommandList.def`, so once the manifest entry
exists you only need to ensure the request struct, command class, and includes are present.

Each `Run*` entrypoint follows this pattern:

1. Construct the concrete command with the manifest profile.
2. Call `ApplyRequest(...)`.
3. Call `PrepareForExecution()`.
4. Call `Execute()` if preparation succeeds.
5. Return an `ExecutionReport`.

## 7. Bind CLI options

In `src/core/command/CommandCatalog.cpp`:

1. Add `Bind<YourCommand>RequestOptions(...)` for command-specific flags.
2. `RegisterCommandSubcommands(...)` will pick it up through the manifest X-macro expansion.
3. Shared flags are bound automatically from the manifest profile.

Shared flags come from the profile:

- `-j,--jobs`
- `-v,--verbose`
- `-d,--database` for `DatabaseWorkflow`
- `-o,--folder`

## 8. Add Python bindings

Update `bindings/CommandApiBindings.cpp`:

1. Bind the new `*Request` type.
2. The `Run*` export list is expanded from the manifest, so no separate binding list needs to be
   maintained.

Shared enums and diagnostics live in `bindings/CommonBindings.cpp`, so that file only needs
changes if the command introduces a new shared enum.

## 9. Tests and documentation

Add or update:

- command unit tests under `tests/core/command/`
- integration tests if the public API or Python surface changes
- developer/user docs if the command changes the documented surface

## 10. Validation checklist

Before merge:

1. The command is implemented in `src/core/command/`.
2. `CommandApi.hpp` contains the new request.
3. `CommandList.def` contains the new manifest entry.
4. `src/CMakeLists.txt` and `tests/CMakeLists.txt` include the new source/test file.
5. CLI binding is added in `CommandCatalog.cpp`.
6. Python bindings are updated if the binding surface changed.
7. Tests and docs are aligned with the final command surface.

Recommended checks:

```bash
cmake --build build --target tests_all -j
ctest --output-on-failure -R "Command|Contract|bindings"
```
