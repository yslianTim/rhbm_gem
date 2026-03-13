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

- `include/rhbm_gem/core/command/CommandApi.hpp`
- `src/core/command/CommandApi.cpp`

CLI runtime binding:

- `src/core/internal/CommandRuntimeRegistry.hpp`
- `src/core/command/CommandRuntimeRegistry.cpp`

Manifest and generated artifacts:

- `src/core/internal/CommandList.def`
- run `python3 scripts/developer/generate_command_artifacts.py`

Python bindings:

- `bindings/CommandApiBindings.cpp`

Tests and docs:

- `tests/core/command/<YourCommand>_test.cpp`
- any integration or contract tests impacted by the new public surface
- docs affected by the new command surface

Optional:

- `src/gui/MainWindow.cpp` and `src/gui/MainWindow.hpp` if the GUI should expose the command
- `docs/developer/commands/<your-command>.md` if you want a command-specific developer note

## 3. Files generated from `CommandList.def`

Do not hand-edit these generated sections/files without rerunning the generator:

- generated `command-id-entries` block in `include/rhbm_gem/core/command/CommandMetadata.hpp`
- generated `command-catalog-entries` block in `src/core/command/CommandCatalog.cpp`
- `src/CommandSources.generated.cmake`
- `tests/cmake/CoreCommandTests.generated.cmake`
- generated blocks in `docs/developer/architecture/command-architecture.md`

Refresh them with:

```bash
python3 scripts/developer/generate_command_artifacts.py
```

## 4. Scaffold (optional)

Create skeleton files:

```bash
python3 scripts/developer/command_scaffold.py --name Example --profile FileWorkflow
```

Create skeletons, append the manifest entry, and regenerate generated artifacts:

```bash
python3 scripts/developer/command_scaffold.py --name Example --profile FileWorkflow --wire
```

Strict wiring mode:

```bash
python3 scripts/developer/command_scaffold.py --name Example --profile FileWorkflow --wire --strict
```

The scaffold does not finish the command for you. You still need to implement the command logic,
update `CommandApi`, add runtime binders, and expose Python bindings.

## 5. Implement the command class

Current command classes are internal implementation types under `src/core/command/`.

Use this shape:

1. Define `Options` deriving from `CommandOptions`.
2. Derive the command from `CommandWithProfileOptions<...>` or `CommandWithOptions<...>`.
3. Implement `ApplyRequest(const XxxRequest&)` as the external configuration entry point.
4. Call `ApplyCommonRequest(request.common)` inside `ApplyRequest(...)`.
5. Keep per-field normalization in setters or `ApplyRequest(...)`.
6. Keep cross-field validation in `ValidateOptions()`.
7. Reset transient runtime state in `ResetRuntimeState()`.
8. Keep `ExecuteImpl()` focused on orchestration.

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

Add the new request struct and run function declaration in
`include/rhbm_gem/core/command/CommandApi.hpp`, then implement the run function in
`src/core/command/CommandApi.cpp`.

Each `Run*` entrypoint follows the same pattern:

1. Construct the concrete command.
2. Call `ApplyRequest(...)`.
3. Call `PrepareForExecution()`.
4. Call `Execute()` if preparation succeeds.
5. Return an `ExecutionReport`.

## 7. Register the command

Add a manifest entry in `src/core/internal/CommandList.def`:

```cpp
RHBM_GEM_COMMAND(
    Example,
    ExampleCommand,
    "example",
    "Run example command",
    FileWorkflow)
```

Then regenerate artifacts:

```bash
python3 scripts/developer/generate_command_artifacts.py
```

## 8. Bind CLI options

Add a runtime binder declaration to `src/core/internal/CommandRuntimeRegistry.hpp`.

In `src/core/command/CommandRuntimeRegistry.cpp`:

1. Add `Bind<YourCommand>RequestOptions(...)` for command-specific CLI flags.
2. Add `Bind<YourCommand>Runtime(...)` that calls `BindRuntime<RequestType>(...)`.
3. Pass the correct `CommonOptionProfile` so shared flags are wired automatically.

Shared flags come from the profile:

- `-j,--jobs`
- `-v,--verbose`
- `-d,--database` for `DatabaseWorkflow`
- `-o,--folder`

Use `command_cli::AddScalarOption(...)`, `AddStringOption(...)`, `AddPathOption(...)`, and
`AddEnumOption(...)` for command-specific flags.

## 9. Add Python bindings

Update `bindings/CommandApiBindings.cpp`:

1. Bind the new `*Request` type.
2. Expose the new `Run*` function.

Shared enums and diagnostics live in `bindings/CommonBindings.cpp`, so that file only needs
changes if the command introduces new shared enum types.

## 10. Tests and documentation

Add or update:

- command unit tests under `tests/core/command/`
- integration tests if the public API or Python surface changes
- developer/user docs if the command changes the documented surface

`tests/core/contract/DocsSync_test.cpp` verifies the generated blocks in
`docs/developer/architecture/command-architecture.md`, so keep the markers intact.

## 11. Validation checklist

Before merge:

1. The command is implemented in `src/core/command/`.
2. `CommandApi.hpp` and `CommandApi.cpp` expose the new request and `Run*` function.
3. CLI binding is added in `CommandRuntimeRegistry.hpp/.cpp`.
4. `src/core/internal/CommandList.def` contains the new manifest entry.
5. Generated artifacts have been refreshed.
6. Python bindings are updated if the binding surface changed.
7. Tests and docs are aligned with the final command surface.

Recommended checks:

```bash
python3 scripts/developer/check_command_sync.py
ctest --output-on-failure -R "Command|DocsSync"
```
