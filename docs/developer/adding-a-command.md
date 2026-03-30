# Adding a Command

Checklist for adding a new top-level command to the current command system.

Related references:

- [`/docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md)
- [`/docs/developer/development-guidelines.md`](/docs/developer/development-guidelines.md)

## 1. Decide whether it should be a new command

Create a new top-level command when the workflow needs its own request type, CLI surface, or
validation/execution path.

Prefer extending an existing command when the change is only a new mode, option, or minor
behavior branch.

## 2. Files that usually need manual changes

Concrete command implementation:

- `/src/core/internal/command/<YourCommand>.hpp`
- `/src/core/command/<YourCommand>.cpp`

Public command API:

- `/include/rhbm_gem/core/command/CommandContract.hpp` if shared diagnostics/defaults or metadata change
- `/include/rhbm_gem/core/command/CommandApi.hpp`

Manifest and internal registration:

- `/include/rhbm_gem/core/command/CommandList.def`
- `/src/core/internal/command/CommandRegistry.hpp`

Tests and docs:

- `/tests/core/command/<YourCommand>_test.cpp`
- `/src/CMakeLists.txt` to add `/src/core/command/<YourCommand>.cpp`
- `/tests/CMakeLists.txt` to add `/tests/core/command/<YourCommand>_test.cpp`
- any integration/contract tests impacted by the public surface
- developer/user docs touched by the command surface

Optional:

- [`/docs/developer/commands/README.md`](/docs/developer/commands/README.md) for the command-doc index and `/docs/developer/commands/<your-command>.md` for command-specific notes

## 3. Manifest shape

Add an entry to `/include/rhbm_gem/core/command/CommandList.def`.

Stable commands are listed directly in the manifest. Commands gated by
`RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` stay in the same file inside the experimental `#ifdef`
section.

Example entry:

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
python3 resources/tools/developer/command_scaffold.py --name Example --profile FileWorkflow
```

Create skeleton files, append the manifest entry, update `CommandRegistry.hpp`, and update the
command/test CMake source lists:

```bash
python3 resources/tools/developer/command_scaffold.py --name Example --profile FileWorkflow --wire
```

## 5. Implement the command class

Current command classes are internal implementation types with headers under
`/src/core/internal/command/` and implementation units under `/src/core/command/`.

Use this shape:

1. Derive the command from `CommandWithRequest<XxxRequest>`.
2. Construct the command with a `CommonOptionProfile`.
3. Keep per-field normalization in `NormalizeRequest()`.
4. Keep cross-field validation in `ValidateOptions()`.
5. Reset transient execution state in `ResetRuntimeState()`.
6. Keep `ExecuteImpl()` focused on orchestration.

`CommandWithRequest<XxxRequest>` already:

- stores the public request as the canonical config
- applies `request.common` through `ApplyCommonRequest(...)`
- synchronizes normalized common fields back into `request.common`
- invalidates prepared state when the request changes

Only keep extra private members for execution-time transient state such as loaded objects, caches,
or derived lookup tables.

Useful base helpers from `CommandBase`:

- `AssignOption(...)` for plain field assignment that should invalidate prepared state
- `MutateOptions(...)` for advanced setters that also manage parse issues, fallback logic, or
  derived state
- `SetRequiredExistingPathOption(...)`
- `SetOptionalExistingPathOption(...)`
- `SetFinitePositiveScalarOption(...)`
- `SetFiniteNonNegativeScalarOption(...)`
- `SetPositiveScalarOption(...)`
- `SetValidatedEnumOption(...)`
- `BuildOutputPath(...)`

## 6. Add the public request and run entrypoint

Add the request struct to `/include/rhbm_gem/core/command/CommandApi.hpp`.

Each public request is also the schema source of truth for CLI and Python, so define
`static void VisitFields(Visitor &&)` on the request and describe every bindable field there.

Shared diagnostics/defaults live in `/include/rhbm_gem/core/command/CommandContract.hpp`, while
`ExecutionReport` lives in `/include/rhbm_gem/core/command/CommandApi.hpp`, so adding a normal
command usually does not require changes to either shared contract layer.

If the command introduces or extends a shared command enum, update the enum declaration and its
CLI / Python / validation mappings together in
`/include/rhbm_gem/core/command/CommandEnumClass.hpp`.

`Run*` declarations/definitions are expanded from `CommandList.def`, so once the manifest entry
exists in the correct fragment you only need to ensure the request struct, command class, and
registry include are present.

Each `Run*` entrypoint follows this pattern:

1. Construct the concrete command with the manifest profile.
2. Call `ApplyRequest(...)`.
3. Call `PrepareForExecution()`.
4. Call `Execute()` if preparation succeeds.
5. Return an `ExecutionReport`.

## 7. CLI registration

CLI registration is generic.

1. `ConfigureCommandCli(...)` expands the manifest in `/src/core/command/CommandOptionSupport.cpp`.
2. `RegisterCommand<...>()` creates the subcommand and binds fields by walking
   `XxxRequest::VisitFields(...)`.
3. Shared flags are filtered from `CommonCommandRequest::VisitFields(...)` using the manifest
   profile.
4. Only add custom parsing code when a field shape cannot be expressed by the built-in field
   specs.

Shared flags come from the profile:

- `-j,--jobs`
- `-v,--verbose`
- `-d,--database` for `DatabaseWorkflow`
- `-o,--folder`

## 8. Python bindings

Python request binding is also generic.

1. `BindCommandApi(...)` expands the manifest in `/src/python/CommandApiBindings.cpp`.
2. `BindRequestType<XxxRequest>(...)` walks `XxxRequest::VisitFields(...)` and emits the
   `def_readwrite(...)` properties.
3. The `Run*` export list is still expanded from the manifest, so no separate handwritten request
   binding list needs to be maintained.

Shared enums and diagnostics live in `/src/python/CommonBindings.cpp`, so that file only needs
changes if the command introduces a new shared enum. The enum token source of truth remains
`/include/rhbm_gem/core/command/CommandEnumClass.hpp`.

## 9. Tests and documentation

Add or update:

- command unit tests under `/tests/core/command/`
- integration tests if the public API or Python surface changes
- developer/user docs if the command changes the documented surface

## 10. Validation checklist

Before merge:

1. The command header is implemented in `/src/core/internal/command/` and the source in `/src/core/command/`.
2. `CommandApi.hpp` contains the new request and its `VisitFields(...)` schema.
3. `CommandList.def` contains the new entry in the correct stable or experimental section.
4. `CommandRegistry.hpp` includes the new command header.
5. `/src/CMakeLists.txt` and `/tests/CMakeLists.txt` include the new source/test file.
6. Tests and docs are aligned with the final command surface.

Recommended checks:

```bash
cmake --build build --target tests_all -j
ctest --output-on-failure -R "Command|Contract|bindings"
```
