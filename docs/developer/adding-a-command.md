# Adding a Command

This guide documents how to add a command to the existing command system.

Related references:

- [`docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md)
- [`docs/developer/development-guidelines.md`](/docs/developer/development-guidelines.md)

## 1. Decide whether it should be a new command

Add a new top-level command when the workflow needs its own:

- request type
- CLI surface
- validation path
- execution path

Extend an existing command when the change is only a new mode, option, or small behavior branch.

## 2. Files you usually need to touch

Concrete command implementation:

- `src/core/internal/command/<YourCommand>.hpp`
- `src/core/command/<YourCommand>.cpp`

Public command API:

- [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp)
- [`include/rhbm_gem/core/command/CommandContract.hpp`](/include/rhbm_gem/core/command/CommandContract.hpp) only if shared metadata, diagnostics, or default paths change
- [`include/rhbm_gem/core/command/CommandEnumClass.hpp`](/include/rhbm_gem/core/command/CommandEnumClass.hpp) only if the command adds or extends a shared enum

Manifest and internal registration:

- [`include/rhbm_gem/core/command/CommandList.def`](/include/rhbm_gem/core/command/CommandList.def)
- [`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp)

Build wiring:

- [`src/CMakeLists.txt`](/src/CMakeLists.txt)
- [`tests/CMakeLists.txt`](/tests/CMakeLists.txt)

Tests and docs:

- existing grouped command tests under `tests/core/command/`
- impacted contract or integration tests
- developer or user docs that describe the command surface

Optional:

- [`docs/developer/commands/README.md`](/docs/developer/commands/README.md)
- `docs/developer/commands/<your-command>.md`

## 3. Add the manifest entry

Add an entry to [`include/rhbm_gem/core/command/CommandList.def`](/include/rhbm_gem/core/command/CommandList.def).

Stable commands go in the main list.
Experimental commands go inside the `#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` block.

Example:

```cpp
RHBM_GEM_COMMAND(
    Example,
    "example",
    "Run example command")
```

Parameters:

1. `COMMAND_ID`: token used by `CommandId`, request/run symbol names, and concrete command names
2. `CLI_NAME`: subcommand token
3. `DESCRIPTION`: CLI help text

## 4. Scaffold command files

Generate skeleton files:

```bash
python3 resources/tools/developer/command_scaffold.py --name Example
```

Generate skeleton files and wire stable registration/build lists:

```bash
python3 resources/tools/developer/command_scaffold.py --name Example --wire
```

`--wire` updates:

- `include/rhbm_gem/core/command/CommandList.def`
- `src/core/command/CommandApi.cpp`
- `src/CMakeLists.txt` stable command source list: `RHBM_GEM_COMMAND_SOURCES`

The scaffold does not:

- add the public request struct to `CommandApi.hpp`
- place entries into the experimental manifest section
- place command sources into `RHBM_GEM_EXPERIMENTAL_SOURCES`
- add enum definitions to `CommandEnumClass.hpp`
- add command-specific validation, execution logic, or documentation beyond the generated stub
- update grouped command test files under `tests/core/command/`

## 5. Implement the command class

Concrete command classes live under [`src/core/internal/command/`](/src/core/internal/command/) and [`src/core/command/`](/src/core/command/).

Use this shape:

1. derive from `CommandWithRequest<XxxRequest>`
2. construct the command with the default `CommandWithRequest<XxxRequest>{}` base
3. keep field normalization in `NormalizeRequest()`
4. keep cross-field validation in `ValidateOptions()`
5. reset transient execution state in `ResetRuntimeState()`
6. keep `ExecuteImpl()` focused on orchestration

`CommandWithRequest<XxxRequest>` already:

- stores the public request on the command
- applies `request.common` through `ApplyCommonRequest(...)`
- synchronizes normalized common values back into `request.common`
- calls `NormalizeRequest()` from `ApplyRequest(...)`
- invalidates prepared state when the request changes

Keep extra private members for execution-time state only, such as loaded objects, caches, or
derived lookup tables.

Useful `CommandBase` helpers:

- `ValidateRequiredPath(...)`
- `ValidateOptionalPath(...)`
- `RequireNonEmptyText(...)`
- `RequireNonEmptyList(...)`
- `RequireCondition(...)`
- `CoerceScalar(...)`
- `CoerceEnum(...)`
- `BuildOutputPath(...)`

## 6. Add the public request and run entrypoint

Add the request struct to [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp).

Each public request must define `VisitFields(Visitor &&)` and describe every bindable field there.
CLI registration and Python bindings both consume that schema.

Update shared contract headers only when needed:

- use `CommandContract.hpp` for shared metadata, diagnostics, or defaults
- use `CommandEnumClass.hpp` for shared enum declarations and their CLI/Python mappings

`Run*` declarations and definitions come from `CommandList.def`, so after the manifest entry is in
place you do not add a separate handwritten `RunXxx(...)` declaration or binding list.

Each `Run*` entrypoint follows this sequence:

1. construct the concrete command
2. call `ApplyRequest(...)`
3. call `PrepareForExecution()`
4. call `Execute()` if preparation succeeds
5. return an `ExecutionReport`

## 7. CLI and Python registration

CLI registration is generic:

1. `ConfigureCommandCli(...)` expands the manifest in [`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp)
2. `RegisterCommand<...>()` creates the subcommand and binds request fields from `VisitFields(...)`
3. common fields from `CommonCommandRequest::VisitFields(...)` are always bound
4. command-specific fields such as `database_path` are bound from the request schema itself

Shared CLI options come from `CommonCommandRequest`:

- `-j,--jobs`
- `-v,--verbose`
- `-o,--folder`

Python registration is also generic:

1. `BindCommandApi(...)` expands the manifest in [`src/python/CommandApiBindings.cpp`](/src/python/CommandApiBindings.cpp)
2. `BindRequestType<XxxRequest>(...)` walks `VisitFields(...)` and emits `def_readwrite(...)`
3. `module.def("RunXxx", ...)` is expanded from the manifest

[`src/python/CommandApiBindings.cpp`](/src/python/CommandApiBindings.cpp) also owns the shared Python
bindings for command enums plus validation/common types. It only needs manual changes when the
command adds a new shared enum or shared validation/common type.

## 8. Build wiring for stable vs experimental commands

[`src/CMakeLists.txt`](/src/CMakeLists.txt) has two command source lists:

- `RHBM_GEM_COMMAND_SOURCES` for stable commands
- `RHBM_GEM_EXPERIMENTAL_SOURCES` for experimental commands

[`tests/core/command/`](/tests/core/command/) is organized by test responsibility rather than one
file per command. In most cases you should extend one of the existing grouped files:

- `CommandBaseLifecycle_test.cpp` for base lifecycle/preflight behavior
- `CommandValidationHelpers_test.cpp` for reusable validation/helper semantics
- `CommandValidationScenarios_test.cpp` for command-specific validation rules
- `CommandWorkflowScenarios_test.cpp` for command workflows and output side effects

Only introduce a new command test file when you are creating a genuinely new test responsibility.
If the command is experimental, keep the experimental-only cases behind the existing feature guards
inside the grouped test files unless a separate responsibility file is warranted.

## 9. Tests and documentation

Add or update:

- command unit tests under [`tests/core/command/`](/tests/core/command/)
- contract tests if command metadata, manifest order, or the public surface changes
- integration tests if the end-to-end command behavior changes
- developer or user docs that describe the command

## 10. Validation checklist

Before merge:

1. the command header exists under [`src/core/internal/command/`](/src/core/internal/command/) and the source exists under [`src/core/command/`](/src/core/command/)
2. `CommandApi.hpp` contains the request struct and its `VisitFields(...)` schema
3. `CommandList.def` contains the manifest entry in the correct stable or experimental section
4. `CommandApi.cpp` includes the new command header
5. `src/CMakeLists.txt` includes the source in the correct stable or experimental list
6. grouped command tests under `tests/core/command/` cover the new validation/workflow behavior
7. tests and docs match the final command surface

Recommended checks:

```bash
cmake --build build --target tests_all -j
ctest --output-on-failure -R "Command|Contract|bindings"
```
