# Adding a Command

This guide describes how to add a command to the current command system.

Related references:

- [`docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md)
- [`docs/developer/development-guidelines.md`](/docs/developer/development-guidelines.md)

## Decide Whether It Should Be a New Command

Add a new top-level command when the workflow needs its own:

- request type
- CLI surface
- validation path
- execution path

Extend an existing command when the change is only a new option, a new mode, or a small behavior
branch on an existing workflow.

## Files You Usually Touch

Concrete implementation:

- `src/core/internal/command/<YourCommand>.hpp`
- `src/core/command/<YourCommand>.cpp`

Public API and contract:

- [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp)
- [`include/rhbm_gem/core/command/CommandContract.hpp`](/include/rhbm_gem/core/command/CommandContract.hpp) only when shared metadata or validation/common contract changes
- [`include/rhbm_gem/core/command/CommandEnumClass.hpp`](/include/rhbm_gem/core/command/CommandEnumClass.hpp) only when the command needs a shared enum

Manifest and build wiring:

- [`include/rhbm_gem/core/command/CommandList.def`](/include/rhbm_gem/core/command/CommandList.def)
- [`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp)
- [`src/CMakeLists.txt`](/src/CMakeLists.txt)

Tests and docs:

- grouped command tests under [`tests/core/command/`](/tests/core/command/)
- relevant contract or integration tests when the public surface changes
- command docs under [`docs/developer/commands/`](/docs/developer/commands/)

You usually do not need to edit [`tests/CMakeLists.txt`](/tests/CMakeLists.txt) unless you add a
new test file with a new test responsibility.

## Add the Manifest Entry

Add an entry to [`include/rhbm_gem/core/command/CommandList.def`](/include/rhbm_gem/core/command/CommandList.def).

Stable commands go in the main list. Experimental commands go inside the
`#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` block.

Example:

```cpp
RHBM_GEM_COMMAND(
    Example,
    "example",
    "Run example command")
```

The manifest drives:

- `CommandId`
- command catalog metadata
- `Run*` declarations in `CommandApi.hpp`
- `Run*` definitions in `CommandApi.cpp`
- CLI registration
- Python bindings

`Run*` declarations are generated from the manifest. Do not add a separate handwritten declaration
list for them.

## Scaffold What You Can

The scaffold script lives at
[`resources/tools/developer/command_scaffold.py`](/resources/tools/developer/command_scaffold.py).

Create skeleton files only:

```bash
python3 resources/tools/developer/command_scaffold.py --name Example
```

Create skeleton files and wire stable registration/build lists:

```bash
python3 resources/tools/developer/command_scaffold.py --name Example --wire
```

`--wire` currently updates:

- `CommandList.def`
- the command include list in `CommandApi.cpp`
- `RHBM_GEM_COMMAND_SOURCES` in `src/CMakeLists.txt`

The scaffold does not currently handle these steps for you:

- adding the public request struct to `CommandApi.hpp`
- defining the request `VisitFields(...)` schema
- placing the command in the experimental manifest block
- adding sources to `RHBM_GEM_EXPERIMENTAL_SOURCES`
- adding new shared enums to `CommandEnumClass.hpp`
- writing command-specific validation or execution logic
- updating grouped tests under `tests/core/command/`

## Add the Public Request Type

Add `<YourCommand>Request` to
[`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp).

Each request must define `VisitFields(Visitor &&)` and list every bindable field there. CLI and
Python registration both consume that schema.

Use the existing field-spec helpers as needed:

- `MakeObjectField(...)`
- `MakeScalarField(...)`
- `MakePathField(...)`
- `MakeEnumField(...)`
- `MakeCsvListField(...)`
- `MakeRefGroupField(...)`

If the command adds a new shared enum, update
[`include/rhbm_gem/core/command/CommandEnumClass.hpp`](/include/rhbm_gem/core/command/CommandEnumClass.hpp)
and the Python/common enum bindings will continue to follow that shared definition.

## Implement the Concrete Command

Concrete command classes live under:

- [`src/core/internal/command/`](/src/core/internal/command/)
- [`src/core/command/`](/src/core/command/)

Use this shape:

1. derive from `CommandWithRequest<XxxRequest>`
2. implement `NormalizeRequest()`
3. implement `ValidateOptions()`
4. implement `ResetRuntimeState()`
5. implement `ExecuteImpl()`

`CommandWithRequest<XxxRequest>` already:

- stores the typed request internally
- binds `request.common`
- invalidates prepared state when a new request is applied
- coerces common options through `CoerceCommonRequest(...)`
- calls `NormalizeRequest()`

The `CommandBase` helper API you will usually use is:

- `ValidateRequiredPath(...)`
- `ValidateOptionalPath(...)`
- `RequireNonEmptyText(...)`
- `RequireNonEmptyList(...)`
- `RequireCondition(...)`
- `CoerceScalar(...)`
- `CoerceEnum(...)`
- `BuildOutputPath(...)`
- `ClearParseIssues(...)`
- `ClearPrepareIssues(...)`
- `AddValidationError(...)`
- `AddNormalizationWarning(...)`

Keep parse-time normalization in `NormalizeRequest()` and prepare-time semantic checks in
`ValidateOptions()`.

## CLI and Python Registration

CLI registration is generic in
[`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp).

`ConfigureCommandCli(...)` expands the manifest and:

1. creates one subcommand per manifest entry
2. binds shared options from `CommonCommandRequest`
3. binds command-specific fields from the request schema
4. routes the callback to the matching `Run*` function

Python registration is generic in
[`src/python/CommandApiBindings.cpp`](/src/python/CommandApiBindings.cpp).

That file expands the manifest to expose:

- request classes
- `ExecutionReport`
- shared validation/common types
- shared enums
- one `Run*` function per command

In most cases you do not need to edit `CommandApi.cpp` or `CommandApiBindings.cpp` beyond adding
the command header include or using the scaffold to do that wiring.

## Build Wiring

[`src/CMakeLists.txt`](/src/CMakeLists.txt) currently keeps command sources in two lists:

- `RHBM_GEM_COMMAND_SOURCES` for stable commands
- `RHBM_GEM_EXPERIMENTAL_SOURCES` for experimental commands

Stable commands belong in `RHBM_GEM_COMMAND_SOURCES`. Experimental commands belong in
`RHBM_GEM_EXPERIMENTAL_SOURCES` and must also stay behind the manifest feature guard.

## Tests and Documentation

Grouped command tests live under [`tests/core/command/`](/tests/core/command/).

In most cases you should extend an existing grouped test file:

- `CommandBaseLifecycle_test.cpp` for base lifecycle and filesystem preflight
- `CommandValidationHelpers_test.cpp` for helper semantics
- `CommandValidationScenarios_test.cpp` for command-specific validation rules
- `CommandWorkflowScenarios_test.cpp` for workflow and output behavior

Also update related tests when needed:

- [`tests/core/contract/CommandCatalog_test.cpp`](/tests/core/contract/CommandCatalog_test.cpp)
- [`tests/core/contract/CommandExecutionContract_test.cpp`](/tests/core/contract/CommandExecutionContract_test.cpp)
- [`tests/integration/CommandApiPipeline_test.cpp`](/tests/integration/CommandApiPipeline_test.cpp)

Optional command-specific developer notes belong in:

- [`docs/developer/commands/README.md`](/docs/developer/commands/README.md)
- `docs/developer/commands/<your-command>.md`

## Validation Checklist

Before merge, verify:

1. the command header exists under `src/core/internal/command/` and the source exists under `src/core/command/`
2. `CommandApi.hpp` contains the request struct and `VisitFields(...)`
3. `CommandList.def` contains the manifest entry in the correct stable or experimental section
4. `CommandApi.cpp` includes the new command header
5. `src/CMakeLists.txt` includes the source in the correct stable or experimental list
6. grouped command tests cover the new validation and workflow behavior
7. command docs match the final surface

Recommended checks:

```bash
cmake --build build --target tests_all -j
ctest --test-dir build --output-on-failure -R "rhbm_tests_core_command|rhbm_tests_core_contract"
```
