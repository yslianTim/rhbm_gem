# Adding a Command

Related references:

- [`docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md)
- [`docs/developer/development-guidelines.md`](/docs/developer/development-guidelines.md)

## Files You Usually Touch

Public API:

- [`include/rhbm_gem/core/command/CommandSystem.hpp`](/include/rhbm_gem/core/command/CommandSystem.hpp)
- [`include/rhbm_gem/core/command/CommandTypes.hpp`](/include/rhbm_gem/core/command/CommandTypes.hpp)

Internal wiring:

- [`src/core/command/detail/CommandRequestSchema.hpp`](/src/core/command/detail/CommandRequestSchema.hpp)
- [`src/core/command/CommandSystem.cpp`](/src/core/command/CommandSystem.cpp) to include the new concrete command header

Concrete implementation:

- [`src/core/command/detail/`](/src/core/command/detail/)
- [`src/core/command/`](/src/core/command/)

## Command List

Top-level command membership is defined in
the `rhbm_gem::command` registry in
[`include/rhbm_gem/core/command/CommandSystem.hpp`](/include/rhbm_gem/core/command/CommandSystem.hpp).

Each entry uses:

```cpp
CommandEntry<ExampleRequest>{
    "example",
    "Run example command",
    "ExampleRequest"
},
```

Stable commands live in the main list. Experimental commands stay inside the
`RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` command registry block.

The command registry drives:

- `ListCommands()`
- CLI registration in `CommandSystem.cpp`
- Python request-type and `RunCommand(...)` overload registration in `CommandSystemBindings.cpp`

## Public Request DTO

Add `<YourCommand>Request` to
[`include/rhbm_gem/core/command/CommandTypes.hpp`](/include/rhbm_gem/core/command/CommandTypes.hpp).

Requests are plain DTOs. Shared options come from `CommandRequestBase`:

- `job_count`
- `verbosity`
- `output_dir`

Command-specific fields live directly on the request type.

Shared default path helpers such as `GetDefaultDatabasePath()` also live in `CommandTypes.hpp`.

Public execution goes through the typed `RunCommand(request)` API in `CommandSystem.hpp`.
Each concrete command still needs a private `CommandSystem.cpp` mapping from `<YourCommand>Request`
to `<YourCommand>Command`, plus an explicit `RunCommand<<YourCommand>Request>(...)` instantiation.

## Internal Request Schema

Add a `CommandRequestSchema<<YourCommand>Request>` specialization to
[`src/core/command/detail/CommandRequestSchema.hpp`](/src/core/command/detail/CommandRequestSchema.hpp).
The request schema helpers live in `rhbm_gem::command_internal`.

That internal schema is the single source for:

- CLI flag registration
- Python field binding

Declare fields directly with the internal field spec `FieldSpec{...}`.
CLI behavior is inferred from the request member type: paths bind as paths, vectors bind as CSV
lists, enum fields use `CommandEnumTraits`, and reference-group maps bind as repeated group
assignments.

## Concrete Command

Implement the concrete command under [`src/core/command/`](/src/core/command/).
Shared command-framework helpers stay under [`src/core/command/detail/`](/src/core/command/detail/).

Use this shape:

1. derive from `CommandBase<XxxRequest>`
2. implement `NormalizeAndValidateRequest()`
3. implement `ValidatePreparedRequest()`
4. implement `ResetRuntimeState()`
5. implement `ExecuteImpl()`

`CommandBase<XxxRequest>` already:

- stores the typed request internally
- coerces `CommandRequestBase` shared options
- uses shared options for lifecycle/preflight
- calls `NormalizeAndValidateRequest()`

## Registration Surfaces

[`src/core/command/CommandSystem.cpp`](/src/core/command/CommandSystem.cpp) owns typed command
execution, `ListCommands()`, internal CLI registration, and `RunCommandCLI(...)` parsing for the
executable entrypoint.

[`src/python/CommandSystemBindings.cpp`](/src/python/CommandSystemBindings.cpp) binds:

- `CommandRequestBase`
- one request type per command
- `CommandResult`
- shared enums
- `ValidationIssue`

The command registry controls which request types and `RunCommand(...)` overloads are registered
there. The request-field list still comes from `CommandRequestSchema`.

## Validation Checklist

Before merge, verify:

1. the command header and source exist together under `src/core/command/`
2. `CommandTypes.hpp` contains the new request DTO
3. `CommandRequestSchema.hpp` contains the internal schema specialization
4. `CommandSystem.hpp` contains the typed entry in the correct stable or experimental section
5. `CommandSystem.cpp` includes the new command header and maps the request type to the command type
6. `src/CMakeLists.txt` includes the source in the correct stable or experimental list
7. grouped command tests cover validation and workflow behavior

Recommended checks:

```bash
cmake --build build --target tests_all -j
ctest --test-dir build --output-on-failure -R "rhbm_tests_core_command|rhbm_tests_core_contract"
```
