# Adding a Command

Related references:

- [`docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md)
- [`docs/developer/development-guidelines.md`](/docs/developer/development-guidelines.md)

## Files You Usually Touch

Public API:

- [`include/rhbm_gem/core/command/CommandTypes.hpp`](/include/rhbm_gem/core/command/CommandTypes.hpp)

Internal wiring:

- [`src/core/command/detail/CommandCatalog.hpp`](/src/core/command/detail/CommandCatalog.hpp)

Concrete implementation:

- [`src/core/command/detail/`](/src/core/command/detail/)
- [`src/core/command/`](/src/core/command/)

## Command List

Top-level command membership is defined in
the internal command catalog in
[`src/core/command/detail/CommandCatalog.hpp`](/src/core/command/detail/CommandCatalog.hpp).

Each entry uses:

```cpp
CommandEntry<ExampleRequest>{
    "example",
    "Run example command",
    "ExampleRequest",
    ExecuteExampleCommand
},
```

Stable commands live in the main list. Experimental commands stay inside the
`RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` command catalog block.

The command catalog drives:

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

Shared default path helper declarations such as `GetDefaultDatabasePath()` also live in
`CommandTypes.hpp`.

Public execution goes through the typed `RunCommand(request)` API in `CommandSystem.hpp`.
Each request type is associated with its executor through the internal
`CommandEntry<XxxRequest>` catalog entry.

## Internal Request Fields

Add a `RequestFieldCatalog<<YourCommand>Request>` specialization to
[`src/core/command/detail/CommandCatalog.hpp`](/src/core/command/detail/CommandCatalog.hpp).
The request field catalog helpers live in `rhbm_gem::command_internal`.

That internal field catalog is the single source for:

- CLI flag registration
- Python field binding

Declare fields directly with the internal request field:

```cpp
VisitFieldList(visitor,
    RequestField{
        "model_file_path",
        "-a,--model",
        "Model file path",
        &Self::model_file_path });
```

The first argument is the public request field name used by Python bindings.
CLI behavior is inferred from the request member type: paths bind as paths, vectors bind as CSV
lists, enum fields use `CommandEnumTraits` from
[`src/core/command/detail/CommandEnumCatalog.hpp`](/src/core/command/detail/CommandEnumCatalog.hpp),
and reference-group maps bind as repeated group assignments.

## Concrete Command

Implement the concrete command under [`src/core/command/`](/src/core/command/).
Shared command-framework helpers stay under [`src/core/command/detail/`](/src/core/command/detail/).
The command class stays local to its `.cpp`; do not add a matching command header.

Use this shape:

1. include `detail/CommandExecutor.hpp`
2. define a local class deriving from `CommandBase<XxxRequest>`
3. implement `NormalizeAndValidateRequest()`
4. implement `ValidatePreparedRequest()`
5. implement `ExecuteImpl()`
6. expose `command_internal::ExecuteXxxCommand(const XxxRequest & request)`

`CommandBase<XxxRequest>` already:

- stores the typed request internally
- coerces `CommandRequestBase` shared options
- uses shared options for lifecycle/preflight
- calls `NormalizeAndValidateRequest()`

## Registration Surfaces

[`src/core/command/CommandSystem.cpp`](/src/core/command/CommandSystem.cpp) owns typed command
execution, internal CLI registration, and `RunCommandCLI(...)` parsing for the executable entrypoint.

[`src/python/CommandSystemBindings.cpp`](/src/python/CommandSystemBindings.cpp) binds:

- `CommandRequestBase`
- one request type per command
- `CommandResult`
- shared enums
- `CommandDiagnostic`

The command catalog controls which request types and `RunCommand(...)` overloads are registered
there. The request-field list comes from `RequestFieldCatalog` in the same internal catalog.

## Validation Checklist

Before merge, verify:

1. the command source exists under `src/core/command/`
2. `CommandTypes.hpp` contains the new request DTO
3. `CommandCatalog.hpp` contains the internal request field catalog specialization
4. `CommandCatalog.hpp` contains the executor declaration and typed request entry in the correct stable or experimental section
5. the command implementation derives from `CommandBase<XxxRequest>` inside the `.cpp`
6. `src/CMakeLists.txt` includes the source in the correct stable or experimental list
7. grouped command tests cover validation and workflow behavior

Recommended checks:

```bash
cmake --build build --target tests_all -j
ctest --test-dir build --output-on-failure -R "rhbm_tests_core_command|rhbm_tests_core_contract"
```
