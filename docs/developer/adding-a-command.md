# Adding a Command

Related references:

- [`docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md)
- [`docs/developer/development-guidelines.md`](/docs/developer/development-guidelines.md)

## Files You Usually Touch

Public API:

- [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp)
- [`include/rhbm_gem/core/command/CommandEnums.hpp`](/include/rhbm_gem/core/command/CommandEnums.hpp) only when the command needs a new shared enum

Internal wiring:

- [`src/core/command/CommandManifest.def`](/src/core/command/CommandManifest.def)
- [`src/core/internal/command/CommandRequestSchema.hpp`](/src/core/internal/command/CommandRequestSchema.hpp)
- [`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp) to include the new concrete command header

Concrete implementation:

- [`src/core/internal/command/`](/src/core/internal/command/)
- [`src/core/command/`](/src/core/command/)

## Manifest

Top-level command membership is defined in
[`src/core/command/CommandManifest.def`](/src/core/command/CommandManifest.def).

Each entry uses:

```cpp
RHBM_GEM_COMMAND(
    Example,
    "example",
    "Run example command")
```

Stable commands live in the main list. Experimental commands stay inside the
`#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` block.

The manifest drives:

- public `Run*` declarations in `CommandApi.hpp`
- `ListCommands()`
- `Run*` definitions in `CommandApi.cpp`
- CLI registration in `CommandCli.cpp`
- Python request-type and `Run*` binding registration in `CommandApiBindings.cpp`

## Public Request DTO

Add `<YourCommand>Request` to
[`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp).

Requests are plain DTOs. Shared options come from `CommandRequestBase`:

- `job_count`
- `verbosity`
- `output_dir`

Command-specific fields live directly on the request type.

Shared default path helpers such as `GetDefaultDatabasePath()` also live in `CommandApi.hpp`.

`Run*` declarations are generated from `CommandManifest.def`, so adding a request DTO does not
require a handwritten `RunYourCommand(...)` declaration.

## Internal Request Schema

Add a `CommandRequestSchema<<YourCommand>Request>` specialization to
[`src/core/internal/command/CommandRequestSchema.hpp`](/src/core/internal/command/CommandRequestSchema.hpp).

That internal schema is the single source for:

- CLI flag registration
- Python field binding

Use the internal helpers there:

- `MakeScalarField(...)`
- `MakePathField(...)`
- `MakeEnumField(...)`
- `MakeCsvListField(...)`
- `MakeRefGroupField(...)`

## Concrete Command

Implement the concrete command under:

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
- binds the request to `CommandBase`
- coerces `CommandRequestBase` shared options through `CoerceBaseRequest(...)`
- calls `NormalizeRequest()`

## Registration Surfaces

[`src/core/command/CommandApi.cpp`](/src/core/command/CommandApi.cpp) owns public `Run*`
definitions and `ListCommands()`.

[`src/core/command/CommandCli.cpp`](/src/core/command/CommandCli.cpp) owns CLI registration
through the internal `ConfigureCommandCli(...)`.

[`src/python/CommandApiBindings.cpp`](/src/python/CommandApiBindings.cpp) binds:

- `CommandRequestBase`
- one request type per command
- `CommandResult`
- shared enums
- `ValidationIssue`

The manifest controls which request types and `Run*` functions are registered there. The
request-field list still comes from `CommandRequestSchema`.

## Validation Checklist

Before merge, verify:

1. the command header exists under `src/core/internal/command/` and the source exists under `src/core/command/`
2. `CommandApi.hpp` contains the new request DTO
3. `CommandRequestSchema.hpp` contains the internal schema specialization
4. `CommandManifest.def` contains the manifest entry in the correct stable or experimental section
5. `CommandApi.cpp` includes the new command header
6. `src/CMakeLists.txt` includes the source in the correct stable or experimental list
7. grouped command tests cover validation and workflow behavior

Recommended checks:

```bash
cmake --build build --target tests_all -j
ctest --test-dir build --output-on-failure -R "rhbm_tests_core_command|rhbm_tests_core_contract"
```
