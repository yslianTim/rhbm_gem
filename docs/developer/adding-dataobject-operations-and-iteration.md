# Adding DataObject Operations and Iteration

This guide explains how to extend behavior on existing `DataObject` types (`ModelObject`, `MapObject`) and manager iteration.

Related references:

- [`architecture/dataobject-typed-dispatch-architecture.md`](architecture/dataobject-typed-dispatch-architecture.md)
- [`architecture/dataobject-io-architecture.md`](architecture/dataobject-io-architecture.md)
- [`adding-a-command.md`](adding-a-command.md)

## 1. Pick the Correct Boundary

Use the narrowest boundary that matches the change:

- `DataObjectWorkflowOps` (`src/core/workflow/DataObjectWorkflowOps.*`)
  - reusable cross-command typed logic on `ModelObject`/`MapObject`
- command-local workflow (`src/core/workflow/*CommandWorkflow*.cpp`)
  - behavior only needed by one command
- `DataObjectManager::ForEachDataObject(...)`
  - iteration policy or callback contract changes

Do not move command-specific branching into shared `DataObjectWorkflowOps` unless multiple commands need the same contract.

## 2. Required File Updates

If you add or change reusable typed operations:

- `src/core/workflow/DataObjectWorkflowOps.hpp`
- `src/core/workflow/DataObjectWorkflowOps.cpp`
- command call sites using the operation
- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp`
- related command tests under `tests/core/command/`
- docs when contract changes

If you change manager iteration behavior, also update:

- `include/rhbm_gem/data/io/DataObjectManager.hpp`
- `src/data/io/DataObjectManager.cpp`
- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp`
- `tests/data/io/DataObjectManager_test.cpp`

## 3. Operation Design Rules

For `DataObjectWorkflowOps`:

- Keep signatures typed (`ModelObject &`, `MapObject &`), not generic base pointers.
- Keep behavior deterministic for the same object state + options.
- Validate preconditions and throw explicit runtime errors when violated.
- Keep command execution side effects out of shared helpers.

Example pattern:

```cpp
struct ModelPruneOptions
{
    bool selected_only{ true };
    bool remove_unknown_atoms{ false };
};

void PruneModelAtoms(
    ModelObject & model_object,
    ModelPruneOptions options = {});
```

## 4. Iteration Extension Rules

Current `ForEachDataObject(...)` contract:

- mutable and const overloads
- optional key filtering via `key_tag_list`
- `IterateOptions::deterministic_order` controls default-order behavior when key list is empty
- explicit key list preserves caller order
- missing keys are skipped with warning logs
- empty callback throws runtime error
- snapshot traversal: callback runs outside map lock

When extending iteration:

1. Preserve default behavior unless there is a deliberate contract update.
2. Add options in `IterateOptions` in header first.
3. Implement option behavior in snapshot-building path in `.cpp`.
4. Keep callback validation and error shape explicit.
5. Keep callbacks outside map lock.

## 5. Typed Dispatch in Callbacks

Use probe helpers for optional branches:

- `AsModelObject(...)`
- `AsMapObject(...)`

Use expect helpers when mismatch must fail fast:

- `ExpectModelObject(...)`
- `ExpectMapObject(...)`

Keep callback body short; move multi-step logic to typed helper functions.

## 6. Command Integration Pattern

Typical flow in command code:

1. load/build typed objects (`GetTypedDataObject<T>()` or command data-loader helpers)
2. apply shared typed operations from `DataObjectWorkflowOps`
3. use manager iteration only when key selection/order belongs to manager layer
4. keep `ExecuteImpl()` orchestration-focused

## 7. Test Checklist

Add or update tests for:

- operation happy path
- operation invalid preconditions
- iteration order (default and explicit key list)
- iteration snapshot behavior under map mutation
- callback-empty rejection behavior

Common suites:

- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp`
- `tests/data/io/DataObjectManager_test.cpp`
- related command suites under `tests/core/command/`

## 8. Documentation Checklist

When contracts change, update:

- `docs/developer/architecture/dataobject-typed-dispatch-architecture.md`
- `docs/developer/architecture/dataobject-io-architecture.md`
- this guide
- command docs if command behavior or options changed
