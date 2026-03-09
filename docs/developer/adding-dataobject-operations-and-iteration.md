# Adding DataObject Operations and Iteration

This guide is for extending behavior on current `DataObject` types (`ModelObject`, `MapObject`)
without adding a new top-level `DataObject` subtype.

Use this when you need to:

- add or refine reusable operations on existing objects
- extend command workflows that process existing objects
- adjust `DataObjectManager::ForEachDataObject(...)` iteration behavior/options

Related references:

- [`architecture/dataobject-typed-dispatch-architecture.md`](architecture/dataobject-typed-dispatch-architecture.md)
- [`architecture/dataobject-io-architecture.md`](architecture/dataobject-io-architecture.md)
- [`adding-a-command.md`](adding-a-command.md)

## 1. Decide the extension boundary first

Choose the narrowest boundary that matches your change:

- `DataObjectWorkflowOps`:
  - cross-command reusable operation on `ModelObject`/`MapObject`
- command-local `*Workflow.cpp`:
  - command-specific behavior that should not become global shared API
- `DataObjectManager::ForEachDataObject(...)`:
  - iteration policy change (ordering, selection, callback contract)

Avoid putting command-specific branching into `DataObjectWorkflowOps`.

## 2. Required change set

For operation extension on existing objects, update these areas together:

- `src/core/workflow/DataObjectWorkflowOps.hpp`
- `src/core/workflow/DataObjectWorkflowOps.cpp`
- command workflow call sites that consume the new operation
- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp`
- relevant command tests under `tests/core/command/`
- architecture/developer docs when operation contract changes

For iteration method extension, additionally update:

- `include/rhbm_gem/core/command/DataObjectManager.hpp`
- `src/data/io/DataObjectManager.cpp`
- `tests/data/io/DataObjectManager_test.cpp`

## 3. Add a reusable operation in `DataObjectWorkflowOps`

Header pattern (`src/core/workflow/DataObjectWorkflowOps.hpp`):

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

Implementation rules (`src/core/workflow/DataObjectWorkflowOps.cpp`):

- keep behavior deterministic for the same model state and options
- keep function contracts typed (`ModelObject &`, `MapObject &`) instead of generic base pointers
- fail with explicit exceptions for invalid preconditions
- keep command-only side effects out of shared operations

## 4. Extend `ForEachDataObject(...)` safely

Current `ForEachDataObject` contract includes:

- mutable and const callback overloads
- optional key filtering by `key_tag_list`
- optional ordering via `IterateOptions::deterministic_order`
- snapshot-based traversal (callbacks run after snapshot capture)

When extending iteration behavior:

1. Keep default behavior backward-compatible.
2. Add new options to `DataObjectManager::IterateOptions` in header.
3. Apply option behavior inside snapshot-building path in `DataObjectManager.cpp`.
4. Preserve callback-empty validation and existing exception shape.
5. Preserve snapshot semantics (do not hold map lock during callback execution).

Example usage pattern at call site:

```cpp
rhbm_gem::DataObjectManager::IterateOptions options;
options.deterministic_order = true;

manager.ForEachDataObject(
    [](rhbm_gem::DataObjectBase & data_object)
    {
        if (auto * model = rhbm_gem::AsModelObject(data_object))
        {
            rhbm_gem::ModelPreparationOptions prep;
            prep.update_model = true;
            rhbm_gem::PrepareModelObject(*model, prep);
        }
    },
    {},
    options);
```

## 5. Add typed operation branches in iteration callbacks

Prefer `As*` helpers for optional branches:

- `AsModelObject(...)`
- `AsMapObject(...)`

Use `Expect*` helpers only when mismatch must terminate the workflow with a hard error:

- `ExpectModelObject(...)`
- `ExpectMapObject(...)`

Keep callback logic short; move multi-step workflows into typed helper functions.

## 6. Command integration pattern

When a command needs the new operation:

1. load/build typed objects through command data-loader helpers or `GetTypedDataObject<T>()`
2. run shared operation from `DataObjectWorkflowOps`
3. use manager iteration only when key-selection/ordering policy is owned by manager-level keys
4. keep command `ExecuteImpl()` orchestration-focused

## 7. Testing checklist

Add or extend tests for:

- operation happy-path behavior
- operation edge cases and invalid preconditions
- iteration order with and without explicit key list
- iteration snapshot behavior under map mutation (`ClearDataObjects()` during callback)
- callback-empty rejection behavior

Common files:

- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp`
- `tests/data/io/DataObjectManager_test.cpp`
- command-specific suites under `tests/core/command/`

## 8. Documentation checklist

When extension changes developer-facing contracts, update:

- `docs/developer/architecture/dataobject-typed-dispatch-architecture.md`
- this guide
- command docs when command behavior or flags changed
