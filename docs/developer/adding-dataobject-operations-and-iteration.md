# Adding DataObject Operations and Iteration

This guide covers how to extend behavior around the current top-level `DataObject` types:

- `ModelObject`
- `MapObject`

It also covers `DataObjectManager` iteration and the shared command helpers that already exist in the repository.

Related references:

- [`architecture/dataobject-io-architecture.md`](architecture/dataobject-io-architecture.md)
- [`adding-a-command.md`](adding-a-command.md)

## 1. Pick the Correct Boundary

Use the narrowest boundary that matches the change.

`command_data_loader` helpers in `src/core/command/CommandDataSupport.hpp`

- use for typed file/database loading with consistent error context
- current helpers are `ProcessModelFile(...)`, `ProcessMapFile(...)`, `OptionalProcessMapFile(...)`, and `LoadModelObject(...)`

`CommandDataSupport` in `src/core/command/CommandDataSupport.*`

- use for reusable typed operations on `ModelObject` or `MapObject`
- current shared operations are `NormalizeMapObject(...)`, `PrepareModelObject(...)`, `ApplyModelSelection(...)`, `CollectModelAtoms(...)`, `PrepareSimulationAtoms(...)`, and `BuildModelAtomBondContext(...)`
- keep this layer focused on logic shared by multiple commands

`MapSampling` in `include/rhbm_gem/core/command/MapSampling.hpp` and `src/core/command/MapSampling.cpp`

- use only for reusable map sampling behavior

command-local code in `src/core/command/*.cpp` or `src/core/workflow/*.cpp`

- use when behavior is specific to one command or one workflow path

`DataObjectManager::ForEachDataObject(...)`

- use only when selection by key, traversal order, or snapshot iteration belongs in the manager layer

`DataObjectDispatch`

- use when generic code needs to probe or enforce `DataObjectBase` runtime type
- keep command policy out of `DataObjectDispatch`

## 2. Required File Updates

If you add or change reusable typed operations:

- `src/core/command/CommandDataSupport.hpp`
- `src/core/command/CommandDataSupport.cpp`
- command call sites using the helper
- `tests/data/DataObjectRuntime_test.cpp`
- related command tests under `tests/core/`
- docs when the contract changes

If you add or change loader helpers in `command_data_loader`:

- `src/core/command/CommandDataSupport.hpp`
- command call sites
- tests that cover the new load path and failure context

If you change manager iteration behavior:

- `include/rhbm_gem/data/io/DataObjectManager.hpp`
- `src/data/io/DataObjectManager.cpp`
- `tests/data/DataObjectRuntime_test.cpp`
- docs when callback or ordering semantics change

## 3. Operation Design Rules

For shared typed helpers:

- keep signatures typed: `ModelObject &`, `MapObject &`, or typed result structs
- keep behavior deterministic for the same object state and options
- validate preconditions and throw explicit `std::runtime_error` when violated
- keep command execution, logging policy, and UI concerns out of shared helpers
- prefer small option structs over long parameter lists when behavior may grow

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

## 4. Iteration Rules

Current `ForEachDataObject(...)` contract:

- mutable and const overloads exist
- `key_tag_list` optionally filters traversal
- `IterateOptions::deterministic_order` only affects the empty-key-list case
- explicit key lists preserve caller order
- missing keys are skipped with warning logs
- empty callbacks throw `std::runtime_error`
- callbacks run after the manager finishes building a `shared_ptr` snapshot

When extending iteration:

1. preserve current defaults unless you are intentionally changing the contract
2. add options to `IterateOptions` in the header first
3. implement the option in the snapshot-building path in `DataObjectManager.cpp`
4. keep callback validation explicit
5. keep callbacks outside the manager map lock

## 5. Typed Dispatch Inside Generic Code

Use probe helpers for optional branches:

- `AsModelObject(...)`
- `AsMapObject(...)`

Use expect helpers when a mismatch is a contract violation:

- `ExpectModelObject(...)`
- `ExpectMapObject(...)`

Use `GetCatalogTypeName(...)` only for top-level persistence routing.

If callback bodies or command-local branches start growing, move the multi-step logic into typed helpers instead of adding more dispatch inline.

## 6. Command Integration Pattern

Typical command flow:

1. load typed objects with `GetTypedDataObject<T>()` or `command_data_loader` helpers
2. apply shared typed operations from `CommandDataSupport` when the behavior is reused
3. use `ForEachDataObject(...)` only when manager-owned key filtering or ordering matters
4. keep `ExecuteImpl()` focused on orchestration

## 7. Test Checklist

Add or update tests for:

- happy-path behavior
- invalid preconditions or type mismatches
- iteration order for default traversal
- iteration order for explicit key lists
- snapshot behavior under map mutation
- empty-callback rejection
- new loader helper failure context when you add a loader wrapper

Common suites:

- `tests/data/DataObjectRuntime_test.cpp`
- related command suites under `tests/core/`

## 8. Documentation Checklist

When contracts change, update:

- `docs/developer/architecture/dataobject-io-architecture.md`
- this guide
- command docs if command behavior or options changed
