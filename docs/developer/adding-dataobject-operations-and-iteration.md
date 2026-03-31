# Adding DataObject Operations and Iteration

This guide covers how to extend behavior around the current top-level `DataObject` types:

- `ModelObject`
- `MapObject`

It also covers `DataObjectManager` iteration and the shared command helpers that already exist in the repository.

Related references:

- [`/docs/developer/architecture/dataobject-io-architecture.md`](/docs/developer/architecture/dataobject-io-architecture.md)
- [`/docs/developer/adding-a-command.md`](/docs/developer/adding-a-command.md)

## 1. Pick the Correct Boundary

Use the narrowest boundary that matches the change.

- use for reusable typed operations on `ModelObject` or `MapObject`
- current shared operations are `NormalizeMapObject(...)`, `PrepareModelForPotentialAnalysis(...)`, `PrepareModelForVisualization(...)`, `ApplyModelSelection(...)`, `CollectSelectedAtoms(...)`, `CollectAtomsWithLocalPotentialEntries(...)`, `PrepareSimulationAtoms(...)`, and `BuildModelAtomBondContext(...)`
- keep this layer focused on logic shared by multiple commands

`DataObjectManager::ForEachDataObject(...)`

- use only when selection by key, traversal order, or snapshot iteration belongs in the manager layer

`DataObjectDispatch`

- use when generic code needs to probe or enforce `DataObjectBase` runtime type
- keep command policy out of `DataObjectDispatch`

## 2. Required File Updates

If you add or change reusable typed operations:

- `/tests/data/DataObjectRuntime_test.cpp`
- related command tests under `/tests/core/`
- docs when the contract changes

If you change manager iteration behavior:

- `/include/rhbm_gem/data/io/DataObjectManager.hpp`
- `/src/data/io/DataObjectManager.cpp`
- `/tests/data/DataObjectRuntime_test.cpp`
- docs when callback or ordering semantics change

## 3. Operation Design Rules

For shared typed helpers:

- keep signatures typed: `ModelObject &`, `MapObject &`, or typed result structs
- keep behavior deterministic for the same object state and options
- validate preconditions and throw explicit `std::runtime_error` when violated
- keep command execution, logging policy, and UI concerns out of shared helpers
- prefer named workflow helpers when only a few combinations are valid in real commands

Example pattern:

```cpp
std::vector<AtomObject *> CollectSelectedAtoms(ModelObject & model_object);
std::vector<AtomObject *> CollectAtomsWithLocalPotentialEntries(ModelObject & model_object);
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

1. load typed objects with `GetTypedDataObject<T>()`
2. use `ForEachDataObject(...)` only when manager-owned key filtering or ordering matters
3. keep `ExecuteImpl()` focused on orchestration

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

- `/tests/data/DataObjectRuntime_test.cpp`

## 8. Documentation Checklist

When contracts change, update:

- `/docs/developer/architecture/dataobject-io-architecture.md`
- this guide
- command docs if command behavior or options changed
