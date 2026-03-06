# DataObject Iteration and Typed Dispatch Manual

This manual documents the current post-refactor architecture for `DataObject`
iteration and runtime type handling.

Use this document when you need to:

- iterate manager-owned `DataObject` batches
- perform runtime type checks for top-level objects (`ModelObject`, `MapObject`)
- implement command workflows and painter ingestion with typed APIs

Read with:

- [`./command-architecture.md`](./command-architecture.md)
- [`./dataobject-io-architecture.md`](./dataobject-io-architecture.md)

## 1. Design Boundary

`DataObjectVisitor` and model traversal mode plumbing were removed.

The architecture now uses:

- manager-owned callback iteration (`DataObjectManager::ForEachDataObject`)
- centralized typed runtime dispatch (`DataObjectDispatch`)
- command-local typed workflow APIs (`DataObjectWorkflowOps`, `MapSampling`, ...)

## 2. Core Contracts

### 2.1 `DataObjectBase`

All `DataObject` types now implement only:

```cpp
virtual std::unique_ptr<DataObjectBase> Clone() const = 0;
virtual void Display() const = 0;
virtual void Update() = 0;
virtual void SetKeyTag(const std::string & label) = 0;
virtual std::string GetKeyTag() const = 0;
```

There is no visitor `Accept(...)` contract in `DataObjectBase`.

### 2.2 `ModelObject`

`ModelObject` no longer exposes framework-level `Traverse(...)` overloads.

Atom/bond operations are handled through typed model APIs:

- `GetAtomList()`
- `GetBondList()`
- `GetSelectedAtomList()`
- `GetSelectedBondList()`

### 2.3 Manager Batch Iteration

`DataObjectManager` provides callback-based traversal:

```cpp
struct IterateOptions
{
    bool deterministic_order{ true };
};

void ForEachDataObject(
    const std::function<void(DataObjectBase &)> & callback,
    const std::vector<std::string> & key_tag_list = {},
    const IterateOptions & options = {});

void ForEachDataObject(
    const std::function<void(const DataObjectBase &)> & callback,
    const std::vector<std::string> & key_tag_list = {},
    const IterateOptions & options = {}) const;
```

Behavior:

- default traversal is deterministic by key sort when `key_tag_list` is empty
- when `key_tag_list` is provided, caller order is preserved
- iteration is snapshot-based, so concurrent map mutation does not invalidate traversal

## 3. Typed Runtime Dispatch

`DataObjectDispatch` provides two dispatch layers:

Non-throwing probe helpers:

- `AsModelObject(...)`
- `AsMapObject(...)`

Throwing expectation helpers:

- `ExpectModelObject(...)`
- `ExpectMapObject(...)`

Catalog naming helper:

- `GetCatalogTypeName(...)`

Error semantics for expectation helpers are preserved in `expected X but got Y` style.

## 4. Typed Workflow Operations

Workflow logic remains typed and visitor-free.

`src/core/DataObjectWorkflowOps.*` exports:

- `NormalizeMapObject`
- `PrepareModelObject`
- `ApplyModelSelection`
- `CollectModelAtoms`
- `PrepareSimulationAtoms`
- `BuildModelAtomBondContext`

Potential-analysis-specific workflow helpers are command-local:

- atom workflow helpers live in `src/core/PotentialAnalysisCommandWorkflow.cpp`
- experimental bond workflow helpers live in `src/core/PotentialAnalysisBondWorkflow.cpp`

`include/core/MapSampling.hpp` provides stateless map sampling through `SampleMapValues(...)`.

## 5. Painter Ingestion Boundary

Built-in painters ingest typed objects via `PainterIngestionInternal.hpp`:

- `AtomPainter`
- `ModelPainter`
- `GausPainter`
- `ComparisonPainter`
- `DemoPainter`

`AddDataObject(...)` and `AddReferenceDataObject(...)` use `dynamic_cast`-based typed checks and preserve existing null/type mismatch error semantics.

## 6. Extension Guidance

### 6.1 When to use manager iteration

Use `DataObjectManager::ForEachDataObject(...)` when ownership of batch selection and ordering belongs to manager-level key policies.

### 6.2 When to use direct typed access

Prefer typed APIs for:

- command-local pipelines
- algorithms with known input types
- atom/bond workflows inside `ModelObject`
- runtime expectation checks with explicit error shaping

### 6.3 Adding a new top-level `DataObject` subtype

1. Add concrete subtype derived from `DataObjectBase`.
2. Extend `DataObjectDispatch` with probe/expect helpers if needed.
3. Extend I/O and persistence layers (`FileFormatRegistry`, DAO/managed store descriptors) where relevant.
4. Add tests for typed dispatch, manager iteration behavior, and persistence round-trip.
5. Update this manual.

## 7. Key Files

Core contracts:

- `include/data/DataObjectBase.hpp`
- `include/data/ModelObject.hpp`

Manager iteration:

- `include/core/DataObjectManager.hpp`
- `src/core/DataObjectManager.cpp`

Typed dispatch and ops:

- `include/data/DataObjectDispatch.hpp`
- `src/data/DataObjectDispatch.cpp`
- `src/core/DataObjectWorkflowOps.hpp`
- `src/core/DataObjectWorkflowOps.cpp`
- `include/core/MapSampling.hpp`
- `src/core/MapSampling.cpp`
- `src/core/PainterIngestionInternal.hpp`

Regression tests:

- `tests/DataObjectDispatchIterationArchitecture_test.cpp`
