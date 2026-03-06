# DataObject Typed Dispatch and Iteration Architecture

This manual defines the current runtime contracts for `DataObject` iteration and typed dispatch.

Use this guide to:

- iterate manager-owned `DataObjectBase` instances safely
- dispatch `DataObjectBase` to concrete top-level types (`ModelObject`, `MapObject`)
- implement typed workflows at command, persistence, and painter boundaries

Related guides:

- [`../development-guidelines.md`](../development-guidelines.md)
- [`./command-architecture.md`](./command-architecture.md)
- [`./dataobject-io-architecture.md`](./dataobject-io-architecture.md)
- [`../adding-dataobject-operations-and-iteration.md`](../adding-dataobject-operations-and-iteration.md)

## 1. Scope

Top-level `DataObject` roots in this architecture are:

- `ModelObject`
- `MapObject`

`AtomObject` and `BondObject` are valid `DataObjectBase` subtypes in runtime dispatch,
but they are treated as internal model-domain objects rather than top-level catalog roots.

## 2. Core Contracts

### 2.1 `DataObjectBase`

All concrete `DataObject` types implement:

```cpp
virtual std::unique_ptr<DataObjectBase> Clone() const = 0;
virtual void Display() const = 0;
virtual void Update() = 0;
virtual void SetKeyTag(const std::string & label) = 0;
virtual std::string GetKeyTag() const = 0;
```

### 2.2 `ModelObject`

Model-domain workflows are expressed through typed APIs such as:

- `GetAtomList()`
- `GetBondList()`
- `GetSelectedAtomList()`
- `GetSelectedBondList()`

### 2.3 `DataObjectManager` key-based typed access

`DataObjectManager::GetTypedDataObject<T>(key_tag)` provides key-based typed retrieval via
`std::dynamic_pointer_cast`, and throws on mismatch (`"Invalid data type for <key>"`).

## 3. Manager Iteration Contract

`DataObjectManager` exposes mutable and const callback traversal with optional ordering options:

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

- Empty `key_tag_list`:
  - `deterministic_order=true` sorts keys lexicographically before callback execution.
  - `deterministic_order=false` iterates current container order.
- Non-empty `key_tag_list`: caller order is preserved; missing keys are skipped with warning logs.
- Empty callback is invalid and throws (`ForEachDataObject(): callback is empty.`).
- Traversal is snapshot-based:
  - manager map access is locked only while building the snapshot vector
  - callbacks run outside the map lock
  - map mutations (for example `ClearDataObjects()`) do not invalidate the active traversal

## 4. Typed Dispatch Contract

`DataObjectDispatch` (`include/data/DataObjectDispatch.hpp`) provides runtime dispatch helpers.

### 4.1 Probe helpers (non-throwing)

- `AsModelObject(...)`
- `AsMapObject(...)`

These return typed pointers or `nullptr`.

### 4.2 Expect helpers (throwing)

- `ExpectModelObject(...)`
- `ExpectMapObject(...)`

These return typed references and throw `std::runtime_error` on mismatch. Error messages include
caller-provided context and resolved runtime type (for example `expected ModelObject but got AtomObject`).

### 4.3 Catalog naming helper

- `GetCatalogTypeName(const DataObjectBase &)`

Contract:

- returns `"model"` for `ModelObject`
- returns `"map"` for `MapObject`
- throws for `AtomObject`/`BondObject` and unresolved types

## 5. Integration Boundaries

### 5.1 File and persistence type checks

`Expect*` helpers are used as strict runtime contracts in I/O and persistence paths:

- `ModelObjectFactory::OutputDataObject(...)` -> `ExpectModelObject(...)`
- `MapObjectFactory::OutputDataObject(...)` -> `ExpectMapObject(...)`
- `ModelObjectDAOv2::Save(...)` -> `ExpectModelObject(...)`
- `MapObjectDAO::Save(...)` -> `ExpectMapObject(...)`

### 5.2 Command workflow helpers

Typed workflow operations live in `DataObjectWorkflowOps`:

- `NormalizeMapObject`
- `PrepareModelObject`
- `ApplyModelSelection`
- `CollectModelAtoms`
- `PrepareSimulationAtoms`
- `BuildModelAtomBondContext`

Map sampling is provided as a stateless helper:

- `SampleMapValues(...)` (`include/core/MapSampling.hpp`)

Potential-analysis-specific typed workflows remain command-local (`PotentialAnalysisCommandWorkflow.*`,
`PotentialAnalysisBondWorkflow.*`).

### 5.3 Painter ingestion checks

Painter ingestion uses `PainterIngestionInternal.hpp` template helpers:

- `AddDataObject(...)`
- `AddReferenceDataObject(...)`

These enforce expected types with `dynamic_cast` and throw mode-specific errors for null or mismatched input.

## 6. Extension Workflow

For existing `ModelObject`/`MapObject` operation and iteration extension, use:

- [`../adding-dataobject-operations-and-iteration.md`](../adding-dataobject-operations-and-iteration.md)

## 7. Key Files

Core interfaces:

- `include/data/DataObjectBase.hpp`
- `include/data/ModelObject.hpp`
- `include/data/MapObject.hpp`

Manager iteration:

- `include/core/DataObjectManager.hpp`
- `src/core/DataObjectManager.cpp`

Typed dispatch:

- `include/data/DataObjectDispatch.hpp`
- `src/data/DataObjectDispatch.cpp`

Typed workflow helpers:

- `src/core/DataObjectWorkflowOps.hpp`
- `src/core/DataObjectWorkflowOps.cpp`
- `include/core/MapSampling.hpp`
- `src/core/MapSampling.cpp`
- `src/core/PainterIngestionInternal.hpp`

Regression tests:

- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp`
