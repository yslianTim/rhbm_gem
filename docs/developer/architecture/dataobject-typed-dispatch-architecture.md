# DataObject Typed Dispatch and Iteration Architecture

This document defines runtime contracts for typed dispatch and manager-level iteration.

Related guides:

- [`../development-guidelines.md`](../development-guidelines.md)
- [`./command-architecture.md`](./command-architecture.md)
- [`./dataobject-io-architecture.md`](./dataobject-io-architecture.md)
- [`../adding-dataobject-operations-and-iteration.md`](../adding-dataobject-operations-and-iteration.md)

## 1. Scope

Top-level dispatch targets in this architecture:

- `ModelObject`
- `MapObject`

`AtomObject` and `BondObject` are valid `DataObjectBase` subtypes but are not top-level catalog roots.

## 2. Dispatch API

`DataObjectDispatch` (`include/rhbm_gem/data/dispatch/DataObjectDispatch.hpp`) provides:

Probe helpers (non-throwing):

- `AsModelObject(...)`
- `AsMapObject(...)`

Expect helpers (throwing):

- `ExpectModelObject(...)`
- `ExpectMapObject(...)`

Catalog type helper:

- `GetCatalogTypeName(const DataObjectBase &)`

Dispatch contract:

- Probe helpers return typed pointer or `nullptr`.
- Expect helpers return typed reference or throw runtime error with caller context + resolved runtime type.
- `GetCatalogTypeName(...)` returns stable root names:
  - `model` for `ModelObject`
  - `map` for `MapObject`
- `GetCatalogTypeName(...)` throws for non-top-level types (`AtomObject`, `BondObject`, unresolved types).

## 3. DataObjectManager Iteration Contract

`DataObjectManager` supports mutable and const traversal:

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
  - `deterministic_order=true`: lexicographic key order.
  - `deterministic_order=false`: container iteration order.
- Non-empty `key_tag_list`: callback order follows input list order.
- Missing keys are skipped and logged as warnings.
- Empty callback throws runtime error.
- Iteration is snapshot-based; callbacks run after snapshot capture, so map mutation does not invalidate active traversal.

## 4. Integration Contracts

File/persistence type safety:

- `WriteDataObject(...)` model branch uses `ExpectModelObject(...)`.
- `WriteDataObject(...)` map branch uses `ExpectMapObject(...)`.
- `ModelObjectDAOSqlite::Save(...)` uses `ExpectModelObject(...)`.
- `MapObjectDAO::Save(...)` uses `ExpectMapObject(...)`.

Typed workflow operations (`DataObjectWorkflowOps`):

- `NormalizeMapObject`
- `PrepareModelObject`
- `ApplyModelSelection`
- `CollectModelAtoms`
- `PrepareSimulationAtoms`
- `BuildModelAtomBondContext`

Map sampling helper:

- `SampleMapValues(...)` (`include/rhbm_gem/core/command/MapSampling.hpp`)

Painter ingestion contract:

- `AddDataObject(...)` / `AddReferenceDataObject(...)` template helpers in `PainterIngestionInternal.hpp`
- Reject null or mismatched type with runtime error.

## 5. Extension Guidance

For extending behavior on existing `ModelObject` / `MapObject` workflows and iteration:

- use [`../adding-dataobject-operations-and-iteration.md`](../adding-dataobject-operations-and-iteration.md)

## 6. Key Files

Core interfaces:

- `include/rhbm_gem/data/object/DataObjectBase.hpp`
- `include/rhbm_gem/data/object/ModelObject.hpp`
- `include/rhbm_gem/data/object/MapObject.hpp`

Manager iteration:

- `include/rhbm_gem/data/io/DataObjectManager.hpp`
- `src/data/io/DataObjectManager.cpp`

Typed dispatch:

- `include/rhbm_gem/data/dispatch/DataObjectDispatch.hpp`
- `src/data/dispatch/DataObjectDispatch.cpp`

Typed workflow helpers:

- `src/core/workflow/DataObjectWorkflowOps.hpp`
- `src/core/workflow/DataObjectWorkflowOps.cpp`
- `include/rhbm_gem/core/command/MapSampling.hpp`
- `src/core/command/MapSampling.cpp`
- `src/core/internal/PainterIngestionInternal.hpp`

Regression tests:

- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp`
