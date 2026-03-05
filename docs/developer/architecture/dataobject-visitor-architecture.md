# DataObject Visitor Developer Manual

This manual describes how the current `DataObject` visitor system works and how to extend it safely.

Use this document when you need to:

- traverse `DataObject` instances through `Accept(...)`
- choose a `ModelVisitMode` for `ModelObject` traversal
- run visitors through `DataObjectManager`
- implement a new visitor or add a new `DataObject` subtype

Read with:

- [`./command-architecture.md`](./command-architecture.md)
- [`./dataobject-io-architecture.md`](./dataobject-io-architecture.md)

## 1. Scope

This document covers runtime dispatch and traversal only:

- `DataObjectBase`
- `DataObjectVisitor` and `ConstDataObjectVisitor`
- concrete objects: `AtomObject`, `BondObject`, `ModelObject`, `MapObject`
- manager traversal entry points: `DataObjectManager::Accept(...)`
- built-in workflow visitors used by commands

Persistence and file/database pipelines are out of scope.

## 2. Quick Start

Direct traversal from a concrete/top-level object:

```cpp
MyVisitor visitor;  // custom DataObjectVisitor implementation
model.Accept(visitor, ModelVisitMode::AtomsAndBondsThenSelf);
```

Manager-level traversal with explicit policy:

```cpp
DataObjectManager manager;
DataObjectManager::VisitOptions options;
options.deterministic_order = true;
options.model_visit_mode = ModelVisitMode::SelfOnly;

MyVisitor visitor;  // custom DataObjectVisitor implementation
manager.Accept(visitor, {}, options);
```

## 3. Core Contracts

### 3.1 `DataObjectBase`

All `DataObject` types must implement:

```cpp
virtual void Accept(DataObjectVisitor & visitor) = 0;
virtual void Accept(ConstDataObjectVisitor & visitor) const = 0;
virtual void Accept(DataObjectVisitor & visitor, ModelVisitMode model_mode) = 0;
virtual void Accept(ConstDataObjectVisitor & visitor, ModelVisitMode model_mode) const = 0;
```

Rules:

- visitor is passed by reference (non-null by type)
- `model_mode` is meaningful for `ModelObject`
- non-model objects accept `model_mode` but ignore it

### 3.2 Visitor Interfaces

Two pure-virtual contracts exist:

- `DataObjectVisitor` (mutable)
- `ConstDataObjectVisitor` (read-only)

Both require full coverage:

- `VisitAtomObject(...)`
- `VisitBondObject(...)`
- `VisitModelObject(...)`
- `VisitMapObject(...)`

There is no default no-op handler.

### 3.3 Dispatch Model

Dispatch is double dispatch:

1. caller holds `DataObjectBase` polymorphically
2. `Accept(...)` resolves by object dynamic type
3. concrete object calls `visitor.VisitXxxObject(*this)`
4. visitor implementation resolves by visitor dynamic type

## 4. Traversal Rules by Object Type

| Object type | `Accept(visitor)` behavior | `Accept(visitor, model_mode)` behavior |
| --- | --- | --- |
| `AtomObject` | `VisitAtomObject(*this)` | same as left (`model_mode` ignored) |
| `BondObject` | `VisitBondObject(*this)` | same as left (`model_mode` ignored) |
| `MapObject` | `VisitMapObject(*this)` | same as left (`model_mode` ignored) |
| `ModelObject` | equivalent to `AtomsThenSelf` | follows `ModelVisitMode` |

`ModelVisitMode` traversal order:

| Mode | Visit order |
| --- | --- |
| `AtomsThenSelf` | all atoms -> model |
| `BondsThenSelf` | all bonds -> model |
| `AtomsAndBondsThenSelf` | all atoms -> all bonds -> model |
| `SelfOnly` | model only |

Compatibility behavior: `model.Accept(visitor)` keeps legacy `AtomsThenSelf` behavior.

## 5. Manager Traversal (`DataObjectManager`)

Overloads exist for both mutable and const visitors:

```cpp
void Accept(DataObjectVisitor & visitor,
            const std::vector<std::string> & key_tag_list = {});

void Accept(DataObjectVisitor & visitor,
            const std::vector<std::string> & key_tag_list,
            const VisitOptions & options);

void Accept(ConstDataObjectVisitor & visitor,
            const std::vector<std::string> & key_tag_list = {}) const;

void Accept(ConstDataObjectVisitor & visitor,
            const std::vector<std::string> & key_tag_list,
            const VisitOptions & options) const;
```

`VisitOptions`:

- `bool deterministic_order = true`
- `ModelVisitMode model_visit_mode = ModelVisitMode::AtomsThenSelf`

Operational behavior:

- traversal uses a snapshot of `shared_ptr` entries built under lock
- actual `Accept(...)` calls run after lock release
- empty `key_tag_list` + `deterministic_order=true`: key-sorted traversal
- empty `key_tag_list` + `deterministic_order=false`: underlying `unordered_map` iteration order
- non-empty `key_tag_list`: preserves caller order
- missing keys emit warning logs and are skipped

## 6. Built-in Workflow Visitors

Current production visitors in `core` are:

- `MapInterpolationVisitor` (`ConstDataObjectVisitor`)
- `MapNormalizeVisitor` (`DataObjectVisitor`)
- `ModelPreparationVisitor` (`DataObjectVisitor`)
- `ModelSelectionVisitor` (`DataObjectVisitor`)
- `ModelSelectedAtomCollectorVisitor` (`DataObjectVisitor`)
- `SimulationAtomPreparationVisitor` (`DataObjectVisitor`)
- `ModelAtomBondContextVisitor` (`DataObjectVisitor`)
- `AtomSamplingVisitor` (`DataObjectVisitor`)
- `AtomGroupingVisitor` (`DataObjectVisitor`)
- `LocalFittingVisitor` (`DataObjectVisitor`)
- `BondSamplingVisitor` (`DataObjectVisitor`)
- `BondGroupingVisitor` (`DataObjectVisitor`)

### 6.1 `MapInterpolationVisitor`

Supported visit target:

- `VisitMapObject(const MapObject&)`

Unsupported targets:

- `VisitAtomObject(...)`, `VisitBondObject(...)`, `VisitModelObject(...)` throw `std::logic_error`

State and output APIs:

- `SetPosition(...)`
- `SetAxisVector(...)`
- `GetSamplingDataList()` for read-only access
- `ConsumeSamplingDataList()` for move-out transfer

Behavior notes:

- every `VisitMapObject(...)` clears previous sampling output first
- null sampler logs warning and returns empty output
- visitor instances are stateful; reuse intentionally and consume/reset output between runs

### 6.2 `MapNormalizeVisitor`

Supported visit target:

- `VisitMapObject(MapObject&)`

Behavior:

- calls `MapObject::MapValueArrayNormalization()`
- unsupported targets throw `std::logic_error`

### 6.3 `ModelPreparationVisitor`

Supported visit target:

- `VisitModelObject(ModelObject&)`

Behavior:

- option-driven preprocessing for model workflows
- supports select-all, optional symmetry filtering, `Update()`, and optional
  local-potential-entry initialization for selected atoms/bonds
- unsupported targets throw `std::logic_error`

### 6.4 `MapSamplingWorkflow`

`MapSamplingWorkflow` is a helper wrapper that reuses an internal
`MapInterpolationVisitor` instance.

Behavior:

- sets position/axis vector
- executes `map.Accept(visitor)`
- returns sampled tuples via `ConsumeSamplingDataList()`

### 6.5 Painter Ingestion Callbacks

Built-in painters (`AtomPainter`, `ModelPainter`, `GausPainter`,
`ComparisonPainter`, `DemoPainter`) implement `DataObjectVisitor` for
ingestion dispatch.

Usage rule:

- call `AddDataObject(...)` / `AddReferenceDataObject(...)` from command workflows
- treat `Visit*Object(...)` overrides as callback-only implementation details

## 7. Extension Playbook

### 7.1 Add a New Visitor

1. Derive from `DataObjectVisitor` or `ConstDataObjectVisitor`.
2. Implement all required `Visit*Object(...)` methods.
3. Keep mutable state explicit and reset/consume per logical run.
4. Add/adjust tests for traversal order and visited-type behavior.

### 7.2 Add a New `DataObject` Type

1. Add new concrete type derived from `DataObjectBase`.
2. Implement all four `Accept(...)` overloads.
3. Add new visit methods to both visitor interfaces.
4. Update all visitor implementations and compile-failing call sites.
5. Add/adjust architecture tests and docs.

## 8. Key Files

Core contracts:

- `include/data/DataObjectBase.hpp`
- `include/data/DataObjectVisitor.hpp`
- `include/data/ModelVisitMode.hpp`

Concrete `Accept(...)` behavior:

- `src/data/AtomObject.cpp`
- `src/data/BondObject.cpp`
- `src/data/ModelObject.cpp`
- `src/data/MapObject.cpp`

Manager traversal:

- `include/core/DataObjectManager.hpp`
- `src/core/DataObjectManager.cpp`

Active visitor and call sites:

- `include/core/MapInterpolationVisitor.hpp`
- `src/core/MapInterpolationVisitor.cpp`
- `src/core/DataObjectWorkflowVisitors.hpp`
- `src/core/DataObjectWorkflowVisitors.cpp`
- `src/core/PotentialAnalysisWorkflowVisitors.hpp`
- `src/core/PotentialAnalysisWorkflowVisitors.cpp`
- `src/core/PotentialAnalysisCommand.cpp`
- `src/core/PotentialAnalysisBondWorkflow.cpp`
- `src/core/MapVisualizationCommand.cpp`

Architecture regression tests:

- `tests/DataObjectVisitorArchitecture_test.cpp`
