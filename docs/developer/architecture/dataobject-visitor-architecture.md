# DataObject Traversal and Typed Operations Manual

This manual documents the current post-refactor architecture for `DataObject` traversal and runtime type handling.

Use this document when you need to:

- traverse `DataObject` instances via visitor
- select `ModelVisitMode` for model traversal
- implement command workflows and painter ingestion
- perform runtime type dispatch without visitor plumbing

Read with:

- [`./command-architecture.md`](./command-architecture.md)
- [`./dataobject-io-architecture.md`](./dataobject-io-architecture.md)

## 1. Design Boundary

`DataObjectVisitor` is now a traversal extension point only.

The following areas are intentionally typed APIs (not visitor-based):

- workflow orchestration
- map sampling/inference utilities
- painter ingestion and reference ingestion
- runtime type dispatch helpers

## 2. Core Contracts

### 2.1 `DataObjectBase`

All `DataObject` types implement only:

```cpp
virtual void Accept(DataObjectVisitor & visitor) = 0;
virtual void Accept(ConstDataObjectVisitor & visitor) const = 0;
```

There is no base-level `Accept(..., ModelVisitMode)` overload.

### 2.2 `ModelObject` Traversal API

`ModelObject` provides explicit mode-aware traversal:

```cpp
void Traverse(DataObjectVisitor & visitor, ModelVisitMode mode);
void Traverse(ConstDataObjectVisitor & visitor, ModelVisitMode mode) const;
```

Compatibility behavior:

- `ModelObject::Accept(visitor)` preserves legacy `AtomsThenSelf`

`ModelVisitMode` order:

| Mode | Visit order |
| --- | --- |
| `AtomsThenSelf` | all atoms -> model |
| `BondsThenSelf` | all bonds -> model |
| `AtomsAndBondsThenSelf` | all atoms -> all bonds -> model |
| `SelfOnly` | model only |

### 2.3 Non-model Objects

`AtomObject`, `BondObject`, and `MapObject` only expose plain `Accept(visitor)` and do not implement mode overloads.

## 3. Manager Traversal Policy

`DataObjectManager::Accept(..., VisitOptions)` keeps traversal policy centralization:

- `VisitOptions.model_visit_mode` is applied only when object dynamic type is `ModelObject` by calling `Traverse(...)`
- all other object types use plain `Accept(visitor)`
- deterministic key ordering and caller-specified key ordering behavior is unchanged

## 4. Typed Runtime Dispatch

`DataObjectDispatch` now uses centralized `dynamic_cast` helpers internally.

Public API remains:

- `ExpectModelObject(...)`
- `ExpectMapObject(...)`
- `GetCatalogTypeName(...)`

Error semantics are preserved, including `expected X but got Y` style messages.

## 5. Typed Workflow Operations

Visitor-based workflow classes were replaced by typed operation functions.

### 5.1 DataObject workflow ops

`src/core/DataObjectWorkflowOps.*` exports:

- `NormalizeMapObject`
- `PrepareModelObject`
- `ApplyModelSelection`
- `CollectModelAtoms`
- `PrepareSimulationAtoms`
- `BuildModelAtomBondContext`

Options/data contracts retained:

- `ModelPreparationOptions`
- `ModelAtomCollectorOptions`
- `SimulationAtomPreparationOptions`
- `SimulationAtomPreparationResult`
- `ModelAtomBondContext`

### 5.2 Potential analysis workflow ops

`src/core/PotentialAnalysisWorkflowOps.*` exports:

- `RunAtomSampling`
- `RunAtomGrouping`
- `RunLocalAtomFitting`
- `RunBondSampling`
- `RunBondGrouping`

Command/bond workflow callers now invoke these functions directly.

## 6. Stateless Map Sampling

Visitor/stateful map interpolation was replaced by:

```cpp
std::vector<SamplingDataTuple> SampleMapValues(
    const MapObject & map,
    const SamplerBase & sampler,
    const std::array<double, 3> & position,
    const std::array<double, 3> & axis_vector);
```

Location:

- `include/core/MapSampling.hpp`
- `src/core/MapSampling.cpp`

`MapInterpolationVisitor` and `MapSamplingWorkflow` are removed.

## 7. Painter Ingestion Boundary

Built-in painters no longer inherit `DataObjectVisitor`.

Affected painters:

- `AtomPainter`
- `ModelPainter`
- `GausPainter`
- `ComparisonPainter`
- `DemoPainter`

`AddDataObject(...)` and `AddReferenceDataObject(...)` now route through typed ingest helpers in `src/core/PainterIngestionInternal.hpp` using `dynamic_cast` and preserving existing null/type mismatch error semantics.

## 8. Extension Guidance

### 8.1 When to use visitor

Use visitor only when you truly need polymorphic traversal over heterogeneous `DataObject` trees.

### 8.2 When not to use visitor

Prefer typed API for:

- command-local pipelines
- algorithms with single supported input type
- painter ingestion and validation
- type expectation checks and error shaping

### 8.3 Adding a new `DataObject` subtype

1. Add concrete subtype derived from `DataObjectBase`.
2. Implement only const/non-const `Accept(visitor)`.
3. Add visit methods to both visitor interfaces if traversal support is needed.
4. Extend typed dispatch helpers/ops where relevant.
5. Update tests and this manual.

## 9. Key Files

Core traversal contracts:

- `include/data/DataObjectBase.hpp`
- `include/data/DataObjectVisitor.hpp`
- `include/data/ModelObject.hpp`
- `src/data/ModelObject.cpp`

Traversal orchestration:

- `include/core/DataObjectManager.hpp`
- `src/core/DataObjectManager.cpp`

Typed dispatch and ops:

- `include/data/DataObjectDispatch.hpp`
- `src/data/DataObjectDispatch.cpp`
- `src/core/DataObjectWorkflowOps.hpp`
- `src/core/DataObjectWorkflowOps.cpp`
- `src/core/PotentialAnalysisWorkflowOps.hpp`
- `src/core/PotentialAnalysisWorkflowOps.cpp`
- `include/core/MapSampling.hpp`
- `src/core/MapSampling.cpp`
- `src/core/PainterIngestionInternal.hpp`

Regression tests:

- `tests/DataObjectVisitorArchitecture_test.cpp`
