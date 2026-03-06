# Extending DataObject Typed Dispatch Features

This guide describes the implementation workflow for extending `DataObject` typed dispatch behavior.
It focuses on concrete code changes and follows current repository contracts.

Related references:

- [`architecture/dataobject-typed-dispatch-architecture.md`](architecture/dataobject-typed-dispatch-architecture.md)
- [`architecture/dataobject-io-architecture.md`](architecture/dataobject-io-architecture.md)
- [`development-guidelines.md`](development-guidelines.md)

## 1. Decide extension scope first

Choose one primary scope before coding:

- extend typed workflow on existing top-level objects (`ModelObject`, `MapObject`)
- add a new top-level `DataObjectBase` subtype with catalog/file/persistence support
- extend painter ingestion checks for additional typed input contracts

Prefer extending existing typed workflows when your feature is model/map-specific and does not require a new catalog root type.

## 2. Required change set by scenario

### 2.1 Extend workflow on existing top-level types

Update these areas:

- `src/core/DataObjectWorkflowOps.hpp`
- `src/core/DataObjectWorkflowOps.cpp`
- command-local workflow files when behavior is command-specific
- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp` and relevant command tests
- architecture/developer docs if public workflow contract changes

### 2.2 Add a new top-level `DataObject` subtype

Update these areas in one change:

- new subtype definition/implementation under `include/data/` and `src/data/`
- `include/data/DataObjectDispatch.hpp` and `src/data/DataObjectDispatch.cpp`
- DAO implementation and registration (`DataObjectDAORegistrar`)
- persistence schema wiring (`ManagedStoreRegistry`, schema validation/migration contracts)
- `DatabaseSchemaManager` object catalog type checks
- file registry/factory wiring if file I/O is required
- tests for dispatch, schema, and I/O behavior
- developer docs and architecture references

## 3. Extend typed dispatch helpers

Dispatch entry points are centralized in `DataObjectDispatch`.

Header updates (`include/data/DataObjectDispatch.hpp`):

```cpp
class VolumeObject;

VolumeObject * AsVolumeObject(DataObjectBase & data_object) noexcept;
const VolumeObject * AsVolumeObject(const DataObjectBase & data_object) noexcept;

const VolumeObject & ExpectVolumeObject(
    const DataObjectBase & data_object,
    std::string_view context);
VolumeObject & ExpectVolumeObject(
    DataObjectBase & data_object,
    std::string_view context);
```

Implementation updates (`src/data/DataObjectDispatch.cpp`):

- add runtime name resolution in the internal type-name helper
- implement `As*` via `dynamic_cast`
- implement `Expect*` with context-rich `std::runtime_error`
- update `GetCatalogTypeName(...)` if this is a new top-level catalog type

Error style should stay consistent with existing helpers:

- `"<context>: expected <Type> but got <RuntimeType>."`

## 4. Add persistence wiring for top-level catalog types

### 4.1 Add DAO and stable catalog name

Create DAO files and register a stable persisted type name:

```cpp
namespace
{
    rhbm_gem::DataObjectDAORegistrar<rhbm_gem::VolumeObject, rhbm_gem::VolumeObjectDAO>
        registrar_volume_dao("volume");
}
```

Where to update:

- `include/data/<NewObject>DAO.hpp`
- `src/data/<NewObject>DAO.cpp`
- `include/data/DataObjectDAOFactoryRegistry.hpp` (usually no API changes, but contract reference)

### 4.2 Add managed store descriptor

Register store-level schema/list-key callbacks in:

- `src/data/persistence/ManagedStoreRegistry.cpp`

Ensure the descriptor defines:

- `object_type` (must match DAO registration stable name)
- `managed_table_names`
- `ensure_schema_v2`
- `validate_schema_v2`
- `list_keys`

### 4.3 Update catalog validation rules

`DatabaseSchemaManager` validates allowed `object_catalog.object_type` values.
When adding a new top-level type, update:

- object catalog table creation check constraint
- v2 validation query that rejects unsupported object types
- migration/bootstrap logic if existing databases must be upgraded safely

File:

- `src/data/DatabaseSchemaManager.cpp`

## 5. Add file I/O wiring when needed

If the new top-level object supports file import/export, update:

- `include/data/FileFormatRegistry.hpp`
- `src/data/FileFormatRegistry.cpp`
- `include/data/FileProcessFactoryBase.hpp`
- `src/data/FileProcessFactoryResolver.cpp`
- object-specific file reader/writer/factory classes

Keep extension-to-kind/backend mapping centralized in `FileFormatRegistry`.
Factory resolver behavior should remain deterministic and exception-driven on unsupported types.

## 6. Extend workflow and ingestion boundaries

Use the narrowest extension point:

- cross-command reusable domain logic -> `DataObjectWorkflowOps`
- command-specific logic -> command-local `*Workflow.cpp`
- painter ingest type contracts -> `src/core/PainterIngestionInternal.hpp`

For manager iteration consumers:

- use `DataObjectManager::ForEachDataObject(...)` for key-driven batch traversal
- use `As*` for optional typed branches
- use `Expect*` when mismatch must fail immediately with explicit context

## 7. Testing updates

At minimum, add/extend tests for:

- dispatch success and mismatch failure (`As*`, `Expect*`, `GetCatalogTypeName`)
- manager iteration ordering and snapshot semantics
- DAO save/load for the new stable catalog name (if top-level)
- schema bootstrap/validation/migration contracts when catalog values change
- file read/write dispatch for new formats (if file I/O is added)

Common suites to update:

- `tests/data/io/DataObjectDispatchIterationArchitecture_test.cpp`
- `tests/data/io/DataObjectManager_test.cpp`
- `tests/data/schema/DataObjectSchemaBootstrap_test.cpp`
- `tests/data/schema/DataObjectSchemaValidation_test.cpp`
- `tests/data/schema/DataObjectSchemaMigration_test.cpp`
- `tests/data/schema/DataObjectIOContract_test.cpp`

## 8. Documentation updates

When extension changes public behavior, update in the same change:

- `docs/developer/architecture/dataobject-typed-dispatch-architecture.md`
- `docs/developer/architecture/dataobject-io-architecture.md` (if I/O or persistence changed)
- this guide when extension workflow rules change
