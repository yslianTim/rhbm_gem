# Adding DataObject Operations

This guide covers how to extend behavior around the current top-level `DataObject` types:

- `ModelObject`
- `MapObject`

It also covers where file loading, database loading, and command-side orchestration belong in the current codebase.

Related references:

- [`/docs/developer/architecture/dataobject-io-architecture.md`](/docs/developer/architecture/dataobject-io-architecture.md)
- [`/docs/developer/adding-a-command.md`](/docs/developer/adding-a-command.md)

## 1. Pick the Correct Boundary

Use the narrowest boundary that matches the change.

- reusable typed object behavior belongs in shared typed helpers or domain-layer code that works on `ModelObject` or `MapObject`
- file format support and extension dispatch belong in `FileIO` and codecs under `/src/data/io/file`
- database persistence behavior belongs in `DataRepository` and `/src/data/io/sqlite`
- command-specific loading, caching, logging, and orchestration belong in `/src/core/command/` and `/src/core/command/detail/`
- do not promote `CommandObjectCache` or other command-private helpers into shared data-layer APIs
- do not add a public type-erased dispatch layer unless the actual public surface is changing
- do not add a shared helper that only forwards to one file I/O or repository call

## 2. Required File Updates

If you add or change reusable typed operations:

- the relevant shared helper or domain file under `/src/data/object/` or `/src/utils/domain/`
- `/tests/data/DataObjectRuntimeBehavior_test.cpp`
- related command tests under `/tests/core/`
- docs when the contract changes

If you change file I/O behavior:

- `/include/rhbm_gem/data/io/ModelMapFileIO.hpp`
- `/src/data/io/file/ModelMapFileIO.cpp`
- format-specific codec files under `/src/data/io/file/`
- `/tests/data/DataObjectFileIO_test.cpp`
- `/tests/data/DataObjectImportRegression_test.cpp` when parser edge cases change

If you change database persistence behavior:

- `/include/rhbm_gem/data/io/DataRepository.hpp`
- `/src/data/io/DataRepository.cpp`
- `/src/data/io/sqlite/SQLitePersistence.*`
- `/src/data/io/sqlite/ModelObjectStorage.*`
- `/tests/data/DataObjectPersistence_test.cpp`
- `/tests/data/SQLitePersistenceTypedApi_test.cpp`
- `/tests/data/DataObjectSchemaBootstrap_test.cpp`
- `/tests/data/DataObjectSchemaCompatibility_test.cpp`
- `/tests/data/DataObjectSchemaValidation_test.cpp`

If you change command-side loading or persistence orchestration:

- `/src/core/command/detail/CommandBase.hpp`
- `/src/core/command/detail/CommandObjectCache.hpp`
- the relevant command implementation under `/src/core/command/`
- related command tests under `/tests/core/command/`

## 3. Operation Design Rules

For shared typed helpers:

- keep signatures typed: `ModelObject &`, `MapObject &`, or typed result structs
- keep behavior deterministic for the same object state and options
- validate preconditions and throw explicit `std::runtime_error` when violated
- keep command execution, logging policy, and UI concerns out of shared helpers
- prefer named workflow helpers when only a few combinations are valid in real commands
- keep owner lookup, model analysis-store access, fit-state clearing, and spatial lookup behind `/src/data/detail/ModelAnalysisAccess.hpp`
- use `ModelAnalysisAccess` directly for local/group entry reads and writes; do not add new `Local*Access` or `Group*Access` forwarding helpers
- do not re-expose internal analysis helpers from `/include/rhbm_gem/data/object/**`
- avoid one-hop shared wrappers that only forward to:
  - `ReadModel(...)`
  - `ReadMap(...)`
  - `DataRepository::LoadModel(...)`
  - `DataRepository::LoadMap(...)`
  - `DataRepository::SaveModel(...)`
  - `DataRepository::SaveMap(...)`

Example pattern:

```cpp
std::vector<AtomObject *> CollectSelectedAtoms(ModelObject & model_object);
std::vector<AtomObject *> CollectAtomsWithLocalPotentialEntries(ModelObject & model_object);
```

## 4. Typed Integration Inside Generic Code

- prefer concrete typed methods over a shared runtime-dispatch layer
- if only `ModelObject` is supported, accept `ModelObject &` directly
- if only `MapObject` is supported, accept `MapObject &` directly
- if a command owns mixed top-level objects, keep type bookkeeping inside `CommandObjectCache`
- keep persistence routing private to `CommandBase` and repository internals; do not add public dispatch helpers just for catalog naming or repository/file routing
- if command-local branches start growing, move the multi-step logic into typed helpers instead of adding another type-erased facade

## 5. Command Integration Pattern

Typical command flow:

1. for file-backed input, call `LoadInputFile<T>(...)` from `CommandBase`, or `ReadModel(...)` / `ReadMap(...)` directly outside command classes
2. for database-backed input, call `AttachDataRepository(...)` and then `LoadPersistedObject<T>(...)`
3. when persisting a command-owned object, call `SaveStoredObject(key_tag, persisted_key)` or use `DataRepository::Save*` directly outside commands
4. let `CommandObjectCache` keep the loaded `shared_ptr` objects keyed by command-local `key_tag`
5. wrap failures with command-specific context near the orchestration boundary
6. keep `ExecuteImpl()` and local workflow helpers focused on typed orchestration

Repository-backed command request structs default `database_path` through `GetDefaultDatabasePath()` in `/include/rhbm_gem/core/command/CommandApi.hpp`, so command changes should preserve that behavior unless the command contract is intentionally changing.

## 6. Test Checklist

Add or update tests for:

- happy-path behavior
- invalid preconditions or type mismatches
- file format matrix and unsupported write behavior
- uppercase extension dispatch
- malformed file input
- database round-trip for model and map objects
- renamed persisted key semantics
- missing database key behavior
- schema bootstrap and invalid schema rejection
- command-cache type mismatch behavior when command-local routing changes
- command-level failure context when file or database loading fails

Common suites:

- `/tests/data/DataPublicSurface_test.cpp`
- `/tests/data/DataObjectFileIO_test.cpp`
- `/tests/data/DataObjectImportRegression_test.cpp`
- `/tests/data/DataObjectRuntimeBehavior_test.cpp`
- `/tests/data/DataObjectDispatchAndIngestion_test.cpp`
- `/tests/data/DataObjectPersistence_test.cpp`
- `/tests/data/SQLitePersistenceTypedApi_test.cpp`
- `/tests/data/DataObjectSchemaBootstrap_test.cpp`
- `/tests/data/DataObjectSchemaCompatibility_test.cpp`
- `/tests/data/DataObjectSchemaValidation_test.cpp`
- `/tests/core/command/CommandObjectCache_test.cpp`
- `/tests/core/command/`

## 7. Documentation Checklist

When contracts change, update:

- `/docs/developer/architecture/dataobject-io-architecture.md`
- this guide
- command docs if command behavior or options changed
