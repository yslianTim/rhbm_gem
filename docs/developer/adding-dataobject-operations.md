# Adding DataObject Operations

This guide covers how to extend behavior around the current top-level `DataObject` types:

- `ModelObject`
- `MapObject`

It also covers where command-side object loading belongs. There is no shared manager-owned iteration layer in the current architecture.

Related references:

- [`/docs/developer/architecture/dataobject-io-architecture.md`](/docs/developer/architecture/dataobject-io-architecture.md)
- [`/docs/developer/adding-a-command.md`](/docs/developer/adding-a-command.md)

## 1. Pick the Correct Boundary

Use the narrowest boundary that matches the change.

- reusable typed object behavior belongs in shared typed helpers or domain-layer code that works on `ModelObject` or `MapObject`
- file format support and extension dispatch belong in `FileIO` and the format codecs under `/src/data/io/file`
- database persistence behavior belongs in `DataRepository` and `/src/data/io/sqlite`
- command-specific loading, caching, logging, and orchestration belong in `/src/core/command/` and `/src/core/command/detail/`
- do not promote command-private cache helpers into shared data-layer APIs
- do not add a shared helper just to forward to a single file or repository call

`DataObjectDispatch`

- use when generic code needs to probe or enforce `DataObjectBase` runtime type
- keep command policy and persistence policy out of `DataObjectDispatch`

## 2. Required File Updates

If you add or change reusable typed operations:

- the relevant shared helper or domain file under `/src/data/object/` or `/src/utils/domain/`
- `/tests/data/DataObjectRuntimeBehavior_test.cpp`
- related command tests under `/tests/core/`
- docs when the contract changes

If you change file I/O behavior:

- `/include/rhbm_gem/data/io/FileIO.hpp`
- `/src/data/io/file/FileIO.cpp`
- format-specific codec files under `/src/data/io/file/`
- `/tests/data/DataObjectFileIO_test.cpp`
- `/tests/data/DataObjectImportRegression_test.cpp` when parser edge cases change

If you change database persistence behavior:

- `/include/rhbm_gem/data/io/DataRepository.hpp`
- `/src/data/io/DataRepository.cpp`
- `/src/data/io/sqlite/SQLitePersistence.*`
- `/src/data/io/sqlite/ModelObjectStorage.*`
- `/tests/data/DataObjectPersistence_test.cpp`
- `/tests/data/DataObjectSchemaBootstrap_test.cpp`
- `/tests/data/DataObjectSchemaCompatibility_test.cpp`
- `/tests/data/DataObjectSchemaValidation_test.cpp`

If you change command-side loading or persistence orchestration:

- `/src/core/command/detail/CommandBase.hpp`
- the relevant command implementation under `/src/core/command/`
- related command tests under `/tests/core/command/`

## 3. Operation Design Rules

For shared typed helpers:

- keep signatures typed: `ModelObject &`, `MapObject &`, or typed result structs
- keep behavior deterministic for the same object state and options
- validate preconditions and throw explicit `std::runtime_error` when violated
- keep command execution, logging policy, and UI concerns out of shared helpers
- prefer named workflow helpers when only a few combinations are valid in real commands
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

## 4. Typed Dispatch Inside Generic Code

Use probe helpers for optional branches:

- `AsModelObject(...)`
- `AsMapObject(...)`

Use expect helpers when a mismatch is a contract violation:

- `ExpectModelObject(...)`
- `ExpectMapObject(...)`

Top-level persistence routing is source-private; do not add new public dispatch helpers for catalog naming or repository/file routing.

If callback bodies or command-local branches start growing, move the multi-step logic into typed helpers instead of adding more dispatch inline.

## 5. Command Integration Pattern

Typical command flow:

1. for file-backed input, call `LoadInputFile<T>(...)` from `CommandBase`, or `ReadModel(...)` / `ReadMap(...)` directly outside command classes
2. for database-backed input, call `AttachDataRepository(...)` and then `LoadPersistedObject<T>(...)`
3. when persisting a command-owned object, call `SaveStoredObject(key_tag, persisted_key)` or use `DataRepository::Save*` directly outside commands
4. wrap failures with command-specific context near the orchestration boundary
5. keep `ExecuteImpl()` and local workflow helpers focused on typed orchestration
6. use ordinary container iteration or typed workflow helpers when traversing loaded objects; there is no shared iteration service to extend

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
- command-level failure context when file or database loading fails

Common suites:

- `/tests/data/DataPublicSurface_test.cpp`
- `/tests/data/DataObjectFileIO_test.cpp`
- `/tests/data/DataObjectImportRegression_test.cpp`
- `/tests/data/DataObjectRuntimeBehavior_test.cpp`
- `/tests/data/DataObjectDispatchAndIngestion_test.cpp`
- `/tests/data/DataObjectPersistence_test.cpp`
- `/tests/data/DataObjectSchemaBootstrap_test.cpp`
- `/tests/data/DataObjectSchemaCompatibility_test.cpp`
- `/tests/data/DataObjectSchemaValidation_test.cpp`
- `/tests/core/command/`

## 7. Documentation Checklist

When contracts change, update:

- `/docs/developer/architecture/dataobject-io-architecture.md`
- this guide
- command docs if command behavior or options changed
