# DataObject I/O Developer Guide

This guide describes the current file and database I/O architecture for top-level `DataObject` instances.

Use it when you need to:

- trace how a file becomes an in-memory object
- trace how a `DataObject` is saved to or loaded from SQLite
- understand the normalized `ModelObject` persistence schema
- add a new file format or a new persistent top-level `DataObject`

Read this together with:

- [`../development-guidelines.md`](../development-guidelines.md)
- [`./command-architecture.md`](./command-architecture.md)

## 1. Scope

This subsystem currently treats only these types as top-level I/O units:

- `ModelObject`
- `MapObject`

`AtomObject` and `BondObject` inherit from `DataObjectBase`, but they are subordinate parts of `ModelObject`. They are not registered as standalone file or database I/O targets.

## 2. Main Responsibilities

- `DataObjectManager`: command-facing orchestration and in-memory ownership
- `FileFormatRegistry`: single source of truth for supported file extensions and backend dispatch
- `FileProcessFactoryRegistry`: optional override layer for choosing the top-level factory (`ModelObjectFactory` or `MapObjectFactory`) from the normalized extension
- `ModelObjectFactory` / `MapObjectFactory`: build or write the concrete object
- `ModelFileReader` / `ModelFileWriter`: select model parser or writer backend from `FileFormatRegistry`
- `MapFileReader` / `MapFileWriter`: select map parser or writer backend from `FileFormatRegistry`
- `DatabaseManager`: own the SQLite connection, schema bootstrap or migration, and cross-table transaction boundaries
- `DatabaseSchemaManager`: schema-version detection, normalized v2 bootstrap, and legacy v1 to v2 migration
- `DataObjectDAOFactoryRegistry`: map runtime types and stored type names to DAO implementations
- `ModelObjectDAO` (`ModelObjectDAOv2`): normalized v2 persistence for `ModelObject`
- `LegacyModelObjectDAO`: legacy v1 read path used only for migration
- `MapObjectDAO`: flat persistence for `MapObject`
- `FileFormatBackendFactory`: shared helper that builds concrete parser or writer backends from `FileFormatRegistry`
- `SQLiteWrapper`: SQL execution, prepared statements, typed bind or read helpers, and RAII transaction handling

## 3. Recommended Call Path

The command layer should normally talk to `DataObjectManager`, not directly to readers, writers, or DAOs.

Typical file import:

```cpp
DataObjectManager manager;
manager.ProcessFile(model_path, "model");
auto model = manager.GetTypedDataObject<ModelObject>("model");
```

Typical database round-trip:

```cpp
DataObjectManager manager;
manager.SetDatabaseManager(database_path);
manager.ProcessFile(model_path, "model");
manager.SaveDataObject("model", "saved_model");

manager.LoadDataObject("saved_model");
auto model = manager.GetTypedDataObject<ModelObject>("saved_model");
```

Practical rules:

- Commands should prefer `DataObjectManager` or helpers in `src/core/CommandDataAccessInternal.hpp`.
- `SetDatabaseManager(...)` must be called before `SaveDataObject(...)` or `LoadDataObject(...)`.
- `ProduceFile(...)` only writes objects already present in memory.
- Reusing an in-memory `key_tag` replaces the previous entry in `DataObjectManager`.
- `SaveDataObject(original, renamed)` changes only the persisted database key. It does not rename the in-memory object or the in-memory map entry.

## 4. Runtime Topology

```mermaid
flowchart LR
    A["Commands / tests"] --> B["DataObjectManager"]

    subgraph F["File I/O"]
        B --> C["FilePathHelper::GetExtension"]
        C --> D["FileProcessFactoryRegistry"]
        D --> E["ModelObjectFactory"]
        D --> F2["MapObjectFactory"]
        E --> G["ModelFileReader / ModelFileWriter"]
        F2 --> H["MapFileReader / MapFileWriter"]
        G --> I["FileFormatRegistry"]
        H --> I
        I --> J["Pdb/Cif or Mrc/CCP4 backend"]
        J --> K["AtomicModelDataBlock or MapObject payload"]
    end

    subgraph P["SQLite persistence"]
        B --> L["DatabaseManager"]
        L --> M["DatabaseSchemaManager"]
        L --> N["DataObjectDAOFactoryRegistry"]
        N --> O["ModelObjectDAOv2"]
        N --> P2["MapObjectDAO"]
        M --> Q["LegacyModelObjectDAO (migration only)"]
        O --> R["SQLiteWrapper"]
        P2 --> R
        Q --> R
    end
```

## 5. File I/O

### 5.1 Dispatch Model

All public file entry points start in `DataObjectManager`:

- `ProcessFile(path, key_tag)`
- `ProduceFile(path, key_tag)`

Extension normalization happens in `FilePathHelper::GetExtension(...)`. The resulting lowercase extension is then used consistently by the rest of the file pipeline.

`FileFormatRegistry` is the single source of truth for supported file formats. Each `FileFormatDescriptor` records:

- normalized extension
- top-level object kind (`Model` or `Map`)
- whether read is supported
- whether write is supported
- concrete backend selection (`ModelFormatBackend` or `MapFormatBackend`)

Current descriptors:

| Kind | Read | Write |
| --- | --- | --- |
| model | `.pdb`, `.cif`, `.mmcif`, `.mcif` | `.pdb`, `.cif` |
| map | `.mrc`, `.map`, `.ccp4` | `.mrc`, `.map`, `.ccp4` |

`FileProcessFactoryRegistry` is no longer the source of default format behavior. Its runtime contract is:

1. if an explicit override was registered with `RegisterFactory(...)`, use that override
2. otherwise ask `FileFormatRegistry` for the descriptor and instantiate the default top-level factory from `descriptor.kind`

This means normal application startup no longer needs `RegisterDefaultFactories()` or global one-time registration just to make built-in formats work.

Readers and writers then use `FileFormatBackendFactory` plus `FileFormatRegistry` to choose the concrete backend implementation.

### 5.2 `ProcessFile(...)`

`DataObjectManager::ProcessFile(...)`:

1. gets the normalized extension
2. asks `FileProcessFactoryRegistry` for a factory
3. calls `CreateDataObject(...)`
4. sets the resulting object's `key_tag`
5. inserts or replaces it in the in-memory map

If parsing fails, `DataObjectManager` wraps the lower-level exception with the file path and requested `key_tag`.

### 5.3 `ProduceFile(...)`

`DataObjectManager::ProduceFile(...)`:

1. looks up the in-memory object by `key_tag`
2. selects a factory from the output extension
3. delegates to `OutputDataObject(...)`

Current behavior to keep in mind:

- missing in-memory keys do not throw; the manager logs a warning and returns
- output format is selected from the target filename, not from the source filename

## 6. Model File Pipeline

### 6.1 Read Path

Model import is split into three layers:

1. `ModelFileReader`
2. `ModelFileFormatBase` implementations (`PdbFormat`, `CifFormat`)
3. `ModelObjectFactory`

`ModelFileReader` asks `FileFormatRegistry` for the backend:

- `.pdb` -> `ModelFormatBackend::Pdb`
- `.cif`, `.mmcif`, `.mcif` -> `ModelFormatBackend::Cif`

The format implementation parses into `AtomicModelDataBlock`, not directly into `ModelObject`.

That intermediate block owns:

- atoms grouped by model number
- bonds
- parsed chemistry dictionaries
- component, atom, and bond key systems
- entity or chain metadata used during parsing
- identifiers such as PDB ID, EMD ID, and resolution

### 6.2 `ModelObjectFactory` Assembly Rules

`ModelObjectFactory::CreateDataObject(...)` converts parsed model data into the supported in-memory form.

Current behavior:

1. read the file into an `AtomicModelDataBlock`
2. choose model number `1` when present
3. otherwise fall back to the first model number in the file
4. move only the selected model's atom list into a new `ModelObject`
5. move the full parsed bond list, then retain only bonds whose endpoints are both in the selected atom set
6. move or copy supported metadata into `ModelObject`

Metadata transferred from the block to `ModelObject`:

- `pdb_id`
- `emd_id`
- `resolution`
- `resolution_method`
- `chain_id_list_map`
- chemical component dictionary
- component key system
- atom key system
- bond key system

Important current limits:

- parser-side entity metadata remains in `AtomicModelDataBlock`; it is not promoted to a public `ModelObject` domain surface
- secondary-structure ranges are applied to atom state during parsing and are not persisted as standalone structures

### 6.3 Write Path

Model export uses:

- `ModelFileWriter`
- `PdbFormat` or `CifFormat`, selected through `FileFormatRegistry`

Write support is intentionally narrower than read support:

| Path | Supported model extensions |
| --- | --- |
| Read | `.pdb`, `.cif`, `.mmcif`, `.mcif` |
| Write | `.pdb`, `.cif` |

`.mmcif` and `.mcif` are read-only aliases today.

## 7. Map File Pipeline

### 7.1 Read and Write Path

Map import and export are simpler than the model path.

Components:

- `MapFileReader` / `MapFileWriter`
- `MapFileFormatBase` implementations
- `MrcFormat`
- `CCP4Format`

`MapFileReader` and `MapFileWriter` also dispatch through `FileFormatRegistry`:

- `.mrc` -> `MapFormatBackend::Mrc`
- `.map`, `.ccp4` -> `MapFormatBackend::Ccp4`

`MapObjectFactory::CreateDataObject(...)` builds `MapObject` directly from:

- grid size
- grid spacing
- origin
- owned voxel array

There is no intermediate object equivalent to `AtomicModelDataBlock` because the map formats already match the in-memory structure closely.

### 7.2 Ownership Detail

`MapFileReader::GetMapValueArray()` transfers ownership of the voxel array to the caller. Treat the returned array as move-only state, not as a reusable buffer.

## 8. SQLite Persistence

### 8.1 Entry Points

Database persistence enters through:

- `DataObjectManager::SaveDataObject(key_tag, renamed_key_tag)`
- `DataObjectManager::LoadDataObject(key_tag)`

`DataObjectManager` itself does not know any schema details. It only:

- ensures a `DatabaseManager` exists
- fetches the in-memory object for save
- delegates to `DatabaseManager`
- inserts the loaded object back into the in-memory map

### 8.2 `DatabaseManager`

`DatabaseManager` owns:

- the SQLite connection
- schema bootstrap and migration via `DatabaseSchemaManager`
- the `object_metadata` dispatch table
- a DAO cache keyed by `std::type_index`
- the cross-table transaction boundary for save and load

`object_metadata` is the polymorphic dispatch table:

```sql
CREATE TABLE IF NOT EXISTS object_metadata (
    key_tag TEXT PRIMARY KEY,
    object_type TEXT
);
```

Save flow:

1. rely on the schema that was already prepared when `DatabaseManager` opened the database
2. open a single transaction
3. create or reuse the DAO for the runtime object type
4. call `dao->Save(...)`
5. upsert `(key_tag, object_type)` into `object_metadata`
6. commit on scope exit

Load flow:

1. rely on the schema that was already prepared when `DatabaseManager` opened the database
2. open a single transaction
3. read `object_type` from `object_metadata`
4. resolve the DAO from `DataObjectDAOFactoryRegistry`
5. call `dao->Load(key_tag)`
6. return the reconstructed object

This is an intentional contract change from the legacy design: DAO implementations are transaction-free. `DatabaseManager` is the only transaction owner for normal save or load operations.

`DatabaseManager` also treats schema initialization as a one-shot concern per instance. Bootstrap, migration, validation, and metadata repair happen when the database handle is opened, not on every hot-path save or load.

### 8.3 Schema Versioning and Migration

`DatabaseSchemaManager` uses SQLite `PRAGMA user_version` as the single schema-version source and makes its decision through `InspectSchemaState()`.

Current versions:

- `DatabaseSchemaVersion::LegacyV1 = 1`
- `DatabaseSchemaVersion::NormalizedV2 = 2`

Schema states:

- `Empty`
- `LegacyV1`
- `NormalizedV2`
- `ManagedButUnversioned`
- `MixedUnknown`

`EnsureSchema()` behavior:

1. if `PRAGMA user_version == 1`, migrate legacy v1 to v2
2. if `PRAGMA user_version == 2`, run `ValidateNormalizedV2Schema()` only
3. if `PRAGMA user_version` is any other non-zero value, reject the database as unsupported
4. if `PRAGMA user_version == 0`, inspect table state:
   - `Empty`: create normalized v2 tables, set version to `2`, validate
   - `LegacyV1`: migrate to v2
   - `ManagedButUnversioned`: create any missing managed tables, run `RepairManagedMetadata()`, set version to `2`, validate
   - `MixedUnknown`: fail fast instead of silently adopting the database

Important intent:

- a healthy v2 database is validated, not rewritten
- metadata repair is reserved for managed-but-unversioned recovery, not ordinary hot-path save/load
- databases with unknown unmanaged tables are not silently claimed by this subsystem

Migration is open-time and in-place:

1. open one transaction
2. create normalized v2 tables
3. read model keys from legacy `model_list` as the authoritative source
4. compare against any existing `object_metadata` model rows and log mismatches
5. load each legacy model through `LegacyModelObjectDAO`
6. save it through `ModelObjectDAOv2`
7. repair `object_metadata`
8. drop only the legacy tables explicitly owned by the migrated keys
9. set `PRAGMA user_version = 2`
10. validate normalized v2

Current migration limits:

- missing legacy metadata that was never persisted, such as historical `chain_id_list_map`, cannot be recovered
- historical sanitize collisions in legacy per-key tables cannot be repaired automatically
- map persistence stays on the existing `map_list` schema; only model persistence is migrated to normalized v2
- partial or stale legacy `object_metadata` does not decide migration scope; `model_list` does
- legacy cleanup is explicit-key based, not broad table-name pattern sweeping

### 8.4 DAO Registration

DAO registration is static and translation-unit driven:

- `ModelObjectDAO` registers as `"model"`
- `MapObjectDAO` registers as `"map"`

The registration mechanism is:

```cpp
DataObjectDAORegistrar<DataObjectType, DAOType>("stable_name")
```

The string name is persisted in `object_metadata.object_type`, so it must remain stable once databases exist on disk.

`LegacyModelObjectDAO` is intentionally not the registered DAO for `"model"`. It exists only as an internal migration reader.

## 9. `ModelObject` Persistence Schema

### 9.1 Current DAO Roles

- `ModelObjectDAO` is the public DAO registered for model persistence
- `ModelObjectDAO` is currently a thin wrapper around `ModelObjectDAOv2`
- `LegacyModelObjectDAO` reads legacy v1 layout only and is used by migration code

### 9.2 Normalized V2 Tables

`ModelObjectDAOv2` persists all model rows into fixed shared tables keyed by `key_tag`.

Implementation is intentionally split:

- `ModelObjectDAOv2.cpp`: orchestration only
- `src/data/model_io/ModelSchemaSql.hpp`: shared SQL constants and scoped table lists
- `src/data/model_io/ModelStructurePersistence.cpp`: structural save or load logic
- `src/data/model_io/ModelAnalysisPersistence.cpp`: local or posterior or group analysis save or load logic
- `src/data/model_io/SQLiteStatementBatch.hpp`: small prepared-statement helper for repeated insert patterns

Core tables:

- `model_object`
- `model_chain_map`
- `model_component`
- `model_component_atom`
- `model_component_bond`
- `model_atom`
- `model_bond`

Analysis-result tables:

- `model_atom_local_potential`
- `model_bond_local_potential`
- `model_atom_posterior`
- `model_bond_posterior`
- `model_atom_group_potential`
- `model_bond_group_potential`

The important design change is that persistence no longer creates per-key tables from sanitized `key_tag` values. This removes schema explosion and prevents collisions such as `a-b` and `a_b` mapping to the same table names.

### 9.3 Save Strategy

`ModelObjectDAOv2::Save(...)` uses scoped replacement by `key_tag`:

- ensure fixed tables exist
- delete only rows for the current `key_tag`
- insert current rows again into each table

It does not clear whole tables for other objects.

### 9.4 What Survives a Database Round-Trip

The normalized model path persists:

- structural metadata from `model_object`
- `chain_id_list_map` through `model_chain_map`
- chemical component dictionary
- atoms and bonds
- local potential entries
- posterior entries
- group potential entries

Selection state is reconstructed indirectly:

- atoms and bonds become selected when corresponding local-potential rows exist
- `ModelObject::Update()` and `SetBondList(...)` rebuild selected lists
- group membership is rebuilt by the existing classifier logic

### 9.5 Round-Trip Limits

A database round-trip still does not attempt to persist every parser-side detail.

Not persisted as standalone public state:

- parser-only entity metadata from `AtomicModelDataBlock`
- derived caches such as KD-tree-like accelerators

## 10. `MapObject` Persistence Schema

`MapObjectDAO` keeps the existing flat schema:

- `map_list`

Stored columns:

- `key_tag`
- `grid_size_x`, `grid_size_y`, `grid_size_z`
- `grid_spacing_x`, `grid_spacing_y`, `grid_spacing_z`
- `origin_x`, `origin_y`, `origin_z`
- `map_value_array` as a BLOB

`MapObjectDAO` now follows the same transaction contract as other DAOs: it does not open its own transaction.

## 11. SQLite Utility Layer

The DAO layer depends on `SQLiteWrapper` for:

- SQL execution via `Execute(...)`
- prepared statements via `Prepare(...)`
- typed binders via `Bind<T>(...)`
- typed column readers via `GetColumn<T>(...)`
- RAII statement cleanup via `StatementGuard`
- RAII transactions via `TransactionGuard`

One practical constraint matters when extending DAO or migration code: `SQLiteWrapper` supports one active prepared statement at a time. Fetch key lists first if a later step needs to prepare another statement on the same connection.

## 12. Supported Surface

| Top-level object | File read | File write | SQLite save/load |
| --- | --- | --- | --- |
| `ModelObject` | `.pdb`, `.cif`, `.mmcif`, `.mcif` | `.pdb`, `.cif` | yes |
| `MapObject` | `.mrc`, `.map`, `.ccp4` | `.mrc`, `.map`, `.ccp4` | yes |

## 13. Key Source Files

Start with these files when debugging or extending this subsystem:

- `include/core/DataObjectManager.hpp`
- `src/core/DataObjectManager.cpp`
- `include/data/FileFormatRegistry.hpp`
- `src/data/FileFormatRegistry.cpp`
- `include/data/FileProcessFactoryRegistry.hpp`
- `src/data/FileProcessFactoryRegistry.cpp`
- `include/data/FileFormatBackendFactory.hpp`
- `src/data/FileFormatBackendFactory.cpp`
- `src/data/ModelObjectFactory.cpp`
- `src/data/MapObjectFactory.cpp`
- `src/data/ModelFileReader.cpp`
- `src/data/MapFileReader.cpp`
- `src/data/ModelFileWriter.cpp`
- `src/data/MapFileWriter.cpp`
- `include/data/DatabaseManager.hpp`
- `src/data/DatabaseManager.cpp`
- `include/data/DatabaseSchemaManager.hpp`
- `src/data/DatabaseSchemaManager.cpp`
- `include/data/DataObjectDAOFactoryRegistry.hpp`
- `src/data/DataObjectDAOFactoryRegistry.cpp`
- `include/data/ModelObjectDAO.hpp`
- `include/data/ModelObjectDAOv2.hpp`
- `include/data/LegacyModelObjectDAO.hpp`
- `src/data/ModelObjectDAO.cpp`
- `src/data/ModelObjectDAOv2.cpp`
- `src/data/LegacyModelObjectDAO.cpp`
- `src/data/model_io/ModelSchemaSql.hpp`
- `src/data/model_io/ModelStructurePersistence.cpp`
- `src/data/model_io/ModelAnalysisPersistence.cpp`
- `src/data/model_io/SQLiteStatementBatch.hpp`
- `include/data/MapObjectDAO.hpp`
- `src/data/MapObjectDAO.cpp`
- `include/data/SQLiteWrapper.hpp`

## 14. Extension Playbooks

### 14.1 Add a New File Format for an Existing Object

For a new model or map format:

1. implement or extend the appropriate `*Format` backend
2. add or update the descriptor in `FileFormatRegistry`
3. update `FileFormatBackendFactory` only if the new descriptor needs a new backend enum or backend branch
4. add read and write tests for the supported matrix
5. update this document's support matrix

Do not add built-in format truth to `FileProcessFactoryRegistry`; that registry is now only the override seam.

If the format parses into a different intermediate shape, decide whether it still fits the current factory contract or whether a new assembly seam is needed.

### 14.2 Add a New Persistent Top-Level `DataObject`

For a new database-persisted top-level object:

1. derive from `DataObjectBase`
2. implement a `DataObjectDAOBase` subclass
3. register it with `DataObjectDAORegistrar<...>("stable_name")`
4. ensure the DAO stays transaction-free
5. extend `DatabaseSchemaManager` if schema bootstrap or migration must create new tables
6. decide whether it also needs file factories and readers or writers
7. add round-trip tests for save and load
8. document the schema and support matrix here

`DatabaseManager` usually does not need new branching logic; it already dispatches through the DAO registry.

## 15. Current Gotchas

- `ProcessFile(...)` and `LoadDataObject(...)` replace an existing in-memory `key_tag`.
- `ClearDataObjects()` only clears the in-memory map; it does not delete database rows.
- `ProduceFile(...)` warns and returns when the key is missing; it does not throw.
- `SaveDataObject(original, renamed)` saves under a new database key but leaves the in-memory object under `original`.
- model write support is intentionally narrower than model read support
- parser-only metadata still does not become full `ModelObject` domain state just because it existed during file parsing
