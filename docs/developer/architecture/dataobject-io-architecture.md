# DataObject I/O Architecture

This document describes file I/O and the current Model-only SQLite
persistence boundary.

## Public surface

File I/O supports both model and map objects:

```cpp
std::unique_ptr<ModelObject> ReadModel(const std::filesystem::path & path);
void WriteModel(const std::filesystem::path & path, const ModelObject & model, int model_parameter = 1);

std::unique_ptr<MapObject> ReadMap(const std::filesystem::path & path);
void WriteMap(const std::filesystem::path & path, const MapObject & map);
```

SQLite persistence is narrower:

```cpp
DataRepository repository{database_path};
repository.SaveModel(model, key_tag);
auto model = repository.LoadModel(key_tag);
```

`DataRepository` does not store maps. Maps use MRC/CCP4 files; models use
PDB/mmCIF files or the SQLite repository.

## File codecs

| Object | Read | Write |
| --- | --- | --- |
| `ModelObject` | `.pdb`, `.cif`, `.mmcif`, `.mcif` | `.pdb`, `.cif` |
| `MapObject` | `.mrc`, `.map`, `.ccp4` | `.mrc`, `.map`, `.ccp4` |

Extension dispatch is case-insensitive. MRC uses its stored floating-point
origin; CCP4 derives origin from integer start indices and grid spacing. Both
codecs decode float32 mode-2 voxels into the runtime's contiguous double map.

## Current persistence contract

| Layer | Supported format |
| --- | --- |
| SQLite | Schema v20 only |
| Stage Result JSON | Version 4 only |
| Joint Result JSON | Schema v19 only |

No layer migrates historical data. A genuinely empty SQLite database is
initialized as v20. An existing database opens only when its `user_version`
is 20 and its exact current schema validates. Unsupported versions or
structures are rejected before schema or row mutation; the repository does
not repair, upgrade, delete, or replace them.

Stage JSON readers require `{"version":4}`. Joint JSON writers always emit
schema 19, and Joint readers require schema 19 with the
`target-only-pure-halo-support-exclusion-v1` parameterization. Version
selection does not depend on convergence, Halo count, or reconstruction
availability. Current documents may represent failed initialization,
unavailable calibration, search failure, an unqualified endpoint, missing
component state, no Halo, unavailable reconstruction, incomplete downstream
results, and nullable diagnostics.

## SQLite schema and lifecycle

Schema v20 contains these tables:

- `model_object`
- `model_chain_map`
- `model_component`
- `model_component_atom`
- `model_component_bond`
- `model_atom`
- `model_bond`
- `model_joint_result`
- `model_stage_result`

Opening an existing database checks the exact table set, ordered columns,
primary keys, required nullability, foreign keys, selection flags, and
referential integrity. Extra schema objects and malformed current schemas are
rejected.

Each repository operation is serialized by its mutex and runs in a SQLite
transaction. Model save replaces the selected key atomically, validates the
Joint/Stage JSON payloads, and relies on foreign-key cascades for child rows.
Joint and Stage rows are loaded as their current JSON values; Stage results are
not reconstructed from other tables.

## Load and save flow

Save writes the current model structure and any Joint and Stage results in one
transaction. A failed write rolls back all rows for that save. A key saved
without a Joint result has no Joint result row.

Load validates the database, loads the current model structure, decodes a
Joint result only as schema 19, and decodes a Stage result only as version 4.
Any malformed JSON or unsupported payload fails the load; it is not rewritten
in place.

## Key files

- `src/data/io/DataRepository.cpp`: schema-version and exact-schema checks.
- `src/data/io/sqlite/ModelObjectStorage.cpp`: current model and result SQL.
- `src/data/io/detail/PersistenceSchema.hpp`: SQLite, Stage, and Joint
  version constants.
- `src/data/io/detail/StageResultJson.cpp`: Stage JSON version 4.
- `src/data/io/detail/JointResultJson.cpp`: Joint JSON schema 19.
- `tests/data/DataObjectSchemaLifecycle_test.cpp` and
  `DataObjectSchemaValidation_test.cpp`: database bootstrap and rejection.
- `tests/data/DataObjectPersistence_test.cpp`: current JSON and model
  round-trips.

See the [Joint runtime guide](../joint-component-runtime.md) for the
target-only parameterization and production solver contract.
