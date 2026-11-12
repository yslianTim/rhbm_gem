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
| SQLite | Schema v25 only |
| Stage Result JSON | Version 9 only |
| Joint Result JSON | Schema v25 only |

SQLite initializes schema v25 when `PRAGMA user_version` is zero and
`sqlite_master` contains no user schema objects. An existing database opens
only when `user_version` is 25 and the exact current schema validates.
SQLite and both JSON codecs require their current formats without migration.
The repository rejects unsupported versions or structures before schema or
row mutation.

Stage JSON readers require `{"version":9}`. Joint JSON writers always emit
schema 25, and Joint readers require schema 25 with the
`target-only-pure-halo-support-exclusion-v1` parameterization. Version
selection does not depend on convergence, Halo count, or reconstruction
availability. Current documents may represent failed initialization,
unavailable calibration, search failure, an unqualified endpoint, missing
component state, no Halo, unavailable reconstruction, incomplete downstream
results, and nullable diagnostics.

Stage JSON encodes the complete `SamplingEntries()` through the `"raw"` array
and restores it through `SetSamplingEntries()` with the recorded geometry
availability. Workflow replacements restore geometry availability. Serialization
preserves ordering, unselected points and the `"y"`, `"distance"`, `"selected"`,
and `"position"` fields.
`GetSelectedSamplingEntries()` is not used to serialize the stored records.
Original-grid peeling retains its `"peeling"` record and stores the
unpeeled voxel `PostJointGridObservation::map_value` under `"raw_response"`,
and its derived peeled response under `"response"`. These keys distinguish
serialization records independently of C++ API spelling. SQLite schema 25,
Stage JSON 9 and Joint JSON 25 retain their existing round-trip contract.
Map Value denotes the currently stored voxel, including after normalization;
Sampling Value denotes an interpolation output. Generic response/observation
containers may also carry derived values without changing their source semantics.

## SQLite schema and lifecycle

Schema v25 contains these tables:

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
Joint result only as schema 25, and decodes a Stage result only as version 9.
Any malformed JSON or unsupported payload fails the load; it is not rewritten
in place.

## Key files

- `src/data/io/DataRepository.cpp`: schema-version and exact-schema checks.
- `src/data/io/sqlite/ModelObjectStorage.cpp`: current model and result SQL.
- `src/data/io/detail/PersistenceSchema.hpp`: SQLite, Stage, and Joint
  version constants.
- `src/data/io/detail/StageResultJson.cpp`: Stage JSON version 9.
- `src/data/io/detail/JointResultJson.cpp`: Joint JSON schema 25.
- `tests/data/DataObjectSchemaLifecycle_test.cpp` and
  `DataObjectSchemaValidation_test.cpp`: database bootstrap and rejection.
- `tests/data/DataObjectPersistence_test.cpp`: current JSON and model
  round-trips.

See [Joint Stage Integration](../joint-stage-integration.md#persistence)
for Joint and downstream payload validation.
