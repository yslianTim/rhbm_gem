# Map Simulation

`map_simulation` prepares per-atom charges once, accumulates an electrostatic potential map, and writes a matching generation record for every map. The CLI, C++ API, and Python command API use the same implementation. No additional option is needed to save the record.

## Responsibility boundaries

- `ComponentHelper::LookupPartialCharge()` is a read-only lookup of residue, spot, structure, and table selection. It returns an optional charge, the selected table when available, and a status. It has no mutable cache and does not log or decide a simulation fallback.
- `GetPartialCharge()` remains the compatibility API: it returns zero for missing data and optionally logs a warning. The new lookup distinguishes an actual table value of zero from missing data.
- `PrepareSimulationAtomList()` snapshots the selected atoms in model order. Coordinates, identity, charge lookup evidence, and the coefficient actually passed to `ElectricPotential` are retained together. Serial IDs do not act as unique storage keys.
- `PopulateMapValueArray()` only evaluates these prepared inputs. Each z plane has one writer and processes its candidate atoms in preparation order, with the existing spherical cutoff. The manifest writer consumes the same prepared records and does not re-query charge tables.

FREE selects buried, SHEET selects sheet, and the existing interval `HELX_P <= structure < TURN_P` selects helix. AMBER selects amber95 regardless of structure. Unsupported structures do not inherit another call's table. Tables and their numerical values are unchanged.

For PARTIAL and AMBER, missing data uses `charge_used=0`, retains its failure status, and emits a warning summarized by reason. NEUTRAL does not look up a table and records `neutral_mode`. In particular, `--charge 1` selects PARTIAL; it does not mean every atom has charge +1.

## Determinism and failures

For the same executable, input, normalized settings, and floating-point environment, changing the requested OpenMP team size does not change voxel summation order. Repeated j1/j2/j4 runs produce identical float64 voxel bits and identical map bytes. OpenMP-disabled builds use the same plane traversal. There is no promise of bitwise identity between different compilers, architectures, math libraries, or floating-point compiler options.

The candidate index stores atom indices by intersected z plane. It does not allocate a full map per worker or scan all atoms at every voxel. Worker exceptions are captured per plane and rethrown on the calling thread. Failed generation leaves the previous in-memory map intact. Non-finite contributions or sums, and values outside the finite float32 range, fail the command before publication.

Grid origin, dimensions, sampling, and potential formulas retain their existing definitions. The outer atom cutoff is separate from the kernel's charge term cutoff: SINGLE_GAUS uses 2.5 Å; FIVE_GAUS_CHARGE uses 3.0 Å for positive blurring widths. Near-zero distance handling and the five-Gaussian charge width lower bound are reported by `ElectricPotential::GetKernelSettings()` using the same constants as evaluation.

`SINGLE_GAUS_USER` still requires amplitude and width arguments that this command does not supply; a non-empty simulation using that option fails on its non-finite potential instead of publishing invalid voxels. This change does not add user-Gaussian parameter support.

## Generation record, schema version 1

Each output `name.map` has an adjacent `name.map.simulation.json`. All floating-point fields are JSON numbers serialized with double round-trip precision. Identifiers are escaped as JSON strings. Atom enum fields preserve the integer codes from `GlobalEnumClass.hpp`; raw component and atom identifiers are saved separately.

| Object or field | Contents |
| --- | --- |
| `schema_version` | Integer `1` |
| `generator` | Project `version`, `source_sha256`, `configuration_sha256`, and combined `build_sha256` |
| `source` | Absolute `model_path`, `pdb_id`, and SHA-256 of the model file bytes |
| `output` | `map_file`, SHA-256 of the saved map bytes, `format=ccp4`, `calculation_precision=float64`, `storage_precision=float32` |
| `settings` | Potential model and charge mode names/codes, this map's exact `blurring_width`, the normalized `blurring_width_list`, grid spacing/size/origin, outer cutoff, coordinate unit, selection flags, and `missing_charge_policy=zero_with_status` |
| `settings` treatment flags | Occupancy, temperature factor, and normalization are currently not applied; each has an explicit false flag |
| `kernel` | `charge_term_cutoff`, `near_zero_distance`, `minimum_charge_width`, and `effective_charge_width`; inapplicable fields are null |
| `execution` | OpenMP availability, normalized requested and actual team sizes, `accumulation=z_planes_in_preparation_order` |
| `charge_semantics` | `argument_passed_to_electric_potential` |
| `atom_count`, `fallback_charge_count` | Selected atom count and count of failed lookups using zero |
| `atoms` | Ordered records described below |

Every atom record contains:

- Zero-based `preparation_index`, `serial_id`, `sequence_id`, `chain_id`, `component_id`, `atom_id`, and `alternate_indicator`.
- `element`, `residue`, `spot`, `structure`, and actual three-dimensional `position` in Å.
- `table`: `buried`, `helix`, `sheet`, `amber95`, or null when no table was selected.
- `lookup_status`: `found`, `unsupported_residue`, `unsupported_structure`, `unsupported_spot`, `table_data_mismatch`, or `neutral_mode`.
- `lookup_charge`: the table value on success, otherwise null.
- `charge_used`: the actual coefficient supplied to potential evaluation, including zero fallback and neutral mode.

The combination of source fingerprint, preparation index, and identity fields identifies an atom. A consumer must not match solely by element or assume serial IDs are globally unique. For SINGLE_GAUS, this coefficient is the physical offset under that model's existing contract; other potential models must be interpreted according to their own formula.

### Build fingerprints

The generated header is refreshed by build dependencies when production source/header/CMake inputs change, including uncommitted edits. A Git checkout is not required. The build directory retains `generated/<configuration>/SimulationSources.txt`, `SimulationConfiguration-CXX.txt`, and `SimulationBuildInfo.hpp`.

The source digest hashes the sorted sequence of relative source paths and each file's SHA-256. The configuration digest hashes the generated configuration text, which includes compiler, selected configuration, flags, definitions, include paths, and dependency information. The combined digest hashes the two digests separated and terminated by newlines. These identify the recorded project source/configuration inputs; they are not an archive or cryptographic fingerprint of every external library or of the executable bytes.

Boost >=1.90 supplies SHA-256 through Boost.Hash2 and the JSON serializer. Boost.JSON is compiled once inside `rhbm_gem`, without requiring a separately deployed Boost.JSON library.

## Publication and compatibility

The command checks all width-derived output names before generating maps. Any duplicate name, including equal widths or distinct widths that round to the same two-decimal filename, is rejected. The record always stores the full width, not a value parsed from that filename.

Each map and record are first written in a unique staging directory next to the destination. The map stream is closed and checked, then the map bytes are hashed and the JSON stream is closed and checked. Existing regular output files are backed up during replacement. The new map is published first and its record last. A caught publication failure removes newly published files and restores the prior pair. If filesystem recovery itself fails, recovery files are retained and the error reports their directory.

Publication is per map, not across an entire width batch: earlier completed pairs remain if a later map fails. Two separate file renames are not crash-atomic. A consumer must require the record and verify its `output.map_sha256` before treating the pair as complete. The command returns success only after both files have been published.

The map format and command request fields remain compatible; the automatic JSON file is an additional output. Existing maps without records are not assigned reconstructed truth. Their original charge choices must be investigated separately. Estimation and convergence do not read these records as hidden inputs.

## Verification

`MapSimulationTest` covers overlapping atoms, exact voxel bits across team sizes, independent serial summation, center/cutoff boundaries, JSON-based map reconstruction through float32 storage, mixed structures/tables, legal and fallback zeros, selection, empty input, width collisions, string escaping, worker exceptions, and output cleanup. `ComponentHelperTest` additionally covers all existing charge tables and concurrent lookup history. Run these with OpenMP ON and OFF; set `OMP_DYNAMIC=FALSE` and permit at least four OpenMP threads for the ON assertions.
