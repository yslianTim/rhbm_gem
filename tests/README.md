# Test Organization Guide

This project uses a two-axis test organization model:

1. **Directory domain**: where tests live.
2. **CTest intent labels**: what behavior they verify.

## Directory Matrix

| Directory | Domain | Typical focus |
| --- | --- | --- |
| `tests/core/command/` | `core` | Command workflows, option handling, command-level validation |
| `tests/core/contract/` | `core` | Command catalog/metadata/surface contracts and docs sync checks |
| `tests/data/` | `data` | Data public-surface guards, file I/O/runtime behavior, and schema/persistence validation |
| `tests/utils/math/` | `utils` | Numeric/statistical/geometry helper algorithms |
| `tests/utils/domain/` | `utils` | Domain helpers (string/logging/file-path/chemistry-related helpers) |
| `tests/utils/hrl/` | `utils` | HRL-specific algorithm and transform tests |
| `tests/integration/` | `integration` | Python binding smoke/validation scripts and end-to-end command pipeline checks |
| `tests/fixtures/` | fixture | Shared fixture files used by C++ and Python tests |

## Label Vocabulary

Required CTest labels:

- `domain:core`
- `domain:data`
- `domain:utils`
- `domain:integration`
- `intent:contract`
- `intent:command`
- `intent:validation`
- `intent:io`
- `intent:schema`
- `intent:migration`
- `intent:algorithm`
- `intent:bindings`

## Common Commands

Build all test targets:

```bash
cmake --build build --target tests_all -j
```

Build output note:

- C++ tests are compiled into a single executable: `build/bin/RHBM-GEM-TEST`
- `ctest` still exposes grouped entries such as `rhbm_tests_core_command`, `rhbm_tests_data_contract`, and `rhbm_tests_data_schema`
- those grouped CTest entries run filtered subsets from the same `RHBM-GEM-TEST` binary

Run all tests:

```bash
ctest --test-dir build --output-on-failure
```

Run by domain:

```bash
ctest --test-dir build -L domain:data --output-on-failure
```

Run by intent:

```bash
ctest --test-dir build -L intent:migration --output-on-failure
```

Run repository guards and install consumer smoke (lint lane):

```bash
cmake --build build --target lint_repo
```

## Adding New Tests

- Place new tests in the matching domain directory.
- For `tests/core/command/`, prefer extending the existing grouped files:
- `CommandBaseLifecycle_test.cpp` for base lifecycle/preflight behavior.
- `CommandValidationHelpers_test.cpp` for reusable helper semantics.
- `CommandValidationScenarios_test.cpp` for command-specific validation rules.
- `CommandWorkflowScenarios_test.cpp` for command workflows and output side effects.
- Only create a new command `*_test.cpp` when you are introducing a new testing responsibility.
- For `tests/data/`, prefer extending the responsibility-based files:
- `DataPublicSurface_test.cpp` for public header surface guards.
- `DataObjectFileIO_test.cpp` for file I/O and map-axis import behavior.
- `DataObjectImportRegression_test.cpp` for CIF/MMCIF parser regression matrices.
- `DataObjectRuntimeBehavior_test.cpp` for `ModelObject` / `MapObject` state behavior.
- `DataObjectDispatchAndIngestion_test.cpp` for painter ingestion contracts.
- `DataObjectPersistence_test.cpp`, `DataObjectSchemaBootstrap_test.cpp`, `DataObjectSchemaCompatibility_test.cpp`, and `DataObjectSchemaValidation_test.cpp` for database persistence/schema behavior.
- Place new fixture files under `tests/fixtures/`.
- Use `tests/support/` for shared test-only seams and reusable helpers/assertions.
- Add the source to the correct grouped target in `tests/CMakeLists.txt`.
- Ensure the target has the correct `domain:*` and `intent:*` labels.
- Prefer searchable suite names (for example `DataObjectSchemaMigrationTest`) over generic names.
