# Test Organization Guide

This project uses a two-axis test organization model:

1. **Directory domain**: where tests live.
2. **CTest intent labels**: what behavior they verify.

## Directory Matrix

| Directory | Domain | Typical focus |
| --- | --- | --- |
| `tests/core/command/` | `core` | Built-in command workflows, option handling, command-level validation |
| `tests/core/contract/` | `core` | Command catalog/metadata/surface contracts and docs sync checks |
| `tests/data/io/` | `data` | DataObjectManager and file I/O integration contracts |
| `tests/data/schema/` | `data` | Schema bootstrap, migration, schema validation, persistence contracts |
| `tests/utils/math/` | `utils` | Numeric/statistical/geometry helper algorithms |
| `tests/utils/domain/` | `utils` | Domain helpers (string/logging/file-path/chemistry-related helpers) |
| `tests/utils/hrl/` | `utils` | HRL-specific algorithm and transform tests |
| `tests/integration/` | `integration` | Python binding smoke/validation scripts and style-guard tests |
| `tests/data/` | fixture | Shared fixture files used by C++ and Python tests |

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
- `intent:style`

Style guard coverage includes:

- repository hygiene guard (`cmake/RepositoryHygieneGuard.cmake`)
- fixture tracking guard (`cmake/TestFixtureTrackingGuard.cmake`)

## Common Commands

Build all test targets:

```bash
cmake --build build --target tests_all -j
```

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

## Adding New Tests

- Place new tests in the matching domain directory.
- Use `tests/support/` for shared test-only seams (for example internal dependency shims).
- Add the source to the correct grouped target in `tests/cmake/*.cmake`.
- Ensure the target has the correct `domain:*` and `intent:*` labels.
- Prefer searchable suite names (for example `DataObjectSchemaMigrationTest`) over generic names.
