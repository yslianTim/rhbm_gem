# Test Organization Guide

This project uses a two-axis test organization model:

1. **Directory domain**: where tests live.
2. **CTest intent labels**: what behavior they verify.

## Directory Matrix

| Directory | Domain | Typical focus |
| --- | --- | --- |
| `tests/core/command/` | `core` | Command workflows, option handling, command-level validation |
| `tests/core/contract/` | `core` | Command catalog/metadata/surface contracts and docs sync checks |
| `tests/core/second_stage/` | `core` | Second-stage fitting state, solvers, acceptance, recovery, finalization, and observation |
| `tests/data/` | `data` | Data public-surface guards, file I/O/runtime behavior, and schema/persistence validation |
| `tests/utils/math/` | `utils` | Numeric/statistical/geometry helper algorithms |
| `tests/utils/domain/` | `utils` | Domain helpers (string/logging/file-path/chemistry-related helpers) |
| `tests/utils/hrl/` | `utils` | HRL-specific algorithm and transform tests |
| `tests/integration/` | `integration` | Python binding smoke/validation scripts and end-to-end command pipeline checks |
| `tests/fixtures/` | fixture | Shared fixture files used by C++ and Python tests |
| `tests/benchmarks/` | benchmark | Small checked-in baselines for opt-in external-data regressions |

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

The 168-atom simulation regression is intentionally excluded unless configured
with `RHBM_GEM_ENABLE_FOLD_168_REGRESSION=ON`. It uses hash-verified external
inputs and can be selected with `-R fold_168_simulation_regression` or
`-L benchmark:external`; see the developer build guide for configuration.
`RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT=ON` builds the developer-only
frozen-IRLS trust-model instrumentation and requires `BUILD_TESTING=ON`.

Run repository guards and install consumer smoke (lint lane):

```bash
cmake --build build --target lint_repo
```

## Adding New Tests

### Second-stage fitting

Keep second-stage defense tests in the following behavior-based files under
`tests/core/second_stage/`:

| File | Responsibility |
| --- | --- |
| `GraphAndBackground_test.cpp` | Graphs, partitions, topology drift, physical halos, uncut component construction, frozen backgrounds, and residual overlays |
| `SeedAndState_test.cpp` | Seed selection and fallback, Gaussian medians, transformed coordinates, damping, and extrapolation |
| `SolverAndProposal_test.cpp` | Conditioning, solver health, joint offsets, Jacobians, refits, and joint polish/correction proposals |
| `CandidateAcceptance_test.cpp` | Objectives and best references, trust radii, backtracking, transaction publication, boundary acceptance, and serial/parallel selection |
| `Recovery_test.cpp` | Suspicious guards, failure masks, quarantine, fallback, ridge guards, and healthy remote clusters |
| `ConvergenceAndFinalization_test.cpp` | Active-coordinate certificates, final dependency polish, persistence, and whole-run intensity scaling |
| `Observation_test.cpp` | Cluster-history lifecycle, performance counters, trust shadow diagnostics, logging, and phase audit neutrality |

Preserve the existing `EstimatorSecondStageDefenseTest` suite and case names.
All seven files belong to `CORE_ESTIMATOR_TEST_SOURCES` and run through the single
`rhbm_tests_core_estimator` CTest group. Do not add per-file CTest groups using
this shared suite filter: each would repeat the entire suite. Pure trust-shadow
policy assertions belong to independent observation cases under
the original trust-experiment guard. Mixed cases that also verify production
acceptance or rollback remain with the production behavior they exercise.

Use `tests/support/SecondStageTestSupport.hpp/.cpp` and its `second_stage_test`
namespace for fixture builders and assertions shared across these files. Keep
single-file helpers in that file's anonymous namespace and support-only helpers
private to the support implementation. Build fresh model/solver state per case;
do not share mutable fixtures. Use `detail` for `rhbm_gem::core::detail` throughout
these tests. `EstimatorTester_test.cpp` retains its workflow tests and fixtures.

Compile the support implementation directly into `rhbm_tests`, outside the
suite-discovery source lists. Keep both observation options consistent with the
library; direct phase collector tests still compile `PhaseAudit.cpp` into the
test target when production phase auditing is disabled.

### General placement

- Place new tests in the matching domain directory.
- For `tests/core/command/`, prefer extending the existing grouped files:
- `CommandRunnerLifecycle_test.cpp` for runner lifecycle/preflight behavior.
- `CommandValidationHelpers_test.cpp` for reusable helper semantics.
- `CommandScenarios_test.cpp` for command-specific validation rules, workflows, and output side effects.
- Only create a new command `*_test.cpp` when you are introducing a new testing responsibility.
- For `tests/data/`, prefer extending the responsibility-based files:
- `DataPublicHeaders_test.cpp` for public header surface guards.
- `DataObjectFileIO_test.cpp` for file I/O and map-axis import behavior.
- `DataObjectImportRegression_test.cpp` for CIF/MMCIF parser regression matrices.
- `DataObjectMapBehavior_test.cpp` for `MapObject` runtime behavior.
- `DataObjectModelAnalysis_test.cpp` for selection, local entry, and group rebuild behavior.
- `DataObjectAssemblySpatialQuery_test.cpp` for model assembly, derived state, and neighbor queries.
- `RHBMTypes_test.cpp` and `LocalPotentialSeries_test.cpp` for analysis value math and local potential series derivations.
- `DataObjectPersistence_test.cpp`, `DataObjectSchemaLifecycle_test.cpp`, and `DataObjectSchemaValidation_test.cpp` for database persistence/schema behavior.
- For `tests/core/painter/`, keep painter-private behavior and ingestion contract tests near the internal painter implementation.
- Place new fixture files under `tests/fixtures/`.
- Use `tests/support/` for shared test-only seams and reusable helpers/assertions.
- Add the source to the correct grouped target in `tests/CMakeLists.txt`.
- Ensure the target has the correct `domain:*` and `intent:*` labels.
- Prefer searchable suite names (for example `DataObjectSchemaMigrationTest`) over generic names.
