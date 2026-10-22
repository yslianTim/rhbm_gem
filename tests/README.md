# Test Organization Guide

This project uses a two-axis test organization model:

1. **Directory domain**: where tests live.
2. **CTest intent labels**: what behavior they verify.

## Directory Matrix

| Directory | Domain | Typical focus |
| --- | --- | --- |
| `tests/core/command/` | `core` | Command workflows, option handling, command-level validation |
| `tests/core/contract/` | `core` | Command catalog/metadata/surface contracts and docs sync checks |
| `tests/core/second_stage/` | `core` | Second-stage fitting state, solvers, acceptance, recovery, convergence, and finalization |
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
- `intent:algorithm`
- `intent:bindings`

## Common Commands

Build all test targets:

```bash
cmake --build build --target tests_all --parallel 4
```

Build output note:

- Default grouped C++ tests are compiled into a single executable: `build/bin/RHBM-GEM-TEST`
- `ctest` exposes grouped entries such as `rhbm_tests_core_command` and `rhbm_tests_data_schema`
- Those grouped CTest entries run filtered subsets from the same `RHBM-GEM-TEST` binary

Run all tests:

```bash
ctest --test-dir build --output-on-failure --parallel 4
```

Run by domain:

```bash
ctest --test-dir build -L domain:data --output-on-failure --parallel 4
```

Run by intent:

```bash
ctest --test-dir build -L intent:algorithm --output-on-failure --parallel 4
```

Run repository guards and install consumer smoke (lint lane):

```bash
cmake --build build --target lint_repo --parallel 4
```

## Joint component regression

The Joint component sources under `tests/core/joint_component/` are split into
`rhbm_tests_joint_component_contract`,
`rhbm_tests_joint_component_workflow`, and
`rhbm_tests_joint_component_numerical`. These run from the ordinary
`RHBM-GEM-TEST` executable. They retain deterministic checks for observation
and support domains, component construction, width calibration, bounded A/C
and Active-set LM solves, FixedNeighbor acceptance and convergence, local
conditioning and stability, geometry-local Halo medians, Original-domain
Peeling, persistence, and workflow behavior.

Run the Joint groups:

```bash
ctest --test-dir build -L 'joint:(contract|workflow|numerical)' --parallel 4 --output-on-failure
```

Research conclusions and measurement limits are summarized in the
[Joint evidence record](../docs/developer/joint-component-evidence.md); the
[runtime guide](../docs/developer/joint-component-runtime.md) describes the
current production contract.

## Adding New Tests

### Second-stage fitting

Keep second-stage defense tests in the following behavior-based files under
`tests/core/second_stage/`:

| File | Responsibility |
| --- | --- |
| `GraphAndBackground_test.cpp` | Graphs, partitions, topology drift, physical halos, uncut component construction, frozen backgrounds, and residual overlays |
| `SeedAndState_test.cpp` | Seed selection and fallback, Gaussian medians, transformed coordinates, damping, and extrapolation |
| `SolverAndProposal_test.cpp` | Conditioning, solver health, joint offsets, frozen-background RHS, and joint polish proposals |
| `CandidateAcceptance_test.cpp` | Objective guards, trust regions, backtracking, transaction publication, and boundary acceptance |
| `Recovery_test.cpp` | Non-finite solver fallback, failure isolation, and quarantine retries |
| `ConvergenceAndFinalization_test.cpp` | Active-coordinate convergence certificates, final dependency polish, and persistence |
| `ProductionFitting_test.cpp` | Nominal endpoint qualification and read-only evaluation |

The 24 `EstimatorSecondStageDefenseTest` cases and two
`ProductionFittingTest` cases run through the single
`rhbm_tests_core_estimator` CTest group. Do not add per-file CTest groups using
these suite filters: each would repeat the same tests. Acceptance and rollback
checks stay with the production behavior they exercise.

Use `tests/support/SecondStageTestSupport.hpp/.cpp` and its `second_stage_test`
namespace for fixture builders and assertions shared across these files. Keep
single-file helpers in that file's anonymous namespace and support-only helpers
private to the support implementation. Build fresh model/solver state per case;
do not share mutable fixtures. Use `detail` for `rhbm_gem::core::detail` throughout
these tests. `EstimatorTester_test.cpp` retains its workflow tests and fixtures.

Compile shared test support directly into `rhbm_tests`, outside the suite source
lists. `SecondStageTestSupport` provides generated inputs for the retained
production behavior tests. Failed-only refinement cases live with the MDPDE
regression tests and use generated endpoints.

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

## Joint validation workflow

The default CTest suite includes current Joint result persistence, unsupported
format rejection, and CLI workflow smoke coverage. Python binding smoke tests
also exercise the production workflow with a synthetic map generated from
`tests/fixtures/test_model.cif`. There is no optional Joint research or
benchmark CTest lane.
