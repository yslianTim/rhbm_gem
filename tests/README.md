# Test Organization Guide

This project uses a two-axis test organization model:

1. **Directory domain**: where tests live.
2. **CTest intent labels**: what behavior they verify.

## Directory Matrix

| Directory | Domain | Typical focus |
| --- | --- | --- |
| `tests/core/command/` | `core` | Command workflows, option handling, command-level validation |
| `tests/core/contract/` | `core` | Command catalog/metadata/surface contracts and docs sync checks |
| `tests/core/second_stage/` | `core` | Second-stage fitting state, solvers, acceptance, recovery, finalization, and performance logging |
| `tests/data/` | `data` | Data public-surface guards, file I/O/runtime behavior, and schema/persistence validation |
| `tests/experiments/` | optional tools | Current Joint measurement, research, and offline diagnosis executables |
| `tests/tools/` | test support | Current runtime and offline validation harness executables |
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
cmake --build build --target tests_all -j
```

Build output note:

- Default grouped C++ tests are compiled into a single executable: `build/bin/RHBM-GEM-TEST`
- `ctest` exposes grouped entries such as `rhbm_tests_core_command`, `rhbm_tests_data_contract`, and `rhbm_tests_data_schema`
- Those grouped CTest entries run filtered subsets from the same `RHBM-GEM-TEST` binary
- Optional Joint offline C++ tests use a separate executable.

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
ctest --test-dir build -L intent:algorithm --output-on-failure
```

Run repository guards and install consumer smoke (lint lane):

```bash
cmake --build build --target lint_repo
```

## Joint component regression

The Joint component sources under `tests/core/joint_component/` are split into
`rhbm_tests_joint_component_contract`,
`rhbm_tests_joint_component_workflow`, and
`rhbm_tests_joint_component_numerical`. These run from the ordinary
`RHBM-GEM-TEST` executable. `joint_component_regression` verifies eight frozen
cases, and `joint_component_physical_smoke` checks double/float32 generation
and CIF/MRC. All required inputs live in `tests/fixtures/joint_component/`; no
historical archive or external model is loaded.

Run the bounded Joint gates with `joint:contract`, `joint:workflow`, and
`joint:numerical` labels. Large Gram and 256-atom diagnostics are separately
listed under `joint:scalability`. `RHBM_GEM_ENABLE_JOINT_EXTENDED_TESTS=ON`
adds the 168-atom float32 case and additional catalog cases to the offline
qualification lane. `RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON` builds the
separate precision/certification targets and tests. Both options default to
OFF; use `joint:offline` for that lane. High precision code must not be linked
into `rhbm_tests` or `joint_test_support`. See the
[runtime guide](../docs/developer/joint-component-runtime.md) for commands and
acceptance policy.

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
| `PerformanceLogging_test.cpp` | Normal performance summary output and quiet-mode behavior |

Preserve the existing `EstimatorSecondStageDefenseTest` suite and case names.
All second-stage files belong to `CORE_ESTIMATOR_TEST_SOURCES` and run through the single
`rhbm_tests_core_estimator` CTest group. Do not add per-file CTest groups using
this shared suite filter: each would repeat the entire suite. Mixed cases that
verify production acceptance or rollback remain with the production behavior
they exercise.

Boundary reference tests cover ordinary/rescue gates, unavailable evidence,
suspicious-correction early exits, strict rejection, and single correction-delta
evaluation. Runtime log checks and debug-level scheduling behavior remain with
the production fitting tests.

Use `tests/support/SecondStageTestSupport.hpp/.cpp` and its `second_stage_test`
namespace for fixture builders and assertions shared across these files. Keep
single-file helpers in that file's anonymous namespace and support-only helpers
private to the support implementation. Build fresh model/solver state per case;
do not share mutable fixtures. Use `detail` for `rhbm_gem::core::detail` throughout
these tests. `EstimatorTester_test.cpp` retains its workflow tests and fixtures.

Compile test support directly into `rhbm_tests`, outside the suite-discovery source
lists. Work-count and frozen-background capture helpers serve permanent numerical
tests; solver failure replay remains in `SolverFailureCapture`.

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

The production outcome codec accepts joint JSON schema 3 inside SQLite v17.
Metadata tests cover actual normalization divisors, input/build fingerprints,
coefficient units and rejection of old documents. CLI, C++ and Python workflow
checks export after deleting their temporary model/map inputs. The installed
consumer smoke target now executes the linked program and checks its exit status,
including capture/export with unknown in-memory input provenance.

Release validation runs all default CTests with Python bindings, then the existing
`joint:extended|joint:offline` lane, plus a testing-disabled installation and
consumer execution. A source copy without Git metadata must also build/install
and run the consumer. Keep the frozen fixture packages and current parity
thresholds unchanged.

Current Joint experiments, build options, validation harnesses and permanent
test owners are listed in the
[Joint experiment inventory](../docs/developer/joint-experiments.md). The sole
benchmark workflow is `tests/integration/joint_benchmark.py`; it includes the
FixedNeighbor search/solve, workflow/postprocess and complete-command scopes.
Partial-selection, weak-halo and noise/mismatch measurements remain current
feature, diagnosis and research tools. Closed campaign drivers and their
machine-readable results are not retained in the current tree.
