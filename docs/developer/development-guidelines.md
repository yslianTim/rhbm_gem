# Development Guidelines

This document defines the engineering principles for new code and for any existing code being modified.
The goal is stable, maintainable, and diagnosable software.

## 1. Rule Levels

- `[Required]` Expected for all new code and modified code.
- `[Recommended]` Should be followed unless there is a clear reason not to.

## 2. Engineering Priorities

- `[Required]` Correctness and safety are higher priority than local convenience.
- `[Required]` Changes should improve or preserve maintainability, testability, and diagnosability.
- `[Required]` Prefer simple designs with clear ownership and clear boundaries.
- `[Recommended]` Favor incremental, reviewable changes over broad coupled rewrites.

## 3. Project Structure and Module Boundaries

- `[Required]` Keep public C++ interfaces in `include/` and implementations in `src/`.
- `[Required]` Keep responsibilities separated: command orchestration in `core`, file/persistence logic in `data`, reusable helpers in `utils`.
- `[Required]` Keep tests in `tests/` and test fixtures in `tests/fixtures/`.
- `[Required]` Keep test sources aligned to domain folders (`tests/core`, `tests/data`, `tests/utils`, `tests/integration`) and intent-specific subfolders where applicable.
- `[Required]` For `src/core`, place command orchestration in `src/core/command/`, workflow implementation in `src/core/workflow/`, painter implementation in `src/core/painter/`, and private helpers in `src/core/internal/`.
- `[Required]` Treat `internal/` path placement as the internal-API marker; do not add redundant `*Internal` filename suffixes inside `internal/`, and prefer moving private headers into an `internal/` subtree over suffix-based labeling.
- `[Required]` For workflow code, keep implementation units (`*.cpp`) in `src/core/workflow/`; place private workflow headers in `src/core/internal/workflow/` when they are not part of shared workflow-facing contracts.
- `[Required]` Keep reusable automation scripts under `scripts/`; do not add scripts at repository root.
- `[Required]` Project-owned public types in `core` and `data` must use namespace `rhbm_gem`.
- `[Required]` Do not use `using namespace` in headers.
- `[Recommended]` Keep header/source pairing clear and predictable for new modules.

## 4. Build and Dependency Policy

- `[Required]` Keep code portable within CMake + C++17 baseline.
- `[Required]` Compile cleanly under project warning settings; do not hide real issues with broad suppressions.
- `[Required]` Keep optional dependencies behind explicit feature gates and provide clean behavior when disabled.
- `[Required]` Public headers must be self-contained.
- `[Recommended]` Reuse existing build patterns (`CompilerFlags`, dependency discovery conventions) instead of creating one-off policies.

## 5. C++ Design and Ownership Rules

- `[Required]` Prefer clear interfaces and explicit ownership.
- `[Required]` Default to `std::unique_ptr` for single ownership.
- `[Required]` Use `std::shared_ptr` only when ownership is truly shared.
- `[Required]` Avoid owning raw pointers in new interfaces unless interop constraints require it and are documented.
- `[Required]` Use references for required non-null inputs; use pointers only when null is meaningful.
- `[Required]` Keep const-correctness and avoid exposing unnecessary implementation detail in public headers.
- `[Recommended]` Use modern C++ features (`enum class`, `explicit`, `override`, `= default`, `nullptr`) where they improve clarity.

## 6. Error Handling and Logging

- `[Required]` Error-handling responsibility must be clear by layer.
- `[Required]` Data-layer and library-like components should throw exceptions with actionable context.
- `[Required]` Command/application boundaries should translate failures into user-facing diagnostics and workflow decisions.
- `[Required]` Avoid duplicate logging of the same failure unless each layer adds new decision-relevant context.
- `[Required]` Error messages should include operation, input/context, and reason.
- `[Required]` Use `Logger` as the main runtime output channel.
- `[Recommended]` Keep debug logging focused on meaningful state transitions and expensive boundaries.

## 7. Command and DataObject I/O Principles

See:

- [`./architecture/command-architecture.md`](./architecture/command-architecture.md)
- [`./adding-a-command.md`](./adding-a-command.md)
- [`./architecture/dataobject-io-architecture.md`](./architecture/dataobject-io-architecture.md)
- [`./adding-dataobject-operations-and-iteration.md`](./adding-dataobject-operations-and-iteration.md)

Command-layer principles:

- `[Required]` Commands orchestrate workflows; they should not own file-format parsing or persistence implementation details.
- `[Required]` Command options should be validated and normalized before main execution.
- `[Required]` Keep command execution phases explicit: validate -> prepare dependencies -> execute workflow -> report result.
- `[Required]` Avoid exposing internal lifecycle hooks, debug-only accessors, or registration internals as public command APIs.

DataObject I/O principles:

- `[Required]` Keep file-format support definitions centralized in a single registry-style source of truth.
- `[Required]` Use explicit fixed routing with centralized registry-driven descriptor lookup.
- `[Required]` Reader/writer and DAO paths should fail fast with clear exceptions when contracts are violated.
- `[Required]` Keep schema versioning, validation, and migration centralized in schema-management components.
- `[Required]` Keep normal save/load hot paths separate from bootstrap/repair logic.
- `[Required]` Keep database transaction boundaries at a manager/service boundary; DAO methods should remain focused on persistence mapping.
- `[Required]` Keep legacy compatibility and migration helpers internal (`src/`) unless explicitly promoted to public API.

Change-integration principles:

- `[Required]` New commands must update CLI registration, Python bindings, tests, and developer docs together unless intentionally internal-only.
- `[Recommended]` Keep manager classes as boundary objects that connect commands, file I/O, and persistence.

## 8. Testing and Regression Policy

- `[Required]` New features, behavior changes, and bug fixes require tests (except documentation-only changes).
- `[Required]` Cover both happy paths and failure paths that are part of supported behavior.
- `[Required]` Keep tests deterministic and independent of machine-local state.
- `[Required]` Add new `*_test.cpp` files to the matching grouped file under `tests/cmake/`.
- `[Required]` Every test target declared through `add_rhbm_gtest_target(...)` must carry both `domain:*` and `intent:*` CTest labels.
- `[Required]` Every file referenced by `TestDataPath(...)` must exist under `tests/fixtures/` and be tracked by git.
- `[Required]` Schema/persistence changes require regression coverage for bootstrap or migration behavior, rollback behavior on failure, and round-trip correctness.
- `[Recommended]` Prefer searchable suite names that encode intent (for example `DataObjectSchemaMigrationTest`) over generic suite names.
- `[Recommended]` Use minimal fixtures and temporary directories/files for generated outputs.

## 9. Public API, CLI, and Python Surface Governance

- `[Required]` Treat headers under `include/`, CLI subcommands/options, Python bindings, and documented workflows as public surfaces.
- `[Required]` Review compatibility impact for all public-surface changes.
- `[Required]` Keep C++ behavior, CLI behavior, Python binding behavior, and documentation aligned.
- `[Recommended]` Add smoke-style tests for key documented user workflows.

## 10. Documentation and Change Hygiene

- `[Required]` Update `README.md` when changes affect build options, dependencies, install/export behavior, public commands, bindings, or examples.
- `[Required]` Update this guideline document when project-wide engineering rules change.
- `[Required]` Keep `CommandList.def`, runtime registration, Python bindings, and public docs aligned.
- `[Recommended]` Prefer concise, principle-driven guidance over temporary implementation details.

## 11. Quality Checks and CI Alignment

- `[Required]` Keep repository guard checks passing through `lint_repo` (structure, hygiene, logging style, fixture tracking, command-manifest sync, and absolute-path checks).
- `[Recommended]` Run targeted formatter checks with `scripts/developer/run_clang_format_check.sh` when touching covered painter/parser paths.
- `[Recommended]` Run targeted static analysis with `scripts/developer/run_clang_tidy_check.sh` (baseline mode when appropriate) when touching covered painter/parser paths.
- `[Recommended]` For broad or risky changes, run `scripts/developer/run_ctest_with_classification.sh` in addition to normal `ctest`.
