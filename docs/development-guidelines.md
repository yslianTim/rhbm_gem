# Development Rules and Existing Conventions

This document describes the development habits that are already established in this repository. 
When adding features, fixing bugs, or refactoring code, prefer extending these existing patterns.
If any rule changes, update both this document and the README.

## 1. Project Structure

- Public interfaces live under `include/core`, `include/data`, and `include/utils`, while implementations live in the corresponding `src/core`, `src/data`, and `src/utils` directories. New classes should generally keep a one-to-one header/source pairing.
- The CLI entry point remains in `src/main.cpp`. Command orchestration belongs in `core`, file formats and database access belong in `data`, and shared algorithms and utilities belong in `utils`.
- Python bindings are centralized in `bindings/`, example scripts in `examples/python/`, unit tests in `tests/`, and test fixtures in `tests/data/`.
- Third-party dependencies are isolated under `third_party/`. Unless you are maintaining a vendored dependency, do not place project logic there.

## 2. Build and Dependency Management

- The project baseline is CMake + C++17 with compiler extensions explicitly disabled. New code should remain portable within C++17.
- Compiler warnings are managed centrally through the top-level `CompilerFlags` interface library. New targets should follow that pattern instead of introducing scattered warning settings.
- Optional features are controlled through `AUTO` / `ON` / `OFF` CMake variables, such as OpenMP, Boost, and ROOT. New optional dependencies should follow the same convention.
- Dependency resolution prefers system packages first, then falls back to `third_party/` when available. If a dependency has no bundled fallback, document that explicitly.
- Install rules, output directories, and exported targets are currently managed consistently through CMake. New libraries, executables, and bindings should include matching install and export rules.

## 3. C++ Style

- Headers use `#pragma once`.
- Types, classes, and functions use PascalCase. Member variables use the `m_` prefix. Command parameters are typically grouped into nested `Options` structs.
- The existing codebase makes heavy use of `enum class`, `explicit`, `override`, `= default`, `nullptr`, brace initialization, and `const auto`. New code should follow the same style.
- Most project types currently live in the global namespace. Unless there is a deliberate repository-wide namespace refactor, follow that existing approach rather than introducing a separate local namespace scheme.
- Parameterless functions are commonly written with `(void)`. When modifying an existing file, preserve consistency within that file.
- Comments are primarily written in English and are used sparingly for ownership, lifetime, format details, or non-obvious logic. Avoid comments that merely restate the code.

## 4. Object Lifetime and Data Flow

- Prefer `std::unique_ptr` by default to represent single ownership. Use `std::shared_ptr` only in genuinely shared contexts such as `DataObjectManager`.
- Data objects move between command and data layers through managers and `key_tag` identifiers rather than through scattered raw-pointer ownership.
- When exposing getters, preserve const-correctness whenever possible. Add mutable overloads only when needed.
- If a shared state already uses mutex-based protection, such as `DataObjectManager` or `Logger`, new related logic should preserve the same locking boundaries and responsibility split.

## 5. Command Layer Pattern

- New CLI commands should derive from `CommandBase` and extend subcommand options through `RegisterCLIOptionsExtend`.
- Each command keeps user input inside `Options : CommandOptions`. Setters are responsible for validation, normalization, and error flag updates where needed.
- Command registration follows the static `CommandRegistrar<...>` pattern inside an anonymous namespace in the `.cpp` file. New commands should keep using that pattern instead of scattering registration logic.
- `Execute()` usually performs data construction or pre-validation first, returns `false` quickly on failure, and only then continues into the main workflow.
- Recoverable input problems generally log a warning and fall back to a safe default. Unrecoverable conditions should mark validation failure or throw an exception.

## 6. Data Processing, Files, and Database Access

- File I/O selects readers, writers, and factories by file extension. When adding support for a new format, update both the factory registry and the corresponding format classes.
- The command layer should cooperate with `DataObjectManager` and `DatabaseManager` rather than embedding file parsing, SQLite logic, or DAO details directly inside commands.
- Error messages should include context such as file name, key tag, and operation stage instead of using vague failure messages.
- If a function can handle failure locally and report it clearly to the user, prefer logging through `Logger`. If the caller must make the decision, throw `std::runtime_error`.

## 7. Logging and Runtime Diagnostics

- Runtime output should go through `Logger`. Outside of `Logger` itself or very limited utility code, avoid writing directly to `std::cout` or `std::cerr`.
- Many core functions log their entry point at `Debug` level. New high-value flows may keep this traceability, but avoid repetitive or low-signal messages.
- Expensive operations can use existing tools such as `ScopeTimer` for diagnostics instead of introducing a separate timing mechanism.
- The log level meanings are already fixed as Error, Warning, Notice, Info, and Debug. New output should match the correct level.

## 8. Testing and Regression Strategy

- Unit tests use GoogleTest, test files follow the `*_test.cpp` naming pattern, and they are listed explicitly in `tests/CMakeLists.txt`.
- New parsers, data transformations, algorithms, and bug fixes should usually come with corresponding tests. The repository already uses many small fixture-driven regression tests.
- Test data should stay small and focused under `tests/data/`, and file names should describe the scenario directly, such as missing fields, malformed formatting, or special tokens.
- I/O-related tests should prefer temporary files or temporary directories and clean up after themselves instead of leaving generated artifacts in the repository.
- Assertion style should follow the existing pattern: use `ASSERT_*` for required preconditions and `EXPECT_*` for detailed verification.

## 9. Python Bindings and Examples

- The Python extension module name is fixed as `rhbm_gem_module`. Do not introduce a second module name.
- Bindings mainly expose command classes, enums, and `Set*` interfaces. When adding new public command capabilities, evaluate whether they should also be exposed to Python.
- Scripts under `examples/python/` are intended to be directly runnable and useful for quick installation verification. New examples should preserve that practical orientation instead of becoming incomplete snippets.

## 10. Documentation Sync Requirements

- If you add CMake parameters, build modes, public modules, Python interfaces, commands, or test entry points, update `README.md` as part of the change.
- If you change any rule documented here, update this file with both the new behavior and the reason for the change so the document does not drift behind the implementation.
