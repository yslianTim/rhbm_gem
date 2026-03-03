# Development Guidelines

This document defines the preferred engineering rules for new code and for any existing code that is being modified. It is not a literal description of every legacy pattern already present in the repository. Legacy code does not need immediate repo-wide cleanup, but touched code should move toward these rules whenever practical.

If a rule changes, update this document and any affected user-facing documentation, especially `README.md`.

## 1. Purpose and Rule Levels

- `[Required]` Rules that new code and modified code are expected to follow.
- `[Recommended]` Rules that should be applied when working in the affected area unless there is a good reason not to.
- `[Legacy note]` Existing historical patterns that may still appear in the codebase but are not preferred for new development.
- The goal of these guidelines is engineering quality: maintainability, correctness, diagnosability, testability, and safe evolution of public interfaces.

## 2. Project Structure and Module Boundaries

- `[Required]` Public C++ interfaces belong under `include/`, implementations belong under the matching `src/` area, Python bindings belong under `bindings/`, runnable examples belong under `examples/python/`, tests belong under `tests/`, and test fixtures belong under `tests/data/`.
- `[Required]` Keep module responsibilities clear: command orchestration in `core`, file formats and persistence in `data`, and reusable utilities or algorithms in `utils`.
- `[Required]` New reusable utility or algorithm modules should prefer a focused namespace rather than expanding the global namespace.
- `[Required]` Project-owned public types in `core` and `data` must live under the `rhbm_gem` namespace. New code must not introduce additional global `core` or `data` APIs.
- `[Required]` Do not use `using namespace` in headers.
- `[Recommended]` Keep a one-to-one header/source pairing for new classes unless there is a strong reason to group multiple small related types together.
- `[Recommended]` Do not place project logic in `third_party/` unless you are intentionally maintaining vendored dependency code.

## 3. Build, Dependencies, and Feature Gating

- `[Required]` New code must remain portable within the project's CMake + C++17 baseline and must not depend on compiler extensions.
- `[Required]` New or modified code must compile cleanly under the existing warning set. Do not add broad warning suppressions to avoid fixing real issues.
- `[Required]` Any target-local compiler flag exception must be narrow in scope and documented with the toolchain-specific reason.
- `[Required]` Public headers should be self-contained and should not rely on indirect includes from unrelated headers.
- `[Required]` Optional dependency code must stay behind its feature gate and must degrade cleanly when the feature is unavailable or explicitly disabled.
- `[Required]` New libraries, executables, and bindings must include matching install and export behavior when they are part of the supported build surface.
- `[Recommended]` Reuse the shared `CompilerFlags` pattern instead of inventing per-target warning policies.
- `[Recommended]` Prefer existing project conventions for dependency discovery: use system packages when available and bundled fallbacks only when intentionally supported.

## 4. C++ Design Rules

- `[Required]` Prefer readability and local consistency over low-value style trivia.
- `[Required]` Headers use `#pragma once`.
- `[Required]` Types, classes, and functions use clear PascalCase names. Member variables use the `m_` prefix.
- `[Required]` New code should preserve const-correctness and narrow interfaces to the minimum mutable surface required.
- `[Required]` Avoid exposing unnecessary implementation details in public headers.
- `[Recommended]` Continue using `enum class`, `explicit`, `override`, `= default`, `nullptr`, brace initialization, and `const auto` where they improve clarity.
- `[Recommended]` Comments should explain why something is done, not restate what the code already makes obvious.
- `[Recommended]` Add comments when they clarify ownership and lifetime, thread-safety invariants, feature-gated behavior, file format assumptions, compatibility workarounds, or non-obvious algorithmic choices.

## 5. Ownership and Lifetime Rules

- `[Required]` Ownership must be explicit in interfaces.
- `[Required]` Default to `std::unique_ptr` for single ownership.
- `[Required]` Use `std::shared_ptr` only when lifetime is genuinely shared across owners, not as a convenience transport type.
- `[Required]` Owning raw pointers are not allowed in new interfaces unless there is a strong interop constraint and the ownership contract is documented.
- `[Required]` Non-owning pointers are acceptable only when null is a meaningful state. Use references for required non-null inputs.
- `[Required]` If a public API returns a raw pointer, document that it is a borrowed view rather than ownership transfer.
- `[Recommended]` Prefer APIs that make ownership transfer obvious at the call site.
- `[Recommended]` Keep lifetime-sensitive comments close to the owning field or factory boundary rather than scattering them across callers.

## 6. Error Handling and Logging

- `[Required]` Error handling responsibilities must be clear by layer.
- `[Required]` Low-level helpers, parsers, and library-like components should prefer throwing exceptions with useful context when the caller must decide how to recover.
- `[Required]` Command and application boundaries should catch failures, convert them into user-facing diagnostics, and decide whether to return `false`, abort the workflow, or continue with a safe fallback.
- `[Required]` Do not both log and rethrow the same failure at multiple layers unless the added log message introduces new decision-relevant context.
- `[Required]` Error messages must identify the operation, the input or context, and the failure reason.
- `[Required]` Recoverable input issues may log a warning and normalize to a safe default. Contract violations and data inconsistency should fail fast.
- `[Required]` `Logger` is the primary channel for application-facing runtime output.
- `[Recommended]` Use `Debug` logging selectively for high-value flow tracing rather than logging every function entry.
- `[Recommended]` Prefer logs that describe state changes, important parameters, expensive-step boundaries, or actionable failure context.
- `[Recommended]` Use direct `std::cout` or `std::cerr` only in low-level helpers, third-party integration boundaries, tests, or short-lived diagnostics where using `Logger` would be inappropriate.
- `[Legacy note]` Some existing helpers still emit direct stream output or overly chatty trace logs. New code should move away from that style.

## 7. Command and Data Layer Responsibilities

- `[Required]` New CLI commands should derive from `CommandBase` and keep command-specific configuration in an `Options : CommandOptions` structure.
- `[Required]` Setters are responsible for validating user input, normalizing values when appropriate, and marking invalid configurations before execution starts.
- `[Required]` `Execute()` should keep a clear phase boundary: validate and construct prerequisites first, then run the main workflow.
- `[Required]` File parsing, database access, and DAO details belong in the data layer, not inside command orchestration logic.
- `[Required]` Support for new file formats must update both the extension-based factory selection path and the corresponding format implementation.
- `[Recommended]` Continue using the existing `CommandRegistrar<...>` registration model for new commands to keep registration localized and predictable.
- `[Recommended]` Managers should remain the boundary objects that move data between commands, files, and persistence layers.

## 8. Testing and Regression Policy

- `[Required]` New features, behavior changes, and bug fixes must include tests unless the change is purely documentation-only.
- `[Required]` Tests must cover not only happy paths, but also invalid input, malformed files, boundary conditions, and fallback behavior where those paths are part of the supported behavior.
- `[Required]` New test files must follow the `*_test.cpp` naming pattern and must be added explicitly to `TEST_SOURCES` in `tests/CMakeLists.txt`.
- `[Required]` Tests must be deterministic and must not depend on machine-specific state or previously generated local artifacts.
- `[Required]` Tests that create files or directories must clean them up.
- `[Required]` If a feature is controlled by optional dependencies or feature modes, add coverage for the enabled or disabled behavior that the change affects.
- `[Recommended]` Use minimal fixtures that reproduce only the behavior under test.
- `[Recommended]` Resolve fixture paths relative to the test source, for example via `std::filesystem::path(__FILE__).parent_path()`.
- `[Recommended]` Use temporary files or temporary directories for generated output.
- `[Recommended]` Use `ASSERT_*` for preconditions and `EXPECT_*` for detailed follow-up verification.
- `[Recommended]` Regression test names should describe the user-visible or bug-specific scenario, not just the internal function being exercised.

## 9. Public API, CLI, and Python Surface Governance

- `[Required]` Treat the following as public interfaces: headers under `include/`, CLI subcommands and options, the `rhbm_gem_module` Python bindings, install/export behavior, and documented example workflows.
- `[Required]` Any change to a public surface must trigger a review of compatibility impact.
- `[Required]` When public behavior changes, update the affected surfaces together: documentation, examples, bindings, install or export rules, and regression coverage.
- `[Required]` If a command is intentionally supported from Python, keep the exposed binding surface aligned with the supported C++ behavior.
- `[Required]` If a command or capability is intentionally CLI-only, document that explicitly rather than leaving the omission ambiguous.
- `[Recommended]` Prefer adding smoke-style coverage for important public workflows that are likely to be used from documentation or examples.

## 10. Documentation Sync

- `[Required]` Update `README.md` whenever a change affects build options, dependency behavior, install or export behavior, public commands, Python bindings, examples, or expected workflow outputs.
- `[Required]` Update this document whenever project-wide engineering rules change.
- `[Recommended]` Keep contributor-facing guidance practical and repository-specific rather than turning it into a generic language style guide.

## 11. Long-Term Engineering Direction

- `[Recommended]` Adopt an explicit formatting policy such as `clang-format` once the team is ready to enforce it consistently.
- `[Recommended]` Add a static-analysis baseline such as `clang-tidy` when the false-positive budget and maintenance cost are acceptable.
- `[Recommended]` Build out feature-matrix CI coverage for representative `AUTO`, `ON`, and `OFF` combinations of ROOT, Boost, and OpenMP.
- `[Recommended]` Evaluate a project-wide namespace strategy if the public API surface continues to grow.

## 12. Legacy Notes

- `[Legacy note]` Some existing files use parameterless `(void)` style. Preserve local consistency when touching those files, but do not treat this as a priority rule for new design.
- `[Legacy note]` Some existing low-level helpers still use direct stream output instead of `Logger`. New code should not copy that pattern without a clear boundary-related reason.
