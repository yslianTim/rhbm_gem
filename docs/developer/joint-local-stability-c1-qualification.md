# Joint Local Stability C1 Qualification

Date: 2026-10-29. This report records the staged C1 split-only diagnostics qualification. A stage is marked complete only after its implementation or experiment has actually run.

## Initial repository state

- Branch: `develop`.
- Initial HEAD and local `origin/develop`: `fddd672f667d680cd88ab05704f7e3786910da76`.
- Working tree: clean; `git diff --check` passed.
- Existing detached worktree: `/Users/yslian/.codex/worktrees/d178/em_map_analysis`; left untouched.
- The local tracking ref is not evidence of the live GitHub branch. The remote will be checked read-only; no push is authorized.
- No Notion documents were accessed or changed.

## S0 — Split diagnostic regression contract

**Status: Complete.** Added an internal split metric accumulator and C++ regression coverage without changing the persisted result type, JSON schema, legacy classifier, thresholds, point-quality, or accepted fitting state.

The regression contract defines prediction normalization as `max(1, ||prediction_0||_2)`, A/C parameter scaling using the N2 reference design column norms and observation scale, and log-width infinity norm. It records the endpoint identity that maximizes each drift and seed metric for every direction. Qualified seed vectors are checked against the vector decomposition; no scalar norm additivity is assumed. An unqualified null reference yields `Unavailable` with no split metrics. An incomplete direction yields `Incomplete` and cannot claim full seed robustness.

The tests cover exact-fit identity; distinct reoptimization and seed differences; representative T0, T4, 300-target, and 1,000-target numeric values; failed and budget-exhausted references; replay failure; incomplete directions; active-bound one-sided probes; zero effective displacement; nonfinite reference values; reference column normalization; and preservation of the legacy `Unstable` classification. These representative-value fixtures test metric arithmetic, not a rerun of the original 300/1,000-target solver workloads.

The 11 new `JointSplitDiagnosticsTest` cases passed in the Release test binary. Production candidate execution, B1/B2 benchmark, difficult-case validation, and promotion have not yet been performed at this stage.
