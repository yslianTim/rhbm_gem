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

## S1 — B1 minimal split-only candidate

**Status: Complete.** Added a transient candidate path alongside the existing local assessment. It executes a baseline-seeded N2 search with the existing local profile solver, constraints, scale, bounds, and precision; requires strict local A/C and width KKT plus Strict Independent Replay; and compares against a separate doubled-budget reference. N2 and perturbation search use the same default profile/update limits. A stopped budget, rank/feasibility issue, failed replay, or high-budget disagreement makes the reference unavailable and prevents Seed Robustness. The accepted state and legacy status remain untouched.

B1 streams qualified perturbation endpoints into the metric accumulator. B2 is a workspace-reuse experiment that uses preallocated prediction and difference buffers. Neither candidate is persisted or passed into Target Evidence.

Release regression: 24 split metric and parameter stability tests passed. A direct B0/B1/B2 `SearchFixedNeighborComponent` regression compared complete `η`, `β`, prediction, residual, objective, sweep count, legacy status/reason, and split qualification/endpoint identities. Arrays and objective were exactly equal for the 2-target fixture.

The first controlled 3-target benchmark (Release; single numeric thread; 1 warmup + 5 measured runs) produced:

| Path | End-to-end wall median | Stability CPU median | N2 CPU median | Peak RSS | Split |
|---|---:|---:|---:|---:|---|
| B0 | 14.196 ms | 10.365 ms | — | 23.69 MB | — |
| B1 | 16.303 ms | 12.561 ms | 2.101 ms | 23.81 MB | 1/1 Qualified |
| B2 | 16.381 ms | 12.617 ms | 2.113 ms | 23.82 MB | 1/1 Qualified |

All three runs had accepted-state fingerprint `a7980d92f9da5a52`; the B1 fit used 6 N2 profile evaluations, 2 accepted updates, and 2 strict replay checks. In the 2-target array regression, B0/B1/B2 were directly equal. Full workload matrix and B2 retention decision are deferred to S2/S3.
