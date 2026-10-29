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

## S2 — B2 workspace reuse experiment

**Status: Complete; B2 not retained for promotion.** B2 reuses one endpoint prediction buffer plus five difference/decomposition buffers while calculating endpoint metrics. It uses the same N2 searches, endpoint qualification, directions, replay policy, and calculations as B1.

Across 3/30/300/1,000 targets, B2 matched B1 split status, N2 qualification, counts, high-budget endpoint, all aggregate Drift/Seed metrics, endpoint fingerprints, and Accepted State fingerprints. A direct 30-target C++ regression also requires equality of each selected block's qualification, reason, all split scalar metrics, endpoint identities, and full accepted arrays.

| Targets | B1 wall | B2 wall | B2 vs B1 | B1 stability CPU | B2 stability CPU | B1 peak RSS | B2 peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 16.410 ms | 16.419 ms | +0.05% | 12.654 ms | 12.655 ms | 23.675 MB | 23.904 MB |
| 30 | 404.536 ms | 403.111 ms | -0.35% | 294.967 ms | 294.804 ms | 35.095 MB | 35.275 MB |
| 300 | 1.23680 s | 1.23659 s | -0.02% | 267.403 ms | 267.388 ms | 87.687 MB | 87.835 MB |
| 1,000 | 3.63787 s | 3.63668 s | -0.03% | 388.818 ms | 387.958 ms | 236.306 MB | 236.175 MB |

The wall differences change sign and remain below 0.5%; B2 does not reduce stability CPU, and its preallocated scratch can raise peak RSS. The implementation adds mode branches, six reusable vectors, and separate buffer lifetime management without measurable production value. Decision: retain B1 semantics only if C1 is later promoted; remove the B2 runtime branch/workspace during final cleanup.

## S3 — Controlled B0/B1/B2 benchmark

**Status: Complete for the required 3/30/300/1,000-target matrix.** Every route used the same frozen chain fixture, initialization (`B=0.55`), global search, KKT/bounds, selected-block policy (cap 8, reserve active-bound claim), Release build, single-thread environment, and 1 warmup + 5 measured runs. The benchmark linked the same Release source with test instrumentation enabled only to select B0/B1/B2. Numeric threads were `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `VECLIB_MAXIMUM_THREADS=1`.

| Targets | B0 wall | B1 wall | B2 wall | B1 cost vs B0 | B0/B1/B2 stability CPU | B1 Null CPU | B1 N2 eval/update/replay | Perturbation replay checks | B0/B1/B2 Peak RSS |
|---:|---:|---:|---:|---:|---|---:|---|---|---|
| 3 | 14.241 ms | 16.410 ms | 16.419 ms | +15.23% | 10.492 / 12.654 / 12.655 ms | 2.122 ms | 6 / 2 / 2 | 6 (each mode) | 23.708 / 23.675 / 23.904 MB |
| 30 | 352.951 ms | 404.536 ms | 403.111 ms | +14.62% | 244.628 / 294.967 / 294.804 ms | 49.234 ms | 18 / 4 / 6 | 18 (each mode) | 35.013 / 35.095 / 35.275 MB |
| 300 | 1.18684 s | 1.23680 s | 1.23659 s | +4.21% | 215.668 / 267.403 / 267.388 ms | 52.446 ms | 24 / 4 / 6 | 18 (each mode) | 87.507 / 87.687 / 87.835 MB |
| 1,000 | 3.54848 s | 3.63787 s | 3.63668 s | +2.52% | 305.936 / 388.818 / 387.958 ms | 84.029 ms | 24 / 10 / 6 | 18 (each mode) | 235.241 / 236.306 / 236.175 MB |

| Targets | Stability workspace B0/B1 | B1 split workspace estimate | Global memory preflight (all modes) | Selected blocks | B1/B2 qualified blocks | Accepted-state FNV-1a64 (B0/B1/B2) |
|---:|---:|---:|---:|---:|---:|---|
| 3 | 5.435 / 5.444 MB | 5.444 MB | 0.447 MB | 1 | 1 / 1 | `a7980d92f9da5a52` |
| 30 | 23.717 / 23.751 MB | 23.751 MB | 3.930 MB | 3 | 3 / 3 | `b18509add52fabf9` |
| 300 | 53.378 / 53.412 MB | 53.412 MB | 36.716 MB | 3 | 3 / 3 | `3535e5d2710ef2ba` |
| 1,000 | 130.208 / 130.242 MB | 130.242 MB | 121.649 MB | 3 | 3 / 3 | `7483b6f0d4edefa8` |

At each scale the B1 and B2 maxima were identical: reoptimization prediction drift `2.8704e-9`, `3.2813e-9`, `1.9387e-7`, `1.9380e-7`; seed prediction robustness `3.0348e-13`, `1.1536e-13`, `3.9139e-12`, `4.0445e-13`. High-budget width, scaled A/C, and prediction differences were zero for all four. The B0/B1/B2 accepted-state fingerprints match at all scales, and each mode's five measured runs repeated its own fingerprint. Direct array tests cover 2 and 30 targets.

B1 is within the exploratory wall gate at 30 targets and above the 15% reference by 0.23 percentage points at 3 targets in this final counter-instrumented run; the earlier 3-target run measured +14.3%. Peak RSS changes remained under 1% in the final matrix. The added workspace estimate was 8.4 KB at 3 targets and 33.6 KB at 30+ targets. Existing memory preflight, global KKT, width KKT, sweep counts, endpoint certification, and independent replay were unchanged. At 3 and 30 targets, Null Search is a visible fixed share of the short end-to-end fit; at 300 and 1,000 targets, global search dominates total wall time.

This matrix covers the four required scale workloads. High-overlap, weak-width, active-bound, and noise-only workload measurements are part of S4 and must not be inferred from the chain fixture.
