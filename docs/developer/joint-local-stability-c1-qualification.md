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

**Status: Complete.** Added a split metric accumulator and C++ regression coverage without changing the persisted result type, JSON schema, legacy classifier, thresholds, point-quality, or accepted fitting state. S5 moved the metric helper to test-only support after C1 was not promoted.

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

**Status: Partial.** The required 3/30/300/1,000-target chain matrix and a separate 3-target noise-only matrix ran with identical B0/B1/B2 fixture, policy, Release build, single-thread environment, and 1 warmup + 5 measured runs. High-overlap, weak-width and active-width did not receive the same full B0/B1/B2 performance matrix. The benchmark linked the same Release source with test instrumentation enabled only to select B0/B1/B2. Numeric threads were `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `VECLIB_MAXIMUM_THREADS=1`.

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

This table covers the four required scale workloads. Noise-only has its own paired matrix below. High-overlap, weak-width, and active-bound received numerical candidate tests but not the same paired B0/B1/B2 performance matrix.

The additional 3-target noise-only paired benchmark used the same 1 warmup + 5 measured runs. B0/B1/B2 accepted-state fingerprint was `60c62facd75f39d2` in all modes; direct B0/B1/B2 arrays, objective, legacy statuses, and Target point-quality were also equal in the regression test.

| Noise-only, 3 targets | B0 | B1 | B2 |
|---|---:|---:|---:|
| Wall median | 29.470 ms | 38.472 ms (+30.5%) | 38.467 ms (+30.5%) |
| Stability CPU | 20.159 ms | 29.187 ms | 29.183 ms |
| Added Null Search CPU | — | 9.016 ms | 8.997 ms |
| Perturbation Replay checks | 6 | 6 | 6 |
| Null eval / accepted updates / replay | — | 22 / 18 / 2 | 22 / 18 / 2 |
| Process peak RSS | 23.806 MB | 24.297 MB | 24.183 MB |
| Qualified / incomplete Split blocks | — | 0 / 1 | 0 / 1 |

The 3-target chain run was +15.23% in the final counter-instrumented repeat (an earlier repeat was +14.3%); noise-only was +30.5%. Both exceed the exploratory 10–15% wall target. RSS stayed below the 10% research gate, but the small-workload wall gate fails. B2 again produced no meaningful improvement over B1. The B0/B1/B2 matrix therefore does not qualify C1 engineering cost for small or noise-only workloads.

## S4 — Difficult-case scientific validation

**Status: Partial qualification; the required synthetic defenses ran, while there is no real-map or complete difficult-case performance matrix.** Added deterministic tests exercised the B1 path directly; legacy status, accepted state, endpoint certification, replay, and Target point-quality were compared against B0 where the case had a qualified fit.

| Case | Actual result |
|---|---|
| Exact / low-noise synthetic | B0/B1 accepted arrays and objective were exactly equal. At the default noise amplitude, maximum absolute A/C/B errors were approximately `6.5e-5 / 1.5e-5 / 3e-6`; prediction RMSE was about `5e-6`. |
| 10× higher noise | B0/B1 arrays remained exactly equal. Maximum absolute A/C/B errors were approximately `6.54e-4 / 1.48e-4 / 3.3e-5`; prediction RMSE was about `5.1e-5`. Both cases met the `1e-2` absolute parameter-accuracy target. |
| Width and A/C scale | Nine exact fixtures with true width `.4/.8/1.0` and A/C amplitude multiplier `.01/1/10` all had qualified baseline, Null reference and Split metrics. This checks finite normalized metrics across the exercised ranges, not all possible scales. |
| Noise-only / weak width | The B0/B1/B2 3-target run converged and certified the global endpoint, but B1 N2 exhausted its 6-update budget (18 accepted updates reported across standard and high-budget searches); split status was `Unavailable`, with no Drift or Seed metric. The old raw evidence independently records a high-budget width difference of `2.398e-5`; its perturbation seeds are excluded. A two-target weak-direction regression qualified its local candidate but was WellConditioned under the conditioning classifier, so it is not counted as a weak-identification result. |
| High-overlap / near-singular | A deliberately near-duplicate support fixture, unlike the earlier cube fixture, measured free A/C condition `5.4683e9`. B1 Null and Split qualification passed; the prediction Drift/Seed values were each below `1e-8`. The legacy classification was unchanged and both atoms retained `Warning` point-quality from IllConditioned conditioning. Small prediction Seed therefore did not remove the warning. |
| Rank-deficient | Two identical support columns classified `Singular`; the B1 Null reference was `Unavailable` with `null-reference-numerical-validity-failed`, no Seed metric, and both point qualities remained `Warning`. Separate low/high width profile attempts failed inner KKT and were not treated as qualified multi-start endpoints. |
| Active width / near boundary | B1 Null reference qualified. It recorded one replay-qualified inward endpoint; the outward endpoint projected to zero displacement and was excluded. Legacy status remained unchanged. |
| Multi-start | On the exact one-target fixture, starts at widths `.4` and `1.1` both passed Strict Independent Replay and converged to width endpoints within `1e-8` in log width and prediction difference within `1e-12`. This is one smooth fixture, not proof of a unique physical solution. Prior raw evidence had 50 replay-qualified starts, maximum prediction difference `3.92e-12`, scaled A/C difference `5.70e-12`, and log-width difference `4.47e-8`. |
| Real cryo-EM map | Not run. No `.mrc`, `.map`, `.ccp4`, or equivalent real-map fixture exists in this checkout. |

The pre-cleanup S4 focused Release run passed 31 Split-metric and Local Parameter Stability tests. An exploratory multi-start attempt on the near-singular fixture failed inner KKT from both starts; those endpoints were excluded. It does not overturn the qualified exact-fixture result or the prior replay-qualified raw starts.

## Preliminary promotion decision

**Do not promote C1. Keep C0 Production.** The synthetic numerical guards support the distinction between common baseline reoptimization drift and seed-to-seed variation, and failed references correctly suppress Seed metrics. Promotion gates remain open because the added Null Search exceeds the small-workload wall target, the full difficult-case performance matrix and real-map evidence are absent. As C1 was not promoted, no persistence schema or consumer path was added. No warning threshold, point-quality, Local/Global KKT, solver, acceptance, endpoint certification, or Independent Replay policy changed.

## S5 — Promotion decision and cleanup

**Result: Partial Qualification; C1 not promoted; C0 retained.** The B1 candidate and B2 workspace experiment were actually executed, compared, and measured in S1–S4. The Release path now contains no Split Diagnostics candidate branch, runtime mode, result field, or persistence field. B2 was rejected for no measurable CPU/wall benefit and its workspace reuse code was removed. The B0 benchmark remains; its B0-only runner no longer has B1/B2 switches.

| Gate | Result | Reason |
|---|---|---|
| Numerical validity | Partial | Qualified synthetic references and endpoint metrics passed KKT/Replay and decomposition regressions. Noise-only N2 budget exhaustion and rank-deficient numerical validity correctly returned Unavailable without Seed metrics. |
| Diagnostic reliability | Partial | Drift/Seed separation held in tested exact, noisy, active-bound, and condition-`5.47e9` cases; conditioning/point-quality warnings remained. Real-map and broader weak-width/multi-start basin coverage are absent. |
| Engineering | Failed for promotion | B1 wall increased 15.23% on the 3-target chain and 30.5% on 3-target noise-only; B2 did not improve B1. RSS stayed within the exploratory gate. |
| Persistence / Target Evidence | Not promoted | No JSON/schema changes. Existing persistence contract and Target Evidence tests pass; candidate diagnostics are not exposed to consumers. |

Cleanup removed the C1 orchestration and B2 reuse branches from `ParameterStability.cpp`, the Split fields/mode and replay counters from the search result/work structures, the B0/B1/B2 benchmark switch, and the temporary noise-only benchmark topology. The pure metric helper moved from `src/` into `tests/support/` and is compiled only into `rhbm_gem_test`; all 11 split formula/qualification regression cases remain. No experimental Production branch, JSON compatibility case, debug output, temporary precision control, extra CMake target, or unused source remains.

Final post-cleanup verification:

| Check | Result |
|---|---|
| Release full build, `--parallel 4` | Passed; clean incremental rebuild reported no work. |
| Full CTest, `--parallel 4` | 19/19 passed. The first attempt set `OMP_NUM_THREADS=1`; an unrelated Eigen thread restoration test then failed because it explicitly sets Eigen to 2. The isolated test passed with that variable unset, and the full suite passed on rerun. |
| Focused Split/Parameter Stability/Conditioning/Target Evidence | 32/32 passed: 11 split metric, 10 parameter stability, 7 conditioning, 4 target evidence. |
| Persistence / JSON contract | `joint_persistence_contract_test` passed; no schema changes. |
| Post-cleanup B0 Release benchmark | 1 warmup + 5 measured, single numeric thread; accepted-state fingerprints matched the B0 fingerprints from S3. |

| Targets | B0 wall median | Stability CPU | Process peak RSS | Accepted-state fingerprint | Endpoint / Independent Replay |
|---:|---:|---:|---:|---|---|
| 3 | 14.165 ms | 10.436 ms | 24.216 MB | `a7980d92f9da5a52` | Qualified / Passed |
| 30 | 340.585 ms | 235.427 ms | 35.619 MB | `b18509add52fabf9` | Qualified / Passed |
| 300 | 1.14346 s | 206.250 ms | 88.179 MB | `3535e5d2710ef2ba` | Qualified / Passed |
| 1,000 | 3.42900 s | 289.526 ms | 236.405 MB | `7483b6f0d4edefa8` | Qualified / Passed |

The B0/B1/B2 measurements in S3 and the noise-only matrix above were completed before cleanup. B1/B2 were then removed, so the final post-cleanup run is B0-only; no post-cleanup B1/B2 runner remains.
