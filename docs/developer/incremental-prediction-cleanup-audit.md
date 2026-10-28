# Incremental Prediction Cleanup Audit

This document records the production baseline and source audit captured before cleanup on `develop`.

## C0 — Baseline

### Git and build

| Item | Baseline |
|---|---|
| Initial HEAD | `661ce0c6ad0f51e824595dd528189e8e85b97a67` |
| Branch / worktree | `develop`, clean |
| Remote tracking | `develop...origin/develop` (equal at baseline) |
| Compiler | AppleClang 21.0.0, arm64 macOS |
| Generator / build type | Ninja / Release |
| CMake options | `BUILD_TESTING=ON`, `BUILD_PYTHON_BINDINGS=OFF`, experimental feature ON, UMAP ON, dependency provider SYSTEM, OpenMP AUTO (OpenMP 5.1 found), ROOT AUTO |
| CTest groups | 19 |
| FixedNeighbor block-coordinate GTest cases | 37 |

The Release `rhbm_gem` shared library compiled without `RHBM_GEM_TEST_INSTRUMENTATION`. Its source list contains no `tests/support` files, and its exported symbols did not include the test-only replay and instrumentation entry points checked at baseline. CMake separately compiles `rhbm_gem_test` from the production source list with test instrumentation for `rhbm_tests`.

The baseline `rhbm_tests_joint_component_numerical` group passed (6.39 s). The Release build completed successfully. `rhbm_fixed_neighbor_closed_loop_benchmark` links `rhbm_gem`, not `rhbm_gem_test`.

### Production E1 benchmark

The checked-in closed-loop benchmark records the starting width, final parameter-vector SHA-256 digests, acceptance counts, search and wall medians, memory preflight, peak RSS, objective and KKT results. Digest input is the vector length followed by the little-endian IEEE-754 bit patterns. Each row used one warm-up and five measured runs in Release; sizes ran sequentially. Initial eta is `log(0.55)` for each target and initial beta is zero.

| Targets | Search median (s) | Wall median (s) | Peak RSS (bytes) | Objective / independent objective | A/C KKT | Projected width KKT | Sweeps; accepted / unchanged / failed | Preflight bytes | Eta SHA-256 | Beta SHA-256 |
|---:|---:|---:|---:|---:|---:|---:|---|---:|---|---|
| 3 | 0.004314583 | 0.005237208 | 23,986,176 | 2.01990781e-10 / 2.01990781e-10 | 3.33066907e-16 | 3.88578059e-15 | 4; 4 / 0 / 0 | 447,445 | `1c07aab49795cbe8e2a75e06e94cd5b779f628330bb0b65a1db8e2071230820c` | `cc02ff34d10aa3e1955399b36d300a492cb8e8872adf60114c8d501be2f77af9` |
| 30 | 0.094831833 | 0.113906250 | 35,454,976 | 1.51707603e-10 / 1.51707603e-10 | 1.60238889e-09 | 1.88066229e-11 | 5; 15 / 0 / 0 | 3,930,461 | `a55389519cefc3d266572140c6463f67efdc938c3ffcedd8f4bb80dd2d1654de` | `1aa3e251e9031e137455c133328b9972ef9ead03e4b1bccd871b0c6a182de9e2` |
| 300 | 1.016224370 | 1.222492710 | 90,161,152 | 1.47904336e-10 / 1.47904336e-10 | 5.06726425e-10 | 1.87827531e-12 | 5; 130 / 0 / 0 | 36,716,466 | `8dc053c3930e234569a40e5ea82316c6759590fe4a2f0ab9d54eb434bd03ebf3` | `54f379bc1c47684ebf58231974a0eba838470a06297adf028c467657ebafbffb` |
| 1,000 | 3.394646710 | 4.083542460 | 236,879,872 | 1.47904359e-10 / 1.47904359e-10 | 2.77527831e-10 | 5.63438185e-13 | 5; 420 / 0 / 0 | 121,648,896 | `7d5b3e38ee28a92a4abccb002a069f8725df11d1f912d2969a46888d89975eb3` | `bf4890cba7bb7d83171c24c12fc5b11159d525db187835ccf0f655e2177b5bc5` |

All four runs returned `block-stationary`, passed endpoint qualification and independent replay, passed memory preflight, used zero full-candidate fallback replays, and reported zero maximum cache-replay error.

### Small and difficult fixture baselines

For these fixtures, the endpoint values below are the E1 branch of the pre-cleanup instrumented diagnostic. Objectives and independent objectives are normalized by the fixed parent scale. `A/U/F` means accepted, unchanged, failed blocks. Every listed case passed memory preflight. A missing independent objective means endpoint assessment was unavailable after the local solve failed.

| Fixture | Initial widths | Objective / independent objective | A/C KKT | Projected width KKT | Sweeps; A/U/F | Search / endpoint / replay |
|---|---|---:|---:|---:|---|---|
| Active upper bound | 0.50 x 3; upper bound 0.50 | 2.0254002767800283e-10 / 2.0254002767799073e-10 | 1.1910380459667635e-09 | 4.4614645311469303e-11 | 5; 10 / 0 / 0 | stationary / verified / passed |
| Initial width 0.35 | 0.35 x 3 | 2.0199080158137056e-10 / 2.0199080158135505e-10 | 4.8809691843842984e-09 | 1.8283163871757324e-10 | 7; 14 / 0 / 0 | stationary / verified / passed |
| Initial width 0.55 | 0.55 x 3 | 2.0199079729737389e-10 / 2.0199079729739116e-10 | 4.3341136235497402e-09 | 1.6235413013987454e-10 | 5; 10 / 0 / 0 | stationary / verified / passed |
| Initial width 0.90 | 0.90 x 3 | 2.0199079727070732e-10 / 2.0199079727069230e-10 | 4.3304900776419686e-09 | 1.6221834986396289e-10 | 23; 46 / 0 / 0 | stationary / verified / passed |
| High-overlap cube-8 | 0.55 x 8 | 1.1487848300043995e-10 / 1.1487848300043186e-10 | 6.9202588104388951e-10 | 1.4783063662093809e-11 | 5; 10 / 0 / 0 | stationary / verified / passed |
| Signed cancellation | 0.55 x 3 | 4.0589644398391108e-10 / 4.0589644398388808e-10 | 2.2204460492503131e-16 | 1.2280135797126945e-09 | 500; 42 / 958 / 0 | sweep budget / not verified / unavailable |
| Rank-deficient constant support | 0.55 x 3 | 0.49999999999999944 / unavailable | 0.48281226720638365 | 0 | 1; 0 / 0 / 1 | recovery exhausted / not verified / unavailable |
| Near-singular overlapping support | 0.55 x 2 | 0.24948812932689915 / unavailable | 0.68976802196285614 | 0 | 1; 0 / 0 / 1 | recovery exhausted / not verified / unavailable |
| Low-signal noise-only | 0.55 x 3 | 3.5558894899614887e-08 / 3.5558894899614893e-08 | 4.4687101427961151e-13 | 1.5473043853131419e-12 | 100; 184 / 16 / 0 | sweep budget / verified / passed |
| Weak-width direction | 0.72, 1.08 | 1.1716558769804921e-32 / 1.0952435371774165e-32 | 1.1102230246251565e-16 | 5.5511151231257827e-17 | 2; 2 / 0 / 0 | stationary / verified / passed |

The exact small-fixture final eta and beta vectors from this baseline are captured here as round-trip decimal values:

| Fixture | Final eta | Final beta (A, C pairs by atom) |
|---|---|---|
| Active upper bound | `[-0.69314718055994529, -0.69315372645241125, -0.69314737955761552]` | `[2.0000092977896125, 0.19999844661139962, 1.9999351609261744, 0.20001477760579911, 1.999990285069869, 0.20000529847809637]` |
| Initial width 0.35 | `[-0.6931438357031382, -0.69315378776784742, -0.69314737115419967]` | `[2.0000242651018496, 0.1999979989737668, 1.9999347020637406, 0.20001481712774341, 1.9999903462814819, 0.20000529340054621]` |
| Initial width 0.55 | `[-0.69314383363060361, -0.69315380395213411, -0.69314736908218821]` | `[2.000024280201524, 0.19999799772082977, 1.9999345841187324, 0.20001482691776665, 1.9999903613786543, 0.20000529214768559]` |
| Initial width 0.90 | `[-0.69314383363141852, -0.69315380394577053, -0.6931473690830029]` | `[2.0000242801955856, 0.19999799772132235, 1.9999345841651093, 0.20001482691391712, 1.999990361372719, 0.20000529214817797]` |
| Signed cancellation | `[-0.69313238734324734, -0.69315371810798543, -0.69314737995827003]` | `[1.8623863311018945e-09, -0.20000155320669047, 1.9999352041743306, 0.20001477549941502, 1.9999905326991614, -0.19999470127506408]` |
| Rank deficient | `[-0.59783700075562041, -0.59783700075562041, -0.59783700075562041]` | `[0, 0, 0, 0, 0, 0]` |
| Near singular | `[-0.59783700075562041, -0.59783700075562041]` | `[0, 0, 0, 0]` |
| Low signal | `[0.78845736036427017, 0.78845736036427017, -1.5879387046886489e-05]` | `[0, 4.2032496951816247e-07, 0, 5.3037448061202324e-07, 0, 1.7721013208787142e-07]` |
| Weak width | `[-0.3285040669720361, 0.076961041136128394]` | `[1.2500000000000064, 0.1899999999999967, 0.82999999999999574, 0.11000000000000372]` |

For high-overlap cube-8, final eta was `[-0.69313871605003374, -0.69315169891460315, -0.69315106161070639, -0.69314739119322577, -0.69314918762653932, -0.69314954990327793, -0.693150352150675, -0.69313826723537253]`; final beta was `[2.0000612641323605, 0.19999496426198482, 1.9999497634050896, 0.200013585769901, 1.9999551482844447, 0.20001272659046143, 1.9999879618205636, 0.2000066598204833, 1.9999712153996854, 0.20001005410874323, 1.999968704757819, 0.20001020790365032, 1.9999614800229633, 0.20001156926991565, 2.0000695623425111, 0.19999152258421249]`.

Input hashes for standard 3/30/300/1,000 target chain fixtures and high-overlap, signed-cancellation, rank-deficient, near-singular, low-signal, and weak-width fixtures remain asserted by the existing numerical tests.

### Test-instrumentation responsibility audit

Counts below are occurrences of `RHBM_GEM_TEST_INSTRUMENTATION` before cleanup.

| Source | Occurrences | Responsibility / action |
|---|---:|---|
| `FixedNeighborBlockCoordinate.cpp` | 35 | D1 objective selector, replay callback, fault injection, observers, memory override, plus telemetry and diagnostic capture. Remove obsolete strategy/callback; retain only hooks needed by transaction and memory regression. |
| `FixedNeighborBlockCoordinate.hpp` | 10 | Test event and policy fields. Move full replay reference out to test support; retain minimal test-only observation/fault controls if integration coverage requires them. |
| `FixedNeighborPolicy.hpp` | 2 | Test configuration callback on the production component policy. Remove or isolate. |
| `ActiveSetBoundedLM.cpp` | 25 | Numerical trace and counters; audit each consumer before removing anything. |
| `TiledDerivative.cpp` | 26 | Numerical correctness/cancellation plus performance counters; preserve cancellation and reference-order behavior. |
| `TiledQR.hpp` | 3 | Scratch materialization/resize counters; retain workspace reuse, remove counters only if no distinct test needs them. |
| `Search.cpp` | 1 | Production `InstrumentedLM<Profile>` integration; retain. |

Classification: **A** `JointBlockCoordinate::ReplaceBlock()` is the only identified Delta-add candidate-state mutator and is used only by one structural/conditional test helper; remove the mutator and its state-update assertions while retaining conditional-response and structural tests. **B** D1 full-objective strategy selection and candidate replay replacement are test-only route selectors with no production role; remove. **C** candidate-fault, sweep-audit perturbation, memory-budget override, retry override, and observers protect accepted-state transaction, replay-audit, memory fallback, and recovery tests; retain only those with direct regression coverage. **D** LM/derivative/QR counters are not prediction experiments; audit them by their test/benchmark call sites. **E** exact affected-row recompute, incremental objective, sweep replay audit, post-audit reconciliation, memory fallback, endpoint certification, and independent replay are required production logic and remain unchanged.

## Cleanup stages

### C1 — Remove the retired Delta-add helper

Removed ReplaceBlock() and renamed its mutable support state to ConditionalReplayState. The conditional-response test now applies candidate parameters to a copied full parameter vector and obtains the next state through the test support's complete replay. It still checks the conditional effective-response identity, local/full objective-delta equality, unchanged parameters outside the core, unchanged predictions outside affected rows, structural core ordering, bridge boundary coupling, affected rows, neighbor atoms, and structural census.

JointBlockCoordinate.cpp/.hpp remain because they still provide the independent full replay and conditional-response calculations used by those assertions. They no longer implement an incremental state update. Their CMake source entry remains necessary.

The C1 source search found no ReplaceBlock, DeltaPrediction, IncrementalPrediction, prediction += change, or residual += change matches in src or tests. rhbm_tests_joint_component_contract passed (1.55 s); the CTest group count remains 19 and the FixedNeighbor block-coordinate case count remains 37. The production numerical source was not changed in C1, so C0 Release E1 results remain the comparison baseline.

Commit: `4e87e273`.

### C2 — Isolate FixedNeighbor test controls from production candidate replay

Removed the test-only D1 full-objective scan, `full_objective_reference_for_testing`, the replaceable candidate replay callback, and the production-side `FullReplayForTesting()` helper. `ReplayAffectedRows()` now has only the E1 incremental-objective calculation and no test-only layout parameter. The component-level `configure_search_for_testing` callback was removed; its former callers now rely on the unchanged fixed production defaults or directly configure `SearchFixedNeighborComponent()` in search-specific numerical tests.

The independent candidate reference now lives in `tests/support/JointBlockCoordinate::Replay()`. When the affected-row path runs, the candidate event references its already-created prediction/residual scratch vectors; it no longer assembles or copies full candidate prediction/residual vectors for observers. Incremental-candidate tests compare affected-row values, unaffected-row preservation, objective, and acceptance decision against the support replay. Full-replay memory fallback remains covered through its integration test. The retired candidate trace and its SHA output were removed from `JointOptimizationBenchmark.cpp`.

The candidate fault injector, sweep-audit perturbation, retry profile override, state/sweep observers, and memory-budget override remain behind `RHBM_GEM_TEST_INSTRUMENTATION` because the numerical suite uses them to verify accepted-state integrity, audit-before-refresh failure, recovery retry, and memory preflight/fallback end to end. Their component policy entry point has no test configuration callback. The Release compile command for the production `FixedNeighborBlockCoordinate.cpp` has no instrumentation definition; the instrumented copy is built only into `rhbm_gem_test`, and the production library exports none of the removed or test-only event symbols. This is source-level isolation for the removed strategy/callback; remaining test hooks are compile-isolated.

Instrumentation occurrence counts changed from 35 to 25 in `FixedNeighborBlockCoordinate.cpp` and from 2 to 0 in `FixedNeighborPolicy.hpp`; `FixedNeighborBlockCoordinate.hpp` remains at 10. Active-set LM, derivative, QR, and Search counts are unchanged in C2. The CTest inventory remains 19 groups and the FixedNeighbor block-coordinate suite remains 37 cases. Release targets `rhbm_tests`, `rhbm_fixed_neighbor_closed_loop_benchmark`, and `rhbm_joint_optimization_benchmark` built successfully. The joint component contract, workflow, and numerical CTest groups passed (3/3; numerical group includes 89 GoogleTest cases across its configured suites).

Commit: `d7182237`.

### C3 — Remove retired numerical-kernel profiling counters

`BoundedLmWorkForTesting` now retains only the damping trace and its capture switch. Removed workspace resize/allocation and peak-byte estimates, dense QR counts/timing, KKT and tolerance evaluation counts, and damping-search timing. Numerical tests still inspect the recorded active-set sequence, damping phase and decision, rank/status, step, and KKT values; bracket and bisection expectations are derived from the trace. The solver calculations and workspace reuse were not changed.

Removed `DerivativeWorkForTesting` and its preparation, tile, reduction, triplet, scratch, cancellation, and timing counters. Cancellation coverage now checks the actual `reference_order` result and compares full and projected derivatives with the independent dense reference. The tile-boundary tests continue to compare derivative and reduction outputs. A new `TiledQRReusesScratchForRepeatedAssembledShape` regression verifies that the existing scratch buffers retain their storage when the assembled shape repeats; the buffers remain in place.

`JointOptimizationBenchmark` no longer reports the removed internal LM, derivative, or QR counters and timings. It retains end-to-end microbenchmark timings and trace/output digests. `Search.cpp` and production `InstrumentedLM<Profile>` are unchanged.

| Source | C0 count | After C3 | Remaining test-only responsibility |
|---|---:|---:|---|
| `ActiveSetBoundedLM.cpp` | 25 | 13 | Record numerical damping/active-set trace when enabled for tests. |
| `ActiveSetBoundedLM.hpp` | — | 1 | Define trace capture storage for tests. |
| `TiledDerivative.cpp` | 26 | 0 | No test-only instrumentation remains; numerical cancellation behavior is checked through output. |
| `TiledDerivative.hpp` | — | 0 | `DerivativeWork` interface removed. |
| `TiledQR.hpp` | 3 | 0 | Scratch buffers remain; only their counters were removed. |
| `Search.cpp` | 1 | 1 | Required production `InstrumentedLM<Profile>` integration. |

The Release targets `rhbm_tests`, `rhbm_fixed_neighbor_closed_loop_benchmark`, and `rhbm_joint_optimization_benchmark` built successfully. The joint component contract, workflow, and numerical CTest groups passed (3/3; the numerical group includes the new scratch-reuse regression). The CTest inventory remains 19 groups and the FixedNeighbor block-coordinate suite remains 37 cases.

Commit: `efc25fbc`.

### C4 — Trim benchmark, support, and CMake residue

Reduced `JointOptimizationBenchmark.cpp` to the bounded-LM solve and tiled-derivative preparation/reduction microbenchmarks. It links the uninstrumented production `rhbm_gem` library. Removed the duplicate synthetic FixedNeighbor diagnostic run, its candidate replay timing split and Amdahl estimate, peak-RSS reporting, and the LM/derivative trace digests. The production `FixedNeighborClosedLoopBenchmark.cpp` remains the end-to-end E1 performance benchmark.

The closed-loop benchmark directly calls `SearchFixedNeighborComponent`, requires a completed sweep and endpoint certification, and checks the independent replay. Its output includes affected-candidate and full-replay-fallback counts; it does not calculate a second full objective for every candidate.

Support call-site review found continuing coverage for each retained helper: `JointBlockCoordinate` supplies a full contributor replay and conditional-response reference to block-coordinate and candidate tests; `JointFixedNeighborDiagnostic` exposes sweep/state evidence to the FixedNeighbor integration tests; `JointComponentChecks`, `JointSyntheticWorkload`, `JointDenseReference`, and `JointTestNumerics` are used by the component, workflow, or numerical suites. The independent dense/reference implementations remain separate from production E1. No remaining support source or function in this audit was unused, and the CTest suite list and regression cases were not reduced.

In CMake, the test instrumentation definition now propagates from `rhbm_gem_test` to `rhbm_tests` through the test library's public interface; it is no longer specified twice. Only `rhbm_tests` uses that instrumented clone. Both benchmarks link the production `rhbm_gem`. The existing configure-time guards still reject test-support sources or `RHBM_GEM_TEST_INSTRUMENTATION` on `rhbm_gem`, and there are no stale source references.

The CTest inventory remains 19 groups and the FixedNeighbor block-coordinate suite remains 37 cases. Contract, workflow, and numerical CTest groups pass after the C4 edits (1.21 s, 2.08 s, and 18.27 s). The Release microbenchmark built and ran against `rhbm_gem`; this run reported medians of 15.6153 ms for bounded LM, 19.3243 ms for derivative preparation, and 5.85488 ms for derivative reduction (13,710 rows by 60 derivative columns). Compile commands confirm the production LM source and both benchmark executables have no test-instrumentation definition, while the instrumented test clone and test sources do.

Commit: `a87efb70`.

### C5 — Final production regression and cleanup audit

The full Release build completed with AppleClang 21.0.0, arm64 macOS, and Ninja. The production library and both benchmarks compiled without `RHBM_GEM_TEST_INSTRUMENTATION`; compile commands show the definition only on `rhbm_gem_test` and its test sources. The production library contains no test-only `ForTesting`, replay, or counter symbols checked in the source audit. Full CTest passed 19/19 in 39.49 s; the test inventory and the 37 FixedNeighbor block-coordinate cases are unchanged.

The Release E1 benchmark runs below each used one warm-up and five measured runs. Every size passed memory preflight, completed endpoint certification, passed independent replay, reported zero full-candidate fallbacks, and produced the C0 eta/beta SHA-256 digests exactly.

| Targets | Candidate evaluation (s) | Sweep replay (s) | Replay audit (s) | Objective reconciliation (s) | Global state (s) | Endpoint certification excluding independent replay (s) | Independent replay (s) | Peak RSS / preflight bytes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0.000462416 | 0.000211417 | 1.917e-6 | 4.2e-8 | 0.000410459 | 0.000104708 | 0.000108833 | 24,133,632 / 447,445 |
| 30 | 0.007064835 | 0.002570751 | 2.6582e-5 | 1.25e-7 | 0.005059292 | 0.001022750 | 0.001071750 | 35,373,056 / 3,930,461 |
| 300 | 0.072989171 | 0.025759001 | 0.000256874 | 4.1e-8 | 0.050417125 | 0.010219125 | 0.010782167 | 88,276,992 / 36,716,466 |
| 1,000 | 0.246198629 | 0.086668084 | 0.000917459 | 8.4e-8 | 0.170182333 | 0.034814501 | 0.037242250 | 240,435,200 / 121,648,896 |

| Targets | C0 search median (s) | C5 search median (s) | Change | C0 wall median (s) | C5 wall median (s) | Change | C0 → C5 peak RSS (bytes) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0.004314583 | 0.004357875 | +1.003% | 0.005237208 | 0.005265958 | +0.549% | 23,986,176 → 24,133,632 |
| 30 | 0.094831833 | 0.094801667 | −0.032% | 0.113906250 | 0.113753375 | −0.134% | 35,454,976 → 35,373,056 |
| 300 | 1.016224370 | 1.014450750 | −0.175% | 1.222492710 | 1.217986000 | −0.369% | 90,161,152 → 88,276,992 |
| 1,000 | 3.394646710 | 3.427572750 | +0.970% | 4.083542460 | 4.121985620 | +0.942% | 236,879,872 → 240,435,200 |

The largest median change is under 1% for search and wall time. The 1,000-target search median is 3.4276 s versus a 3.3996–3.5564 s range across the five current runs; the C0 median is 3.3946 s. This is not a clear or reproducible regression. State digests, objective, KKT, sweep and acceptance results match the baseline.

Hard-fixture post-cleanup values were collected from the same production FixedNeighbor search through the test-only diagnostic adapter; the adapter did not install observers or alter candidate decisions. Final eta/beta vectors printed at 17-digit round-trip precision match the C0 vectors exactly. Values below also match the C0 objective, KKT, and acceptance rows above.

| Fixture | Final / independent objective | Global A/C KKT | Projected width KKT | Sweeps; accepted / unchanged / failed | C5 stop and qualification |
|---|---:|---:|---:|---|---|
| Initial width 0.35 | 2.0199080158137056e-10 / 2.0199080158135505e-10 | 4.8809691843842984e-09 | 1.8283163871757324e-10 | 7; 14 / 0 / 0 | stationary; endpoint and replay passed |
| Initial width 0.55 | 2.0199079729737389e-10 / 2.0199079729739116e-10 | 4.3341136235497402e-09 | 1.6235413013987454e-10 | 5; 10 / 0 / 0 | stationary; endpoint and replay passed |
| Initial width 0.90 | 2.0199079727070732e-10 / 2.0199079727069230e-10 | 4.3304900776419686e-09 | 1.6221834986396289e-10 | 23; 46 / 0 / 0 | stationary; endpoint and replay passed |
| Active width upper bound | 2.0254002767800283e-10 / 2.0254002767799073e-10 | 1.1910380459667635e-09 | 4.4614645311469303e-11 | 5; 10 / 0 / 0 | stationary; endpoint and replay passed |
| High-overlap cube-8 | 1.1487848300043995e-10 / 1.1487848300043186e-10 | 6.9202588104388951e-10 | 1.4783063662093809e-11 | 5; 10 / 0 / 0 | stationary; endpoint and replay passed |
| Signed cancellation | 4.0589644398391108e-10 / 4.0589644398388808e-10 | 2.2204460492503131e-16 | 1.2280135797126945e-09 | 500; 42 / 958 / 0 | sweep budget; endpoint not certified; independent replay did not qualify |
| Rank-deficient constant support | 0.49999999999999944 / unavailable | 0.48281226720638365 | 0 | 1; 0 / 0 / 1 | recovery exhausted after 3 local attempts; endpoint/replay unavailable |
| Near-singular overlapping support | 0.24948812932689915 / unavailable | 0.68976802196285614 | 0 | 1; 0 / 0 / 1 | recovery exhausted after 2 local attempts; endpoint/replay unavailable |
| Low-signal noise-only | 3.5558894899614887e-08 / 3.5558894899614893e-08 | 4.4687101427961151e-13 | 1.5473043853131419e-12 | 100; 184 / 16 / 0 | sweep budget; endpoint and independent replay passed |
| Weak-width direction | 1.1716558769804921e-32 / 1.0952435371774165e-32 | 1.1102230246251565e-16 | 5.5511151231257827e-17 | 2; 2 / 0 / 0 | stationary; endpoint and replay passed |

The rank-deficient and near-singular cases retain their original recovery failures. Signed cancellation and low-signal retain their budget-stop limitations. Active-width and weak-width qualification still pass. No KKT, width, convergence, objective, or certification thresholds changed.

The transaction and qualification regressions remain end-to-end tests: exact affected-row values and incremental objective versus full replay; accepted-state integrity for objective increase and nonfinite candidates; replay-audit failure before cache reconciliation; unchanged-candidate replay skip; recovery retry; memory preflight and full-replay fallback; sweep objective reconciliation; endpoint certification; and independent replay. The full numerical, contract, and workflow groups all passed as part of the 19-group CTest run.

Final source search for `DeltaPrediction`, `IncrementalPrediction`, `ReplaceBlock`, `full_objective_reference_for_testing`, and `candidate_replay_for_testing` returned no hits in `src/core/detail/joint_component` or `tests`. `candidate_fault_for_testing` and `RHBM_GEM_TEST_INSTRUMENTATION` remain only as required, compile-isolated test support; their responsibilities are summarized in the final audit table below. There are zero retired experimental references and zero unexpected residuals.

| Source | C0 occurrences | After C4 | Production binary | Remaining responsibility |
|---|---:|---:|---|---|
| `FixedNeighborBlockCoordinate.cpp` | 35 | 25 | No test branches | Test fault injection, candidate/sweep/state observation, memory override, retry and audit perturbation. |
| `FixedNeighborBlockCoordinate.hpp` | 10 | 10 | No test fields | Declarations for the test-only candidate event, diagnostic controls, and observers. |
| `FixedNeighborPolicy.hpp` | 2 | 0 | No test fields | Component policy test configurator removed. |
| `ActiveSetBoundedLM.cpp` | 25 | 13 | No trace branch | Capture active-set, damping, rank, step, and KKT trace for numerical tests. |
| `TiledDerivative.cpp` | 26 | 0 | No test instrumentation | Cancellation correctness checked through `reference_order` and dense numerical comparison. |
| `TiledQR.hpp` | 3 | 0 | No test instrumentation | QR scratch buffers remain; only their counters were removed. |
| `Search.cpp` | 1 | 1 | No test metadata branch | Test-only capture of bounded subproblem details; production `InstrumentedLM<Profile>` remains. |

The retained FixedNeighbor fault injectors are still in the production solver source under `RHBM_GEM_TEST_INSTRUMENTATION`; they are absent from the production build and binary. Candidate and sweep observers are used by integration tests to compare E1 candidate rows/objectives to the independent replay and to verify audit-before-reconciliation and accepted-state transactions. Memory-budget/retry overrides and fault injection cover fallback and recovery end to end. Active-set tracing supports numerical damping, rank, release, and KKT tests. These hooks remain necessary test-only source coupling; they are not described as source-decoupled.

### Final decision

1. Retired Delta-add prediction and `ReplaceBlock()` are removed; `JointBlockCoordinate` remains for independent full replay and conditional-response reference tests.
2. Production has one E1 candidate objective; the D1 selector and candidate replay replacement callback are removed.
3. Fault injection and observers remain behind the test-only compile definition for transaction, replay, fallback, and recovery regression coverage.
4. The obsolete LM/derivative/QR counters and timers are removed; numerical LM traces and actual QR scratch reuse remain covered.
5. `JointOptimizationBenchmark` is reduced to bounded-LM and derivative/QR timings on the production library. `FixedNeighborClosedLoopBenchmark` remains the production E1 performance benchmark.
6. Production `rhbm_gem` has no test instrumentation or test-support source. Endpoint certification, sweep replay audit, objective reconciliation, memory fallback, and independent replay remain in the production path.
7. Numerical results are bitwise identical on all standard and hard fixtures measured. Median benchmark changes were below 1%; at 1,000 targets the C0 search median is only 0.15% below the minimum of the five C5 runs, so the measurements do not show a clear regression.

Commit: recorded in the Git history for this stage.
