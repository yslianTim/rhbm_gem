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

Commit: recorded in the Git history for this stage.
