# FirstStage cleanup dependency audit

Audited at initial develop `503a6607674e068e4d3ef3db52aa770abcf3aea4` before removing production paths. Initial tree clean; fetched origin/develop identical. Existing other detached worktree left untouched. No new branch/worktree or Notion changes.

All APIs below are private `src/core/detail` interfaces unless an `include/` path is specified. Research readers include the C++ executable, Python aggregators, three research test files, and historical Markdown/JSON; no production persistence or serialization encodes a sampling/interpolation policy or call timer. The raw samples and OLS/MDPDE estimator records are official persisted data and stay unchanged.

| Symbol | Definition | Production callers | Research/test callers at initial HEAD | Action |
| --- | --- | --- | --- | --- |
| `JointSeedSamplingMethod` | `src/core/detail/JointSeedSampling.hpp` | MapSampler sampling dispatch | `tests/core/joint_component/SeedSampling_test.cpp`; `tests/core/joint_component/SeedSensitivity_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | remove enum and both dispatch branches; use the shared Fibonacci workflow |
| `JointSeedSamplingPolicy` | `src/core/detail/JointSeedSampling.hpp` | FirstStageInitialization / PotentialFittingWorkflow defaults | `tests/core/joint_component/SeedSampling_test.cpp`; `tests/core/joint_component/SeedSensitivity_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | remove policy and forwarding parameters |
| `JointSeedSamplingReport` | `src/core/detail/JointSeedSampling.hpp` | optional output on FirstStage/workflow call chain (unused by public entrypoints) | `tests/core/joint_component/SeedSampling_test.cpp`; `tests/core/joint_component/SeedSensitivity_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | remove run-owned report and pointer parameters |
| `JointSeedAtomSamplingDiagnostics` | `src/core/detail/JointSeedSampling.hpp` | optional sampling report | none | remove counts/shell diagnostics |
| `RunJointSeedSamplingWorkflow` | `src/core/MapSampler.cpp` | PrepareJointWidthInitialization | `tests/core/joint_component/SeedSampling_test.cpp`; `tests/core/joint_component/SeedSensitivity_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | replace with RunPotentialSamplingWorkflow, retaining clamped worker count |
| `JointSeedSamplingName` | `src/core/MapSampler.cpp` | FirstStage sampling log | `tests/research/FirstStageSamplingStudy.cpp` | remove dispatch/name API; fixed algorithm description |
| `MapInterpolationMethod` | `src/core/detail/MapInterpolation.hpp` | MapSampler / QScoreHelper | `tests/core/MapInterpolation_test.cpp`; `tests/core/joint_component/SeedSampling_test.cpp`; `tests/core/joint_component/SeedSensitivity_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | remove enum; direct Tricubic and Trilinear calls |
| `MapInterpolationName` | `src/core/MapSampler.cpp` | FirstStage sampling log | `tests/research/FirstStageSamplingStudy.cpp` | remove name dispatch |
| `InterpolateMapValue` | `src/core/detail/MapInterpolation.hpp` | QScoreHelper / BuildLocalPotentialSampleList | `tests/core/MapInterpolation_test.cpp`; `tests/core/joint_component/SeedSampling_test.cpp` | remove dispatcher; preserve kernels verbatim |
| `BuildQScoreSeedSamples` | `src/core/detail/QScoreSampling.cpp` | retired alternative MapSampler branch only | `tests/core/QScoreHelper_test.cpp` | remove FirstStage seed builder |
| `QScoreSeedSamples` | `src/core/detail/QScoreSampling.hpp` | retired seed builder | `tests/core/QScoreHelper_test.cpp` | remove wrapper |
| `QScoreShellDiagnostics` | `src/core/detail/QScoreSampling.hpp` | populated in QScore shell; not read by QScore correlation | none | remove counters, preserve retry loop and acceptance |
| `QScoreShellSamples` | `src/core/detail/QScoreSampling.hpp` | QScoreHelper reads points only | none | return SamplingPointList |
| `GenerateSpiralSpherePoints` | `src/core/detail/QScoreSampling.cpp` | SampleQScoreShell | `tests/core/QScoreHelper_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | make private; exercise geometry via public QScore radial API |
| `SampleQScoreShell` | `src/core/detail/QScoreSampling.cpp` | GetRadialPointsForQScore | `tests/core/QScoreHelper_test.cpp` | retain with plain point-list return |
| `IsInsideNeighborExclusion` | `include/rhbm_gem/utils/domain/SampleFilter.hpp + src/utils/domain/SampleFilter.cpp` | QScoreSampling only | `tests/support/FibonacciNeighborSelection.hpp`; `tests/utils/domain/SampleFilter_test.cpp` | move into QScore detail module with validation and accumulation order; move boundary tests |
| `FirstStageTimingReport` | `src/core/detail/FirstStageInitialization.hpp` | optional argument; no public workflow enables it | `tests/core/joint_component/WidthSource_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | remove struct, pointer, chrono and per-worker arrays |
| `CalibrateQualifiedWidthDonors` | `src/core/detail/joint_component/WidthCalibration.cpp` | CalibrateGlobalWidth | `tests/support/ResearchWidthInitialization.hpp` | make private; preserve one calibration implementation |
| `ResearchWidthInitialization` | `tests/support/ResearchWidthInitialization.hpp` | none | `tests/core/joint_component/WidthSource_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | remove OLS source, fallback, calibration/cohort and projection bridge |
| `FibonacciNeighborSelection` | `tests/support/FibonacciNeighborSelection.hpp` | none | `tests/core/joint_component/SeedSampling_test.cpp`; `tests/core/joint_component/SeedSensitivity_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp`; `tests/support/FibonacciNeighborSelection.hpp` | remove FO/FQ masks/controls |
| `FirstStageSamplingFixture` | `tests/support/FirstStageSamplingFixture.hpp` | none | `tests/core/joint_component/SeedSensitivity_test.cpp`; `tests/core/joint_component/WidthSource_test.cpp`; `tests/research/FirstStageSamplingStudy.cpp` | extract shared Gaussian map/model to FirstStageSynthetic.hpp; no duplicated generator |
| `FirstStageSamplingStudy` | `tests/research/FirstStageSamplingStudy.cpp` | none | `tests/CMakeLists.txt`; `tests/research/WidthInitializationStudy.hpp` | remove executable source |
| `WidthInitializationStudy` | `tests/research/WidthInitializationStudy.hpp` | none | `tests/research/FirstStageSamplingStudy.cpp` | remove research include implementation |
| `first_stage_sampling_study` | `tests/CMakeLists.txt` | none; EXCLUDE_FROM_ALL, not installed or CTest registered | `tests/CMakeLists.txt`; `tests/research/FirstStageSamplingStudy.cpp`; `tests/research/run_first_stage_sampling_study.py`; `tests/research/run_neighbor_selection_study.py`; `tests/research/run_width_initialization_study.py` | remove target/link/include configuration |

## Candidate controls

| Controls | Definitions/callers | Decision |
| --- | --- | --- |
| F1, Q0, Q1 | MapSampler private policy, SeedSampling/SeedSensitivity tests, FirstStageSamplingStudy and sampling runner | Remove FirstStage alternatives; retain independent QScore feature. |
| FO, FQ | FibonacciNeighborSelection.hpp, research executable/neighbor runner, SeedSampling/SeedSensitivity tests | Remove ownership-only and QScore-distance-only FirstStage masks. |
| Q8, S0, S1 | Research executable, sampling runner and reports | Remove repeated-center and matched-spiral controls. S0/S1 section labels in reports are history. |
| MM, OM, MO, OO | ResearchWidthInitialization.hpp, WidthInitializationStudy.hpp, WidthSource tests, width runner | Retain production MDPDE initialization/calibration; remove all source-switch/factorial execution. |
| Common bound | Research executable, support bridge and SeedSensitivity/WidthSource tests | Remove research execution, preserve production bounded projection. |

The exact-word search of FO/FQ/Q8/S0/S1/MM/OM/MO/OO found no unrelated production policies. JSON evidence and reports use these names historically. OLS names in Gaussian fitting, runtime tests and persistence are production, not retirement targets.

## Build, data and documentation boundaries

- CMake: one noninstalled research target; normal rhbm_tests groups use source-extracted suite filters with a registration guard. No public FitOptions or CLI switches exist for these studies.
- WidthSource_test directly read docs/research/first-stage-sampling/baseline.json. C0 extracts the necessary original numeric fields into tests/fixtures/joint_component/first_stage.json and removes that dependency.
- Existing baseline.json contains aggregate sample counts, alpha/parameters, calibration and full final JCE state/check identities. C0 supplements it with exact per-column IEEE-754 fingerprints over 48 production raw vectors (isolated and crowded) before production changes. Timing/RSS are excluded.
- Shared synthetic generator moves once to tests/support/FirstStageSynthetic.hpp. Production truth, endpoint, persistence and invalidation checks survive in FirstStageInitialization_test. WidthCalibration_test already covers native/refined donors, 20/95% limits, Type-7, fallback exclusion and immutable projection.
- QScoreHelper_test covers geometry, retry/empty shells, hydrogen/crowding, center weighting, interpolation boundary, correlation and invalid inputs. Only its FirstStage seed-builder tests retire; predicate boundary tests move from SampleFilter to QScore.
- JointAnalysis persistence schema, original OLS/MDPDE fields and failed-only recovery are unchanged. Existing runtime, numerical, data schema, downstream and command tests remain.
- Developer joint-stage-integration.md describes the retired current policies and needs correction. Study reports retain observations/limits with an archive notice and no current runner instructions. All 9 JSON files will be indexed by Git object/size/SHA-256 before pruning; their unique evidence remains recoverable from the initial commit.

## Verification sequence

1. C0 baseline/test extraction -> full Release CTest before removing production code.
2. C1 fixed Fibonacci/Tricubic -> sampling, FirstStage, calibration and QScore regressions.
3. C2 QScore API cleanup -> full QScore/predicate/kernel regressions.
4. C3 width/timer cleanup -> FirstStage/calibration/JCE exact regression.
5. C4 tools/test/CMake cleanup -> all registered GoogleTests, no missing/disabled/orphan sources.
6. C5 evidence/doc archive -> no test/runtime dependencies on study JSON.
7. C6 Release/Debug all targets/CTest, lint/install, static audit, fixed-workload timing/RSS, clean commits and ordinary push.

## Final token disambiguation

The final broad short-token scan also finds `MO` in ChemicalDataHelper.cpp's
molybdenum element lookup and `S1` in MapSimulation_test.cpp's synthetic mmCIF
entity fixture. Both are unrelated production/test identities and are retained.
All fully named retired policy/tool symbols are absent from src/include/tests/cmake.
`CalibrateQualifiedWidthDonors` and `GenerateSpiralSpherePoints` remain private
production implementations; the QScore distance predicate is in the private
QScore module, outside SampleFilter's public API.
