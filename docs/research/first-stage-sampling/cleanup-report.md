# FirstStage production cleanup delivery

## Initial state and constraints

- Initial HEAD: `503a6607674e068e4d3ef3db52aa770abcf3aea4`.
- Existing branch: `develop`; working tree clean, no user edits.
- `git fetch origin` succeeded; HEAD/origin/develop divergence 0/0.
- Existing secondary detached worktree untouched. No branch, worktree or detached checkout created; no reset or force push. Notion is not accessed.
- Existing Release/Debug configurations have BUILD_TESTING=ON. Build parallelism <=4, CTest parallel=1, OMP_NUM_THREADS=2, OMP_THREAD_LIMIT=4, OPENBLAS_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1; existing group overrides retained.
- No dependencies, numerical policies, qualification thresholds, convergence criteria, fitting radii or persistence schemas change.

See [symbol inventory](cleanup-inventory.md) for definitions, callers and decisions.

## Stage log

### C0 — production regression protection

Extracted the original committed coarse F0 atom/alpha/OLS/MDPDE, calibration,
JCE A/C/Width/gradient and certification identities into
`tests/fixtures/joint_component/first_stage.json` (27,453 bytes). Deduplicated
identical component checks; no alternate-policy records or runtime telemetry.
Captured exact per-column raw sample fingerprints for 24 isolated and 24
crowded atoms, plus peeled original-grid responses and Local OLS/MDPDE, before
production edits. New tests use precise JSON parsing and exact comparisons.

Extracted the shared synthetic generator once to `FirstStageSynthetic.hpp`.
`FirstStageInitialization_test.cpp` retains truth at two spacings, native
qualification, MDPDE seeds, distinct fallback/calibration pools, calibration,
projection, full endpoint state, KKT/Replay/Endpoint status, QScore, downstream,
codec/SQLite round-trip and raw-edit/grid-invalidation independence.
`WidthCalibration_test.cpp` already protects refined donors, minimum/coverage,
Type-7 and immutable nonzero projection. Retired duplicate production baseline
from WidthSource; its research-only cases remain until C3/C4.

Validation: initial Release all-target build and CTest 21/21 (42.01 s).
C0 Release build and CTest 21/21 (41.82 s), including registration coverage.
All original retained numerical fields match exactly. No ordinary GoogleTest
reads the research directory. Standalone temporary F0 timing harness uses
existing production calls; three pre-cleanup measurements saved for C6.

### C1 — one production sampling/interpolation path

Removed JointSeedSampling.hpp, all sampling/interpolation enums, reports,
name/dispatch functions and policy parameters in FirstStage/workflow. Joint
now calls the existing Fibonacci/Tricubic sampling workflow with the same
clamped worker count. Its parallel exception collection, deterministic commit
order, editor contract, progress and quiet behavior stay intact. Removed the
single-use generic sampling dispatcher. Both interpolation kernels are
unchanged; QScore calls Trilinear directly. Kernel tests retain all fields,
small grids and distinct boundary/nonfinite behavior.

Retired dependent factorial/neighbor SeedSampling cases and SeedSensitivity
file now that C0 protects their production contracts. Retained exact F0 parity,
serial/four-worker/repeated initialization, QScore isolation and native worker
failure atomicity for every prior raw vector. The old Q1-only boundary failure
is replaced by an unattached-atom sampler exception; finite out-of-grid cubic
sentinel behavior remains explicitly tested, rather than imposing Trilinear's
exception contract on Tricubic.

Validation: Release all-target build and CTest 21/21 (39.45 s), exact C0
sampling/FirstStage/calibration/JCE/QScore/peeled/Local fixtures unchanged.
The EXCLUDE_FROM_ALL retired research source is intentionally not built and
will be removed with its target in C4; no ordinary test uses the removed policy.

### C2 — QScore owns its production sampling

Removed BuildQScoreSeedSamples, QScoreSeedSamples, QScoreShellDiagnostics and
QScoreShellSamples. The shell helper now returns SamplingPointList, retaining
50 attempts, +2 candidates, first-success minimum/all accepted ordering and
empty failure behavior. Spiral geometry is private. The validated strict
squared-distance predicate moved verbatim to the private QScore detail module;
it is no longer a general SampleFilter public API. Its internal declaration
remains for QScore's numerical boundary/invalid-geometry tests, not as a
research-only interface. Ownership/cone filtering is untouched.

Removed two FirstStage seed-builder tests; rewrote diagnostic assertions as
public radial-output checks. Retained deterministic geometry, first-success
ordering, empty shells, hydrogen/crowding, center weighting, correlation,
Trilinear boundaries and invalid arguments. Three predicate tests moved intact
from SampleFilter to QScore. No new exception contracts or thresholds.

Validation: Release build and CTest 21/21 (40.34 s). Exact F0/QScore and all
C0 fixtures still pass; git diff --check passes.

### C3 — production MDPDE initialization only

Removed ResearchWidthInitialization.hpp and its six remaining OLS selection,
same-source fallback, calibration/cohort, projection-bridge and factorial tests
in WidthSource_test.cpp. Production protection was moved and passed in C0;
WidthCalibration_test remains intact. Removed FirstStageTimingReport, optional
pointer, chrono calls, timing arrays/fills and conditional accumulation only.
The same worker batches, alpha/estimator order and exception handling remain.
Moved CalibrateQualifiedWidthDonors into WidthCalibration.cpp's anonymous
namespace without modifying donor selection, thresholds, quantiles or bound
math. OLS/MDPDE results, estimator implementation and persistence are untouched.

Validation: Release build and CTest 21/21 (39.39 s); exact MDPDE seeds,
original/used Width, qualification, fallback, calibration, projected starts,
JCE state/objective and downstream regression remain passing.

### C4 — retired tooling and consolidated test registration

Removed FirstStageSamplingStudy.cpp, WidthInitializationStudy.hpp, all three
study Python runners and FibonacciNeighborSelection.hpp; ResearchWidthInitialization
was removed in C3 and the old FirstStageSamplingFixture was extracted once in C0.
All eight requested original paths are absent. Removed the research executable,
link and include CMake configuration. Renamed the retained sampling test file to
FirstStageSampling_test.cpp. Permanent synthetic support now keeps only the
separation and spacing needed by production regressions; retired phase/noise
controls were removed with exact fixture parity preserved.

Validation: Release all-target build and CTest 21/21 (38.76 s). Registration
covers all 1,092 enabled GoogleTests (initial registry verified at 1,105), zero
disabled, missing or unexpected duplicate coverage. All 83 remaining *_test.cpp
files are in compile_commands.json; no uncompiled test source. The generated
Ninja target inventory contains no study target. No retired symbols remain in
src/include/tests/cmake. Numerical baseline remains exact; no tolerances or
CTest timeout/qualification settings changed.

### C5 — historical conclusions and evidence archive

Retained all three scientific reports, their numerical tables, qualifications
and limitations. Added explicit archive notices, removed current executable
commands and dated superseded retention recommendations. Updated only the
FirstStage/QScore description in joint-stage-integration.md to the single
production path and existing calibration rules. All original table numeric
tokens were checked unchanged against the initial commit.

Pruned nine study JSON files (full factorial results, per-atom/endpoint audits
and validation ledgers). Before deletion, each was parsed and verified byte
identical to its immutable initial-commit Git object. evidence-archive.md records
all sizes, SHA-256 values, blob IDs and immutable recovery links. No unique
scientific evidence is destroyed: the ancestor commit is the archive; no
independent external archive was created. Ignored detailed build outputs remain
supplemental, not the sole evidence source. The production fixture is independent.

Validation: no CMake, source, ordinary tests, Python scripts or developer docs
refer to the retired research paths/policies; git diff --check and lint_repo
pass. No persistence, CLI/help schema or Notion document changes.

### C6 — final coverage review and build-mode reference audit

The final review retains the original crowded F0 Width RMSE <0.03 Å contract
alongside the unchanged isolated 0.002/0.02 Å bounds. The atomicity test now uses
a different replacement map so an early partial write cannot pass by coinciding
with prior responses. Removed the copied single-policy test's unnecessary model
copy/scope. No production code or existing qualification/tolerance changes.

The first full Debug run found only a one-ULP QScore/Release-reference mismatch:
0.93585969694861493 versus 0.93585969694861482. Rebuilding the initial commit's
QScoreHelper/QScoreSampling/SampleFilter sources and interpolation header with
the same Debug compiler flags reproduced 0.93585969694861493 exactly. This is
an existing optimization-mode difference, not cleanup drift. The fixture now
stores the original exact Debug and Release references; EXPECT_EQ remains in
both modes. All other FirstStage/calibration/JCE/raw/peeled/Local baseline
assertions passed in that first run. Updated Debug FirstStage suites pass 9/9;
complete Release/Debug acceptance is recorded below after the final runs.

Broad short-token matches MO (molybdenum lookup) and S1 (mmCIF fixture identity)
are unrelated to the retired factorial controls and remain unchanged. Fully
named retired symbols have zero src/include/tests/cmake matches.

## Final numerical regression

All exact claims below concern the same configuration/input. Elapsed time and
RSS are measurements. No existing tolerance, KKT threshold, convergence rule,
calibration threshold, fitting domain or persistence schema was changed.

| Numerical contract | Result and evidence |
| --- | --- |
| F0 positions / distances / ordering | Exact IEEE-754 fingerprints, 24 isolated + 24 crowded targets; direct serial/four-worker equality |
| F0 selection / counts | Exact fingerprints/counts, including rejected crowded samples; strict Ownership/Cone tests retained |
| F0 Tricubic responses | Exact fingerprints for every raw row, selected or rejected |
| FirstStage Alpha / OLS / MDPDE | Exact original 24-atom fixture; two real 120-target F0 runs match committed atom records |
| MDPDE seeds / original_tau / used_tau / seed source | Exact coarse baseline; fallback and unavailable/target/Halo cases protected without OLS source controls |
| Donor eligibility / quantiles / bound | Exact coarse and real references; native/refined/fallback, 20 minimum, 95% coverage and Type-7 tests retained |
| Bounded projection | Original immutable projection tests unchanged; real initial violation counts remain 1/0 |
| JCE objective / A/C / Width / gradient | Exact coarse fixture and both committed real MDPDE endpoints |
| KKT / Independent Replay / Endpoint | Same check identities/status/reasons; real global numerical check values/thresholds also exact; runtime and target qualification pass |
| QScore | Exact within Release and Debug; the original modes differ by one ULP, independently reproduced from the initial source; no tolerance added |
| Post-Joint | Exact coarse peeled-response fingerprints and Local OLS/MDPDE; real coverage unchanged |
| Persistence / invalidation | Codec equality, SQLite raw geometry round-trip and raw-edit/grid-state independence pass; existing schema/Group/invalidation groups pass |

Only F0 production was rerun on real data; no factorial/promotion study was
reopened. Both input CIF/map pairs match all four committed byte counts and
SHA-256 digests. Current atom fields and full final state were compared with
the committed width-study MDPDE production records as well as prior detailed
F0 snapshots; 6Z6U also matches the original committed full baseline.

| Workload | Width bound Å | Final objective | Qualified donors | Peeling available/total | Local / Group available |
| --- | ---: | ---: | ---: | ---: | ---: |
| coarse phase 0 (24 targets) | 1.8782695995625727 | 1.631044078500785e-7 | 24/24 | Exact fixture | 24 / fixture eligibility |
| 6Z6U / EMD-11103 (120 targets) | 3.1823057120925764 | 0.012784319782223002 | 120/120 | 33868/33868 | 120 / 45 |
| 8RQB / EMD-19436 (120 targets) | 2.80705057158481 | 0.016970359134848675 | 119/120 | 141457/141457 | 120 / 42 |

## Final verification

| Check | Result |
| --- | --- |
| Release all configured targets | Pass; library, CLI, test library/executable and resources |
| Debug all configured targets | Pass; same target set |
| Release CTest | 21/21, 38.66 s |
| Debug CTest | 21/21, 467.71 s, after exact build-mode reference audit |
| Updated Debug FirstStage suites | 9/9; complete final CTest also passes |
| GoogleTest registry | 1092 enabled, 1092 covered, 0 missing, 0 disabled |
| Duplicate registration | 20 intentional persistence overlaps, 0 unexpected |
| Test source inventory | All 83 remaining *_test.cpp files compiled |
| Production source inventory | All enabled production sources compiled; Python bindings remain OFF in the existing configurations |
| Research target/source inventory | No retired executable/runner/support path or generated study target |
| Static symbol scan | Zero fully named retired policy/tool symbols in src/include/tests/cmake |
| Numerical function audit | Interpolation, Ownership/Cone, distance predicate, calibration/projection and postprocessing bodies unchanged; Gaussian design, Fibonacci and persistence files unchanged |
| git diff --check | Pass |
| lint_repo | Pass |
| lint_install_smoke | Pass; installed library/CLI/resources consumed by CMake example |
| Build/test limits | Build jobs <=4, CTest parallel=1; OMP_THREAD_LIMIT=4, BLAS/VECLIB=1; existing group environments unchanged |

Only the three pre-existing missing-field-initializer warnings in
DataObjectPersistence_test.cpp remain (lines 422, 876, 915). That file is
unchanged. No timeout was increased. Optional Python bindings were not enabled
or tested; this retains the original build configuration.

## Fixed production workload measurements

A standalone temporary harness called the existing production sampling,
prepared FirstStage and initialized JCE functions. Three sequential one-worker
Release processes before and after cleanup use the same 24-target coarse map.
Map generation/problem construction are outside the reported phase times; RSS
is whole-process high-water through JCE. No production instrumentation remains.

| Metric | Before median (range) | After median (range) |
| --- | ---: | ---: |
| Sampling ms | 4.092 (3.820–4.095) | 4.016 (3.897–4.136) |
| Prepared FirstStage ms | 54.060 (52.475–54.091) | 52.243 (52.183–53.283) |
| Sampling + FirstStage ms | 58.152 (56.295–58.186) | 56.258 (56.080–57.419) |
| Initialized JCE ms (including calibration/projection) | 48.613 (46.851–48.918) | 47.115 (46.785–47.239) |
| Peak RSS decimal MB | 34.488 (34.423–34.603) | 34.341 (34.144–34.374) |

The objective remains exactly 1.631044078500785e-7 in every measurement.
These few short process measurements do not establish a speedup; timing/RSS
variation is not numerical parity evidence. No material memory increase is
observed in this workload.

## Removed and retained interfaces

- Removed production enums: JointSeedSamplingMethod, MapInterpolationMethod
  (four values); removed test-only WidthSource (two values).
- Removed production structs: JointSeedSamplingPolicy, JointSeedSamplingReport,
  JointSeedAtomSamplingDiagnostics, QScoreSeedSamples, QScoreShellDiagnostics,
  QScoreShellSamples, FirstStageTimingReport. Removed test structs SelectedWidth,
  WidthSelection, NeighborSelectionRow, NeighborSelectionShell, NeighborSelectionMasks.
- Removed callable dispatch/report interfaces: RunJointSeedSamplingWorkflow,
  JointSeedSamplingName, MapInterpolationName, InterpolateMapValue,
  BuildQScoreSeedSamples and the single-use RunPotentialSamplingWorkflowImpl
  template. Removed policy/report arguments from PrepareJointWidthInitialization
  and RunJointComponentWorkflow, interpolation arguments from sample conversion,
  and the optional call-timing pointer. No public MapSampler overload was removed.
- Removed research diagnostics, timer arrays/branches, geometry/mask controls,
  source-specific fallback/calibration/cohort/projection bridges and common-bound
  execution. No unpromoted selector survives as an unused enum or hidden branch.
- Retained deterministic Fibonacci (20 × 10, 0.05..1.95 Å), original search,
  hydrogen exclusion, strict Ownership + 30-degree cone and Tricubic.
- Retained QScore Spiral, distance exclusion/validation, adaptive retry, original
  center weighting, Trilinear and correlation.
- Retained alpha training, OLS/MDPDE estimator and records, failed-only recovery,
  MDPDE seeds, qualification/median fallback, 20/95% Type-7 calibration with
  4*p95 and bounded projection.
- Retained BCD, constrained A/C least squares, bounded LM, KKT, Independent
  Replay, Endpoint Certification and original-grid peeling/Local/Group workflow.

## Actual reduction

Metrics are measured from `git diff 503a6607` after the final report, with Git's
default rename detection. JSON bytes are reported separately from line counts.

| Metric | Actual |
| --- | ---: |
| Modified / deleted / added / renamed paths | 21 / 20 / 7 / 1 |
| Added lines, all files | 1261 |
| Removed lines, all files | 3096 |
| Net line reduction, all files | 1835 |
| Code/build added / removed lines (src/include/tests, excluding JSON/Markdown) | 587 / 2681 |
| Code/build net line reduction | 2094 |
| Production code net line reduction | 207 |
| Removed production / test-only enums | 2 / 1 |
| Removed production / test-support structs | 7 / 5 |
| Production if sites, in changed source set | 119 -> 106 (13 removed) |
| Production switch dispatches / case arms | 4 -> 1 / 8 -> 2 (3 dispatches, 6 arms removed) |
| GoogleTests, enabled | 1105 -> 1092 (net 13 fewer) |
| Registry names removed / added | 26 / 13, including retained tests moved or renamed |
| Removed CMake executable targets | 1 |
| Retired original tool/support paths | 8, including the fixture extracted once |
| Archived JSON files / bytes removed from active checkout | 9 / 2718293 |

The registry-name difference includes renamed/moved production tests and is
not a claim that 26 numerical contracts were discarded. The symbol inventory,
stage log and retained regression mapping explain the replacements.

## Git delivery and remaining concerns

| Commit | Verified change |
| --- | --- |
| a351d792 | C0 independent production fixture and regression protection |
| ee7aaf02 | C1 fixed production sampling/interpolation |
| 89727b60 | C2 QScore research-support removal |
| 38f94b04 | C3 width-source/timer cleanup and private calibration helper |
| 826e567c | C4 tooling/target removal and test consolidation |
| 9b9d1d75 | C5 historical report/evidence archive |
| e83c3bbc | C6 final regression contracts and exact build-mode reference |
| Final documentation commit | Complete verification, metrics and delivery report |

Initial HEAD is 503a6607674e068e4d3ef3db52aa770abcf3aea4. All changes are committed
on the existing develop branch. The final documentation commit's hash, final
working-tree state and ordinary push result are reported in the delivery reply.
No branch/worktree, detached checkout, reset, force push or Notion modification
was performed. The pre-existing other worktree was not used.

No cleanup item remains pending. Remaining related identifiers are required
production implementations: private Spiral geometry and calibration math,
private-module QScore distance predicate, formal OLS/MDPDE data and fitting,
and the regression synthetic fixture. MO/S1 short-token matches are unrelated
chemical/fixture identities. Retired names elsewhere are explicitly historical
reports/audit/archive records and have no production caller. No separate
follow-up algorithm change is required by this cleanup.

Scope limits: real validation uses the historical first 120 non-hydrogen targets
per structure, not full-model truth; only production F0 was rerun. Historical
research numbers were not regenerated. The immutable reachable Git ancestor is
the evidence archive, without an independent external backup. Timing uses only
three short repetitions and is not a performance claim.

## Stage file inventory

The following lists are taken from each committed stage diff. M/A/R denote modified, added or renamed; D denotes deleted. C6 acceptance adds only this report after the verified test commit.

| Stage | Commit | M/A/R files | D files |
| --- | --- | --- | --- |
| C0 | a351d792 | `docs/research/first-stage-sampling/cleanup-inventory.md`; `docs/research/first-stage-sampling/cleanup-report.md`; `tests/CMakeLists.txt`; `tests/core/joint_component/FirstStageInitialization_test.cpp`; `tests/core/joint_component/SeedSensitivity_test.cpp`; `tests/core/joint_component/WidthSource_test.cpp`; `tests/fixtures/joint_component/README.md`; `tests/fixtures/joint_component/first_stage.json`; `tests/research/FirstStageSamplingStudy.cpp`; `tests/support/FirstStageSamplingFixture.hpp → tests/support/FirstStageSynthetic.hpp` | — |
| C1 | ee7aaf02 | `docs/research/first-stage-sampling/cleanup-report.md`; `src/core/MapSampler.cpp`; `src/core/QScoreHelper.cpp`; `src/core/detail/FirstStageInitialization.cpp`; `src/core/detail/FirstStageInitialization.hpp`; `src/core/detail/MapInterpolation.hpp`; `src/core/detail/PotentialFittingWorkflow.cpp`; `src/core/detail/PotentialFittingWorkflow.hpp`; `tests/CMakeLists.txt`; `tests/core/MapInterpolation_test.cpp`; `tests/core/joint_component/SeedSampling_test.cpp` | `src/core/detail/JointSeedSampling.hpp`; `tests/core/joint_component/SeedSensitivity_test.cpp` |
| C2 | 89727b60 | `docs/research/first-stage-sampling/cleanup-report.md`; `include/rhbm_gem/utils/domain/SampleFilter.hpp`; `src/core/QScoreHelper.cpp`; `src/core/detail/QScoreSampling.cpp`; `src/core/detail/QScoreSampling.hpp`; `src/utils/domain/SampleFilter.cpp`; `tests/core/QScoreHelper_test.cpp`; `tests/utils/domain/SampleFilter_test.cpp` | — |
| C3 | 38f94b04 | `docs/research/first-stage-sampling/cleanup-report.md`; `src/core/detail/FirstStageInitialization.cpp`; `src/core/detail/FirstStageInitialization.hpp`; `src/core/detail/joint_component/WidthCalibration.cpp`; `src/core/detail/joint_component/WidthCalibration.hpp`; `tests/CMakeLists.txt` | `tests/core/joint_component/WidthSource_test.cpp`; `tests/support/ResearchWidthInitialization.hpp` |
| C4 | 826e567c | `docs/research/first-stage-sampling/cleanup-report.md`; `tests/CMakeLists.txt`; `tests/core/joint_component/FirstStageInitialization_test.cpp`; `tests/core/joint_component/SeedSampling_test.cpp → tests/core/joint_component/FirstStageSampling_test.cpp`; `tests/support/FirstStageSynthetic.hpp` | `tests/research/FirstStageSamplingStudy.cpp`; `tests/research/WidthInitializationStudy.hpp`; `tests/research/run_first_stage_sampling_study.py`; `tests/research/run_neighbor_selection_study.py`; `tests/research/run_width_initialization_study.py`; `tests/support/FibonacciNeighborSelection.hpp` |
| C5 | 9b9d1d75 | `docs/developer/joint-stage-integration.md`; `docs/research/first-stage-sampling/README.md`; `docs/research/first-stage-sampling/cleanup-report.md`; `docs/research/first-stage-sampling/evidence-archive.md`; `docs/research/first-stage-sampling/neighbor-selection-study.md`; `docs/research/first-stage-sampling/ols-vs-mdpde-initialization-study.md` | `docs/research/first-stage-sampling/6z6u-endpoint-audit.json`; `docs/research/first-stage-sampling/8rqb-endpoint-audit.json`; `docs/research/first-stage-sampling/baseline-6z6u.json`; `docs/research/first-stage-sampling/baseline.json`; `docs/research/first-stage-sampling/neighbor-selection-results.json`; `docs/research/first-stage-sampling/neighbor-selection-validation.json`; `docs/research/first-stage-sampling/ols-vs-mdpde-initialization-results.json`; `docs/research/first-stage-sampling/ols-vs-mdpde-validation.json`; `docs/research/first-stage-sampling/results.json` |
| C6 regression | e83c3bbc | `docs/research/first-stage-sampling/cleanup-inventory.md`; `docs/research/first-stage-sampling/cleanup-report.md`; `tests/core/joint_component/FirstStageInitialization_test.cpp`; `tests/core/joint_component/FirstStageSampling_test.cpp`; `tests/fixtures/joint_component/README.md`; `tests/fixtures/joint_component/first_stage.json` | — |

## Verification log digests

Commands were executed locally with the limits above. This index records the observed execution logs; ignored/tmp logs are supplemental rather than permanent research fixtures.

| Log | SHA-256 |
| --- | --- |
| firststage-final-release-build.log | `6eb3a3dc1ac841848f3e75a07b4dda34a021df315a9922508d579463967fdf93` |
| firststage-final-release-ctest.log | `6b872677b06722a988673ba1daa7fc1ce64556dd4355e5194d477b2967d137b3` |
| firststage-final-debug-build.log | `2554bec3635c1229132c7bbbb7e5810f1230f47c2e5024fb3bb5675f1997c3bb` |
| firststage-final-debug-ctest.log | `c8e1b7a264f119e40049766e89572f175305b30fc0ed786b3b7c8767e6cdacf1` |
| firststage-final-debug-focused.log | `bdc4256040da7f8d49a211f3b86360cb06e56dca6ed761f8d04f2ac06b195cba` |
| firststage-final-registration.log | `790437bd3fa99966f1022f9545bccf52be69d0f5f4175b960176be5b039ff16c` |
| firststage-final-lint.log | `8e41f6bba64e78a149b46222f84b6ce1f03086cf6d9cd1a28cdefbb0c343ad9c` |
| firststage-final-install.log | `d02ba4513ee26ad631c268223d7b05b9d471a79ad3172a9c2aad091018aedcd7` |
| firststage-original-qscore.log | `bc5523c11d260f56b74477cc06c7828837e026bd0501dc38ec2512f708cc0acb` |
| firststage-perf-before.csv | `3a031151525b924b9fc76a5524d64b69f1d4b9596a6d7c4d9386b07e81d9e83e` |
| firststage-perf-after.csv | `de58aacb8b69e0fae6be082c0953c4ff1eab1a011a25153695383049955d1ce5` |

Real one-off results match committed MDPDE baseline vectors/checks exactly. Their captured output SHA-256 values are 6Z6U `2b4948bd1c3dff90be796ea4b0ab70e79d692fd60e68caba0bdddfe9437336ff` and 8RQB `7e139723716bd1120e857a9495be271b2c8ff17c8b5fb4514a14448d2381d456`.
