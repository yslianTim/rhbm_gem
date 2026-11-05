# Gauss–Seidel BCD naming and contract validation

## P0 baseline

- Initial HEAD: `0b99d5f893ea547813679e9cfa8ed4237eba657f` on `develop`.
- Working tree: clean. After `git fetch origin develop`, HEAD and
  `origin/develop` were identical (0 ahead, 0 behind).
- Existing configurations: `build/release` and `build/debug`, AppleClang,
  system dependencies, UMAP/ROOT/OpenMP enabled, Python bindings disabled.
- The existing Release cache had testing disabled and stale CTest entries.
  Reconfiguring with `cmake --preset release -DBUILD_TESTING=ON` restored the
  current 21-test registration. The stale run's missing research script and
  two failures caused by an initial two-thread cap were configuration artifacts.
- Current baseline: Release build passed; **21/21 CTest tests passed**, including
  GoogleTest registration coverage, Joint contract/workflow/numerical groups,
  persistence, CLI integration, bounded LM/recovery, replay, and halo exclusion.
- Resource limits: builds use `--parallel 4`; CTest uses `-j 1` with
  `OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 OPENBLAS_NUM_THREADS=1
  VECLIB_MAXIMUM_THREADS=1`. Joint test groups retain their one-thread setting.

The inventory covered `src`, `include`, `tests`, `docs`, CMake, Python bindings,
CLI, progress, memory and persistence. It found 652 matching lines and 65
unique old identifiers/strings. Each occurrence was classified as internal
symbol, source path, logger/progress, serialization, diagnostic, test or
documentation. Actual neighbor/support/halo relationships are outside this
rename. Search settings and component execution settings are separate types:
`BcdSearchPolicy` and `BcdExecutionPolicy`.

Before the change, Joint JSON was schema 22 with solver provenance v5. Its
validator accepted both the current 300-sweep policy and a historical
200-sweep policy, and permitted null solver metadata. SQLite v20 and Stage
JSON v4 are independent formats. The new Joint structure requires schema 23,
complete solver provenance v6 and only the current BCD policy.

A temporary test probe captured eight deterministic cases using existing
synthetic workloads and the coupled-pair adaptive-budget fixture: single
block, multiple blocks, active width upper bound, nominal budget exhaustion,
failed-block recovery, adaptive extension, memory rejection, and hard cap.
Both direct search and production component orchestration/assembly were run.
Floating-point values were captured as their `uint64_t` bit patterns, including
eta/tau, beta/A/C, prediction, residual, objective, sweep KKT, replay discrepancies
and width bounds. The probe also captured block acceptance, recovery,
certification/qualification, runtime convergence, progress order/counts,
profile/LM/damping-QR/replay counts and deterministic memory estimates.
Timing and measured RSS were excluded. No probe is compiled into production.

Baseline capture: 764 lines, 858826 bytes, SHA-256
`aa45a5e4256345c1410b087fd44a9844140e0629497ac5c4a4e103e498c973dc`.

## P1 core naming

Renamed the core source/header, budget, policy, search/component/sweep records,
endpoint certification, conditioning entrypoint and test instrumentation.
Search and execution policy responsibilities are unchanged. Release production
and test targets build successfully; the eight-case capture is byte-for-byte
identical to P0. Joint contract, workflow and numerical regression groups pass.
No compatibility wrappers/aliases were added, and no numerical expressions,
constants, thresholds or traversal logic changed. `git diff --check` passes.

## P2 progress, logger and memory

Progress uses `JointBcdProgress`, `JointProgressEvent.bcd` and
`JointSolverConfiguration.bcd_core_atoms`. The logger identifies Gauss-Seidel
BCD / Variable Projection / Bounded Active-set LM, with a method identifier in
debug detail. Preparation uses `MemoryStage::BcdPreparation` and
`bcd-preparation`. Memory formulas and checks, observer cadence and event order
are unchanged. Release build and all three Joint regression groups pass; new
logger/stage-label assertions pass. The complete capture remains identical to
P0, including progress values, LM/QR/replay counts and memory estimates.

## P3 persisted contract

Joint Result JSON is schema 23. Solver provenance requires
`search_method=gauss-seidel-bcd`, `contract_version=joint-solver-provenance-v6`,
`sparse_backend=EIGEN`, `bcd_core_atoms=12`,
`bcd_policy_version=gauss-seidel-bcd-production-v1`,
`bcd_maximum_sweeps=300`, `bcd_order=forward`, and
`bcd_local_search=active-set-lm`. The historical policy acceptance and nullable
solver reader paths are removed. Initialization/calibration failures now carry
the same versioned solver identity and retain their original unstarted-search
restrictions. SQLite and Stage versions are unchanged. Preparation/mapping
resource-limit reasons use the BCD prefix.

Release build, 21/21 full CTest tests, provenance round-trip, CLI output, legacy
method/policy/key and mixed-key rejection, unknown/missing/null field rejection
and invalid core/sweep/order tests pass. Pre-search failure round-trip retains
BCD provenance. The eight-case numerical/work/progress capture remains exactly
identical to P0. No migration, reader fallback or old-policy alias exists.

## P4 tests and documentation

Tests now use `GaussSeidelBcd_test.cpp`, `ComponentSweepBudget_test.cpp`,
`JointBcdDiagnostic`, `JointBcdSearchTest` and
`JointComponentSweepBudgetTest`. CMake source lists and persistence filters
match. A source-level declaration audit preserved all 950 original tests under
the naming map and added one provenance rejection test (951 total). Full
Release CTest passes 21/21, including the registration coverage audit.

Developer guides explain sequential accepted-state updates, Effective Response,
Variable Projection, Constrained A/C Active-set QR, Bounded Active-set LM and
schema 23. Source links point to the renamed files. The only remaining old
method/key strings are seven lines of explicit rejection data in the persistence
test. Production code and developer documentation contain no old algorithm
names. Legitimate neighbor/support/peeling concepts retain their names. A
case-insensitive scan also finds the second-stage offset-fit diagnostic at
`src/core/detail/second_stage/JointFitting.cpp:316`; it describes subtracting an
actual fixed neighbor's response and is intentionally unchanged.

An audit of all 256 production C++ files found that, after exact declared
identifier/string substitutions and comment removal, every file except the
provenance assembly, Joint JSON validator and schema constants is unchanged.
Those three differences implement only the new persisted contract.

## Numerical baseline comparison

The final Release capture with Python bindings enabled is byte-for-byte equal
to P0 (same SHA-256, length and line count). Both direct search and component
assembly/progress were captured. The following values are identical before
and after; KKT columns are the final sweep's A/C and projected Width values.

| Case | Objective | A/C KKT | Width KKT | Sweeps | Accepted updates | Profiles | Recovery success/attempt |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Single block | 2.0199078132865863e-10 | 3.8857805861880479e-16 | 4.0945025148175773e-13 | 3 | 3 | 6 | 0/0 |
| Three forward blocks | 2.0199079006559097e-10 | 3.1542578549093037e-09 | 1.1782208542143735e-10 | 6 | 18 | 36 | 0/0 |
| Active width upper bound | 2.0254002647560009e-10 | 3.3306690738754696e-16 | 1.1102230246251565e-16 | 2 | 2 | 4 | 0/0 |
| Nominal sweep limit | 2.0958526701391608e-10 | 2.9403302566066714e-06 | 1.0939575745272379e-07 | 4 | 12 | 24 | 0/0 |
| Failed-block recovery | 1.7934999880322473e-10 | 1.7828503395822537e-08 | 8.1889350855846033e-10 | 5 | 9 | 20 | 1/1 |
| Adaptive extension | 8.9386328233049828e-11 | 7.2592128419479707e-09 | 6.9365724275627372e-10 | 34 | 68 | 136 | 0/0 |
| Memory rejection | unavailable | unavailable | unavailable | 0 | 0 | 0 | 0/0 |
| Hard sweep cap | 8.938746938111263e-11 | 3.0175387077946425e-08 | 2.883422611787978e-09 | 31 | 62 | 124 | 0/0 |

All A/tau/C, prediction/residual and objective bit patterns are unchanged, as are
block acceptance decisions, termination reasons, qualification, independent
replay, runtime convergence, width bounds, observer events and deterministic
workspace estimates. The active-bound case retains one active width bound.
Subproblem/profile/damping-QR/replay counts are unchanged; no extra block
solve, replay or dense workspace was introduced. No numerical tolerance or
convergence threshold changed. Timing was not used as an equivalence criterion.

## P5 final verification

| Check | Result |
| --- | --- |
| Release build, existing preset, Python bindings enabled | Passed |
| Release full CTest | 23/23 passed, 0 failed (36.10 s) |
| Debug build, existing preset | Passed |
| Debug full CTest with local runner time allowance | 21/21 passed, 0 failed (392.83 s) |
| GoogleTest registration coverage, both configurations | 1061 enabled / 1061 uniquely covered, 0 missing |
| Repository lint | Passed |
| Installed C++ consumer build, public API, JSON/SQLite round-trip | Passed |
| Eight-case numerical/work/progress capture | Exact P0 byte equality |
| Renamed includes / GitHub source links | No unresolved dependencies |
| Diff whitespace check | Passed |

The first Debug CTest run passed 20/21 groups and hit the existing 180-second
runner timeout in the numerical group, without an assertion failure. For the
complete rerun, only that group's generated `build/debug/tests/CTestTestfile.cmake`
TIMEOUT was temporarily raised from 180 to 1800 seconds. All 107 numerical
GoogleTests passed in 233.47 seconds; the 300-atom full Replay audit alone took
199.258 seconds. The generated file was restored afterward. No repository
timeout, compiler flag, test data, numerical expectation or tolerance changed.

The installed consumer fixture also required complete BCD provenance and had
a stale schema-21 assertion. It now uses schema 23, verifies method/contract
identity, and round-trips the current policy through the installed public API
and SQLite. The strict reader was retained.

All checks ran on macOS with AppleClang 21.0.0.21000334. No Linux/Windows or
external real-map data run is claimed. Builds were limited to four jobs; tests
were sequential with at most four runtime threads. No new dependency, branch
or worktree was introduced, and Notion was not accessed. External Notion links
to renamed source files may need a separate update by their owners.
