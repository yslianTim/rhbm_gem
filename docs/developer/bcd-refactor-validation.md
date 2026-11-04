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
eta/B, beta/A/C, prediction, residual, objective, sweep KKT, replay discrepancies
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
