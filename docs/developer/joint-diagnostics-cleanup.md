# Joint diagnostics cleanup acceptance record

## C0 inventory and baseline

Initial develop HEAD and local origin/develop: `23ff0457ddb9c6004d1d14d344bea2e3696d5fa3`.
The checkout was clean. GitHub DNS lookup failed; origin/develop is a local
tracking reference, not a verified live GitHub HEAD. The existing detached
worktree is excluded from this task.

Release uses AppleClang, Ninja, system dependencies and OpenMP. Build and
CTest parallelism <=4. Benchmarks and snapshots use OMP_NUM_THREADS=1,
OMP_THREAD_LIMIT=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1,
VECLIB_MAXIMUM_THREADS=1. Full CTest uses OMP_NUM_THREADS=2 because the
Eigen thread-count scope tests exercise a two-thread setting; Joint groups
set OMP_NUM_THREADS=1. Baseline Release build and full CTest: 19/19 passed.

| Classification | Surface | Responsibility |
|---|---|---|
| Preserve | Numerics.cpp, Search.cpp, Replay, endpoint assessment | Strict local A/C certificate 1e-10; strict independent local replay; rank and feasibility |
| Preserve | ConvergenceCriteria, bounded LM, FixedNeighbor search/certification, Components assembly | practical-v2 global A/C 1e-7, projected width 1e-9; candidate acceptance, rollback/recovery, complete sweep confirmation |
| Preserve | LocalConditioning, MemoryStage::LocalConditioning | Separate local rank, spectrum and active-face evidence |
| Preserve | TargetRuntimeConvergence, IsQualifiedJointAnalysisResult, row masks, Halo, Peeling, Group | Formal numerical qualification and downstream behavior |
| Remove | ParameterStability.cpp/.hpp, SplitDiagnostics.cpp/.hpp, their CMake entries and dedicated tests | Perturbation search/classification, thresholds, seeds, representative selection, conditional reprofile |
| Remove | FixedNeighbor policy/work/result, Components and estimator remapping/budgets, ResourceWork | Diagnostic-only claims, counters, copies, workspace stage |
| Rewrite | JointEstimationTypes, TargetEvidence, JointResultJson, persistence tests, benchmark | Remove diagnostic types and fields; preserve solver state |
| Requires contract decision | JointAtomPointValidation and consumers | v2 quality uses convergence, validity and conditioning only; ActiveFace remains Warning; explicit NotAssessed stays NotAssessed |
| Remove/rewrite | Retired local stability and width-promotion studies, tracked data, indexes and current guides | Keep current global width contract in runtime documentation |

The audit covered public headers, src, tests/support, tests/benchmarks, JSON,
SQLite consumers, result remapping, CLI exports, and developer documentation.
No SQLite column depends on the retired members: it stores a versioned JSON
payload. Runtime target qualification merges component convergence, parameter
validity, selected target checks, row availability and assembly checks; it
does not use point-quality as a substitute for numerical qualification.
The stability module runs after endpoint qualification and receives copies
of the accepted local parameters; removing it must leave accepted state intact.
Dedicated stability tests assert diagnostic classification/reprofiling; retained
solver tests already cover independent replay veto, strict local KKT failure,
bounds, rank, recovery and unchanged accepted state.

Full round-trip-safe JSON numbers for eta, beta (interleaved A/C), B,
prediction, residual, objective, KKT, sweep/acceptance/recovery counts,
certification/replay, all persisted conditioning fields and per-point quality
are saved in `build/cleanup-baseline/states.jsonl` (also gzip). These ignored
local artifacts are preserved, not committed as a permanent research dataset.
Reproduce using the C0 commit's `rhbm_joint_cleanup_snapshot` executable.
Baseline file SHA-256: `a05dbeca73752e0dd65e25b8dfd16d7d18f195a9e25fffeb8b37f10d556d8f66`.
The snapshot fixture inputs also carry SHA-256 digests. Five measured runs
following one warmup are in `build/cleanup-baseline/performance.txt`.

| Fixture | Qualified | Stop reason | Old quality counts |
|---|---|---|---|
| high-overlap-cube-8 | True | block-stationary | {'Warning': 4, 'Available': 4} |
| rank-deficient-constant-support | False | block-recovery-exhausted | {'Unavailable': 3} |
| near-singular-overlapping-support | False | block-recovery-exhausted | {'Unavailable': 2} |
| low-signal-noise-only | True | block-stationary | {'Warning': 3} |
| weak-width-direction | True | block-stationary | {'Available': 2} |
| chain-3 | True | block-stationary | {'Available': 3} |
| chain-30 | True | block-stationary | {'Warning': 24, 'Available': 6} |
| chain-300 | True | block-stationary | {'Warning': 24, 'NotAssessed': 274, 'Available': 2} |
| chain-1000 | True | block-stationary | {'Warning': 36, 'NotAssessed': 964} |
| active-width | False | block-recovery-exhausted | {'Unavailable': 3} |

The active-width fixture starts below a 0.45 upper bound and subsequently
fails recovery with bounded-initial-width-infeasible. This existing failure
is recorded rather than changing the solver during cleanup. Dedicated width
bound, active LM and endpoint/replay tests remain the safety authority.
