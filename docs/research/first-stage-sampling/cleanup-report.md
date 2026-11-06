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
