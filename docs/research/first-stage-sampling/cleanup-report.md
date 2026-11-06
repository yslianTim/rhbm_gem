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
