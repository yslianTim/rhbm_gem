# Production FirstStage regression

`first_stage.json` keeps the 24-target coarse-phase-0 F0 + MDPDE baseline from
commit `503a6607674e068e4d3ef3db52aa770abcf3aea4`. Atom identity, alpha, OLS/MDPDE,
calibration and complete final A/C/Width/gradient values are copied without
numerical changes from the original committed study baseline. The identical
component check contract is stored once. No time/RSS or alternative policies
are reference data.

C0 additionally captures isolated and crowded raw sampling fingerprints and
original-grid peeling/Local results on the same unmodified production code.
Each fingerprint is FNV-1a over low-byte-first IEEE-754 words for positions,
distances (squared distance for grid observations), selection and responses,
with separate raw/selected counts. Comparisons are exact, without tolerance.
They protect ordering and rejected responses without retaining thousands of
raw rows. Direct serial/parallel comparisons retain per-sample diagnostics.

The shared normalized Gaussian generator lives in
`tests/support/FirstStageSynthetic.hpp`; it is independent of retired study
runners. `InitialStage_test.cpp` also covers truth at two spacings,
median fallback, certified endpoints, persistence and grid invalidation.
`WidthCalibration_test.cpp` owns donor gating, Type-7 quantiles, 20/95% limits
and immutable bounded projection. Existing workflow/schema tests cover Group
eligibility and failed-only recovery. The fixture uses precise JSON parsing;
platform changes must be investigated rather than silently regenerating it.

QScore keeps separate exact Release and Debug references. Its original Debug
implementation from commit `503a6607`, rebuilt with the current Debug compiler
flags, returns 0.93585969694861493; the original Release reference is
0.93585969694861482. This pre-existing one-ULP optimization-mode difference is
not a cleanup regression. Both modes use exact equality, without loosening the
correlation tolerance or changing production code.
