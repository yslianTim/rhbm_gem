# Production InitialStage regression

`initial_stage.json` records the 24-target coarse-phase-0 F0 + MDPDE baseline.
It covers atom identity, alpha, InitialStage OLS/MDPDE, qualified donor count and
coverage, the donor median `q50`, and the global upper bound `4 × q50`. It also
covers final Joint A/C/Width, objective, gradient and certification evidence.
The identical component check contract is stored once.

The fixture also captures isolated and crowded raw sampling fingerprints and
original-grid peeling/Local results.
Each fingerprint is FNV-1a over low-byte-first IEEE-754 words for positions,
distances (squared distance for grid observations), selection and responses,
with separate raw/selected counts. Comparisons are exact, without tolerance.
They protect ordering and rejected responses without retaining thousands of
raw rows. Direct serial/parallel comparisons retain per-sample diagnostics.

The shared normalized Gaussian generator lives in
`tests/support/InitialStageSynthetic.hpp`. `InitialStage_test.cpp` also covers
truth at two spacings, median fallback, certified endpoints, persistence and
grid invalidation.
`WidthCalibration_test.cpp` covers the 20-donor and 95%-coverage gates, the
median-based bound, and projection of the working initialization while preserving
the original estimates. Workflow and schema tests cover Group eligibility and
failed-only recovery. The fixture uses precise JSON parsing and exact comparisons.

QScore has separate exact Release and Debug references for the one-ULP
optimization-mode difference. Release expects 0.93585969694861482 and Debug
expects 0.93585969694861493.
