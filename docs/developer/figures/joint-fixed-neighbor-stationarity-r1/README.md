# FixedNeighbor stationarity-confirmation qualification

This S1 campaign records complete-sweep coordinate changes and compares three offline stop candidates against per-sweep endpoint assessments. The six controls are chain/cube at 32, 128, and 256 atoms. Search runs use 128-atom cores, forward serial Gauss-Seidel, local `LegacyCompact`, one Eigen thread, and the existing `stop_after_stationarity=false` experiment setting.

The sweep observer records the global A/C KKT, width gradient, objective, eta change, scale-aware beta change, accepted and unchanged blocks, and whether a prior complete sweep exists. Full endpoint assessment runs only in the benchmark observer. It does not change block-step acceptance or update the search state.

`analysis.json` compares eta-only, scaled-beta-only, and eta-plus-beta confirmation at the existing `1e-10` estimator-level thresholds. It reports each rule's first confirmation sweep, the first sweep passing all five endpoint checks and `RuntimeConvergence`, whether confirmation could stop too early, delay after certification, and `RuntimeConvergence` by sweep. The first coordinate-change observation is excluded from confirmation because no previous complete sweep exists.

The per-case JSON files preserve the complete search telemetry and per-sweep ground truth. `summary.csv` contains the sweep-level comparison used to choose the smallest safe confirmation rule. This qualification does not change numerical search behavior or production defaults.
