# FixedNeighbor confirmed-stationarity 512 certification

This S5 campaign uses 128-atom cores, forward serial Gauss-Seidel, the local `LegacyCompact` solver, one Eigen thread, and `CertifiedLocal` disabled. The eta-only confirmation rule uses the existing KKT `1e-10`, width-gradient `1e-12`, and eta-change `1e-10` thresholds. No global profile coefficients are written back to the FixedNeighbor state.

Both chain-512 and cube-512 passed the full endpoint and `RuntimeConvergence` gates. Their global `LegacyCompact` and `OperatorPcg` baselines also passed `RuntimeConvergence`; objective and scaled A/C parity passed against `LegacyCompact`. The full runs stop at confirmed sweeps 7 and 17, respectively. Cube-512 first met the cheap criteria at sweep 13 and continued until coordinate confirmation at sweep 17.

The former full case JSON files contained sweep and block telemetry, endpoint checks, both global baselines, parameter/objective comparisons, total elapsed time, and process peak RSS. The per-run telemetry is historical execution data; the current tree keeps the qualification result and summary below. All three ranks are 512 for both topologies.

The gate result and per-case metrics are in `analysis.json` and `summary.csv`. Peak RSS is process-wide after endpoint assessment and global baselines; it is not a claim of constant memory use.
