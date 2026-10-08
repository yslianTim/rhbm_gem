# FixedNeighbor confirmed-stationarity requalification

This S3 campaign reruns the six nonlinear controls using the qualified eta-only outer stopping rule. All runs use 128-atom cores, forward serial Gauss-Seidel, local `LegacyCompact`, one Eigen thread, and `CertifiedLocal` disabled. The endpoint is assessed after the search and is independent of the confirmation condition.

The former full result files included first-order and confirmed stationarity sweeps, sweep telemetry, endpoint checks, `RuntimeConvergence`, global `LegacyCompact` and `OperatorPcg` references, and objective/parameter comparisons. The current tree retains the F2 endpoint and parity gate in `analysis.json`; `summary.csv` compares the prior cheap-stop sweep with the confirmed sweep and records the final endpoint metrics.

All six cases passed endpoint certification and `RuntimeConvergence`; the global objective and scaled A/C parity gate also passed. The single-block 32/128 cases use one additional sweep. Chain-256 confirms at sweep 7 after the cheap candidate at 6. Cube-256 confirms at sweep 14 after the cheap candidate at 11; sweep 13 was the first fully certified state.
