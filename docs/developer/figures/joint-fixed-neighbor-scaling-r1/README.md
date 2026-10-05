# FixedNeighbor confirmed-stationarity nonlinear scaling

This S6 campaign measured chain and cube workloads at 128, 256, 512, 768, 1024, and 2048 atoms. It used SPQR with one Eigen thread, 128-atom cores, forward serial Gauss-Seidel, local `LegacyCompact`, confirmed stationarity, and `CertifiedLocal` disabled. The 512 points reuse the S5 endpoint-audit trajectories because the settings match exactly.

The 128–768 runs include full endpoint assessment. The 1024 and 2048 runs are search-only; their output scope is `fixed-neighbor-search-only`. No global endpoint claim is made for those four results.

The formal four-point gate uses 256, 512, 768, and 1024 atoms separately by topology. Confirmed sweeps have log-log slopes of 0.000 (chain) and 0.142 (cube), with endpoint ratios 1.000 and 1.214. Both remain below the predeclared growth gate (slope 0.25 and endpoint ratio 1.5), so the classification is **stable**. Cheap-sweep slopes and confirmation overhead are included in `analysis.json`; overhead remains bounded across the formal points.

Both 2048 frontier searches completed and confirmed stationarity. Chain used 7 sweeps in 4496.4 seconds; cube used 18 sweeps in 11314.9 seconds. Endpoint assessment was skipped as planned. Local factor columns remained bounded at 256 across all sizes; maximum local rows were 58300 (chain) and 52464 (cube). Process peak RSS varies by workload and scope and is not evidence of constant memory.

The formal scaling gate supports proceeding to converged order-sensitivity measurements. It does not by itself promote FixedNeighbor; LegacyCompact remains the production default. The raw measurements and per-case status are in `individual-results/`, `analysis.json`, and `summary.csv`.
