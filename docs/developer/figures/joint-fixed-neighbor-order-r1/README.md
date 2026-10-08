# FixedNeighbor converged order sensitivity

This S7 campaign compared forward and reverse serial Gauss-Seidel sweeps with SPQR, one Eigen thread, 128-atom cores, local `LegacyCompact`, and `CertifiedLocal` disabled. Both orders used the confirmed-stationarity rule: A/C KKT `<=1e-10`, width gradient `<=1e-12`, and eta change `<=1e-10` from the previous complete sweep.

Only search-converged endpoints were compared. The scale-aware A/C difference is the largest absolute coefficient difference multiplied by the larger support-derived scale weight from the two endpoints. The `1e-10` eta and scaled A/C interpretation bounds reuse existing confirmation and coefficient-comparison semantics. Objective differences are reported without adding a new objective threshold.

| Case | Cheap sweep F/R | Confirmed sweep F/R | Objective difference | Eta infinity difference | Scaled A/C difference | Endpoint and runtime |
|---|---:|---:|---:|---:|---:|---|
| chain-256 | 6 / 6 | 7 / 7 | `1.10e-17` | `9.62e-12` | `1.12e-12` | both passed |
| chain-512 | 6 / 6 | 7 / 7 | `1.89e-17` | `9.62e-12` | `5.27e-13` | both passed |
| cube-256 | 11 / 11 | 14 / 14 | `1.69e-15` | `5.86e-12` | `6.56e-13` | both passed |
| cube-512 | 13 / 13 | 17 / 16 | `4.64e-16` | `2.96e-11` | `1.85e-12` | both passed |
| cube-1024 | 14 / 14 | 17 / 18 | `5.91e-16` | `3.75e-11` | `2.09e-12` | search-only |

All five pairs converged in both orders. The largest parameter differences remain below the existing `1e-10` interpretation bounds; the maximum confirmed-sweep difference is one. For all fully assessed 256/512 cases, both endpoint certification and `RuntimeConvergence` passed in both orders. The 1024 case is explicitly search-only and carries no endpoint or runtime-convergence claim.

Global KKT and width-gradient values are preserved for each direction in `order-analysis.json`. The former raw case JSONs are historical execution telemetry. The largest observed local factor width was 256 columns. This is bounded by the 128-atom core, while process RSS remains workload-wide and is not constant-memory evidence.

The order-sensitivity gate passes, and coordinated/shared-parameter blocks remain deferred. FixedNeighbor is not promoted to a selectable estimator route: 1024/2048 scaling runs did not include full endpoint assessment, and the campaign has no matched large-size baseline establishing an end-to-end resource benefit. `LegacyCompact` remains the production default. `CertifiedLocal` remains diagnostic support.

The current tree retains the machine-readable comparison and CSV in `order-analysis.json` and `summary.csv`; raw per-run measurements remain recoverable from Git history.
