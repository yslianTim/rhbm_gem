# Current FixedNeighbor backend qualification

This is the bounded promotion evidence for the current Joint Component
estimator at commit `4f17f189fc90c04cd19cfffc3647c60a6f1d4cfd`.

Every run used the same production contract: FixedNeighbor, core 12, Forward
order, 30 maximum sweeps, at most one trusted accepted local update per block
visit, local LegacyCompact profile search, Eigen threads 1, and the existing
KKT/width-gradient/eta confirmation thresholds. The only numerical variable
was the configure-time sparse factor backend.

Both backends passed the small chain/cube endpoint cases, the observable/
nuisance and multi-component postprocessing workflows, and the partial-
selection workflow. A/C and log-width parity on chain-8 and cube-8 is at
floating-point roundoff level.

Completed resource probes show a reproducible SPQR wall-time advantage:

| case | EIGEN total (s) | SPQR total (s) | EIGEN RSS (MiB) | SPQR RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| chain-128 | 24.2025 | 10.7096 | 68.48 | 82.61 |
| cube-128 | 42.8090 | 22.5991 | 65.23 | 81.98 |
| chain-256 | 52.3671 | 24.7471 | 114.33 | 142.75 |
| cube-256 | 99.7612 | 55.7053 | 107.64 | 142.89 |

The 512 search-only probes were bounded at 60 seconds and timed out for both
backends; neither backend receives a completion claim from those probes.
SPQR is selected because its approximately 1.8–2.3x current FixedNeighbor
wall-time advantage is reproducible on both topologies and sizes that
completed, despite a moderate RSS increase. The dependency and deployment
cost is accepted for that production-level time benefit. Historical operator,
ordering, and benchmark campaigns are not part of this decision.

`analysis.json`, `summary.csv`, and `campaign-manifest.json` are the complete
compact evidence. Individual result files and process logs are intentionally
not committed.
