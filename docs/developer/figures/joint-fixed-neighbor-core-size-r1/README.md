# FixedNeighbor core-size performance study

This P5 study keeps the P4-selected OneAcceptedLocalUpdate policy, forward serial Gauss-Seidel order, SPQR backend, one Eigen thread, and frozen stationarity checks. It varies only core_atoms over 64, 128, and 256 for chain-512, cube-512, and cube-1024. All runs are search-only; no full endpoint claim is made for 1024.

The selected core is the minimum aggregate search time among core sizes that are search-correct on every study case. Local work, profile factor telemetry, RSS, and maximum local columns remain in the machine-readable artifact.

Correctness gate: **passed**. Selected core: **64**.
