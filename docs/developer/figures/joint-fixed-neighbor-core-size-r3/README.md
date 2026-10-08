# FixedNeighbor outer-core performance study

This search-only study keeps the production-compatible OneAcceptedLocalUpdate policy, forward serial Gauss-Seidel order, SPQR backend, one Eigen thread, and frozen stationarity checks. It varies only outer core_atoms over 8, 12, 16, 24, 32 for chain-512, cube-512.

The selected core is the minimum aggregate search time among core sizes that are search-correct on every study case. Local work, profile factor telemetry, RSS, and maximum local columns remain in the machine-readable artifact.

Correctness gate: **passed**. Selected core: **12**.
