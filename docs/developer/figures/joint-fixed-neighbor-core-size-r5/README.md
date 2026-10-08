# FixedNeighbor outer-core repeated frontier

This P4 frontier keeps LegacyCompact, OneAcceptedUpdate, Forward serial Gauss-Seidel, SPQR, one Eigen thread, and the existing search gates. It uses one warmup and 3 interleaved measured repetitions. The 512-atom cases compare finalists 12, 16 against control 64; the fastest correct finalist by aggregate median search time is then compared with control on both 1024-atom cases.

Gate C4 requires at least 10% median search improvement for the selected core on chain and cube at both sizes, no material RSS/total-time/local-geometry regression, and correctness on every measured run. Sweeps above 24 are reported as warnings.

Qualification gate: **passed**. Selected core: **12**.
