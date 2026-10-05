# FixedNeighbor endpoint decomposition

This F2.5 diagnostic compares the raw FixedNeighbor state with primary and independent-reference global profiles at the exact same width vector. The full snapshots, including each state’s assessment and spectrum, are in `individual-results/`. `analysis.json` contains the machine classification.

The analyzer classifies the cube-256 trajectory as `mixed`. At sweep 11, same-width reprofiling resolves the inner coefficient disagreement (`1.047e-9` to `7.4e-15`) but does not pass endpoint certification: the profiled width gradient is `1.15e-12` and the local correction is `1.24e-9`. At sweep 13, without a profile being written back, the raw block state and both same-width profile controls pass all endpoint checks. Sweeps 14 and 15 remain passed with essentially unchanged search metrics.

This supports measuring the width-update contribution before selecting a single repair. It does not support A/C polish as a sufficient repair at sweep 11.
