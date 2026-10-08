# FixedNeighbor local certification comparison

This benchmark compares SearchOnlyLocal at its existing cheap-stationarity point with a CertifiedLocal strategy that assesses each core conditional candidate before committing it. The local assessment reuses the existing inner, width-gradient, local-correction, identifiability, and trust checks. It does not assemble a global projected-width matrix or write reference coefficients into the block state.

CertifiedLocal rejected both candidates in its first sweep for both 256-atom cases. On chain-256, local inner, gradient, identifiability, and trust checks passed, but local corrections were `3.878e-10` and `3.879e-10`, above `1e-10`. On cube-256, both local width gradients and corrections failed: gradients were `2.70e-11` and `2.78e-11`; corrections were `4.75e-8` and `4.91e-8`. No local candidate was committed, both searches stopped as uncertified, and both global RuntimeConvergence results failed.

The maximum local design was 58,300 × 256 on chain and 50,039 × 256 on cube, within the 128-atom core bound. Local assessments took `29.3s` and `28.7s` total per case. This does not satisfy the target gate and the CertifiedLocal path is not retained as a repair. The compact result is in `analysis.json`; raw measurements remain recoverable from Git history.
