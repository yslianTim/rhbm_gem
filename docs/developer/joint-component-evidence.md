# Joint component evidence and historical findings

This document separates the current production contract from retained research
conclusions. Current implementation and test owners are listed in [Joint
experiments](joint-experiments.md) and [Joint component runtime](joint-component-runtime.md).

## Current estimator

Production uses FixedNeighbor with core size 12, forward serial Gauss-Seidel,
at most 200 sweeps, and at most one trusted accepted Profile LM update per
block visit. Selected targets are the only free FullABC parameters. Supported
halos contribute validated FirstStage `(A, B, 0)` values as fixed responses;
the parent-normalized objective uses original observations for its scale.
Target-only projected A/C KKT, width-gradient stationarity, complete sweep
confirmation, monotone global-objective acceptance, trust checks, and
independent full-prediction endpoint replay define numerical convergence.

A converged endpoint becomes the component state even when local conditioning
is poor or parameter probes show sensitivity. Local conditioning describes
target free A/C, projected-width, and full profile-derivative behavior within
each FixedNeighbor block while outside targets and fixed halos are held fixed.
Parameter stability perturbs targets only. Neither provides a global
identifiability certificate.

Conditioning classes distinguish full-rank/ill-conditioned matrices from
well-conditioned matrices using the provisional column-normalized
`kappa * epsilon_double >= 1e-8` warning. The condition number describes local
roundoff amplification; it is not an estimate of statistical parameter error.
Stability selection sorts eligible blocks by conditioning class, available
local condition score, and block identity. Selected blocks probe a weak
projected-width direction, a related coordinate direction, and common mode,
subject to the fixed block, seed, evaluation, and CPU budgets. Parameter and
prediction changes are reported separately, and an incomplete probe is never
reported as stable.

Target point validation uses atom-to-block indexes with work proportional to
points, blocks, and core memberships rather than scanning all blocks for every
point. Invalid diagnostic mappings are separately reported and do not change
runtime convergence.

Target records use `target-point-validation-v1` to report component convergence,
parameter validity, local conditioning, parameter stability, point quality, and
reasons conditional on the fixed FirstStage halo contribution. This numerical
qualification does not certify halo statistical accuracy. Historical
`target-estimability-v1` through `v4` and their rank records keep their original
meanings in persistence and are not converted to new evidence. Post-fit
peeling, local MDPDE covariance, and group fitting keep their own validation
and qualification steps.

## Current memory controls

`ResourceWork` preflights problem construction, structural construction,
FixedNeighbor preparation, global streaming evaluation, endpoint verification,
reference solves, local conditioning, and parameter stability. It uses checked
byte arithmetic and a 3,200,000,000-byte estimated peak budget below the
4,000,000,000-byte hard limit. The estimate is not observed RSS.

The search and endpoint checks stream predictions, residuals, KKT reductions,
and width gradients over support entries. Production does not allocate global
free-design factors, projection coefficients, or projected-width factors.
Local conditioning and perturbation stability use bounded block workspaces.
Construction and FixedNeighbor preparation remain independent memory risks.

A deterministic synthetic shape with 5,938,087 rows, 37,406 contributors,
and 3,118 blocks previously measured 2,892,736 prepared-mapping bytes and an
876,119,666-byte global-evaluation estimate. A no-sweep preparation test
process recorded 474,546,176 bytes peak RSS. These results cover the synthetic
fixture only; they do not measure a physical-data fitting run or its unknown
support graph. See the [memory census](joint-component-memory-census.md) for
which historical measurements have been retired.

## Historical FixedNeighbor decisions

- Fixed-B conditional updates were feasible but were not promoted; FixedNeighbor
  superseded them.
- Recorded chain/cube experiments selected core size 12. Against the historic
  core-64 control, median search time improved by 33.0%–64.3% at 512 and 1024
  atoms. These are selection measurements, not a current benchmark guarantee.
- Forward was adopted as the serial Gauss-Seidel order. Reverse was an offline
  comparison only.
- OneAccepted was selected over Full and TwoAccepted local update policies and
  remains intrinsic to each block visit.
- Prepared mappings remain covered by deterministic tests. They do not reuse
  numerical factors or accepted states across sweeps.

Further historical choices are in [Joint FixedNeighbor](joint-fixed-neighbor-experimental.md).

## Retired global-rank research

The former global rank qualification pipeline and its memory planner have been
removed. Earlier dense, bounded, projected-width, and weak-direction studies
remain available as research records in [bounded global rank research](joint-bounded-rank-research.md)
and [Joint component memory census](joint-component-memory-census.md). Their
rank thresholds, allocation estimates, and unavailable outcomes describe the
retired pipeline. They are not current runtime checks, memory estimates, or
production prerequisites.

Historical profile-derivative and high-precision tests still cover local
numerical primitives. Local SparseQR active-set checks, CompactSvd, KKT, and
profile derivative utilities remain in the production solver where they are
needed.

## Scientific and workflow limits

Local full rank does not rule out cross-block degeneracy. A bounded local
perturbation study does not replace a global multi-start or block-order study.
Joint point-estimate error can propagate through peeling. The Joint fit does
not produce a global parameter covariance; local MDPDE covariance and group
fitting retain their separate workflows.

Historical weak-halo observations showed that a numerically full-rank design
can still have poor practical sensitivity. Small partial-selection controls
showed that omitted halo contributors can shift target estimates. The retained
noise and mismatch research is documented in
[joint-noise-mismatch-validation.md](joint-noise-mismatch-validation.md), with
its current receipts in `figures/joint-validation/`.

The 6Z6U Local Profile LM replay study and FullABC width study are bounded to
their recorded inputs and host. They document a sweep 4 / block 21 local A/C
KKT failure and an accepted large-width step followed by severe cancellation.
Neither study established a full-data converged endpoint or target recovery.
See [local replay validation](../research/joint-local-replay-validation.md)
and the [FullABC width study](../research/joint-fullabc-width-degeneracy.md).

Closed campaign code and raw receipts are archived in Git history. Historical
results apply to their saved source, host, and cases; they are not current
performance claims.
