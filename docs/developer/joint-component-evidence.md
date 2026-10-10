# Joint component historical evidence

This document keeps the conclusions needed to understand the current Joint
estimator and its numerical limits. Current tools, permanent test owners and
production behavior are listed in [Joint experiments](joint-experiments.md)
and [Joint FixedNeighbor](joint-fixed-neighbor-experimental.md).

## Current estimator

The production route is Joint → FixedNeighbor → local LegacyCompact →
LinearSolve → EIGEN. Its policy is core 12, Forward serial Gauss-Seidel,
at most 30 sweeps, and at most one trusted accepted local update per block
visit. The frozen KKT, width-gradient, eta-confirmation, replay and endpoint
contracts are recorded in the [production contract](joint-fixed-neighbor-experimental.md).

## Historical algorithm decisions

- Fixed-B conditional block updates and replay accounting were numerically
  feasible, but Fixed-B was not promoted. FixedNeighbor superseded it.
- Repeated qualification selected core 12 over tested larger cores. Across the
  recorded chain/cube 512 and 1024 cases, core 12 reduced median search time
  by 33.0%–64.3% against the historical core-64 control. These are decision
  measurements, not a current benchmark option.
- Forward became the production serial Gauss-Seidel order. Reverse was a
  qualification diagnostic only.
- OneAccepted was selected over Full and TwoAccepted local updates. The
  at-most-one trusted accepted update is now an intrinsic production rule.
- Prepared mappings remain covered by deterministic tests; preparation does
  not add a convergence rule or reuse numeric factors across sweeps.

The complete concise record of those FixedNeighbor choices is in
[joint-fixed-neighbor-experimental.md](joint-fixed-neighbor-experimental.md).

## Numerical and resource findings

- Historical fixed-state compact-SVD checks passed their recorded numerical
  gates. The campaign showed substantial reductions in selected compact
  reference and derivative calculations, but did not establish a current
  end-to-end speedup.
- A tiled derivative study reduced sampled peak RSS from 1,252 MiB to 359 MiB
  on one heterogeneous-168 input. That is a host- and case-specific storage
  result, not a general scalability guarantee.
- Fixed normal actions used two Q/Q′ calls versus six in the composed
  comparison. This fixed-state result did not justify a different global
  search route.
- Same-width endpoint reprofiling could remove an inner coefficient
  disagreement at one checkpoint without passing the width-gradient and
  local-correction endpoint gates. Later checkpoints passed with the raw
  state and no profile write-back. A/C polish alone was therefore not a
  sufficient repair.
- Per-block CertifiedLocal rejected both first-sweep candidates in the
  recorded chain-256 and cube-256 studies; it did not pass the unchanged
  endpoint and RuntimeConvergence gates. It was not promoted.
- Persistent workspace residency did not produce a material wall-time gain in
  the repeated chain/cube 256–1024 study and increased sampled RSS. Exact
  numeric-factor reuse opportunities were limited to initial-profile
  evaluations; no factor-cache policy was adopted.
- Tested derivative tile sizes did not reach the 10% improvement gate. The
  production local solver remains LegacyCompact.
- The endpoint-reference comparison passed its recorded fixed-state audits.
  On its historical host, single-512 analysis/export changed from 307.805 s
  to 300.564 s and sampled peak RSS from 2.500 GiB to 1.872 GiB. This
  incremental comparison is not a current performance guarantee.
- Bounded rank and projected-width/tail prototypes matched completed small
  controls, but large reconstructions or resource gates remained incomplete.
  Unavailable and timed-out results do not establish deficient rank or
  numerical failure.

The sparse-backend and global OperatorPcg/Schwarz conclusions, including the
bounded EIGEN/SPQR timing and memory comparison, are summarized in
[joint-operator-search.md](joint-operator-search.md). SPQR passed bounded
parity and was faster on completed 128/256 probes; EIGEN used less peak RSS
there. Both bounded 512 probes timed out. The original wall-time-first
qualification selected SPQR; later production policy prioritized memory,
dependency simplicity and maintenance, so EIGEN became the sole backend.

## Structural, workflow and scientific limits

The historical full-ABC 6Z6U study found 20 singleton halo atoms whose A and C
columns were proportional, bounding the design rank at 4,314 of 4,334 for
that fixed support. Observable-halo profiling removed that diagnosed initial
obstruction without changing the observation rows, but a complete 6Z6U
solve and qualification were not established.

Weak-halo evidence shows that a numerically full-rank design can still have
poor practical sensitivity. The active weak-halo diagnosis remains available
under [joint-validation](figures/joint-validation/weak-halo.json); its
interpretation and limits are described in the
[weak-halo guide](joint-weak-halo-attribution.md).

Small partial-selection controls showed that omitting halo contributors could
cause target compensation, while complete contributor closure recovered
near-zero objective and target errors. They did not establish large-scale
performance. Current partial-selection behavior has a current measurement
driver and permanent test coverage.

Historical stage and persistence comparisons preserved explicit unavailable
states and downstream invalidation. They did not establish a general
postprocessing speedup. Current persistence compatibility is owned by the
ordinary workflow and schema tests, not by historical experiment receipts.

The active noise/mismatch research is documented in
[joint-noise-mismatch-validation.md](joint-noise-mismatch-validation.md), with
its current machine-readable results under
[figures/joint-validation](figures/joint-validation/). The fixed-seed pilot is
unbalanced, covers one noise-correlation length and one position-mismatch
type, and is not a population success-rate, coverage, or real-data validity
study.

Historical connected-component probes did not establish a completed 6Z6U
workflow or a general large-component scalability bound. A timeout or resource
stop is neither a numerical failure nor a convergence result. Process
completion alone is not scientific qualification, and lower objective does
not prove parameter recovery or a global optimum.

## Evidence retention

Current machine-readable research and diagnosis are retained only where they
answer an active question: weak-halo diagnosis and noise/mismatch research in
figures/joint-validation. Permanent fixtures and regression tests remain in
their ordinary test locations.

Closed experiment implementations, campaign drivers, analysis JSON, CSV,
manifests and raw outputs are intentionally absent from the current tree.
Git history is the archive; retrieve a historical file with git show
COMMIT:path. These historical results are bounded to their recorded source,
host and cases and are not current performance claims.
