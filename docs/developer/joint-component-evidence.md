# Joint component historical evidence

This document keeps the conclusions needed to understand the current Joint
estimator and its numerical limits. Current tools, permanent test owners and
production behavior are listed in [Joint experiments](joint-experiments.md)
and [Joint FixedNeighbor](joint-fixed-neighbor-experimental.md).

## Current estimator

The production route is Joint → FixedNeighbor → local Profile LM → conditional
A/C profile solve → EIGEN SparseQR. Its policy is core 12, Forward serial
Gauss-Seidel, at most 200 sweeps, and at most one trusted accepted local update
per block visit. Endpoint verification uses an independent streaming replay at
the accepted parameters and constrained A/C KKT; it does not run a second
global A/C solve. Full numerical qualification separately requires design
rank, local correction and numerical identifiability. The current thresholds
and versioned checks are recorded in the
[production contract](joint-fixed-neighbor-experimental.md).

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

- Current `ResourceWork` preflights use actual support memberships and a
  3,200,000,000-byte (3.200 GB, 2.980 GiB) estimated-peak budget, leaving 20%
  below the nominal 4 GB limit. The estimate is not process RSS. A resource
  rejection is reported as `Unavailable` with a stage-specific reason.
- Compact prepared row mappings remove the old `rows × blocks` allocation.
  At 5,938,087 rows and 3,118 blocks, that old index map alone is
  148,119,642,128 bytes (148.120 GB, about 137.9 GiB). The replacement's exact
  bytes require the actual affected-row and support-membership counts; those
  counts cannot be derived from rows and atoms alone.
- Global sweep stationarity and endpoint verification stream prediction,
  residual, A/C and width reductions from support entries. The local
  conditional solve remains EIGEN SparseQR. The independent global reference
  solver remains only an offline/test oracle; production does not build its
  dense reference `R`.
- For 37,406 atoms, a 74,812-column DesignSpectrum `ReferenceQR` factor alone
  is estimated at 44.775 GB (41.7 GiB). Production no longer runs a global
  independent reference solve, but DesignSpectrum remains the dense rank
  primitive when preflight permits it. Large rank and target assessments
  report unavailable when their stage estimates exceed budget; no approximate
  spectral result is treated as a certificate.
- A deterministic 37,406-atom, 5,938,087-row single-chain fixture has 74,811
  support memberships and 3,118 blocks. It records 2,892,736 prepared-mapping
  bytes and an 876,119,666-byte global-evaluation estimated peak. A direct
  no-sweep test process recorded 474,546,176 bytes peak RSS (0.475 GB,
  0.442 GiB). The fixture disables sweeps and endpoint certification; it does
  not represent the unknown physical support graph or qualify the full
  physical search and certification.
- A separate serial census adds endpoint streaming and both memory preflights.
  On its synthetic single-chain shape, the actual full A/C design has 74,812
  columns but only 37,406 unique supported rows. The production structural
  row-count check therefore proves full-design rank deficiency without QR or
  SVD; it does not give an exact rank or classify the 37,406-column free face.
  The dense DesignSpectrum estimate remains 189,484,131,163 bytes, projected
  width and local correction remain 61,453,268,315 bytes each, and target
  identifiability remains 317,422,079,259 bytes. Those dense and downstream
  stages are not run. The earlier process peak of 622,624,768 bytes covered
  input, preparation, verification and the rejected preflights only. Per-stage
  shape, timing and RSS details are in the
  [component memory census](joint-component-memory-census.md).
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
  production local solver remains Profile LM.
- The legacy endpoint-reference comparison passed its recorded fixed-state audits.
  On its historical host, single-512 analysis/export changed from 307.805 s
  to 300.564 s and sampled peak RSS from 2.500 GiB to 1.872 GiB. This
  incremental comparison is not a current performance guarantee.
- The historical 168-atom reference diagnostics remain useful evidence that
  KKT and coefficient agreement are distinct numerical guarantees. The
  sweep-107 coefficient difference was 1.2573e-9 while primary KKT was
  4.5219e-11; the sweep-150 difference was 2.8179e-12. These values are not
  results from replaying the new production contract and do not prove parameter
  accuracy from KKT alone.
- Earlier bounded-rank and projected-width/tail prototypes matched completed
  small controls, but large reconstructions or resource gates remained
  incomplete. The current bounded-rank study proves the 37,406-row versus
  74,812-column synthetic full-design obstruction by structural row count and
  keeps that one-way failure check in production. Its general streaming
  floating-point bounds remain experimental. Eigen SparseQR disagreed with
  dense SVD at the strict near-threshold equality and its 1,024-atom controls
  were stopped by memory preflight; it was not promoted. Full results,
  per-ordering fill and RSS measurements, and the promotion decision are in
  [bounded global rank research](joint-bounded-rank-research.md).
- An opt-in affected-row candidate replay experiment compared the incremental
  objective/residual delta with the existing full replay on one-sweep synthetic
  chain/cube workloads from 128 through 1024 atoms. Across 339 candidate
  comparisons, all acceptance decisions agreed; the maximum absolute objective
  difference was `1.85e-14`, the maximum residual difference was `1.12e-15`,
  and the maximum sweep-end cache replay error was `1.12e-15`. At 1024 atoms,
  full candidate replay took 1.463 s (chain) and 1.618 s (cube), compared with
  0.00243 s and 0.00245 s for the incremental calculation, including local
  core prediction. Those are candidate-phase timings from a Release synthetic
  benchmark, not measured total-search savings: every run still performed the
  full replay oracle. The experiment did not execute a separate no-oracle
  trajectory or cover borderline acceptance, weak-signal, rank-deficient and
  active-amplitude-boundary searches. Incremental replay remains experimental;
  production acceptance continues to use full replay.

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
