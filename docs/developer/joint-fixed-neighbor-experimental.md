# Fixed-neighbor joint search: qualification evidence

This document records the historical qualification and resource evidence for
FixedNeighbor. It is now a supported production Joint search method selected
explicitly through `FitOptions::joint_search_method` or
`--joint-search fixed-neighbor`; `LegacyCompact` remains the default and no
automatic route selection is added. The broader Full/128-core and diagnostic
controls described below remain historical reproduction controls, not the
production policy.

## Fixed-B qualification

The cube-256 fixed-B plateau was classified as a global roundoff lock. The local conditional solves continued to return valid certificates and nonzero updates, but the global replay could increase by about `9.91e-16`, within the existing objective replay enclosure of about `1e-12`. The existing strict global acceptance test rejected those updates, leaving the sweep KKT at `7.0516e-9`.

The selected repair makes block acceptance roundoff-aware only when the local conditional objective does not increase, the local/global delta identity is enclosed by the replay error, and any positive global replay delta is within the existing replay enclosure. The global KKT target remains `1e-10`.

After that repair, fixed-B chain/cube controls passed at 256 atoms in both orders. Forward scaling remained stable-looking through 1024 atoms: chain required `5, 5, 5` sweeps at 256/512/1024; cube required `9, 11, 11`. The cube trend has only three multi-block sizes, so the formal four-point growth gate was not evaluated. The 1024 cube run reached about 1.48 GB peak RSS.

Detailed evidence is in:

- [fixed-B floor attribution](figures/joint-fixed-b-floor-n1/analysis.json)
- [fixed-B coupling diagnostics](figures/joint-fixed-b-coupling-n1/coupling.json)
- [fixed-B requalification](figures/joint-fixed-b-requal-r1/analysis.json)
- [fixed-B scaling](figures/joint-fixed-b-scaling-r1/analysis.json)

## Initial nonlinear endpoint diagnosis

The initial nonlinear experiment stopped at the first cheap stationarity sweep. It used disjoint 128-atom cores, serial Gauss-Seidel updates, local `LegacyCompact` width search, and the existing global replay and endpoint certification. Global `LegacyCompact` and `OperatorPcg` returned runtime-converged references for all six cases.

| Case | Sweeps | FixedNeighbor objective | Global A/C KKT | Width gradient infinity norm | Scaled A/C difference vs LegacyCompact | Endpoint | Runtime |
|---|---:|---:|---:|---:|---:|---|---|
| chain-32 | 1 | `1.5342583e-10` | `1.80e-16` | `2.11e-17` | `2.18e-16` | passed | passed |
| cube-32 | 1 | `9.8793029e-11` | `1.25e-16` | `2.19e-17` | `1.51e-14` | passed | passed |
| chain-128 | 1 | `1.4868379e-10` | `9.71e-17` | `7.46e-18` | `1.21e-16` | passed | passed |
| cube-128 | 1 | `8.7674368e-11` | `1.11e-16` | `1.56e-16` | `1.47e-16` | passed | passed |
| chain-256 | 6 | `1.4847246e-10` | `3.99e-12` | `1.60e-14` | `7.62e-12` | passed | passed |
| cube-256 | 11 | `8.2223286e-11` | `5.66e-11` | `2.19e-13` | `1.38e-10` | **failed** | **failed** |

All six FixedNeighbor objectives were within the experiment's objective parity tolerance. Cube-256 also met the global A/C KKT target, but the existing endpoint checks failed: inner coefficient difference was `1.05e-9` against `1e-10`, width stationarity was `1.15e-12` against `1e-12`, and local correction was `1.24e-9` against `1e-10`. Its endpoint trust result was `reference-disagreement`. Objective parity and KKT alone do not certify this endpoint.

An exploratory continuation past the initial stationary point reduced cube-256 KKT to `1.97e-13` by sweep 14, then produced the same KKT at sweep 15. Endpoint certification had still not passed through sweep 14. This did not identify a supported repair, so no sweep-budget or certification semantics were changed.

That initial F2 result identified premature stopping, rather than a failed block update. Stationarity qualification selected eta-only coordinate confirmation at the existing `1e-10` change threshold. With that stopping rule, all six nonlinear 32/128/256 cases passed endpoint certification and `RuntimeConvergence`; chain-512 and cube-512 also passed. The qualification, six-case requalification, and 512 gate are documented in [stationarity confirmation](figures/joint-fixed-neighbor-stationarity-r1/README.md), [nonlinear requalification](figures/joint-fixed-neighbor-confirmed-r1/README.md), and [512 certification](figures/joint-fixed-neighbor-confirmed-512-r1/README.md).

Formal 256/512/768/1024 nonlinear sweep scaling was stable for both topologies, with confirmation overhead bounded at one sweep for chain and at most four for cube. Chain and cube 2048 search-only frontier runs completed with confirmed stationarity. Local factor columns remained bounded at 256; process peak RSS varies and is not constant-memory evidence. See the [scaling campaign](figures/joint-fixed-neighbor-scaling-r1/README.md).

The [converged order campaign](figures/joint-fixed-neighbor-order-r1/README.md) passed: all five forward/reverse pairs confirmed convergence, maximum eta infinity difference was `3.75e-11`, maximum scaled A/C difference was `2.09e-12`, and fully assessed endpoint/runtime results agreed and passed through 512. Cube-1024 remains search-only. Coordinated/shared-parameter blocks remain deferred because the measured endpoints show no material order divergence, sweep growth, or observed persistent oscillation.

The historical qualification established bounded local factor width and
confirmed search, but 1024/2048 did not receive full endpoint assessment.
Those results remain search-only evidence; they do not imply full estimator
certification. `LegacyCompact` remains the default, and `CertifiedLocal`
remains available only for diagnostics and reproduction.

## Local-work and core-size qualification

The matched local-work attribution found substantial repeated local factor work in the
Full policy. The bounded prototype therefore counted trusted accepted updates, not rejected
trials, while leaving the global replay and every endpoint threshold unchanged.

- `OneAcceptedUpdate` was selected over `Full` and `TwoAcceptedUpdates`: all four 256/512-atom
  chain/cube full-endpoint cases passed endpoint certification and `RuntimeConvergence`, and
  `OneAcceptedUpdate` was the fastest policy in each case.
- The 1024 qualification was intentionally search-only. OneAccepted reduced local work by
  67.607% (chain) and 75.414% (cube), with search reductions of 66.798% and 74.722%.
  The final A/C KKT and width-gradient checks remained at the existing `1e-10` and `1e-12`
  thresholds.
- The core-size study selected 64 atoms from 64/128/256 on chain-512, cube-512, and cube-1024.
  Every point passed the search correctness gate; aggregate search seconds were 1226.860,
  2598.111, and 6267.089 for cores 64, 128, and 256.

The machine-readable records are in the [inexact qualification](figures/joint-fixed-neighbor-inexact-qualification-r2/analysis.json)
and [core-size study](figures/joint-fixed-neighbor-core-size-r1/analysis.json) artifacts.

## Optimized matched frontier

P6 compared the P1 LegacyCompact, OperatorPcg, and original Full/128-core FixedNeighbor
records with OneAccepted/64-core FixedNeighbor on the same six search-only workloads. The
formal envelope was 600 seconds and 4 GiB RSS; non-passing formal runs received a 7200-second
diagnostic retry.

| Case | OperatorPcg search / RSS | Original FixedNeighbor search / RSS | Optimized FixedNeighbor search / RSS | Formal result |
|---|---:|---:|---:|---|
| chain-512 | 15.508 s / 306 MiB | 1161.205 s / 198 MiB | 169.145 s / 288 MiB | pass |
| cube-512 | 65.720 s / 631 MiB | 2762.299 s / 477 MiB | 323.885 s / 324 MiB | pass |
| chain-768 | 22.690 s / 392 MiB | formal timeout | 262.766 s / 283 MiB | pass |
| cube-768 | 135.293 s / 1376 MiB | formal timeout | 541.565 s / 515 MiB | pass |
| chain-1024 | 30.270 s / 581 MiB | formal timeout | 368.137 s / 310 MiB | pass |
| cube-1024 | 220.797 s / 2681 MiB | formal timeout | 769.116 s / 276 MiB | diagnostic pass |

OperatorPcg was fastest in all six groups. The optimized FixedNeighbor route used less RSS
than OperatorPcg in all six groups and materially improved the original FixedNeighbor search,
but OperatorPcg also completed every point inside the formal envelope. This is a resource
observation, not a demonstrated formal resource-survival case.

The complete reports, route ranking, promotion gate, and provenance are in the
[optimized frontier artifact](figures/joint-fixed-neighbor-optimized-frontier-r1/README.md).

## FixedNeighbor local-search experiment

The P6 frontier above compares the global `OperatorPcg` route with the outer
FixedNeighbor route. It must not be confused with the later benchmark-only
experiment that changed only the local width-search engine inside a
FixedNeighbor block. The outer algorithm remained serial Forward
Gauss-Seidel, with 64-atom cores, at most 30 sweeps, `OneAcceptedUpdate`,
global replay, and the existing convergence and endpoint checks.

The local candidates were `LegacyCompact` and `OperatorPcg` with Identity,
Diagonal, or Schwarz preconditioning. Production remained `LegacyCompact`.
The chain-256 comparison was:

| Local candidate | Search (s) | Local search (s) | RSS (MiB) | Sweeps / block solves | Profile evaluations | PCG mean / max | Assessment local | Runtime |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| LegacyCompact | 82.474 | 77.021 | 232.1 | 8 / 40 | 130 | — | `3.01e-13` | Passed |
| OperatorPcg + Identity | 24.952 | 20.169 | 207.1 | 7 / 35 | 109 | 6.04 / 7 | `1.458662526e-10` | **Failed** |
| OperatorPcg + Diagonal | 24.880 | 20.108 | 207.0 | 7 / 35 | 109 | 6.02 / 7 | `1.458662500e-10` | **Failed** |
| OperatorPcg + Schwarz | 20.277 | 15.505 | 207.8 | 7 / 35 | 109 | 2.31 / 4 | `1.458662500e-10` | **Failed** |

Identity and Diagonal operator setup took about `2.08 s` with no separate
preconditioner setup. Schwarz setup took `2.08 s` for the operator and
`1.42 s` for the preconditioner. All three OperatorPcg runs reached
`block-stationary`, global A/C KKT `1.01e-13`, width-gradient
`1.61e-13`, trusted endpoint state, and endpoint certification. They still
failed the existing local assessment threshold because
`assessment_local = 1.4586625e-10 > 1e-10`, so `RuntimeConvergence` failed.
Objective parity was about `1.0e-16` versus LegacyCompact, but objective
parity and global KKT do not replace endpoint certification.

Therefore the original, pre-geometry Decision Gate A was **failed**. The
OperatorPcg local route was not promoted, no automatic fallback or threshold
relaxation was added, and the 512/1024 local-operator matrix plus repeated
large-case promotion campaign were not run in that historical experiment.
The revised independent-local-Schwarz campaign and its bounded-polish
follow-up are recorded below; neither changes the production route.

Because Gate A failed, the exact LegacyCompact tile screen was tested next on
chain-512. The 8192-row control was the only full-endpoint run; the other
rows were search-only:

| Tile rows | Search (s) | Change vs 8192 | Derivative reduction (s) | Jacobian TiledQR (s) | Allocation/copy (s) | Householder (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1024 | 172.693 | +3.78% | 112.520 | 88.416 | 1.349 | 82.526 |
| 2048 | 169.900 | +2.10% | 109.311 | 85.112 | 1.553 | 79.228 |
| 4096 | 168.263 | +1.12% | 107.844 | 83.473 | 1.295 | 77.932 |
| 8192 | 166.401 | 0.00% | 106.244 | 82.110 | 1.296 | 76.667 |
| 16384 | 167.432 | +0.62% | 106.406 | 82.010 | 1.282 | 76.621 |

All five chain-512 searches converged, and the 8192 control passed endpoint
certification. Householder arithmetic was `93.1%`–`93.3%` of Jacobian
TiledQR, while allocation/copy was only `1.3%`–`1.6%`; no tile reached the
10% improvement gate. The tile direction therefore stopped without a
production QR/scratch change.

## Current supported route

The P6 promotion decision is complete: FixedNeighbor is a supported production
Joint search method, without a performance gate and without automatic routing.
The selected production policy is core 64, forward serial Gauss-Seidel order,
at most 30 sweeps, and one trusted accepted local `LegacyCompact` update per
block visit. Existing replay, KKT, width-gradient, eta-confirmation and
endpoint certification semantics remain unchanged.

- `OperatorPcg` is the general fast route and was fastest on the matched frontier.
- `LegacyCompact` remains the compatibility/reference route and the production default.
- `FixedNeighbor` is the explicitly selected bounded-memory block-coordinate route; its
  lower RSS trades against longer wall time on the measured frontier.

The benchmark-only local `OperatorPcg` experiment did not change this policy:
its correctness gate failed at the existing local assessment threshold, so the
production FixedNeighbor local solver remains `LegacyCompact`. The global
`OperatorPcg` bullet above refers to the separate global route, not to a
promoted local solver inside FixedNeighbor.

Atom-count auto-routing is not introduced. Coordinated or shared-parameter
blocks remain outside this policy, and `CertifiedLocal` remains diagnostic-only
and is not part of production semantics.

## Prepared-block hardening

The current production kernel separates the prepared component from the
historical diagnostic wrapper. An immutable prepared component owns the local
snapshot, parent/local mappings, identities and support domain. Each structural
FixedNeighbor block is prepared once with its core mapping, affected rows,
local domain and context/rank template. Sweeps update only the state-dependent
effective response, widths, coefficients and numerical basis values.

Each prepared block owns one persistent `LinearWorkspace`. This enables SPQR
symbolic reuse for a stable block pattern and scope, while numeric
factorization, active-set decisions and coefficient solves remain fresh for
every width state. The optimization does not change objective, effective
response, serial Gauss-Seidel order, accepted-step semantics, replay, or any
KKT, width-gradient, eta-confirmation or endpoint threshold.

The repeated workspace-residency campaign is now closed as Route C. It covered
chain/cube at 256, 512 and 1024 atoms with one warmup and three measured
repetitions per case; all numerical comparisons passed, while the persistent
wall-time gate was `not-material`. Search improvement ranged from `-0.32%` to
`+0.16%`, and persistent RSS was higher in every case. Exact numeric-factor
reuse opportunities were `1.17%`–`15.00%`, all at initial-profile evaluation,
but numeric factorization remains fresh. There is therefore no production
workspace-policy change, exact numeric-factor reuse, or candidate-replay
optimization. The large-case bottleneck is instead local derivative
preparation plus reduction: `291.99 s`/`83.7%` of chain-1024 search and
`600.93 s`/`81.5%` of cube-1024 search. See the [workspace attribution](joint-fixed-neighbor-workspace-attribution.md).

Production output keeps minimal aggregate search telemetry. Full local trial
trajectories, reverse order, alternate core sizes and local certification are
experiment/test diagnostics. The reproducible matched qualification driver is
`tests/integration/joint_fixed_neighbor_prepared_block.py`; its large cases are
search-only and its symbolic attribution is separate from any wall-time gate.

## Historical verification record

- `cmake --build build/joint-spqr -j4`: passed.
- P1 route, P3/P4/P5 policy, and P6 analyzer tests: 7/7 passed.
- `joint_benchmark_test.py`: 45/45 passed (the expected argparse diagnostic is printed by one
  negative test).
- `ctest --test-dir build/joint-spqr -L joint:runtime --output-on-failure`: 10/11 tests
  passed; the 213-test `rhbm_tests_joint_component` process reached the ctest 1500-second
  timeout without an assertion failure, while the remaining ten runtime tests passed.

This is retained as a historical regression record, not as the current runtime
status. No production default or numerical threshold was changed by the
performance experiments.

## Current implementation verification

The current production integration is covered by the runtime and persistence
checks documented in the main [joint runtime contract](joint-component-runtime.md).
The current `joint:runtime` lane is 11/11 green after the long frontier and
diagnostic suites were separated into `joint:scalability`.

## Revised local Schwarz geometry qualification (2026-10-08)

This campaign corrected the experimental confounder in the historical local
`OperatorPcg + Schwarz` result. The three sizes are now recorded separately:

- `outer_core_atoms = 64`: atoms permitted in the historical production
  control used by this benchmark-only local-Schwarz qualification.
- `local_atoms <= 64`: atoms in the particular local Operator problem; the
  chain/cube 256 workloads had local min/mean/max of `4/51.20/64` and
  `1/42.67/64`, respectively.
- `local_schwarz_core_atoms`: the explicit structural core used only to build
  the local Schwarz preconditioner.

The benchmark-only local policy is explicit and independent of the global
`SchwarzPolicy`. For these experiments its resource envelope was fixed at
`local_schwarz_max_block_atoms = 64`; the global default
`core_atoms=128, overlap_hops=1, max_block_atoms=512` was not changed. The
JSON records the requested policy, partition count, min/mean/max core and
overlap atoms, min/mean/max realized block atoms, and mean/max coverage ratio
without retaining per-action block matrices. The coverage ratio is
`realized_block_atoms / local_atoms`.

### 128 versus 64 sanity

The deterministic structural test compared the historical `128/1/512`
policy with `64/1/64` and `64/0/64` on a 64-atom local domain. The result was
**equivalent** in all three cases: one partition, one realized block of 64
atoms, zero realized overlap atoms, identical core/overlap assignments,
identical PCG iteration/state trajectory, and identical objective, eta/beta
and replay results. The tail blocks also realized the complete local domain,
not a 128-atom block. Thus the historical `Schwarz-128` result was a
full-local preconditioner result, not a qualified 128-atom local additive
Schwarz experiment.

### 256 search screen

All candidates converged in the search-only screen and passed the existing
search correctness checks. `assessment_local` and `RuntimeConvergence` were
intentionally `NotRun` here; endpoint assessment was done only after the
geometry screen. `partition_count` was 35 for chain-256 and 84 for cube-256
for every Schwarz candidate.

The following values are `chain/cube`; setup is aggregate preconditioner
setup, and PCG is `mean/max` iterations.

| Local candidate | Search seconds | Schwarz setup seconds | PCG mean/max |
| --- | ---: | ---: | ---: |
| Identity | 22.790 / 98.325 | — / — | 6.037/7 / 7.338/13 |
| Diagonal | 22.979 / 97.196 | — / — | 6.023/7 / 7.248/9 |
| Schwarz 128/1 | 19.285 / 59.682 | 1.422 / 2.894 | 2.312/4 / 2.193/8 |
| Schwarz 64/0 | 19.368 / 59.949 | 1.410 / 2.854 | 2.312/4 / 2.193/8 |
| Schwarz 64/1 | 19.370 / 60.286 | 1.446 / 2.922 | 2.312/4 / 2.193/8 |
| Schwarz 32/0 | 20.095 / 89.484 | 1.200 / 2.469 | 2.951/5 / 5.423/9 |
| Schwarz 32/1 | 19.899 / 82.124 | 1.284 / 5.051 | 2.764/4 / 4.342/8 |
| Schwarz 16/0 | 20.721 / 93.158 | 1.232 / 2.510 | 3.069/5 / 5.788/9 |
| Schwarz 16/1 | 20.014 / 87.130 | 1.317 / 6.369 | 2.852/5 / 4.642/8 |
| Schwarz 8/0 | 21.305 / 96.527 | 1.390 / 2.997 | 3.375/5 / 6.562/9 |
| Schwarz 8/1 | 21.314 / 91.426 | 1.616 / 9.163 | 3.229/5 / 5.184/8 |

| Local candidate | Realized block min/mean/max (chain) | Realized block min/mean/max (cube) | Mean coverage ratio (chain/cube) |
| --- | ---: | ---: | ---: |
| Identity / Diagonal | — | — | — / — |
| Schwarz 128/1 | 4/51.20/64 | 1/42.67/64 | 1.000 / 1.000 |
| Schwarz 64/0 | 4/51.20/64 | 1/42.67/64 | 1.000 / 1.000 |
| Schwarz 64/1 | 4/51.20/64 | 1/42.67/64 | 1.000 / 1.000 |
| Schwarz 32/0 | 4/25.60/32 | 1/21.33/32 | 0.500 / 0.500 |
| Schwarz 32/1 | 4/26.60/34 | 1/38.83/64 | 0.516 / 0.785 |
| Schwarz 16/0 | 4/14.22/16 | 1/11.13/16 | 0.278 / 0.261 |
| Schwarz 16/1 | 4/15.67/18 | 1/29.65/61 | 0.301 / 0.559 |
| Schwarz 8/0 | 4/7.53/8 | 1/5.22/8 | 0.147 / 0.122 |
| Schwarz 8/1 | 4/9.24/10 | 1/21.63/45 | 0.174 / 0.403 |

The overlap result is topology-dependent. On the chain, one-hop overlap is
small and remains close to the nominal sub-block geometry. On the cube, the
same overlap can inflate a nominal 32/16/8 core to mean realized blocks of
38.83/29.65/21.63 atoms, with maximum blocks of 64/61/45. This is why the
realized coverage ratio, not the requested core alone, is the selection
metric.

### Gate S and endpoint qualification

The full endpoint screen used chain-256 and cube-256 with the existing
assessment, trust, certification, replay and RuntimeConvergence gates. No
Schwarz candidate passed both topologies. `Identity` and `Diagonal` were
search controls and did not receive this endpoint screen.

| Local candidate | `assessment_local` chain/cube | Runtime chain/cube |
| --- | ---: | --- |
| Identity | NotRun / NotRun | NotRun / NotRun |
| Diagonal | NotRun / NotRun | NotRun / NotRun |
| Schwarz 128/1 | 1.458662500e-10 / 1.002286348e-09 | Failed / Failed |
| Schwarz 64/0 | 1.458662500e-10 / 1.002286348e-09 | Failed / Failed |
| Schwarz 64/1 | 1.458662500e-10 / 1.002286348e-09 | Failed / Failed |
| Schwarz 32/0 | 1.458662500e-10 / 1.002286291e-09 | Failed / Failed |
| Schwarz 32/1 | 1.458662496e-10 / 1.002286381e-09 | Failed / Failed |
| Schwarz 16/0 | 1.458662500e-10 / 1.002286377e-09 | Failed / Failed |
| Schwarz 16/1 | 1.458662500e-10 / 1.002286426e-09 | Failed / Failed |
| Schwarz 8/0 | 1.458662502e-10 / 1.002286293e-09 | Failed / Failed |
| Schwarz 8/1 | 1.458662500e-10 / 1.002286380e-09 | Failed / Failed |

The unchanged local correction threshold is `1e-10`; values such as
`1.4587e-10` remain failures. Therefore Gate S selected **no correctness-
qualified geometry**. `Schwarz 32/0/64` was retained only for attribution
because it had the smallest cube correction, nearly the smallest chain
correction, lower setup cost than 32/1, and a measured coverage near 0.5.
It is not a production selection.

### Correction attribution

With 32/0/64, the per-sweep diagnostic classified both chain-256 and
cube-256 as **plateaued**, not decreasing, oscillatory, rank-boundary-driven,
or active-face-driven. Rank stayed 256. The largest chain correction was
`1.4586625e-10` at atom 195; the largest cube correction was
`1.0022863e-9` at atom 73. The corresponding LegacyCompact controls reached
`3.0066468e-13` and `2.8718686e-11`, respectively. Operator and Legacy
trajectories therefore differ at a small number of coordinates, but the
existing endpoint gate still rejects the pure local Operator route.

### Local-work qualification

The selected geometry was held fixed at 32/0/64. Every local-work budget
failed the unchanged endpoint gate:

| Work policy | Chain search / PCG / correction | Cube search / PCG / correction | Runtime |
| --- | --- | --- | --- |
| OneAccepted | 20.006 s / 2.951/5 / 1.459e-10 | 89.328 s / 5.423/9 / 1.002e-9 | Failed / Failed |
| TwoAccepted | 24.193 s / 2.680/5 / 5.182e-10 | 104.280 s / 5.407/9 / 1.113e-9 | Failed / Failed |
| Full | 30.387 s / 2.507/5 / 5.028e-10 | 133.751 s / 5.454/9 / 1.067e-9 | Failed / Failed |

Gate A therefore failed for pure tuned OperatorPcg. No pure Operator local-
work policy was selected.

### Bounded Legacy polish

Because pure OperatorPcg failed, the bounded hybrid prototype was tested with
the same outer FixedNeighbor route and 32/0/64 local geometry. Confirmation
was recomputed after every solver switch; a one-polish run was not accepted
when its eta confirmation still exceeded `1e-10`.

| Policy | Topology | Operator s | Polish s | Total search s | `assessment_local` | Runtime |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Operator only | chain | 20.006 | 0.000 | 20.006 | 1.459e-10 | Failed |
| +1 Legacy polish | chain | 19.773 | 129.627 | 29.749 | 1.899e-11 | Failed (eta 1.451e-10) |
| +2 Legacy polish | chain | 19.954 | 138.822 | 39.968 | 1.899e-11 | Passed |
| Legacy control | chain | 0.000 | 81.484 | 81.484 | 3.007e-13 | Passed |
| Operator only | cube | 89.328 | 0.000 | 89.328 | 1.002e-9 | Failed |
| +1 Legacy polish | cube | 89.114 | 112.561 | 98.268 | 9.116e-11 | Failed (eta 9.111e-10) |
| +2 Legacy polish | cube | 89.092 | 122.256 | 107.589 | 2.325e-11 | Passed |
| Legacy control | cube | 0.000 | 149.227 | 149.227 | 2.872e-11 | Passed |

The strict hybrid analyzer selected `+2 Legacy polish` for the 256 cases;
`+1` was correctly rejected for confirmation even when its endpoint
assessment passed.

### 256/512 route gate and promotion decision

The selected route was then qualified with full endpoint checks on chain/cube
256/512. Three cases passed every required gate. Cube-512 passed endpoint
assessment, trust, certification and the RuntimeConvergence assessment, but
failed the separate eta confirmation/search-convergence gate:

| Case | Total search s | `assessment_local` | Eta change | Search convergence | Runtime |
| --- | ---: | ---: | ---: | --- | --- |
| chain-256 | 39.968 | 1.899e-11 | 9.99e-16 | Passed | Passed |
| cube-256 | 107.589 | 2.325e-11 | 6.82e-11 | Passed | Passed |
| chain-512 | 82.964 | 2.457e-11 | 2.47e-11 | Passed | Passed |
| cube-512 | 223.119 | 5.526e-11 | 1.547e-10 | **Failed** | Passed assessment, route rejected |

Since one 512 case failed the required confirmation gate, 1024 repeated
timing/frontier work was **not run** and Gate B is not passed. The production
decision is **DO NOT PROMOTE**. There is no `SearchMethod::Hybrid`, no
production local Operator policy, no public auto-routing, and no changed
threshold. At that historical checkpoint production remained outer
FixedNeighbor core 64 with local `LegacyCompact`, serial Forward Gauss-Seidel,
`OneAcceptedUpdate`, and the existing endpoint/replay semantics. The later
small outer-core qualification superseded that default with core 12; see the
[outer-core qualification](joint-fixed-neighbor-outer-core-qualification.md).

The campaign drivers are `tests/integration/joint_fixed_neighbor_local_schwarz.py`,
`joint_fixed_neighbor_operator_diagnostic.py`,
`joint_fixed_neighbor_local_schwarz_local_work.py`,
`joint_fixed_neighbor_hybrid.py`, and
`joint_fixed_neighbor_operator_route.py`. Their search-only outputs label
endpoint assessment and RuntimeConvergence as `NotRun`; search evidence must
not be described as full endpoint certification.
