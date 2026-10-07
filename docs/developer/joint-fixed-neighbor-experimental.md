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
