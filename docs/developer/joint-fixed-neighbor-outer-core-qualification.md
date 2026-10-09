# FixedNeighbor small outer-core qualification

Status: **qualified and adopted** on 2026-10-08. The production
`FixedNeighborSearchPolicy` default outer core is now **12 atoms**. FixedNeighbor
is the default Joint production route; during Checkpoint A, the remaining
top-level LegacyCompact route is only a temporary reference path.

## Scope and contract

The study was restricted to the production-compatible route:

- `FixedNeighbor` outer blocks
- `LegacyCompact` local search
- at most one trusted accepted local update per block visit
- forward serial Gauss-Seidel order
- SPQR and one Eigen thread
- maximum 30 sweeps

OperatorPcg, local Schwarz, hybrid routing, replay optimization and threshold
changes were excluded. The existing A/C KKT (`1e-10`), width-gradient
(`1e-12`), eta confirmation, endpoint, RuntimeConvergence and replay
semantics were preserved.

## Qualification gates

### C2: coarse and refinement sweeps

The coarse 512-atom sweep compared outer cores `8, 16, 32, 48, 64` on chain
and cube. Every run passed the existing search/replay correctness gate; core
16 was the fastest coarse finalist. The refinement sweep compared
`8, 12, 16, 24, 32`; core 12 was fastest on both topologies, while core 16
was close enough to remain a finalist.

### C3: full endpoint qualification

Finalists 12 and 16 were compared with the historical/control core 64 on
chain/cube × 256/512. All 12 runs completed with `block-stationary` search,
endpoint certification, `RuntimeConvergence`, inner/gradient/local/identified
assessment, endpoint trust and replay limits. Final beta/eta/scaled A/C
parity and projected/corrected/normalized rank evidence matched the control;
the endpoint campaign did not require bitwise-identical trajectories.

See the [C3 README](figures/joint-fixed-neighbor-core-size-r4/README.md),
[C3 analysis](figures/joint-fixed-neighbor-core-size-r4/analysis.json), and
the endpoint phase of the [outer-core qualification driver](../../tests/integration/joint_fixed_neighbor_outer_core_qualification.py).

### C4: repeated frontier

The selected finalists and control were measured with one warmup and three
interleaved measured repetitions on chain/cube-512. Core 12 won by aggregate
median search time, so only core 12 and control 64 were measured on
chain/cube-1024.

| Case | Core 12 median search | Core 64 median search | Improvement | Core 12 median RSS | Core 64 median RSS | Sweeps |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| chain-512 | 61.18 s | 171.26 s | 64.3% | 278.6 MiB | 302.2 MiB | 8 / 8 |
| cube-512 | 151.02 s | 322.81 s | 53.2% | 289.0 MiB | 624.5 MiB | 20 / 16 |
| chain-1024 | 167.73 s | 346.28 s | 51.6% | 668.8 MiB | 565.3 MiB | 8 / 8 |
| cube-1024 | 490.85 s | 732.64 s | 33.0% | 646.4 MiB | 806.4 MiB | 23 / 18 |

Every measured candidate and control run passed search convergence and replay
correctness. The C4 gate required at least 10% median search improvement on
both topologies, no material RSS/total-time/local-geometry regression, and
correctness on every measured run. No sweep warning was emitted: the largest
candidate sweep count was 23.

The complete matched records are in the [C4 README](figures/joint-fixed-neighbor-core-size-r5/README.md),
[C4 analysis](figures/joint-fixed-neighbor-core-size-r5/analysis.json), and
[C4 summary CSV](figures/joint-fixed-neighbor-core-size-r5/summary.csv). The
reproducible frontier phase is exposed by the
[outer-core qualification driver](../../tests/integration/joint_fixed_neighbor_outer_core_qualification.py).

## Production decision

The qualified default was adopted in
`src/core/detail/joint_component/FixedNeighborPolicy.hpp` and locked by the
`FixedNeighborUsesComponentRouteAndQualifiedDefaults` regression test. A
caller or benchmark may still construct an explicit 64-atom policy; core 64
is retained as the historical/control value in the qualification artifacts.

The adoption changed only the outer-core default. It did not change the
production route, local solver, block order, maximum sweeps, numerical
thresholds, endpoint semantics or persistence contracts.

The earlier campaign artifacts remain useful as historical controls:
[C2 coarse sweep](figures/joint-fixed-neighbor-core-size-r2/README.md),
[C2 refinement](figures/joint-fixed-neighbor-core-size-r3/README.md), and
[C3 endpoint qualification](figures/joint-fixed-neighbor-core-size-r4/README.md).
