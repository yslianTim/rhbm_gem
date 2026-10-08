# Fixed-neighbor workspace attribution

Status: **Route C — study closed**. The repeated campaign does not justify a
production workspace-policy change, exact numeric-factor reuse, or candidate
replay optimization. Production solver semantics, ordering, tolerances and
certification gates are unchanged.

## Scope

This measurement compares the existing persistent per-block workspace with a
matched benchmark control that creates a fresh `LinearWorkspace` for every
block visit.

- topology: chain and cube
- core size: 64 atoms
- block order: Forward
- local policy: `OneAcceptedUpdate`
- maximum sweeps: 30
- measurement scope: fixed-neighbor search only
- treatment: `reuse_block_workspace=true`
- control: `reuse_block_workspace=false`

The aggregate work counters and timers are opt-in telemetry. The production
component route does not enable telemetry, and the default workspace policy
remains the treatment.

## Results

The canonical repeated campaign covers six cases (chain/cube × 256/512/1024)
with one warmup and three matched measured repetitions per case. The numerical
gate passed for every measured pair. The persistent-workspace wall-time gate is
`not-material`.

| Case | Fresh search (s) | Persistent search (s) | Improvement | Fresh / persistent RSS (MB) | Exact numeric reuse |
| --- | ---: | ---: | ---: | ---: | ---: |
| chain-256 | 82.167 | 82.434 | -0.32% | 146.1 / 177.8 | 15.00% |
| chain-512 | 168.892 | 168.630 | +0.16% | 263.9 / 355.7 | 13.19% |
| chain-1024 | 348.867 | 348.998 | -0.04% | 491.2 / 568.0 | 12.50% |
| cube-256 | 150.360 | 150.565 | -0.14% | 295.4 / 365.0 | 12.24% |
| cube-512 | 325.039 | 324.943 | +0.03% | 543.7 / 632.4 | 1.17% |
| cube-1024 | 735.772 | 737.640 | -0.25% | 645.4 / 747.9 | 7.07% |

Exact reuse opportunities occur only during initial-profile evaluation; trial,
reference and accepted-endpoint opportunities are zero. The treatment still
performs fresh numeric factorization: the reuse numbers are a census of
matching inputs, not a production reuse policy.

## Attribution

The repeated campaign identifies local derivative preparation plus reduction as
the dominant path in the production `LegacyCompact` local search:

| Case | Search (s) | Derivative preparation (s) | Derivative reduction (s) | Combined share |
| --- | ---: | ---: | ---: | ---: |
| chain-256 | 82.434 | 18.188 | 53.168 | 86.3% |
| chain-512 | 168.630 | 36.960 | 107.552 | 85.7% |
| chain-1024 | 348.998 | 74.541 | 217.447 | 83.7% |
| cube-256 | 150.565 | 38.593 | 88.802 | 84.6% |
| cube-512 | 324.943 | 83.451 | 193.199 | 85.1% |
| cube-1024 | 737.640 | 181.767 | 419.167 | 81.5% |

For the two large cases, chain-1024 spends 291.99 s (83.7% of search) and
cube-1024 spends 600.93 s (81.5%) in these two derivative stages. This is the
new optimization target; candidate replay and vector copy are not material.

Persistent workspace reuse substantially reduced symbolic factorizations in the
two 1024 cases, while numeric factorizations stayed equal:

| Case | Variant | Symbolic factorizations | Symbolic reuses | Numeric factorizations | Peak RSS (MB) |
| --- | --- | ---: | ---: | ---: | ---: |
| chain-1024 | control / treatment | 226 / 107 | 46 / 165 | 272 / 272 | 327.7 / 472.6 |
| cube-1024 | control / treatment | 690 / 360 | 31 / 361 | 721 / 721 | 327.9 / 659.7 |

The lower symbolic count did not produce a material wall-time gain, and the
persistent treatment used more peak memory in these repeated runs. The
dominant cost is therefore the local derivative path, not the full candidate
replay, vector copy, or numeric factorization.

## Decision and next direction

Route C is closed with no production workspace-policy change, no exact
numeric-factor reuse, and no candidate replay optimization. Performance work
then moved to FixedNeighbor local `SearchProfile` derivative preparation and
reduction. The search-specific attribution and local-solver decision are now
closed: the local OperatorPcg candidate failed the existing correctness gate,
and the chain-512 tile screen found no 10% LegacyCompact improvement. The
production local route remains `LegacyCompact`; see the [qualification
evidence](joint-fixed-neighbor-experimental.md) for the measured gate failure
and tile results.

## Reproduction and validation boundaries

The matched timing harness was a closed campaign driver and is retired from the
current tree. The compact analysis above is the retained evidence; the
historical implementation remains recoverable from Git history if reproduction
is needed.

The canonical campaign uses one warmup and three measured repetitions for each
of the six cases above. Its machine-readable results are in
`figures/joint-fixed-neighbor-workspace-residency-r2/analysis.json`.

Endpoint certification in the timing harness remains `not-run-search-only`;
the KKT and width-gradient values above are the search-produced certification
fields, not an independent endpoint rerun.
