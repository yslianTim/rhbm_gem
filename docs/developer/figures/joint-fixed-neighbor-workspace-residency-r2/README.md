# Matched FixedNeighbor workspace timing

Control uses a fresh block-local `LinearWorkspace` on every block visit. Treatment uses the persistent prepared-block workspace and SPQR symbolic reuse. Both are search-only and use one warmup run plus three measured runs for each topology/size pair.

The matched policy is `core_atoms=64`, 30 maximum sweeps, forward block order, `OneAcceptedUpdate`, `LegacyCompact`, SPQR, and one Eigen thread. The r2 phase timers are non-overlapping: basis and linear-certificate timings are measured explicitly around their calls.

Numerical qualification: **passed**. Wall-time gate: **not-material**.

| Case | Search control / treatment (s) | Persistent RSS increase | Exact numeric reuse |
| --- | ---: | ---: | ---: |
| chain-256 | 82.167 / 82.434 | 21.6% | 15.00% |
| chain-512 | 168.892 / 168.630 | 34.8% | 13.19% |
| chain-1024 | 348.867 / 348.998 | 15.6% | 12.50% |
| cube-256 | 150.360 / 150.565 | 23.6% | 12.24% |
| cube-512 | 325.039 / 324.943 | 16.3% | 1.17% |
| cube-1024 | 735.772 / 737.640 | 15.9% | 7.07% |

The exact reuse opportunities are all initial-profile opportunities; trial, reference, and accepted-endpoint opportunities are zero. The treatment performs no numeric-factor skip: the census only measures whether all factor inputs match.

Gate 1 decision: **Route C**. Persistent residency has no material wall-time benefit and its RSS increase is not uniformly above the strict 25% threshold on the large cases. Route B is also not justified: exact reuse is low and linear numeric factorization is small relative to derivative preparation/reduction (the dominant phase). No production workspace policy change is made.
