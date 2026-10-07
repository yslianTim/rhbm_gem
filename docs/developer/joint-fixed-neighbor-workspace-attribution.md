# Fixed-neighbor workspace attribution

Status: **Gate B — stop after the measurement commits**. The measurements do
not justify localized candidate-replay changes. Production solver semantics,
ordering, tolerances and certification gates are unchanged.

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

The chain-256 row is a qualification smoke campaign with one warmup and three
matched measured pairs. The 1024 rows are single matched pairs run after the
smoke campaign; they are exploratory measurements, not a repeated statistical
campaign.

| Case | Control search (s) | Treatment search (s) | Search improvement | Treatment candidate replay + copy | Result |
| --- | ---: | ---: | ---: | ---: | --- |
| chain-256, 1 warmup + 3 measurements | 82.842 median | 83.213 median | -0.45% | 1.76% | not major |
| chain-1024, single pair | 358.609 | 349.155 | +2.64% | 5.51% | not major |
| cube-1024, single pair | 745.278 | 754.583 | -1.25% | 6.41% | regression in this pair |

The chain-256 qualification numerical gate passed for all three matched
comparisons. Its wall-time gate was `not-material`. For the 1024 pairs, both
variants converged with the same sweep count, reason, block trajectory,
objective, final eta, final beta, AC scaling weights, and certification values.

| Case | Sweeps | Search reason | Final global AC KKT | Final width-gradient norm |
| --- | ---: | --- | ---: | ---: |
| chain-1024, both variants | 8 | block-stationary | 7.54e-14 | 6.83e-16 |
| cube-1024, both variants | 18 | block-stationary | 4.04e-13 | 2.34e-15 |

## Attribution

The cube-1024 treatment is the larger measured search and provides the clearest
breakdown:

- local search: 670.47 s, 88.85% of search time;
- candidate replay: 48.37 s;
- candidate state copy: 0.016 s;
- candidate replay plus copy: 6.41% of search time;
- 360 block solves and 356 full candidate replays;
- 7,339,991 affected-row updates.

Persistent workspace reuse substantially reduced symbolic factorizations in the
two 1024 cases, while numeric factorizations stayed equal:

| Case | Variant | Symbolic factorizations | Symbolic reuses | Numeric factorizations | Peak RSS (MB) |
| --- | --- | ---: | ---: | ---: | ---: |
| chain-1024 | control / treatment | 226 / 107 | 46 / 165 | 272 / 272 | 327.7 / 472.6 |
| cube-1024 | control / treatment | 690 / 360 | 31 / 361 | 721 / 721 | 327.9 / 659.7 |

The lower symbolic count did not produce a material wall-time gain, and the
persistent treatment used more peak memory in these direct runs. The dominant
cost is therefore the local block search/profile path, not the full candidate
replay or vector copy.

## Decision and next direction

This is Gate B. Commits 4–7 (localized replay, partial state restore, and
related production-path changes) are intentionally not implemented. The
current evidence does not meet the roadmap's material replay threshold or show
a scaling-driven need for that risk.

If performance work resumes, the next measurement should attribute the local
`SearchProfile` path itself: local factorization work, profile evaluations,
local state evaluation, and the effect of workspace residency on memory. A
repeated 512/1024 campaign should be used before changing production behavior.

## Reproduction and validation boundaries

The matched campaign harness is:

```sh
PYTHONPATH=tests/integration python3 \
  tests/integration/joint_fixed_neighbor_workspace_timing.py \
  --build-dir build/qualification \
  --output-dir <output-directory>
```

The harness requires at least three measured repetitions. The 1024 direct
measurements above were run one control/treatment pair at a time because of
their runtime and must not be reported as a full campaign.

Endpoint certification in the timing harness remains `not-run-search-only`;
the KKT and width-gradient values above are the search-produced certification
fields, not an independent endpoint rerun.
