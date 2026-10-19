# Bounded Active-set LM Production Promotion — 2026-10-20

## Decision

The production Pure Halo-Support Exclusion Joint path now defaults to the
bounded Active-set LM candidate. Calibration uses
`Bmax = 4 * Q95(qualified FirstStage target widths)` and requires at least 20
qualified donors with at least 95% target coverage. A failed calibration is
persisted as Unavailable, skips downstream processing, and never retries the
unbounded solver. Fixed-Halo compatibility retains its separate Profile LM
path and does not inherit the Pure Exclusion calibration.

This promotion was carried out on `develop` at the user's explicit request.
The complete 37,406-target fitting RSS gate remains open: the exact-shape
preflight below performs zero search sweeps and zero Active-set LM subproblem
solves. This promotion record does not claim that a full-scale fit met the
4 GB target. Real 6Z6U post-promotion reruns were also unavailable because the
required model and map inputs were not present in the accessible workspace;
prior candidate runs are not presented as post-promotion results.

## Production route and persisted identity

- The ordinary initialized entry selects bounded Active-set LM for the
  `pure-halo-support-exclusion-v1` policy. Calibration failure returns an
  unavailable result; there is no silent unbounded fallback.
- The `first-stage-fixed-halo-v1` compatibility route remains unbounded and
  keeps its prior parameterization.
- Active-set LM retains boundary feasibility, projected Width KKT, the
  existing constrained A/C solve, endpoint certification, Independent Replay,
  global practical convergence limits, and downstream qualification rules.
  No correctness or replay tolerance was loosened.
- Solver provenance and CLI progress identify the bounded route as
  `active-set-lm`; Fixed-Halo compatibility continues to identify as
  `profile-lm`. JSON persistence accepts and reloads the new method value.
- No public solver or width-bound CLI option was added. No Notion content was
  changed.

## Verification

| Check | Result |
| --- | --- |
| Release build | Passed: `cmake --build build --parallel 4` |
| Full CTest | Passed: 45/45, `--parallel 4`, OpenMP 2, BLAS/VECLIB 1 |
| Production CLI bounded fit, persistence, reload and downstream | Passed in the 25-target integration fixture |
| Small selection | Calibration Unavailable is persisted; no Joint endpoint or downstream results are published |
| Failure lifecycle | Post-fit search failure and qualified-fit/no-donor outcome persist through API/CLI paths |
| Python bindings | Passed; stored solver provenance and calibration show the production bounded method |
| Bounded synthetic | Passed in the full test suite |
| 37,406-target resource preflight | Passed; zero sweeps and zero Active-set LM solves; not a full-fit measurement |
| Real 6Z6U simulation/unsharpened rerun | Not run: source model/map files unavailable in this workspace |
| `git diff --check` | Passed before commit |

The post-promotion exact-shape preflight used 5,938,087 observations, 37,406
targets/contributors, one connected component and 3,118 blocks. It reported:

| Measure | Bytes |
| --- | ---: |
| FixedNeighbor preparation peak estimate | 912,690,261 |
| Global state evaluation peak estimate | 975,208,628 |
| Endpoint verification peak estimate | 627,146,451 |
| Preflight budget | 3,200,000,000 |
| Baseline RSS | 24,985,600 |
| Measured process peak RSS | 1,721,008,128 |

RSS uses decimal bytes (`1 GB = 1,000,000,000 bytes`). The measured peak is for
the preflight test process only. FirstStage, calibration, search sweeps,
local A/C solves, Local Conditioning, conditional reprofiling, Parameter
Stability, assembly, persistence preparation, and downstream were not run at
this scale. The full-workflow RSS gate therefore remains open; see the
[giant-component resource validation](joint-bounded-giant-component-resource-validation-20261020.md).

## Earlier promotion stages

| Stage | Commit | Evidence |
| --- | --- | --- |
| P0 — Local criteria root-cause audit | `32cb6362a3d988d7fa77f715515351a9071d58dd` | [Root-cause audit](joint-local-stability-root-cause-audit-20261020.md) |
| P1 — Local stability practicalization | `ddfabe470f37886086c800c4e657b4b6395aea6a` | [Practicalization](joint-local-stability-practicalization-20261020.md) |
| P2 — Representative probe coverage | `08151497f68f6c49b250c97f94f4dbe562dc092f` | [Probe coverage](joint-local-stability-probe-coverage-20261020.md) |
| P3 — Production downstream workflow | `975ed2c72bfb35df695a387770f6ef6d796bba44` | [Downstream validation](joint-bounded-production-downstream-validation-20261020.md) |
| P4 — Giant-component resource validation | `d9fcbd9ec219ddab04d383cf24ce72b0950a49a7` | [Resource validation](joint-bounded-giant-component-resource-validation-20261020.md) |
| P5 — Production default promotion | recorded in this commit | This report |
