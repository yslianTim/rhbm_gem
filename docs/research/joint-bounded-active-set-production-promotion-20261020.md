# Bounded Active-set LM Production Promotion — 2026-10-20

> Historical promotion report. FixedHalo estimation was retired in a later
> change; the current repository keeps only legacy result compatibility.

## Decision

The production Pure Halo-Support Exclusion Joint path now defaults to the
bounded Active-set LM candidate. Calibration uses
`Bmax = 4 * Q95(qualified FirstStage target widths)` and requires at least 20
qualified donors with at least 95% target coverage. A failed calibration is
persisted as Unavailable, skips downstream processing, and never retries the
unbounded solver. At the time of this promotion, the Fixed-Halo Profile LM
compatibility path still existed; it was retired afterward.

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
- The `first-stage-fixed-halo-v1` compatibility route existed at the time of
  this report and kept its prior parameterization; it has since been removed.
- Active-set LM retains boundary feasibility, projected Width KKT, the
  existing constrained A/C solve, endpoint certification, Independent Replay,
  global practical convergence limits, and downstream qualification rules.
  No correctness or replay tolerance was loosened.
- Solver provenance and CLI progress identify the bounded route as
  `active-set-lm`; JSON persistence accepts and reloads the new method value.
  Historical `profile-lm` solver metadata remains preserved by the decoder.
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

## Follow-up Release reproduction after research cleanup

The 6Z6U simulation inputs became available for a later Release run. The
no-instrumentation production reproducibility tool completed on chain A with
1,418 targets and 180,478 observations. The model, map and simulation-truth
SHA-256 values match the saved B4 record. This follow-up does not change the
historical availability statement above or close the giant-component RSS gate.

| Measure | Saved B4 | Follow-up Release run |
| --- | ---: | ---: |
| Qualified FirstStage donors | 1,415 / 1,418 | 1,415 / 1,418 |
| Global P95 (Å) | 0.5689788338510722 | 0.5689788338510722 |
| Bmax = 4×P95 (Å) | 2.2759153354042887 | 2.2759153354042887 |
| Candidate objective | 1.2330854370810358e-13 | 1.2330854370810358e-13 |
| Prediction RMSE | 9.206498889018314e-7 | 9.206498889018622e-7 |
| A/C KKT | 9.434574186999167e-8 | 9.434574186999167e-8 |
| Projected Width KKT | 8.776510069988752e-10 | 8.776509519137221e-10 |
| Search sweeps / profile evaluations | 83 / 22,908 | 83 / 22,908 |
| Final active-bound atoms | 0 | 0 |
| Independent Replay / RuntimeConvergence | Passed / Passed | Passed / Passed |
| Peak RSS | 586,285,056 B | 616,595,456 B |
| Wall time | 297.87 s | 305.00 s |

The reported B RMSE also matches within floating-point variation. The saved
B4 A/C truth RMSE values (9.50818e-6 / 2.92104e-6) differ from this run's
4.18660e-6 / 1.28618e-6 by the map normalization divisor: multiplying the
saved values by 0.44031545732594374 reproduces the follow-up values within
floating-point variation. The saved record labels its A/C truth scale as
original-domain, but does not retain per-target final A/C estimates; aggregate
RMSE parity is established after this scale conversion, while per-target
endpoint equality is not claimed. The follow-up reports A/B/C RMSE of
4.18660e-6 / 1.46475e-7 Å / 1.28618e-6, with truth coverage 1,418/1,418 and
`truth_used_for_calibration=false`.

The Release run used `joint_width_repro_no_instrumentation` with
`instrumentation_enabled=false`. Its machine-readable summary and fingerprints
were written to `/private/tmp/joint-post-promotion-cleanup-parity-20261020`;
they are not committed. The 37,406-target / 5,938,087-observation complete-fit
4 GB RSS gate remains open; this 1,418-target run is not that measurement.

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
