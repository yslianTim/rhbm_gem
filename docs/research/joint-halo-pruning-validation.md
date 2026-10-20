# Joint halo-pruning validation

> Historical validation: the threshold-sweep executable was retired after this
> study closed. Its results below remain the research record; current Halo
> eligibility and Pure Exclusion behavior are covered by
> `HaloEligibility_test.cpp` and `PureHaloSupportExclusion_test.cpp`.

## Production rule under study

Production uses `m_min = 3` and requires both at least three halo support rows
on the selected observation domain and at least three exactly distinct squared
distances. The distance comparison is exact; the gate does not implement a
floating-point tolerance or claim full-component identifiability. This is the
smallest necessary geometric gate for the three FullABC parameters.

The deterministic validation executable
`tests/experiments/joint_halo_pruning_validation.cpp` fits three known-truth
datasets at `m_min = 3, 4, 6, 8, 12`. It reports halo retention and exclusion
reasons, pre/post component counts, solver completion and qualification,
objective and residual RMS, target A/C/B errors, ranks, construction preflight,
elapsed time, and process peak RSS. It creates each observation vector once
from the full synthetic model, then keeps that observation domain fixed across
the sweep. Each fitted input uses the same policy validator as the production
constructor. Peak RSS is process-wide; it is not attributed to an individual
threshold. Construction preflight is not a total solver peak estimate.

The shared candidate set contains halos with 3, 4, 6, 8, and 12 support rows,
including a four-row halo with only two exact radii. The weak case adds a weak
two-row halo; the adversarial case replaces it with a strong two-row halo; the
overlap case adds three strong two-row halos whose supports overlap the target
and one another. These are focused deterministic examples, not a population
simulation or an estimate of real-map error rates.

## Threshold sweep results

Numerical success and scientific accuracy are separate. `numerical_success`
means the production solver passed its current runtime qualification for the
pruned model. Target deltas compare the fitted target A/C/B to the generating
target parameters. A qualified result can still have appreciable target error
because omitted halo signal remains in the observations.

The tables below are from the final Release executable. `block failures / no
state` counts failed block searches and components without an assembled state.
All three synthetic inputs have one component before pruning and one afterward.
The weak and strong cases each have 7 candidate halos (2 with q below 3); the
overlap case has 9 candidates (4 with q below 3). At each threshold, the weak
case excludes 1, 2, 4, 5, and 6 halos for insufficient rows and 1, 1, 0, 0,
and 0 for insufficient radial information, respectively. The strong case has
the same split. The overlap case excludes 3, 4, 6, 7, and 8 for insufficient
rows and 1, 1, 0, 0, and 0 for insufficient radial information.

### Weak omitted halo

| m_min | Retained / excluded | Retention | Block failures / no state | Search complete / numerical success | ΔA / ΔC / ΔB | RMS | Objective | Preflight bytes | Time (s) |
| ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |
| 3 | 5 / 2 | 0.7143 | 1 / 1 | no / no | NA / NA / NA | NA | NA | 28,916 | 0.008385 |
| 4 | 4 / 3 | 0.5714 | 0 / 0 | yes / yes | -0.0011355 / +0.0007074 / +0.0002577 | 0.001005 | 6.62518e-6 | 28,718 | 0.001791 |
| 6 | 3 / 4 | 0.4286 | 1 / 1 | no / no | NA / NA / NA | NA | NA | 28,152 | 0.001794 |
| 8 | 2 / 5 | 0.2857 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 27,938 | 0.064447 |
| 12 | 1 / 6 | 0.1429 | 0 / 0 | yes / yes | +0.0135723 / +0.0032302 / +0.0037456 | 0.005863 | 2.25526e-4 | 27,276 | 0.001579 |

### Strong low-support halo

| m_min | Retained / excluded | Retention | Block failures / no state | Search complete / numerical success | ΔA / ΔC / ΔB | RMS | Objective | Preflight bytes | Time (s) |
| ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |
| 3 | 5 / 2 | 0.7143 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 28,916 | 0.048046 |
| 4 | 4 / 3 | 0.5714 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 28,718 | 0.036111 |
| 6 | 3 / 4 | 0.4286 | 1 / 1 | no / no | NA / NA / NA | NA | NA | 28,152 | 0.002927 |
| 8 | 2 / 5 | 0.2857 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 27,938 | 0.088067 |
| 12 | 1 / 6 | 0.1429 | 1 / 1 | no / no | NA / NA / NA | NA | NA | 27,276 | 0.001969 |

### Multiple overlapping omitted halos

| m_min | Retained / excluded | Retention | Block failures / no state | Search complete / numerical success | ΔA / ΔC / ΔB | RMS | Objective | Preflight bytes | Time (s) |
| ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |
| 3 | 5 / 4 | 0.5556 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 29,154 | 0.047957 |
| 4 | 4 / 5 | 0.4444 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 29,084 | 0.068948 |
| 6 | 3 / 6 | 0.3333 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 28,518 | 0.028408 |
| 8 | 2 / 7 | 0.2222 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 28,048 | 0.028133 |
| 12 | 1 / 8 | 0.1111 | 0 / 1 | no / no | NA / NA / NA | NA | NA | 27,386 | 0.016997 |

The executable's process peak RSS was 23,543,808 bytes. The preflight values
above cover construction only, not total solver memory. Only the weak case
produced qualified target estimates: `m_min=4` had smaller target error and
residual RMS than `m_min=12`; the latter removed more contributors. The strong
and overlap cases produced no assembled state, so their target errors and
objective values are unavailable rather than zero. These examples show that
pruning can leave omitted signal that a retained model cannot absorb, and that
numerical solver success does not establish scientific accuracy. They do not
support increasing production `m_min` above 3 or establish 3 as an optimal
accuracy threshold. The production value remains the minimum geometric gate.

## Full-data validation: 6Z6U

The full workflow used source HEAD `b53c1d8bf90bf73b2b3920d4d6a4a94beafae5fc`
on `develop`, the supplied CIF and map, and a fresh SQLite backup at
`/private/tmp/rhbm-6z6u-width-stability-measure2-b53c1d8b.sqlite` with unique
key `6z6u_width_measure2_b53c1d8b`. The source database was opened read-only
for the backup and was not modified. The executable was an AppleClang 21.0.0
Release build on arm64 with `-O3 -DNDEBUG` and
`RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON`.

| Measure | Current run |
| --- | --- |
| Observations / selected targets | 229,783 / 1,418 |
| Halo candidates / retained / excluded | 35,988 / 626 / 35,362 |
| Exclusion reason | 35,362 had fewer than 3 support rows; 0 excluded only for radial information; production remains `m_min=3`, `q_min=3` |
| Connected components after pruning | 1 component, 2,044 contributors, 229,783 rows; pre-pruning component count unavailable |
| Search | 5,056 profile evaluations; 620 accepted updates across 3 completed sweeps; stopped in one-based sweep 4 / block 21; `search_completed=false` |
| First blocking failure | `inner-kkt-failed`; local A/C design 930 × 24, 2,894 nonzeros, rank 24, free rank 23, condition `6.6234727e10` |
| Local constrained A/C qualification | Feasible, no active-set releases, solver `solved`; KKT `4.656612873e-10` (Fail against unchanged `1e-10`) |
| Accepted-state width stationarity | Fail; gradient infinity norm `1.756e2` against `1e-12`; search stopped before stationarity |
| Endpoint verification / global Full Design Rank / projected-width rank | NotRun / NotRun / NotRun |
| Target estimability / full numerical qualification | NotRun / Unavailable |
| Global assembled state / selected-target A/C/B recovery | Unavailable / Unavailable |
| Estimated peak memory / measured child maximum RSS | Not recorded / 1,340,866,560 bytes (1.341 GB, 1.249 GiB) |
| CLI result / numerical result | Exit 0 and printed `Completed`; search incomplete, no qualified state |
| Wall time | 33.645 s wrapper; CLI callback `33.479 s` |

The bounded accepted-update trace recorded all 620 updates from the 2,044-atom
component. For reporting, “extreme width” here means an observed accepted
width of at least 10,000 Å; it is only a report filter, not a solver bound.
Each row gives first-stage width, maximum accepted width, and its one-based
sweep/block:

| Atom ID | Initial B (Å) | Maximum accepted B (Å) | Sweep / block |
| ---: | ---: | ---: | ---: |
| 11483 | 0.4690891689 | 44,291.5791 | 3 / 19 |
| 11721 | 0.5107510641 | 85,164.2930 | 3 / 21 |
| 11724 | 0.5268537435 | 113,195.6873 | 3 / 25 |
| 16409 | 0.4576126073 | 369,272.3596 | 3 / 58 |
| 16610 | 0.4422681954 | 18,768.1336 | 3 / 63 |
| 17007 | 0.5425051850 | 160,916.1571 | 3 / 70 |
| 17247 | 0.4715698279 | 239,833.7129 | 2 / 78 |
| 55975 | 0.5015657953 | 23,358.7370 | 3 / 82 |
| 63256 | 0.4518990389 | 20,102.6638 | 3 / 177 |
| 70821 | 0.4208963272 | 68,150.5152 | 3 / 184 |
| 73371 | 0.4159698614 | 406,905.5114 | 3 / 120 |

Across the accepted trace, the largest observed local condition indicator was
about `1.0503e12` at atom 17247's block 78. The largest normalized A/C
cancellation ratio, defined as
\(\frac{\max(|A|\|\mathbf g\|_2,|C|\|\mathbf k\|_2)}{\left|A\|\mathbf g\|_2+C\|\mathbf k\|_2\right|}\),
was `8.3764e10` at atom 17247, sweep 2 / block 78, with normalized terms
`+1.15177471e10` and `-1.15177471e10`. The largest individual normalized
coefficient term was `1.55850578e10` at the same atom in sweep 3 / block 78.
For atom 11721, the accepted sweep 3 step moved B from `0.4942385633` to
`85,164.2929583` Å; at that accepted point A/C cancellation was not yet
present. At the later failing sweep 4 / block 21 state, the terms were
`+9.23067447e8` and `-9.23067446e8`. Atom 11650 moved from `0.4101723088` to
`0.1232364573` Å and its observed Gaussian norm was `2.10279e-38`.

The run logs are
`/private/tmp/rhbm-6z6u-width-stability-measure2-b53c1d8b.log` and the
failure fixture is `tests/fixtures/joint_local_ac_kkt/block-21-sweep-4.json`.
The failure occurs before replay trust for that local profile, endpoint
verification, global rank qualification, target estimability, or trusted
assembly. The CLI exit code is not a numerical success signal.
