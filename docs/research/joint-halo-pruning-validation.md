# Joint halo-pruning validation

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

The final production run used the supplied CIF and map with a fresh SQLite
backup at `/private/tmp/rhbm-6z6u-corrected-14d8c2cd9b.sqlite` and unique key
`6z6u_replay_final_14d8c2cd9b`. The source database was not modified. The
production source was at `c7e6dd5d` (`fix(joint): keep multiprecision out of
production`); the later test-only block 21 fixture commit did not change solver
code. The executable was an AppleClang 21.0.0 Release build on arm64.

| Measure | Final result |
| --- | --- |
| Observation / selected target count | 229,783 / 1,418 |
| Halo candidates / retained / excluded | 35,988 / 626 / 35,362 |
| Excluded support / radial counts | 35,362 had fewer than 3 rows; 0 were excluded only for radial information |
| Components after pruning | 1 component (2,044 atoms, 229,783 rows); pre-pruning component count unavailable |
| Original sweep 2 / block 82 trial | Current replay passes: prediction difference 0, width-gradient difference `7.5871e-15`, KKT difference 0; primary scaled KKT `5.7795e-12` |
| Final failure | Component 1, one-based sweep 4, block 21; `inner-kkt-failed` before replay trust; no rejected LM trial capture |
| Failure local design | 930 × 24, 2,894 nonzeros, rank 24, free rank 23, condition indicator `6.6235e10` |
| Failure local A/C certificate | Feasible; 0 active-set releases; solver reason `solved`; scaled projected KKT `4.656612873e-10`, above unchanged `1e-10` |
| Search progress | 5,056 evaluations; 620 accepted updates over 3 completed sweeps; sweep 4 stopped at block 21; `search_completed=false` |
| Accepted-state diagnostics | Objective `2.567e-3`, width-gradient infinity norm `1.756e2`; unqualified and not a scientific result |
| Final local profile objective | `2.0159238843`; local state did not pass KKT |
| Global A/C KKT / log-B gradient | Not qualified / unavailable |
| Endpoint verification / Full Design Rank / free-design rank / projected-width rank | Not run |
| Target estimability / numerical qualification / trusted assembly | Not reached / failed to establish / unavailable |
| Process result | CLI exit 0 and printed `Completed`; numerical search incomplete |
| Stage times: construction / initialization / search / assessment / assembly | 1.970 / 5.891 / 19.126 / 0.918 / 0.000010 s |
| Estimated peak memory / measured child maximum RSS | Not recorded / 1,447,378,944 bytes (1.447 GB, 1.348 GiB) |
| Elapsed wall time | 33.442 s |

The run log is `/private/tmp/rhbm-6z6u-corrected-14d8c2cd9b.log`; the
failure-only local fixture was committed at
`tests/fixtures/joint_local_ac_kkt/block-21-sweep-4.json`. The 6Z6U failure
moved beyond the original replay disagreement. It is now blocked by an
ill-conditioned but full-rank local A/C problem that does not meet the existing
KKT gate. There is no endpoint certification, qualified numerical state, or
selected-target A/C/B estimate. The capture records `sweep_index=4` directly.
The earlier baseline run stopped at sweep 2 / block 82; the block 82 fixture
now replays as trusted, while the final full-data run stops later at block 21.
These are distinct decisions and must not be conflated.
