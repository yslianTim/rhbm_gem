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

The validation ran the requested Release workflow against the supplied 0.50 Å
map. To preserve the user's SQLite database, it used a fresh SQLite backup and
the unique key `6z6u_gaus_charge1_bw050_kkt_refined_20261014`; the original
database was not written. The CLI exited 0 and printed `Completed`, but the
numerical search did not complete. Unavailable measurements are marked
explicitly; the selected observation domain and pruned additive model are not
directly comparable to the historical all-contributor objective.

| Measure | Result |
| --- | --- |
| Full 6Z6U run | Executed; CLI exit status 0. Numerical completion was not achieved. |
| Observation / target / candidate / retained / excluded counts | 229,783 / 1,418 / 35,988 / 626 / 35,362. |
| Excluded support / radial counts | 35,329 with 0 support rows; 20 singletons; 13 with 2 rows. Their q counts are 35,329 / 20 / 13 for q=0 / 1 / 2. All were excluded as `insufficient-observations`; no record was reason-coded only for insufficient radial information. |
| Components before / after pruning | Before: unavailable in this run's diagnostics. After: 1 component (2,044 atoms, 229,783 rows). |
| Search status and stop reason | `block-search-failed`, `search_completed=false`; component 1, one-based sweep index 2, block 82 stopped at `untrusted-trial` (`replay-disagreement`). It has 12 core atoms and 881 affected rows. |
| Local A/C problem | Design 881 × 24, 3,368 nonzeros, numerical rank 24, free rank 22; rank-deficient halo count 0. Active A atoms are 8 and 11. The initial local profile is now `qualified-inner`: feasible, zero releases, scaled projected KKT `5.82077e-11` (passes the unchanged `1e-10` gate). |
| Independent replay discrepancies | Prediction `1.28373e-08`; width-gradient infinity norm `6.29113e-05`; KKT `1.16415e-10`. Existing replay trust checks reject the trial; no replay or qualification threshold was changed. |
| Objective / final A/C KKT / log-B gradient | Final qualified metrics unavailable. Candidate-only diagnostics are `uncertified` and are not scientific estimates. |
| Endpoint verification / global ranks / target estimability | Not run or unavailable because no trusted component state exists. No selected-target recovery is available. |
| Runtime convergence / global assembly | Runtime convergence and target runtime convergence are unavailable; regular certificate is `not-run`; assembled state is null. |
| Construction / initialization / search / assessment / assembly | 1.969437 / 5.902256 / 6.276690 / 2.754835 / 0.000008 seconds. Callback time about 21.68 seconds; shell wall time about 21.9 seconds. Measured RSS is unavailable. |
| Local KKT failure / remaining blocker | The original local KKT failure is resolved. Full estimator qualification remains blocked by the independent replay disagreement, before endpoint verification or target recovery. |

The result was saved only in `/private/tmp/rhbm-gem-joint-local-ac-kkt-final-20261014.sqlite`. This run is not a numerical success despite the CLI exit status. The local A/C state passes its conditional KKT gate, but independent prediction, width-gradient and KKT comparisons reject its trial; the estimator therefore has no trusted state to verify, qualify or assemble. The diagnostic code records one-based sweep indexes; this capture reports sweep 2, while an earlier summary labeled the block Sweep 1. That prior label discrepancy remains. The failure is localized by component and block identity.
