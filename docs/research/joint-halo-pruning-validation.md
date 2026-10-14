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

The validation uses the requested local command and writes its saved result to
the supplied SQLite test database. Report fields come from the command output
and saved schema-11 result. Unavailable measurements are marked explicitly;
the selected observation domain and pruned additive model are not directly
comparable to the historical all-contributor objective.

| Measure | Result |
| --- | --- |
| Full 6Z6U run | Completed with the requested Release command; exit status 0. |
| Observation / target / candidate / retained / excluded counts | 229,783 / 1,418 / 35,988 / 626 / 35,362. |
| Excluded support / radial counts | 35,329 with 0 support rows; 20 singletons; 13 with 2 rows. Their q counts are 35,329 / 20 / 13 for q=0 / 1 / 2. All were excluded as `insufficient-observations`; no record was reason-coded only for insufficient radial information. |
| Components before / after pruning | Before: unavailable in this run's diagnostics. After: 1 component (2,044 atoms, 229,783 rows). |
| Updated search status | `block-search-failed`, `search_completed=false`; the assembled state is unavailable. The failure occurred at sweep 1, block 82, with 12 core atoms, 881 affected rows, and 626 retained halos. |
| Failed block ranks | Local design: 24 rows, 24 columns, rank 24; rank-deficient halo count 0. The inner KKT check failed (`kkt-failed`). |
| Objective / A-C KKT / log-B gradient | Final qualified metrics unavailable. Provisional diagnostics showed objective 0.01804, KKT 0.01368, and gradient infinity norm 0.0007427; these do not qualify the result. |
| Endpoint verification / global design ranks / target estimability | Unavailable; no assembled state or global rank evidence was produced. |
| Runtime convergence / global assembly | Runtime convergence failed. Global assembly completed as untrusted; regular certificate was `not-run`. |
| Estimated peak / measured RSS / elapsed time | Estimated total peak unavailable; measured peak RSS 1,227,636,736 bytes (~1.14 GiB); 23.67 seconds real time (23.102 seconds callback time). |
| Construction / initialization / search / assessment | 2.085 / 6.064 / 6.486 / 2.818 seconds; assembly 0.000014 seconds. |
| Rank failure resolved / remaining failures | The failed block contained no rank-deficient halo, so the removed low-support halo rank cause was absent there. End-to-end success was not achieved: the inner KKT search failed, leaving global rank and target qualification unavailable. |

The existing database row under this key was snapshotted before the run. It
held a schema-10 result built from a 0.30 Å map, while the requested validation
used the supplied 0.50 Å map. The prior row is therefore not a same-domain
baseline, and its component count or failure details are not treated as the
pre-pruning state for this run. The requested run overwrote that key after the
snapshot was saved under `/private/tmp/joint_6z6u_database_before.sqlite`.

For the specified map, the pruning diagnostics show that the contributor
population was sharply reduced before problem construction. The observed
search failure remains unresolved; the available diagnostics isolate an inner
KKT failure in a locally full-rank block, while final global convergence and
target recovery cannot be assessed without an assembled qualified state.
