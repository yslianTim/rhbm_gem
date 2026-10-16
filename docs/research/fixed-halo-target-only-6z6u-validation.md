# Fixed-Halo / Target-Only 6Z6U Validation

Date: 2026-10-15  
Production input: `/Users/yslian/data/6z6u.cif` and
`/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map`  
Simulation truth: the map's matching `.simulation.json` sidecar.

## Formal production run

The Release CLI ran `potential_analysis --estimator joint-components
--select-chain A` and exited 0 after saving its result. The normal CLI map
normalization was applied with divisor `0.4403154573`. The selected observation
domain contained 229,783 rows and 1,418 targets. Contributor selection found
35,988 non-target candidates; it retained all 659 with nonempty target-domain
support and excluded 35,329 zero-support atoms. The retained set includes 20
one-row halos and 13 two-row halos.

The run stopped during FirstStage validation, before Joint optimization:

| Measure | Formal result |
| --- | ---: |
| Fixed halos retained | 659 |
| FirstStage-qualified fixed halos | 655 / 659 |
| Target components fitted | 0 |
| Largest fitted component | Unavailable; solver did not start |
| Target-only blocks / profile evaluations / accepted updates / sweeps | 0 |
| Converged components / available target estimates | 0 / 0 |
| Objective, target A/C KKT, width gradient | Unavailable |
| Independent endpoint replay | Not run |
| CLI wall time / process peak RSS | 12.261 s / 1,098,203,136 bytes |

The result records
`fixed-halo-initialization-unavailable:first-stage-fitting-failure`. These
four halos have finite candidate A/B values, but their native FirstStage status
is `MAX_ITERATIONS_REACHED`; the estimator correctly refuses to use them as
fixed parameters.
No fixed-halo snapshot was accepted; the exported CSV labels targets
`TargetOnlyUnestimated` and supported halos `FixedHaloUnavailable`.

| Atom ID | Joint support rows | Distinct radii | FirstStage samples | CLI A | CLI B |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 16664 | 73 | 73 | 63 | 14.37693 | 0.478284 |
| 57190 | 1 | 1 | 140 | 16.62917 | 0.445106 |
| 62463 | 81 | 81 | 99 | 12.37041 | 0.508381 |
| 63216 | 50 | 50 | 71 | 7.88796 | 0.416487 |

CLI A values use its normalized map units. The raw-map research study below
uses raw-map units and preserves the original observation scale across each
controlled fit. No oracle values or failed FirstStage values enter production.

## FirstStage halo error and low-support checks

The truth comparison is research-only. It uses finite FirstStage A/B values
even for the four failed records, so these summary errors describe candidate
parameters rather than a production-eligible fixed-halo field.

| Measure | Raw-map candidate result |
| --- | ---: |
| FirstStage halo A RMSE | 1.63277 |
| FirstStage halo B RMSE | 0.0315442 |
| Fixed-halo prediction RMSE | 0.0783596 |
| Fixed-halo prediction L2 error | 37.5622 |
| Truth halo prediction L2 norm | 77.4963 |
| Relative fixed-halo prediction error, `||H_first - H_truth||₂ / max(1, ||H_truth||₂)` | 0.484697 |
| Maximum absolute local halo prediction error | 1.18324 |

| Joint support rows | Halo count | FirstStage qualified / failed | FirstStage sample range | Per-row halo-field error RMSE |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0 | 0 / 0 | — | — |
| 1 | 20 | 19 / 1 | 40–200 | 0.00178955 |
| 2 | 13 | 13 / 0 | 34–200 | 0.00191819 |
| 3+ | 626 | 623 / 3 | 20–200 | 0.0783305 |

FirstStage sample counts are independent of Joint-domain support counts: the
single-row halo 57190, for example, used 140 FirstStage samples. All 33 halos
with one or two Joint support rows are now retained even though their number
of distinct radii cannot meet the retired three-radius FullABC gate. Four
FirstStage failures are rejected regardless of their support count.

## Controlled target-only candidate fits

These are research-only endpoint candidates, not qualified estimates. To
measure controlled scenarios despite the four invalid production records, the
study injects finite candidate values and marks the resulting scenarios as
research-only. The oracle cases use simulation truth only in this experiment.
All scenarios use the same 229,783 raw observation rows and original scale
`s = 346.78395`.

| Scenario | Halo settings | Target A RMSE | Target B RMSE | Target C RMSE | Target prediction RMSE | Full prediction RMSE | Objective |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| FirstStage candidate | FirstStage A/B, C=0; all 659 | 0.516450 | 0.0121137 | 0.0712743 | 0.0447331 | 0.0667794 | 0.00426045 |
| Oracle fixed halo | Truth A/B/C; all 659 | 0.116850 | 0.00402451 | 0.0326313 | 0.0125070 | 0.0125070 | 0.000149443 |
| FirstStage A/B, oracle charge | FirstStage A/B, truth C; all 659 | 0.134749 | 0.00434847 | 0.0349543 | 0.0143068 | 0.0302132 | 0.000872093 |
| Reduced legacy set | FirstStage A/B, C=0; 626 legacy-eligible halos | 0.516449 | 0.0121137 | 0.0712742 | 0.0447331 | 0.0667794 | 0.00426045 |

Each case formed one target-only component of 1,418 targets and used 142
blocks per sweep, 2,272 profile evaluations, 1,136 accepted updates, and the
eight-sweep research cap. None converged. At the final candidate endpoint,
target A/C KKT ranged from `1.6056e-3` to `1.6319e-3` (required at most
`1e-10`), and target width-gradient infinity norm ranged from `8.9172e-6` to
`9.1469e-6` (required at most `1e-12`). Independent replay, target conditioning,
and target parameter stability were unavailable because no candidate passed
convergence. Per-case time was 33.40–33.68 s; process peak RSS across the study
was 809,484,288 bytes.

The reduced control keeps the old support/radial eligibility set of 626 halos;
the other 33 supported halos have zero amplitude in that research scenario.
Three of the 626 have failed FirstStage statuses, so their finite values are
also candidates used only for this control. The close FirstStage and reduced
candidate metrics suggest that the 33 low-support contributions have little
effect at this unfinished endpoint; they do not establish that conclusion for
a converged fit. Oracle values substantially improve candidate target errors,
which indicates sensitivity to fixed-halo model error. Giving halos their
truth charge while retaining FirstStage A/B also improves the candidate, so
the forced `C=0` model can matter in this simulation. These comparisons are
directional only because all four target fits stopped before the unchanged
qualification gates.

## Comparison with the previous FullABC run

The prior run documented in
[`joint-halo-pruning-validation.md`](joint-halo-pruning-validation.md) retained
626 free halos and 1,418 targets in one 2,044-contributor component. It stopped
after 3 completed sweeps, with 5,056 profile evaluations and 620 accepted
updates, then failed at sweep 4 / block 21. Its local A/C KKT was
`4.6566e-10` against the unchanged `1e-10` gate; accepted widths exceeded
400,000 Å. Wall time was 33.645 s and process peak RSS was 1,340,866,560 bytes.

The new formal run never entered the target optimizer, so it cannot show that
width runaway or the old block-21 failure is fixed. The old failing region was
reported around atom 11721; none of the four new FirstStage failures is that
atom. Block numbers are not comparable across the target-only partition. The
eight-sweep research candidates are incomplete and do not record a trusted
width trajectory, so width-runaway status remains unavailable. Their resource
numbers are not a completed-fit comparison with the old FullABC run.

## Large synthetic resource regression

`ProductionScaleSingleComponentPreparationStaysBelowBudget` passed with the
requested 5,938,087 observations and 37,406 contributors:

| Measure | Result |
| --- | ---: |
| Targets / fixed halos | 27,406 / 10,000 |
| Target-only components / blocks | 1 / 2,284 |
| Target A/C design columns / total A,B,C target parameters | 54,812 / 82,218 |
| Total / target support entries | 74,812 / 54,812 |
| Fixed-halo prediction vector | 47,504,696 bytes |
| Effective observation vector | 47,504,696 bytes |
| Construction preflight estimate / budget | 969,939,078 / 3,200,000,000 bytes (Pass) |
| Endpoint verification peak estimate | 624,026,451 bytes |
| Observed process peak RSS | 1,585,790,976 bytes |

The test also verified that the target parameter vector has three scalars per
target, that component preparation includes targets only, and that no global
dense matrix is allocated.

## Implementation and validation conclusion

The production builder now selects every non-hydrogen halo with at least one
target-domain support row. The production estimator validates and snapshots
FirstStage A/B with C=0, then optimizes only target A/B/C. Target-only
components, block cores, KKT, width gradients, conditioning, and stability use
target parameters; full replay and post-fit peeling retain fixed-halo signal.
Result export now writes fixed halo A/B/C values as `FixedFromFirstStage` in
the CSV parameterization column; target rows distinguish `JointEstimated`
from `TargetOnlyUnestimated`.
The parameterization and contributor-policy contracts distinguish this result
from historical FullABC output. Numerical gates remain KKT `1e-10`, width
gradient `1e-12`, and complete-sweep log-width change `1e-10`.

This replaces the Free-Halo FullABC route in the production workflow. The
6Z6U formal fit did **not** converge and produced no usable target estimates:
initialization correctly stopped on four FirstStage maximum-iteration
failures. The smallest next step is to diagnose and resolve those four local
FirstStage failures under FirstStage's existing success criteria, then rerun
the formal integration study. Their finite endpoint values must not be
silently accepted as fixed parameters.

## Final verification

The final `Release` build completed with `cmake --build build --parallel 4`.
The complete CTest suite passed 40/40 in 75.20 s with parallelism 4. The final
6Z6U CLI and result export both exited 0; the large synthetic resource
regression passed separately.

## P0 follow-up: failed FirstStage endpoint diagnostics

The four failures were reproduced with the formal chain-A input and the CLI's
map normalization (`MapValueArrayNormalization()` before sampling). Source
SHA-256 values:

| Input | SHA-256 |
| --- | --- |
| `/Users/yslian/data/6z6u.cif` | `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa` |
| `sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map` | `6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a` |

The normalization divisor was `0.44031545732594374`; the sampled Joint domain
again had 229,783 observations, 1,418 targets, and 659 supported halos. The
fixture in `tests/fixtures/mdpde/6z6u-fixed-halo-first-stage.json` records all
raw sample distances, responses and selection flags for each failed halo, the
resulting positive-response log dataset and design matrix, local alpha, OLS
start, solver options, fixed-point trace, and both recovery comparisons. Its
replay test is independent of the local source files.

| Atom | Joint rows | Samples / dataset rows | Alpha | Native 100: beta change² / relative variance change | Native 300 / 1,000 | Existing failed-only refinement |
| ---: | ---: | ---: | ---: | ---: | --- | --- |
| 16664 | 73 | 63 / 45 | 0.7 | `1.0213e-7` / `1.2630e-4` | SUCCESS at 268; same endpoint at 1,000 | Rejected: root residual `1.2620e-4` exceeds `1e-8`; root Jacobian condition `6.996e6` |
| 57190 | 1 | 140 / 68 | 0.4 | `1.0506e-9` / `5.9137e-5` | SUCCESS at 122; same endpoint at 1,000 | RefinedSuccess; fresh residual `2.08e-16`; reference residual `9.54e-11`; branch and weight-floor checks pass |
| 62463 | 81 | 99 / 72 | 0.3 | `8.4141e-7` / `3.6159e-4` | SUCCESS at 237; same endpoint at 1,000 | RefinedSuccess; fresh residual `5.47e-16`; reference residual `9.91e-11`; branch and weight-floor checks pass |
| 63216 | 50 | 71 / 53 | 0.7 | `1.7638e-9` / `3.9099e-5` | SUCCESS at 129; same endpoint at 1,000 | RefinedSuccess; fresh residual `8.77e-16`; reference residual `9.83e-11`; branch and weight-floor checks pass |

All four 100-iteration endpoints meet the beta-change criterion but fail the
relative-variance-change criterion (`1e-5`). Extending the unchanged native
fixed-point algorithm to a hard cap of 300 iterations qualified all four;
the 1,000-iteration control stopped at the same iteration and endpoint for
each. The amplitude/width design matrices have rank 2 and condition numbers
12.6–16.3. The MDPDE estimating-equation effective rank is 2, and all
variances remain above `1.9e-3` (weight floor `1e-8`). The OLS residuals are
not at the roundoff-exact-fit boundary. No non-finite intermediate values
occurred. The 16664 variance path had two direction reversals but no amplitude
or width reversals and then converged; the other three variance paths had no
direction reversals. This is consistent with a slow variance fixed point, not
a persistent oscillation, exact-fit boundary, or singular local design.

For 57190, 62463, and 63216, the endpoint refiner's root residual, independently
continued fixed-point reference, positive covariance, branch coordinate check,
weight consistency and floor mask all pass. For 16664, the root solve stopped
without a reference residual and failed its equation residual gate; its large
root-Jacobian condition and remaining-coordinate error reject that candidate.
The qualified 300-iteration native endpoint remains available for this atom.
Therefore no finite unqualified endpoint is needed: the evidence supports a
failed-only, bounded 300-iteration retry, with existing endpoint refinement
attempted first and the bounded retry used only if refinement is rejected.
Successful FirstStage fits retain their current 100-iteration path and do no
recovery work. The fixed-halo model and `C=0` rule are unchanged.

## P1 failed-only FirstStage recovery

The implementation runs recovery only for non-target fixed halos whose native
FirstStage solver did not return `SUCCESS`. It attempts the existing endpoint
refinement first. Only when that result is rejected and native status is
`MAX_ITERATIONS_REACHED` does it retry the same native fixed-point solver with a
300-iteration cap and unchanged tolerance and weight floor. A retry must return
native `SUCCESS`. In either recovery path, the persisted native status remains
`MAX_ITERATIONS_REACHED`; the fixed parameters remain the accepted FirstStage
A/B and C=0.

The formal 6Z6U chain-A initialization reported `valid-widths`, retained all 659
supported halos, and entered the target-only solver with a 1,418-target
component. The saved schema-15 snapshot contains 655 `NativeSuccess` and 4
`RefinedSuccess` fixed halos, all qualified, each with C=0. Recovery was tried
for 4/4 previously failed halos: endpoint refinement accepted 3; the remaining
atom 16664 used the bounded native retry and succeeded at iteration 268. Its
native beta change squared was `2.058e-10`, relative variance change was
`9.822e-6`, and diagnostic estimating-equation residual was `8.378e-6`.
The latter is retained as diagnostic evidence and is not tested against the
stricter endpoint-refinement root residual gate. The 300-iteration native retry
returned `SUCCESS` under FirstStage's existing stopping rule; the 1,000-iteration
control had already established that it reaches the same endpoint.

The four recoveries took 3.04 ms total; the largest per-halo recovery time was
1.50 ms. Fixed-halo JSON round-trip and historical schema-14 decoding pass.
Targeted recovery, fixed-halo validation, Joint snapshot and SQLite persistence
tests pass. When an initialization record remains unqualified, its per-atom
reason now records the recovery failure or native status along with its atom
identity and recovery diagnostics; it is not silently promoted or removed.

## P2 formal target-only validation

The supported CLI was rerun with the 6Z6U CIF, normalized map, saved database,
chain-A selection, and one worker. Initialization passed with 659/659 qualified
fixed halos (655 `NativeSuccess`, 4 `RefinedSuccess`), all C=0. The solver
entered one target-only component with 1,418 targets and 229,783 observations.
This is an initialization pass only; the target fit did not qualify.

| Measure | Formal result |
| --- | ---: |
| Completed sweeps / budget | 200 / 200 |
| Target-only blocks per sweep | 142 |
| Block visits: accepted / unchanged / failed | 15,310 / 13,090 / 0 |
| Accepted local updates | 13,123 |
| Profile evaluations | 56,490 |
| Recovery attempts / successes | 0 / 0 |
| Final candidate objective | 0.0041020821175608135 |
| Target A/C KKT | `2.682679252e-10` (required `<=1e-10`) |
| Target width-gradient infinity norm | `7.214524806e-12` (required `<=1e-12`) |
| Maximum target B | `0.6392831793` |
| Independent replay | Failed; no trusted endpoint |
| Converged components / available target estimates | 0 / 0 |
| Solver search time / whole CLI wall time | 519.0 s / 530.266 s |
| Peak RSS | 1,249,853,440 bytes |
| CLI exit / saved API outcome | 3 / `JointNotConverged` |

The sweep-by-sweep values are preserved in
[`fixed-halo-target-only-6z6u-sweep-trajectory.csv`](fixed-halo-target-only-6z6u-sweep-trajectory.csv).
The run used the existing 12-atom forward Gauss-Seidel cores, one trusted local
update per block visit, existing constrained A/C solve and Profile LM, width
safeguards, bounded block recovery, and independent verification. The gates
were unchanged. No block failed or needed recovery. B settled near 0.6393 with
no runaway. The objective flattened near 0.004102; in the last 25 sweeps KKT
varied from `2.214e-10` to `2.683e-10` and the reported width gradient stayed
at `7.215e-12`. The search therefore stopped on a stationarity
plateau at the sweep budget, not on an initialization, local-block, memory, or
width-runaway failure. Endpoint evidence records primary KKT failure and
independent width-gradient consistency failure (`per-coordinate-replay-enclosure`).
The candidate contains finite values for the selected targets, but none is a
qualified Joint estimate. Increasing the sweep cap alone is unsupported by
this trajectory; further target-stationarity diagnostics are needed before a
solver change can be justified.

The CLI now separates process completion from estimator qualification:
`CommandResult.succeeded` is true when the diagnostic is saved, while
`CommandResult.outcome` reports `JointNotConverged`; the CLI exits 3. A failed
initialization exits 2 after saving its diagnostic, and a qualified Joint fit
exits 0. The Python and C++ API and the CLI integration test cover these
outcomes. The reported exit 3 was captured from the CLI process itself; the
macOS `time -l` wrapper could not query `kern.clockrate` in this sandbox, so
wall time and RSS were measured by a Python child-process wrapper.

## P3 controlled fixed-halo model-error study

The four controlled scenarios were run sequentially in isolated research
processes against the same raw 229,783-row observation domain, original scale
`s = 346.7839502619515`, 1,418 targets, target initial widths, target-only
FixedNeighbor solver, unchanged gates, and 200-sweep maximum. The production
and Charge Oracle cases use the qualified 659-halo FirstStage A/B snapshot;
the two Gaussian Oracle cases replace halo A/B with simulation truth. All
Joint fixed-halo entries still have C=0. The Charge Oracle and Full Oracle
truth charge signals are represented only by research-side effective
observations; no production halo charge parameter or optimization is
introduced. Per-sweep trajectories and
isolated process wall time and peak RSS are recorded in
[`fixed-halo-target-only-6z6u-model-error.json`](fixed-halo-target-only-6z6u-model-error.json).

The halo-field error decomposition is:

| Signal comparison | L2 error | Per-row RMSE |
| --- | ---: | ---: |
| FirstStage A/B, C=0 vs truth halo signal | 37.56479 | 0.07836498 |
| Truth A/B, C=0 vs truth halo signal (missing charge) | 38.81258 | 0.08096802 |
| FirstStage A/B, truth C vs truth halo signal (Gaussian A/B error) | 13.59347 | 0.02835772 |
| Truth A/B/C vs truth halo signal | 0 | 0 |

Across the 659 halos, FirstStage parameter RMSE was 1.63253 for A and 0.031577
for B. With widths held at truth, the isolated Gaussian amplitude-error
contribution had L2/RMSE `14.06620 / 0.0293439`; the subsequent width-error
contribution had L2/RMSE `5.95145 / 0.0124155`. The missing charge signal
(`38.81258 / 0.0809680`) is larger than the combined FirstStage Gaussian
signal error (`13.59347 / 0.0283577`). These vector errors are not orthogonal:
the production total (`37.56479 / 0.0783650`) includes their interaction and
partial cancellation. The decomposition identity's maximum absolute residual
was `9.99e-16`.

| Halo scenario | Qualification | Sweeps | Candidate/verified target RMSE A / B / C | Target prediction RMSE | Full prediction RMSE | Objective | Scenario halo-signal RMSE | Process wall / fit seconds | Peak RSS |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Production: FirstStage A/B, C=0 | Finite candidate, unqualified | 200 | candidate `0.509995 / 0.011678 / 0.065912` | candidate `0.043347` | candidate `0.065527` | candidate `0.00410208` | 0.07836498 | 532.08 / 522.24 | 804,798,464 B |
| Gaussian Oracle A/B, C=0 | Finite candidate, unqualified | 200 | candidate `0.548196 / 0.012362 / 0.071264` | candidate `0.046432` | candidate `0.066625` | candidate `0.00424079` | 0.08096802 | 530.24 / 521.29 | 807,845,888 B |
| FirstStage A/B, Oracle C | Finite candidate, unqualified | 200 | candidate `0.066842 / 0.001588 / 0.011621` | candidate `0.006845` | candidate `0.027509` | candidate `0.00072299` | 0.02835772 | 551.97 / 543.00 | 805,224,448 B |
| Full Oracle A/B/C | Converged and verified | 143 | verified `4.74e-8 / 2.47e-9 / 2.65e-9` | verified `6.41e-9` | verified `1.72e-8` | `2.81935e-16` | 0 | 580.28 / 571.31 | 809,385,984 B |

The three unqualified cases had zero available formal target estimates and
1,418 finite candidate values each. Their candidate KKT values were
`2.174e-10`, `4.857e-10`, and `1.837e-10`; width-gradient values were
`2.184e-12`, `4.991e-12`, and `1.163e-12`, respectively. These miss one or
both unchanged gates (`1e-10` and `1e-12`); all three exhausted the sweep
budget, and no trusted full-observation replay was available. Their target
RMSE, target/full prediction RMSE, and objective entries above are explicitly
candidate-only and cannot quantify final target bias. All cases had zero
failed blocks and zero block-recovery attempts.

Full Oracle completed one target component and passed full-sweep confirmation
and independent replay. Its final target KKT was `1.494e-11`, width-gradient
infinity norm `6.927e-14`, and maximum accepted log-width step on the final
sweep `9.129e-11`; all are within their existing gates. It recovered the
simulation target parameters to near floating-point accuracy. Thus the
target-only solver can converge on this observation domain when the fixed
halo model matches the generating signal. This does not validate the other
three models' unqualified target candidates.

For the model-error attribution, the signal-level evidence shows that omitting
halo charge causes the larger fixed-halo error in this simulation. The
finite-candidate comparison also moves substantially when truth charge is
supplied, but it is not a qualified target-bias estimate. FirstStage A/B's
target-accuracy effect therefore remains unquantified; the production versus
Gaussian-Oracle candidates are both unqualified and do not establish which
target estimates are more accurate. The experiment does not establish that
the production C=0 assumption is scientifically acceptable. The required
production C=0 contract remains unchanged pending a converged comparison under
that model.

Compared with the historical Free-Halo FullABC run, the new target-only
trajectory has no block failures, no recovery attempts, and no width runaway;
the Full Oracle case converges with the same target solver. This is better
observed search stability on this workload, but the production C=0 fit itself
still stops at a stationarity plateau. No evidence from these runs justifies
changing FixedNeighbor, SparseQR, or Profile LM convergence settings.

## P4 API and integration closeout

The public API comments and developer guide now distinguish a raw
`BuildJointProblem` from a problem carrying a qualified FirstStage fixed-halo
snapshot. They state that `FitJointComponents` widths belong to selected
targets only, `FitJointComponentsInitialized` attaches qualified FirstStage
halos before fitting, and initialization validity is separate from target
convergence. `EstimateJointComponents` is documented as the construction,
initialization, and target-only fitting convenience path. Historical Free-Halo
JSON remains readable under its stored contract and is not treated as a new
fixed-halo snapshot.

The final Release build passed with `cmake --build build --parallel 4`. The
first complete test run exposed stale contract assertions: current-schema
writers now emit schema 15, the CLI returns 3 for a saved but unqualified
Joint fit, and the experiment cleanup test needed the FirstStage diagnostic
source in its allowlist. Those expectations were updated; the five affected
CTest entries then passed, followed by the final full suite:

| Verification | Result |
| --- | ---: |
| Release build, max parallelism 4 | Passed |
| Affected workflow, persistence/schema, cleanup, and CLI smoke tests | 5 / 5 passed |
| Complete CTest, `--parallel 4` | 40 / 40 passed in 73.84 s |
| Fixed-halo initialization and schema round-trip | Passed; current schema 15, historical schema 14 remains readable |
| PostFitPeeling / Local / Group integration coverage | Passed in Joint workflow and persistence tests |
| Large component memory preparation regression | Passed (`joint_component_mapping_memory_scalability`) |

The large memory test is a preparation/preflight regression, not a full
27,406-target fit. Its documented synthetic input remains 5,938,087
observations, 37,406 contributors, 27,406 targets, and 10,000 fixed halos;
the earlier recorded process peak was 1,585,790,976 bytes. The persistence
schema 15 assertions verify new writes while the historical schema 14 decode
paths continue to be exercised separately.

## P2 Practical global convergence validation (2026-10-16)

The same chain-A CIF, mode-2 float32 map, simulation sidecar, CLI map
normalization, FirstStage initialization, support selection, component
partition, initial target widths, 12-atom core, forward Gauss-Seidel order,
one trusted update per block visit, and 200-sweep budget were used for all four
re-runs. The normalized observation scale was `787.5806867375935`; the map
normalization divisor was `0.44031545732594374`. Outputs de-normalize A/C back
to the original map units before parameter comparisons. The production
scenario used 659/659 qualified FirstStage A/B halos with immutable C=0 and
1,418 target atoms in one component across 229,783 observations. Both
production Strict and Practical runs preserved the original observation
scale; the controlled-scenario scale factor was `1.0000000000000164` in the
Practical run.

The Strict production endpoint is an extended-search **candidate**. It
exhausted 200 sweeps and failed strict global KKT and width-gradient
qualification; its 1,418 parameter values are not a qualified reference. The
Practical production run completed full-sweep confirmation at sweep 79 and
passed endpoint verification, independent replay, full-observation replay,
and target availability checks. All 1,418 target estimates were qualified in
this controlled production scenario.

| Production FirstStage A/B, fixed C=0 | Strict candidate | Practical qualified |
| --- | ---: | ---: |
| Completed / confirmed / certified sweep | 200 / none / 200 | 79 / 79 / 79 |
| Candidate objective | 0.004102082117560809 | 0.004102082117643850 |
| Global A/C KKT | 3.568481e-10 | 9.197902e-8 |
| Global width-gradient infinity norm | 2.580471e-12 | 4.231091e-10 |
| Complete-sweep `max |Δlog B|` | Not saved for this unqualified candidate; the 200-sweep trajectory recorded a max accepted width step of 7.78516e-11 | 5.666148e-7 |
| Available qualified target estimates | 0 | 1,418 |
| Independent component replay | Ran; stationarity gates failed, consistency gates passed | Passed all stationarity and consistency gates |
| Full-observation / fixed-halo replay | Unavailable for unqualified endpoint | Passed |
| Failed blocks / maximum final target B | 0 / 0.639283 | 0 / 0.639283 |
| Profile evaluations / accepted updates | 56,745 / 13,093 | 22,439 / 10,785 |
| Search wall time / process peak RSS | 529.65 s / 682,672,128 B | 343.88 s / 665,731,072 B |

The Practical run used 121 fewer sweeps and saved 185.77 seconds (35.1%) of
search wall time versus the Strict 200-sweep candidate. Both runs reached the
same objective plateau to about `8.3e-14` absolute difference. Search time
here is measured inside the controlled fitting process; the Strict candidate
is not formal ground truth.

### A/B/C stability against the extended-search candidate

Differences below compare all 1,418 Practical parameters with the 200-sweep
Strict candidate after conversion to original-map A/C units. B remains in its
physical width units. This measures numerical stability only; it does not
measure scientific estimation error for the misspecified fixed-halo C=0
model.

| Parameter | Mean absolute difference | Median | P95 | P99 | Maximum | Atoms >0.0025 / >0.005 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A | 7.224e-7 | 5.600e-8 | 2.788e-6 | 1.288e-5 | 5.453e-5 | 0 / 0 |
| B | 2.307e-8 | 1.448e-9 | 8.038e-8 | 4.413e-7 | 1.945e-6 | 0 / 0 |
| C | 1.488e-7 | 8.490e-9 | 4.703e-7 | 2.878e-6 | 1.756e-5 | 0 / 0 |
| Combined A/B/C values | 2.981e-7 | 9.461e-9 | 9.371e-7 | 5.530e-6 | 5.453e-5 | 0 atoms / 0 atoms |

No atom exceeded either parameter-difference limit. The largest per-atom
difference was atom 400: `|ΔA|=5.453e-5`, `|ΔB|=1.870e-6`,
`|ΔC|=1.183e-5`. Its block 36 was classified WellConditioned with full-rank
free A/C design condition 6.91 and projected-width/profile condition about
1.24. Across all 142 Practical blocks, every block was WellConditioned and
full rank; the maximum free A/C condition was 7.014 and maximum projected
width/profile condition was 1.264. There are no problematic atoms above the
requested 0.0025 threshold.

Target and full prediction differences between the Practical endpoint and
Strict candidate were identical up to their shared fixed-halo signal: maximum
absolute difference `2.6794e-5`, mean absolute difference `1.3475e-7`, and
RMSE difference `6.6924e-7`, in CLI-normalized map units. The candidate
objective difference was `8.3040e-14`.

### Independent replay and the previous coordinate-enclosure report

At the Practical production endpoint, primary and independent values were:

| Gate / diagnostic | Primary | Independent | Result |
| --- | ---: | ---: | --- |
| Global A/C KKT | 9.197902e-8 | 9.197902e-8 | Both ≤ 1e-7; difference 0 |
| Width-gradient infinity norm | 4.231091e-10 | 4.231091e-10 | Both ≤ 1e-9; coordinate agreement passed |
| Prediction replay max difference | — | 1.066e-14 | Passed |
| Residual replay max difference | — | 1.066e-14 | Passed |
| Objective replay difference | 2544.4531891830266 | 2544.4531891830156 | Difference 1.0914e-11; original objective enclosure passed |
| Parameter feasibility | — | — | Passed |
| Full observation assembly including fixed halos | — | — | Passed |

The largest Practical per-coordinate width-gradient difference was
`1.1411e-18` at atom 2497, against its unchanged `1.0000e-13` enclosure
(ratio `1.1411e-5`). The 200-sweep Strict candidate's primary and independent
KKT differed by `1.7347e-18`; its maximum coordinate width-gradient
difference was `1.2602e-18` at atom 2656. KKT agreement, gradient agreement,
prediction, residual, objective, and parameter-feasibility consistency all
passed on that candidate, but both Strict stationarity gates failed, so its
overall replay qualification and full-observation assembly remained
unavailable. The earlier `per-coordinate-replay-enclosure` failure was not
reproduced in these re-runs; no failed coordinate or failure ratio exists in
the current Strict/Practical evidence. The P1 regression still demonstrates
that a deliberately perturbed coordinate fails the unchanged enclosure.
This change did not loosen any replay-consistency bound.

### Full Oracle control

Full Oracle replaces the research scenario's fixed halo A/B with simulation
truth while fixed halo C remains zero; truth halo charge is represented only
in the research effective observations. Practical and Strict both passed
full-sweep confirmation and independent/full-observation replay, with all
1,418 target states available.

| Full Oracle A/B/C | Strict | Practical |
| --- | ---: | ---: |
| Sweeps / wall time | 143 / 599.79 s | 79 / 343.34 s |
| Objective | 2.81935e-16 | 7.95680e-14 |
| Candidate A/B/C RMSE to truth | 4.740e-8 / 2.466e-9 / 2.652e-9 | 3.273e-6 / 1.158e-7 / 8.201e-7 |
| Independent replay | Passed | Passed |
| Failed blocks / full-observation replay | 0 / passed | 0 / passed |

Practical Full Oracle parameters differed from the Strict qualified endpoint
by at most `5.320e-5` in A, `1.900e-6` in B, and `1.714e-5` in C; no target
exceeded 0.0025. The Practical objective is larger by `7.929e-14`, while
target A/B/C truth errors remain below `3.3e-6`. This is no material loss of
near-truth recovery at the requested parameter precision.

The fixed C=0 model error remains a separate scientific concern. The
Practical convergence result does not change that model contract or establish
that production target bias is below 0.01.

The complete four-run sweep trajectories, raw-unit parameter checkpoints,
compact replay/fit summary, and local-conditioning diagnostics are saved in:

- [`fixed-halo-target-only-6z6u-practical-trajectories.csv`](fixed-halo-target-only-6z6u-practical-trajectories.csv)
- [`fixed-halo-target-only-6z6u-practical-checkpoints.csv`](fixed-halo-target-only-6z6u-practical-checkpoints.csv)
- [`fixed-halo-target-only-6z6u-practical-summary.json`](fixed-halo-target-only-6z6u-practical-summary.json)

**P2 decision:** the practical global criteria passed the 6Z6U production and
Full Oracle controls, preserved independent replay checks, and met the
0.0025 A/B/C early-stop stability target. The Strict production 200-sweep
candidate remains unqualified. This evidence permits the conditional P3
production promotion; it does not claim statistical accuracy for the
misspecified C=0 model.
