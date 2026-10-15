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
