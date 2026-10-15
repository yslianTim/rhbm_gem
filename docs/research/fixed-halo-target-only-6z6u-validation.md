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
