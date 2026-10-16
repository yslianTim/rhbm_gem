# Dynamic Halo Mean Joint Estimator Research

## Scope and current status

This research keeps the production estimator unchanged. Production continues to use qualified FirstStage fixed Halo parameters `(A, B, C) = (A_FirstStage, B_FirstStage, 0)`, 12-target FixedNeighbor cores, and the `practical-v1` convergence criteria. The dynamic implementation is research-only.

| Phase | Status | Evidence |
|---|---|---|
| P0: Production baseline and ownership | Completed | Baseline rerun, deterministic ownership helper and contract tests, saved checkpoint |
| P1: Post-fit AB/ABC means | Completed | Fixed-target halo and full-field comparison in checkpoint |
| P2: Coupled Dynamic Halo optimization | Not completed | No coupled optimizer or coupled stationarity certification yet |
| P3: Robustness and sensitivity | Not completed | Requires a qualified P2 solver |

## Production contract and baseline

The Joint problem selects 1,418 Chain A heavy-atom Targets from 2,077 supported contributors. Production fits all Targets and subtracts 659 qualified fixed Halos. The baseline uses 229,783 observation rows, one Target component, 12-target cores, and original observation scale `s = 787.5806867375935`. Map normalization divisor was `0.44031545732594374`.

The baseline was rerun on `develop` at source commit `c84eeb42c0abaac314fce20309114313ada5af40` using `/Users/yslian/data/6z6u.cif`, `/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map`, and its `.simulation.json` truth sidecar. The production CLI exited successfully. Model SHA-256: `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa`; map SHA-256: `6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a`.

| Measure | Rerun |
|---|---:|
| Target states qualified | 1,418 / 1,418 |
| Fixed Halos qualified | 659 / 659 |
| Sweeps | 79 |
| Profile evaluations | 22,436 |
| Accepted updates | 10,796 |
| Objective | 0.004102082117643847 |
| Global A/C KKT | 9.197884806924628e-8 |
| Global log-width gradient | 4.2311477169481625e-10 |
| Complete-sweep log-width change | 5.666472405740919e-7 |
| Runtime convergence | passed (`practical-v1`) |
| Independent replay | passed |
| Full-observation replay | passed |
| Search time | 334.417 s (reported by production) |

This agrees with the prior reference values to their reported precision. The saved checkpoint preserves row and contributor identities, support edges and squared distances, Target initial/final states, Halo FirstStage states and qualification, simulation truth, block partition, and ownership. It avoids dense observation-by-contributor matrices.

## Halo ownership

For each supported Halo, ownership uses the count of unique observation rows shared with each block's union of Target support rows. The maximum overlap wins; equal counts are resolved by lexicographically stable block ID. The mapping is constructed before loading simulation truth. Every one of the 659 Halos has exactly one owner; 410 intersect more than one block, and none cross Target components in this dataset.

Across all blocks, Target affected-row memberships total 405,716; adding owned-Halo support yields 423,469 memberships. The largest dynamic block has 4,929 rows and the largest per-block expansion is 559 rows. The helper stores sparse row memberships and does not create one full-length observation mask per block. Unit tests cover deterministic tie breaking under reversed traversal, unique ownership, affected-row expansion, cross-component detection, physical arithmetic means, and the width chain-rule factor.

## P1: post-fit halo means

P1 leaves every qualified Production Target parameter fixed and recomputes only the Halo field. It evaluates the production Gaussian and charge kernels. Raw-map RMSE converts normalized fit amplitudes by the measured map normalization divisor. Every objective uses the original fixed observation scale.

| Halo model | Halo RMSE (raw map units) | Full prediction RMSE (raw map units) | Objective |
|---|---:|---:|---:|
| Production FirstStage A/B, C=0 | 0.0783650 | 0.0655265 | 0.0041020821 |
| Post-fit block mean AB, C=0 | 0.0950808 | 0.0826364 | 0.0065239863 |
| Post-fit block mean ABC | 0.1139574 | 0.1235819 | 0.01459081598 |

AB increases Halo-field RMSE by 21.33%; ABC increases it by 45.42%. The ABC field is worse than AB, so adding mean Halo C does not explain or fix the degradation. The P1 full-prediction error and observation objective also worsen. Target truth RMSE is unchanged by construction: A `0.509995` raw units, B `0.0116784 Å`, and C `0.0659125` raw units. These are diagnostic candidates, not Dynamic Joint results.

By Halo support size, the 626 Halos with three or more observation rows dominate the field and degrade under both post-fit means. ABC improves the small one-row and two-row Halo strata slightly, while AB is approximately neutral there; those strata are too small to offset the larger-support degradation.

## P2 and P3

P1 does not update Targets and cannot establish whether coupled re-optimization would produce a different result. P2 still requires the augmented A/C design, exact arithmetic-mean width chain rule in the profile Jacobian, block transaction rollback/commit, coupled global projected KKT and width gradient, and independent full-objective replay. Until these checks are implemented and pass, no Dynamic Halo candidate is qualified.

P3's E0/E1/E2/E3 full optimization comparison, forward/reverse traversal, alternate mean ranges, synthetic robustness runs, and large-component resource checks remain pending. No recommendation to promote Dynamic Halo follows from P1; the current evidence supports retaining the production fixed-Halo behavior while P2 is investigated.

## Reproduction and checkpoint

The research-only executable is `build/bin/joint_dynamic_halo_mean`. It consumes the model, map, a qualified Production result dump, an output directory, and the source baseline commit. The compressed checkpoint contains `summary.json`, row identities and normalized observations, sparse contributor support edges, contributor states and truth, plus the production result dump and atom CSV. Archive SHA-256: `b36bac6a26bf5b1d3beec9b3a9c3ef5229190c85a4b3b25f1be0298acd729084`.

Focused verification: `rhbm_tests_joint_component_contract` and `joint_experiment_cleanup_test` passed after building with parallelism 4. The production baseline was separately rerun before this research-only code was added.
