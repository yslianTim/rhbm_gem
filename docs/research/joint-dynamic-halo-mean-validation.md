# Dynamic Halo Mean Joint Estimator Research

> Historical study: its Dynamic Halo implementation and runners have been
> retired. The results below are retained as research evidence. Current
> Production behavior is documented in
> [`joint-bounded-active-set-production-promotion-20261020.md`](joint-bounded-active-set-production-promotion-20261020.md).

## Scope and current status

This study evaluated an earlier Production baseline using qualified FirstStage
fixed Halo parameters `(A, B, C) = (A_FirstStage, B_FirstStage, 0)`, 12-target
FixedNeighbor cores, and the `practical-v1` convergence criteria. That baseline
predates Pure Halo-Support Exclusion and Active-set LM; the dynamic
implementation and its runners are retired.

| Phase | Status | Evidence |
|---|---|---|
| P0: Production baseline and ownership | Completed | Baseline rerun, deterministic ownership helper and contract tests, saved checkpoint |
| P1: Post-fit AB/ABC means | Completed | Fixed-target halo and full-field comparison in checkpoint |
| P2: Coupled Dynamic Halo optimization | Completed | E0–E3 comparison; E2/E3 pass coupled stationarity and replay; E1 retained as finite-candidate-unqualified |
| P3: Robustness and sensitivity | Completed | Forward/reverse, overlapping-target mean, runner-up owner, six fixed-seed synthetic cases, and large-scale arithmetic preflight |

## Research hypothesis and experimental groups

The hypothesis is that a block's fitted Target mean may be a better local proxy for nearby non-target contributors than independent FirstStage fits. The competing risk is that a block mean smooths real atom-to-atom variation, and that a shared mean C can introduce charge-kernel error. Target truth is the main scientific criterion; a smaller observation objective alone is not evidence of better Target estimates.

| Group | Halo parameters | Halo updates | Interpretation |
|---|---|---|---|
| E0 Production | qualified FirstStage A/B, C=0 | fixed | Formal baseline |
| E1 Dynamic AB | current owner-core mean A/B, C=0 | after each accepted block transaction | Tests dynamic A/B coupling |
| E2 Dynamic ABC | current owner-core mean A/B/C | after each accepted block transaction | Primary research group |
| E3 Static Mean ABC | whole-core mean from common initial Target state | fixed after initialization | Separates initial mean from dynamic feedback |

The common Target start is A=C=0 with the qualified FirstStage B widths, matching the inspected production Target solver initialization. Simulation truth is opened only after ownership is fixed and is used only for evaluation.

## Coupled model and numerical checks

For each supported Halo `h`, the owner is the block maximizing the count of unique shared observation rows between `V_h` and the union of that block's Target supports. Stable block ID breaks ties. If `H_b` is the owned Halo set and `n_b` the core size, the dynamic models use

`A_h = mean_i A_i`, `B_h = mean_i B_i`, and, for E2, `C_h = mean_i C_i` over Targets in the owner core. E1 fixes Halo C to zero; E3 takes all means from the common initial state once.

The coupled objective is `F = ||T(theta_T) + H(M(theta_T)) - y||^2 / (2s^2)` with the original fixed scale `s`. Dynamic block rows are the union of Target-core and owned-Halo support rows. For fixed widths, each Target A/C design column contains its own kernel plus `1/n_b` of every owned Halo kernel on those Halo rows. Halo charge columns are omitted for E1. The constrained linear solve keeps A nonnegative and C signed and is accepted locally only with projected KKT at most `1e-10`.

Widths use `eta_i = log(B_i)` while Halo width remains the physical arithmetic mean. The coupled Halo width derivative is evaluated as `d log(B_h)/d eta_i = B_i/(n_b B_h)` times the kernel's log-width derivative. The profile Jacobian projects these full coupled width derivatives against the free A/C design columns; active A columns remain excluded from that projection. Unit tests compare A/C and log-width derivatives against finite differences, including an active A boundary and shared Halo supports.

Every block proposes a Target state, derives its owned Halo state, recomputes the local/full objective delta, verifies it against the Production block replay enclosure, and commits all parameters together only if the full objective is monotone. Rejections leave the accepted Target/Halo state intact. A separate direct contributor replay rebuilds all Halo means and the full prediction. Global stationarity uses effective per-Target A/C columns after summing overlapping Target/Halo contributions row by row, plus the exact coupled log-width gradient. The practical comparison thresholds are the unchanged Production values: A/C KKT `1e-7`, width gradient `1e-9`, and complete-sweep change `1e-6`; two consecutive stationary sweeps are required.

## Production contract and baseline

The Joint problem selects 1,418 Chain A heavy-atom Targets from 2,077 supported contributors. Production fits all Targets and subtracts 659 qualified fixed Halos. The baseline uses 229,783 observation rows, one Target component, 12-target cores, and original observation scale `s = 787.5806867375935`. Map normalization divisor was `0.44031545732594374`.

The baseline was rerun on `develop` at source commit `c84eeb42c0abaac314fce20309114313ada5af40` using `/Users/yslian/data/6z6u.cif`, `/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map`, and its `.simulation.json` truth sidecar. The production CLI exited successfully. Model SHA-256: `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa`; map SHA-256: `6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a`.

Simulation truth sidecar SHA-256: `33a3c221ec3f2fff18dcd2a0ea6484f68334971fc688f7cd100e220ad5848931`.

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

| Halo support rows | Halo count | Production RMSE | Post-fit AB RMSE | Post-fit ABC RMSE |
|---|---:|---:|---:|---:|
| 1 | 20 | 0.191817 | 0.191837 | 0.168199 |
| 2 | 13 | 0.180328 | 0.180357 | 0.171487 |
| 3+ | 626 | 0.169119 | 0.205195 | 0.245995 |

Halo parameter RMSEs further separate mean components. Production FirstStage Halo A/B/C RMSEs are `1.63253` raw units / `0.0315770 Å` / `0.408623` raw units. Post-fit AB means give `1.22116` / `0.0552455` / `0.408623`; post-fit ABC gives `1.22116` / `0.0552455` / `0.452192`. A error improves, B error worsens, and adding mean C worsens Halo C parameter error and overall field reconstruction.

Across the 1,418 Target estimates, within-owner-core deviations from the block mean have RMSE `0.941424` raw units for A, `0.0389478 Å` for B, and `0.445602` raw units for C. The optional Halo-field error projection onto the Target Jacobian was not computed; the saved checkpoint does retain owner overlaps, support counts, and parameter states needed for follow-up.

## P2: coupled Dynamic Halo optimization

The research-only solver implements the coupled objective described above. It derives each Halo from its owner Targets at every candidate, adds owned-Halo supports to block rows, uses the augmented A/C design and exact log-width chain rule, and commits or rejects each Target/Halo state atomically after full objective replay. Endpoint checks independently rebuild every Halo and the full prediction. No production source or default was changed.

All groups used the same `A=C=0`, FirstStage-width Target start, observations, block partition, local solver, trust checks, fixed observation scale, and 79-sweep cap. The table reports Target parameter errors in truth units. Improvements use `(E0 RMSE - candidate RMSE) / E0 RMSE`; negative values mean worse accuracy.

| Group | Status | Target A RMSE (improvement) | Target B RMSE Å (improvement) | Target C RMSE (improvement) |
|---|---|---:|---:|---:|
| E0 Production | qualified | 0.509995 (—) | 0.0116784 (—) | 0.0659125 (—) |
| E1 Dynamic AB | finite-candidate-unqualified | 0.494875 (+2.96%) | 0.0117066 (-0.24%) | 0.0712300 (-8.07%) |
| E2 Dynamic ABC | qualified | 0.337742 (+33.78%) | 0.0097210 (+16.76%) | 0.0751745 (-14.05%) |
| E3 Static Mean ABC | qualified | 0.498010 (+2.35%) | 0.0115611 (+1.00%) | 0.0723486 (-9.77%) |

E2 improves Target A and B, but worsens C. E1 narrowly misses the unchanged global A/C KKT threshold (`1.1619e-7` versus `1e-7`) at the sweep cap. It remains `finite-candidate-unqualified`; its passing width-gradient or replay checks do not qualify it.

| Group | Halo prediction RMSE | Target prediction RMSE | Full prediction RMSE | Objective | Sweeps / profiles | Accepted / rejected / failed | Global A/C KKT | Width gradient | Sweep change | Search seconds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| E0 Production | 0.0783650 | 0.0433466 | 0.0655265 | 0.00410208 | 79 / 22,436 | 10,796 / — / — | 9.1979e-8 | 4.2311e-10 | 5.6665e-7 | 334.417 |
| E1 Dynamic AB | 0.473410 | 0.242885 | 0.416737 | 0.00623662 | 79 / 44,311 | 11,218 / 0 / 0 | 1.1619e-7 | 4.2547e-10 | 9.1825e-7 | 425.297 |
| E2 Dynamic ABC | 0.449233 | 0.218320 | 0.425364 | 0.00649751 | 79 / 44,481 | 11,218 / 0 / 0 | 7.6910e-8 | 3.2515e-10 | 6.0316e-7 | 426.824 |
| E3 Static Mean ABC | 0.833862 | 0.222413 | 0.803996 | 0.0232131 | 79 / 30,531 | 11,218 / 0 / 0 | 8.4268e-8 | 3.8380e-10 | 5.1939e-7 | 364.478 |

RMSE columns are raw map units. E1, E2, and E3 all have exact Halo-mean consistency and zero independent replay prediction/objective discrepancy. E2 and E3 satisfy the unchanged coupled thresholds (`1e-7` A/C KKT, `1e-9` width gradient, `1e-6` complete-sweep change) for two confirmed sweeps. E1 is the only unqualified group. The E2 observation objective is 58.4% above E0; its Halo, Target, and Full prediction RMSE are respectively 5.73×, 5.04×, and 6.49× the E0 values. Therefore the Target A/B improvements are not accompanied by better field reconstruction or objective fit. The experiment does not establish that these Target parameter changes are identifiable from this observation model.

The full forward comparison peaked at 652,247,040 bytes RSS. The largest 6Z6U local dense design is 946,368 bytes; no global dense design is constructed. For the requested large-component arithmetic preflight (5,938,087 rows; 37,406 contributors), the forbidden global design would be 3,553,921,317,152 bytes, two prediction/residual buffers total 95,009,392 bytes, and row-to-block vector headers total 142,514,088 bytes. This is an arithmetic estimate, excludes input/support/allocator storage, and is not a measured large-map RSS result.

## P3: robustness and sensitivity

### S1: block traversal order

The E2 forward and reverse runs started from the same initial state and both qualified. Reverse converged after 68 sweeps (forward used 79); both had zero rejected transactions, zero failed blocks, and exact independent replay. Target parameter differences, computed per Target across the two final tables, were:

| Parameter | RMSE of forward-minus-reverse | Mean absolute difference | Maximum absolute difference |
|---|---:|---:|---:|
| A, raw units | 3.48e-6 | 1.36e-6 | 5.08e-5 |
| B, Å | 1.20e-7 | 4.52e-8 | 1.74e-6 |
| C, raw units | 7.20e-7 | 2.21e-7 | 8.66e-6 |

Forward/reverse objectives were `0.006497513583164643` and `0.006497513583172617`. The final Target errors were effectively unchanged, so this dataset shows little traversal-order sensitivity after qualification.

### S2: mean scope

The secondary overlapping-target mean averages over the owner-block Targets whose observation supports intersect the Halo. Its owner-mean sizes ranged from 1 to 12; 197 Halos use one overlapping Target. It also qualified after 79 sweeps with zero rejected/failed blocks and exact independent replay.

| E2 Halo mean | Target A RMSE | Target B RMSE Å | Target C RMSE | Halo prediction RMSE | Full prediction RMSE | Objective |
|---|---:|---:|---:|---:|---:|---:|
| Whole owner core | 0.337742 | 0.0097210 | 0.0751745 | 0.449233 | 0.425364 | 0.00649751 |
| Overlapping Targets | 0.368841 | 0.0105393 | 0.0848194 | 0.442291 | 0.429314 | 0.00661874 |

Restricting the mean to overlapping Targets slightly improves Halo RMSE but worsens Target A/B/C and Full prediction relative to the primary whole-core E2. Both mean scopes retain the same scientific tradeoff against E0: better Target A/B but worse Target C and field prediction.

### S3: owner sensitivity

The primary owner is the maximum-overlap block with stable block-ID tie-breaking. Of 659 Halos, 249 touch one block, 165 touch two, 81 touch three, 63 touch four, 47 touch five, 26 touch six, 16 touch seven, 11 touch eight, and 1 touches nine. Thus 410 Halos (62.2%) are shared. Owner overlap counts have quartiles 25, 55, and 115 rows (range 1–370). No Halo crosses Target components.

The sensitivity run assigns a shared Halo to its second-ranked block under the same overlap-count and stable-ID ordering; single-block Halos retain their only owner. It reassigns all 410 shared Halos. Runner-up owner overlap counts have quartiles 11, 33, and 64 rows (range 1–231). This is an intentionally deterministic perturbation, not a proposed production ownership rule.

| E2 owner policy | Status | Target A/B/C RMSE | Halo / Target / Full prediction RMSE | Objective | Sweeps / profiles | A/C KKT | Width gradient | Search seconds |
|---|---|---|---|---:|---:|---:|---:|---:|
| Maximum overlap | qualified | 0.337742 / 0.0097210 / 0.0751745 | 0.449233 / 0.218320 / 0.425364 | 0.00649751 | 79 / 44,481 | 7.6910e-8 | 3.2515e-10 | 426.824 |
| Runner-up overlap | qualified | 0.334524 / 0.0093620 / 0.0719280 | 0.439302 / 0.208837 / 0.420721 | 0.00635645 | 79 / 43,663 | 8.4467e-8 | 4.0396e-10 | 429.017 |

The runner-up improves the three Target RMSEs by about 1.0%, 3.7%, and 4.3% relative to the maximum-overlap E2, and reduces Halo/Full prediction RMSE by about 2.2%/1.1%. It does not approach E0 field error. The reassigned-owner row expansion is 441,488 total affected memberships, 4,748 maximum rows in one block, and 1,051 maximum added rows; its maximum local dense design estimate is 911,616 bytes. The primary owner run had 423,469 affected memberships, 4,929 maximum rows, and 946,368 local-design bytes. Both runs had zero failed/rejected transactions, zero mean-consistency error, and exact replay.

### S4: synthetic perturbations

Six deterministic one-block synthetic cases use the production Gaussian/charge kernels, constrained A/C fits, fixed Target widths set to truth, and fixed observation-noise seeds (11, 17, 19, 23, 29, 31). They vary homogeneous/heterogeneous Target A/B/C, matched and mismatched Halo truth, C/O/N labels, mixed-sign charge, Halo support counts 1/2/8/16, and noise standard deviations 0.005/0.02/0.05/0.1. Matched cases generally improve Halo/full prediction; the mixed-sign, one-row mismatch worsens Target A/C and Halo/full prediction. This isolates A/C mean coupling and does not test nonlinear width estimation or multi-block synthetic convergence.

| Synthetic case (seed) | Halo rows / noise SD | Target A RMSE fixed → mean | Target C RMSE fixed → mean | Halo RMSE fixed → mean | Full RMSE fixed → mean |
|---|---:|---:|---:|---:|---:|
| Homogeneous matched (11) | 8 / 0.005 | 0.04385 → 0.00466 | 0.03847 → 0.00285 | 0.1089 → 0.00010 | 0.08798 → 0.00148 |
| Heterogeneous matched charge (17) | 16 / 0.005 | 1.136 → 0.01729 | 0.06063 → 0.00430 | 0.1150 → 0.00045 | 0.05168 → 0.00242 |
| Heterogeneous Halo mismatch (19) | 16 / 0.02 | 1.360 → 1.360 | 0.2021 → 0.1444 | 0.7241 → 0.5130 | 0.4238 → 0.2833 |
| Mixed elements, short support, matched charge (23) | 2 / 0.02 | 0.02133 → 0.05460 | 0.01812 → 0.01090 | 0.1203 → 0.00185 | 0.1151 → 0.00652 |
| Mixed-sign charge, one-row mismatch (29) | 1 / 0.05 | 0.09726 → 0.6423 | 0.03740 → 0.06897 | 0.1527 → 0.1856 | 0.1497 → 0.1916 |
| Higher noise, matched (31) | 8 / 0.1 | 0.04761 → 0.05555 | 0.02050 → 0.00941 | 0.07648 → 0.00459 | 0.06613 → 0.01949 |

Target B was held at its truth width in this controlled A/C solve, so these synthetic cases do not measure B estimation error. The results show that a mean can work well when Halo truth matches the Target mean and can fail sharply when a weakly supported Halo has conflicting signed charge.

### Recommendation

Retain the current Production fixed-Halo estimator. E2 is numerically qualified and materially improves Target A/B truth RMSE, but it worsens C truth RMSE, all three prediction RMSEs, and the observation objective. The runner-up order and mean-scope experiments do not remove that tradeoff. P3's large real input was unavailable; the deterministic preflight and synthetic cases are the documented substitute, with no claim of measured RSS at 5,938,087 rows.

## Reproduction, artifacts, and verification

The research-only executables are `build/bin/joint_dynamic_halo_mean`, `build/bin/joint_dynamic_halo_joint`, and `build/bin/joint_dynamic_halo_synthetic_robustness`. The coupled runner accepts group, traversal order, mean scope, and owner policy; the sensitivity outputs preserve `summary.json` and `targets.csv`. Synthetic output records all six cases and seeds.

| Artifact | Contents | SHA-256 |
|---|---|---|
| `docs/research/artifacts/dynamic-halo-p0-p1-6z6u.tar.gz` | P0 checkpoint, row identities, observations, supports, contributor truth/states, owners, production result and P1 metrics | `70418b2d130aeb55170f4f545f7f6ecd75d9ba09155efb5367e97465d368a376` |
| `docs/research/artifacts/dynamic-halo-p2-joint-6z6u.tar.gz` | Forward E0–E3 joint comparison summaries and Target tables | `8fc3ca83c33776bd71963faa17530a9cc19888445f972ada6696126c2ddf0055` |
| `docs/research/artifacts/dynamic-halo-p3-sensitivity-6z6u.tar.gz` | Reverse-order, overlapping-Target, runner-up-owner, and synthetic summaries/tables | `c58145a1b6c817fc924b8d9c50cffb70b9f5eab5ec23eee4f583d0500494b174` |
| `docs/research/artifacts/dynamic-halo-synthetic-robustness-6z6u.json` | Six deterministic fixed-width synthetic perturbations | `955b2748cfb4177e523f64d4e1e97b18e989a3d445f5fb5d65cea8d7d1a5c192` |

The focused `rhbm_tests_joint_component_contract` and `joint_experiment_cleanup_test` passed. The full `ctest --test-dir build --output-on-failure -j 4` run passed all 42 tests. During final verification, the statistical smoke harness was corrected to summarize the selected Target plus any fixed Halo actually returned by the estimator instead of assuming two Target parameter rows. The production CLI baseline was separately rerun at source commit `c84eeb42c0abaac314fce20309114313ada5af40` before adding research-only code. Build configuration is Release with Unix Makefiles; every build used parallelism 4, and all optimization experiments used one compute thread. No production estimator source, default, or persistence contract changed.
