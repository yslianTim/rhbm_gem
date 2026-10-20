# Fixed-target Observation-driven Halo FullABC — P3 validation

Date: 2026-10-21
Branch: `develop`
Status: research implementation; current Production reconstruction remains `geometry-local-target-median-v1`.

## Contract and experiment

The runner `tests/experiments/joint_halo_fullabc_validation.cpp` creates each synthetic map on a 51×63×16 grid at 0.5 Å spacing, selects 20 Targets, and uses deterministic sinusoidal noise. It executes the actual FirstStage initialization and Production Pure Halo-Support Exclusion Target Joint estimator. A scenario proceeds to Halo reconstruction only if that first Target endpoint qualifies. The same captured Target A/B/C values are then used by geometry kNN-4, sparse observation A/C at fixed kNN widths, FullABC, and Oracle signal scoring. Postprocessing checks the Target endpoint bit-for-bit and rejects any drift.

All Halo-Halo centers in the normal scenarios are strictly farther than 1.5 Å apart. The requested distances 6.2, 4, 3, 2, and 1.6 Å are present. Noise amplitudes are 0.003 and 0.005. Similar and contrasting A/B/C profiles, positive and negative charges, width contrast, a map-edge Halo, and reconstruction-only support are covered. Separate cases put exactly 3 or 4 distinct OriginalInput observations on a free Halo.

The raw artifacts are [results.json](results/halo_fullabc_p3/results.json), [summary.csv](results/halo_fullabc_p3/summary.csv), and [halos.csv](results/halo_fullabc_p3/halos.csv). The runner accepts an optional scenario name for focused reproduction. Its CMake target is intentionally an experiment executable, not part of the default production flow.

## Physical geometry and observation coverage

| Scenario | Noise | Minimum Halo-Halo center distance | Pairs sharing OriginalInput voxels / pairwise shared voxels | `n_h` histogram | Fixed / free | Target / FullABC qualified |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| separated 6.2 Å, similar | 0.003 | 6.2 Å | 0 / 0 | 0:1, 1:1, 2:1, 102:1, 109:1 | 3 / 2 | yes / yes |
| separated 6.2 Å, higher noise | 0.005 | 6.2 Å | 0 / 0 | 0:1, 1:1, 2:1, 102:1, 109:1 | 3 / 2 | yes / yes |
| overlap 4 Å, similar | 0.003 | 4.0 Å | 1 / 16 | 0:1, 1:1, 2:1, 191:1, 109:1 | 3 / 2 | yes / yes |
| overlap 3 Å, contrasting | 0.003 | 3.0 Å | 1 / 40 | 0:1, 1:1, 2:1, 146:1, 109:1 | 3 / 2 | yes / yes |
| overlap 2 Å, width/charge contrast | 0.003 | 2.0 Å | 1 / 61 | 0:1, 1:1, 2:1, 105:1, 109:1 | 3 / 2 | yes / yes |
| overlap 1.6 Å, contrasting | 0.003 | 1.6 Å | 1 / 70 | 0:1, 1:1, 2:1, 103:1, 109:1 | 3 / 2 | yes / yes |
| overlap 4 Å, width/charge contrast | 0.003 | 4.0 Å | 1 / 16 | 0:1, 1:1, 2:1, 191:1, 109:1 | 3 / 2 | yes / yes |
| map-edge reconstruction-only Halo | 0.003 | 3.0 Å | 1 / 40 | 0:2, 1:1, 2:1, 146:1, 109:1 | 4 / 2 | yes / yes |
| exactly 3 observation free-Halo stress | 0.003 | 3.0 Å | 1 / 40 | 0:1, 1:1, 2:1, 3:1, 146:1, 109:1 | 3 / 3 | yes / no |
| exactly 4 observation free-Halo stress | 0.003 | 4.0 Å | 1 / 16 | 0:1, 1:1, 2:1, 4:1, 191:1, 109:1 | 3 / 3 | yes / no |

The row counts are distinct `OriginalInput` voxel IDs. Both the 3-row and 4-row Halos enter the free set as specified; neither is pre-pruned by rank, conditioning, or labels. Both then fail numerical qualification. The map-edge Halo has zero observations, is fixed to the qualified Target median before solving, and is retained in the reconstruction contributor set.

The synthetic support measurements distinguish atom-center separation from sampled support overlap. At 4, 3, 2, and 1.6 Å, supports shared 16, 40, 61, and 70 discrete voxel rows respectively. The 6.2 Å pair shared none. The synthetic Halo-to-nearest-Target distances range from 3.0 to 4.84 Å in the primary scenes; the additional reconstruction-only edge Halo is 6.04 Å from the nearest selected Target.

The real 6Z6U chain A census is recorded in [the P0 report](joint-halo-observation-census-p0-20261020.md). It has 35,988 Halo candidates, of which 35,362 (98.2605%) have `n_h < 3`; among 659 supported Halos, 33 are fixed and 626 are free. The real-coordinate minimum Halo-Halo distance is 1.1782 Å for all candidates and 1.1790 Å among supported Halos, so the >1.5 Å assumption is false for this structure. Those close pairs are covalent neighbors, not duplicate centers. Discrete support overlap remains separately measured.

## Estimator used

The first Target Joint uses the unchanged Production Pure Exclusion path. P2 accepts only that qualified endpoint and checks `JointProblem::OriginalInput()` identity. It subtracts every frozen Target A/B/C prediction from OriginalInput responses, then subtracts every pre-fixed Halo contribution. Qualified Target componentwise medians supply complete fixed A/B/C profiles for all `n_h < 3` Halos. They remain reconstruction contributors.

Every `n_h >= 3` Halo enters the free FullABC layout in deterministic Atom-ID order. The second optimization uses the saved OriginalInput rows and sparse Halo support incidence; it does not allocate a dense observation-by-Halo design. The tested scenes retain 10,346 OriginalInput rows, 211–300 free support entries, and 10,062–10,204 zero-sensitivity rows in the objective/replay. Those rows preserve the residual constant and original normalization.

The free width bound reuses qualified FirstStage Target widths: `Bmax = 4 Q95`. The eight qualified mixed scenes have 20/20 qualified width donors, Q95 from 1.5721 Å to 1.6608 Å, and Bmax from 6.2885 Å to 6.6431 Å. The estimator uses the existing FixedNeighbor block-coordinate and constrained A/C solve with Active-set LM log-width profiling, KKT/convergence checks, endpoint verification, and independent replay. It does not synthesize Halo FirstStage estimates.

## Parameter and signal results

The following paired aggregates use only the eight cases where both the qualified Target endpoint and FullABC endpoint are valid. Metrics are means across those eight scenarios; per-scenario values, worst cases, and unavailable cases remain in the raw JSON.

| Method | Halo A RMSE | Halo B RMSE | Halo C RMSE | Halo-only signal RMSE | Shared-support signal RMSE | Original-domain total prediction RMSE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Production geometry kNN-4 | 0.2570 | 0.1099 | 0.1233 | 0.01149 | 0.05741 (6 cases) | 0.01242 |
| Observation A/C, kNN width fixed | 0.2866 | 0.1021 | 0.01787 | 0.0003152 | 0.003081 (6 cases) | 0.002592 |
| Observation-driven FullABC | 0.3722 | 0.05373 | 0.01661 | 0.00009410 | 0.0005809 (6 cases) | 0.002563 |
| Oracle Halo truth, diagnostic only | 0 | 0 | 0 | 0 | 0 (6 cases) | 0.002565 mean |

FullABC lowers Halo-only signal RMSE in all eight paired cases, and lowers Original-domain total prediction RMSE to the noise/Target-endpoint floor. Parameter accuracy is mixed: A RMSE improves in six cases but its mean worsens because the separated high-noise case has A RMSE 1.413 and maximum absolute A error 3.150. B RMSE improves in six cases and C RMSE in seven. FullABC mean A/B/C biases are +0.1967, +0.0315, and +0.00286. Mean per-scenario median absolute errors are 0.0368, 0.00733, and 0.000592; mean per-scenario P95 absolute errors are 0.6735, 0.0972, and 0.0305. The largest observed absolute A/B/C errors are 3.150, 0.438, and 0.249. The large A outlier is retained in the result and is a material qualification risk.

The geometry kNN-4 mean biases are -0.0635, +0.0211, and +0.0160; mean per-scenario median absolute errors are 0.1460, 0.0318, and 0.00337; mean per-scenario P95 errors are 0.4132, 0.1923, and 0.2125. Its largest absolute A/B/C errors are 0.638, 0.428, and 0.381. Observation A/C has mean biases -0.1352, +0.0309, and +0.00379; mean per-scenario median absolute errors 0.0368, 0.00733, and 0.000592; mean per-scenario P95 errors 0.5068, 0.1862, and 0.0329; and largest absolute errors 1.102, 0.428, and 0.249. Observation A/C is scored for parameters and signal; it is not passed through a separate Peeling reconstruction contract.

## Peeling and downstream results

| Metric, eight paired cases | Geometry kNN-4 | FullABC | Paired outcome |
| --- | ---: | ---: | --- |
| Peeling response RMSE | 0.008194 mean | 0.002300 mean | FullABC lower in 8/8 |
| Tail Peeling Ratio MAE | 0.01373 mean | 0.00710 mean | FullABC lower in 6/8 |
| Peeling sample coverage | 32,000 / 32,000 | 32,000 / 32,000 | 100% both |
| Post-Joint Local A RMSE | 0.2807 mean | 0.2800 mean | FullABC lower in 4/8 |
| Post-Joint Local B RMSE | 0.09388 mean | 0.09295 mean | FullABC lower in 4/8 |
| Local qualification | 20/20 each case | 20/20 each case | complete both |
| Local effective samples | mean 99.8 per Target; minimum case minimum 90 | 100 per Target | FullABC consistent in these scenes |
| Local information | 20/20 available | 20/20 available | mean precision trace recorded in JSON |
| Group MDPDE A/B | 1.00173 / 1.10335 mean | 1.00053 / 1.10332 mean | estimates are close to truth (1 / 1.10) |
| Group posterior A/B RMSE | 0.26537 / 0.08830 mean | 0.26478 / 0.08749 mean | each lower in 4/8 |
| Group members / outliers | 20 / 0 each case | 20 / 0 each case | equal |
| Group alpha | 0 in all cases | 0 in all cases | equal |

Group prior and member uncertainty values and Group posterior differences are in `results.json`. Their mean differences are small and not consistently favorable. The synthetic evidence supports improved Peeling residuals; it does not establish a reliable Local or Group improvement. The fixed n=0 boundary Halo remains included through the FullABC Peeling context, with no loss of the 4,000 Peeling samples.

## Solver and resource evidence

Eight FullABC endpoints qualified. Runtime was 0.099–0.554 seconds per successful Halo solve, with estimated block/endpoint workspace 2,518,051–4,163,379 bytes. The two exact-3/4-observation stress cases did not return usable estimates. Their reports show failed A/C feasibility, conditional and independent A/C optimality, KKT consistency, independent replay, and endpoint verification; width stationarity and numerical qualification are unavailable because no trusted state exists. Active-set subproblem qualification is absent when that endpoint is missing. These statuses were not converted to fallback profiles. Rank deficiency was not separately measured in these scenarios.

The full ten-scene experiment took 112.20 seconds summed over per-case timers. Process peak RSS was 50,151,424 bytes (about 47.8 MiB); this is a small synthetic run and does not validate the 4 GB large-component limit. The independently estimated workspace is reported separately from process RSS. One repeated 4 Å scenario produced identical Halo profiles and all reported accuracy values across executions; elapsed-time fields naturally vary.

All ten First Target Joint endpoints qualified, and every downstream postprocessing replay measured Target A/B/C drift of exactly zero. The eight qualified FullABC scenes had no reported Halo solver evidence failure. The two free-Halo stress failures are retained separately from observation-count fixing: `n_h=3` and `n_h=4` were both correctly classified free before the solver ran.

Not measured here: a 5.94M-row/37,406-contributor component, actual large-component RSS/time, production-data FullABC Peeling/Local/Group comparison, SQL save/reload, UMAP, or Two Stage. The experiment does not establish SQLite schema compatibility because no new Halo persistence contract was promoted. P2 targeted tests separately cover OriginalInput identity, fixed Target and Halo subtraction, replay, no-Halo/all-fixed paths, distinct-row counts, chemical-label invariance, and reconstruction-only contributors.

## Reproduction

```sh
cmake --build build/debug --target joint_halo_fullabc_validation --parallel 4
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  build/debug/bin/joint_halo_fullabc_validation docs/research/results/halo_fullabc_p3
```

The run wrote 10 scenario records, 40 method rows, and 53 per-Halo rows. P3 validates synthetic physics and the research estimator only. It does not authorize Production promotion; the feasibility decision is recorded in the P4 report.
