# Global FirstStage Width Bounds: Real-map Phase (P3b)

## Result

The paired audit passes for both the 6Z6U simulation map and the unsharpened EMDB map. The simulation pair passed in the initial run; after a test-only evidence-context correction, the unsharpened pair also passed. The simulation comparison had no active width bounds, so it confirms compatibility but does not show a bounded-search effect. On the unsharpened map, B2 and B4 stop the serial-1032 width runaway at their calibrated caps and pass independent replay and projected-width stationarity. The unsharpened A baseline remains an unqualified candidate with a failed replay.

No simulation truth was used to calibrate the bounds. The unsharpened map has no truth; its residual and objective values are descriptive fit metrics, not parameter-accuracy claims.

## Calibration and provenance

Calibration uses qualified widths returned by `FirstStage::InitializeContributors` for the same Chain A non-hydrogen target selection used in Joint. The P95 convention is Hyndman–Fan Type 7. Fallback and unqualified donor values are excluded.

| Dataset | Targets | Qualified donors | P95 (Å) | B2 = 2×P95 (Å) | B4 = 4×P95 (Å) | Initial violations B2/B4 |
|---|---:|---:|---:|---:|---:|---:|
| 6Z6U simulation | 1,418 | 1,415 | 0.568979 | 1.137958 | 2.275915 | 0 / 0 |
| 6Z6U unsharpened | 1,418 | 1,405 | 0.802109 | 1.604218 | 3.208436 | 2 / 1 |

The unsharpened FirstStage initialization had 1,414 native successes, 12 invalid-width fallbacks, four unqualified outputs, and no initialization failures. None of the 12 fallback widths entered calibration. The simulation initialization had 1,415 native successes, three unqualified outputs, and no fallback or invalid widths.

| Input | Model SHA-256 | Map SHA-256 | Problem fingerprint | FirstStage initialization fingerprint |
|---|---|---|---|---|
| 6Z6U simulation | `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa` | `6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a` | `b9c6816ebdb612d14a57dbee06c3b07680ad4ecfd369b4f9265db5487c46eae4` | `a4260ff5097fb6651ec780f8f19661cd28843b36b3332b05a4f6a8510b3463d3` |
| 6Z6U unsharpened | `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa` | `fa9a6b15ad86d471130fc83495754cbf0a2e4084adc5327bf09f21d7e5ec60d1` | `b111feb20669e8d703d546f81ea2ad313341aed604f4fb271c0396150c7ff344` | `46cb8a262102775ff6bfcc1a050061164cd18cce47a60ea59dcb7f65f06dc5c6` |

Each formal A/B2/B4 group used the same input fingerprint, FirstStage initialization fingerprint, calibration record, observation scale, chain selection, and model/map hashes. Simulation truth was read only after calibration and used only to score completed estimates (`truth_leakage=false`).

## Simulation-map comparison

All three groups passed runtime convergence, qualification, independent replay, and width stationarity. B2 and B4 did not make a bound active; all groups returned serial-1032 width 0.500000 Å. Truth coverage under both bounds was 100% (1,418/1,418).

| Group | Objective | Prediction RMSE | A RMSE | B RMSE (Å) | C RMSE | Active bound atoms |
|---|---:|---:|---:|---:|---:|---:|
| A | 1.2333160e-13 | 9.20736e-7 | 9.50872e-6 | 1.46484e-7 | 2.92110e-6 | 0 |
| B2 | 1.2330854e-13 | 9.20650e-7 | 9.50818e-6 | 1.46475e-7 | 2.92104e-6 | 0 |
| B4 | 1.2330854e-13 | 9.20650e-7 | 9.50818e-6 | 1.46475e-7 | 2.92104e-6 | 0 |

The complete per-seed and paired audit records are in [`6z6u-simulation/`](6z6u-simulation/). The small metric differences are descriptive; because neither bound was active, this map does not establish an effect from the cap.

## Unsharpened-map comparison

The initial unbounded A run stopped after 70 sweeps with `block-recovery-exhausted`, serial-1032 width 6,055.801 Å, and failed independent replay. It is an unqualified candidate. The bounded groups completed 124 sweeps and passed global runtime qualification, independent replay, bound feasibility, and projected-width stationarity.

| Group | Serial-1032 width (Å) | Bound (Å) | Active atoms | Objective | Prediction RMSE | Projected width KKT | Runtime / replay |
|---|---:|---:|---:|---:|---:|---:|---|
| A | 6,055.801 | unbounded | 0 | 0.018417515 | 0.9260385 | unavailable | unavailable / failed |
| B2 | 1.604218 | 1.604218 | 3 | 0.018417630 | 0.9260414 | 9.4008e-10 | passed / passed |
| B4 | 3.208436 | 3.208436 | 2 | 0.018417539 | 0.9260391 | 9.4027e-10 | passed / passed |

Both bounded groups had zero bound violation, zero boundary-gradient sign violation, zero complementarity residual, and zero active-set subproblem failures. The maximum active-set subproblem KKT residual was 1.6264e-19. Their raw width-gradient norms (1.9381e-7 for B2 and 4.7538e-8 for B4) exceed the 1e-9 threshold, while independent projected KKT values are below it. The evidence therefore uses the bound-aware projected criterion; the threshold was not changed.

The bounded fits did not materially change the observed-data objective or residual RMSE from the unqualified A candidate. There is no truth for this map, so this comparison supports the width-cap behavior and numerical qualification only; it cannot establish parameter accuracy.

The formal corrected outputs and all passed audit checks are in [`recheck-unsharpened/`](recheck-unsharpened/). The first bounded attempt is preserved in [`6z6u-unsharpened/`](6z6u-unsharpened/); its audit failed only the B2/B4 bound-aware stationarity checks. That exposed a test-only propagation gap: search and component certification had the research cap, but final assembled evidence replayed with the unbounded context. The test-only policy now carries the same cap into the assembled assessment. A regression test covers this path. The A/B2/B4 numerical solutions and acceptance thresholds were not changed.

### Projected-initialization controls

The auxiliary Aproj2/Aproj4 runs apply the initial B2/B4 projection and then use the unbounded A search. They are controls, not formal bounded groups. Aproj2 ended at 4,588.729 Å after recovery exhaustion with failed replay; Aproj4 ended at 2,494.618 Å after the 200-sweep budget with failed replay. Initial projection alone did not prevent subsequent width runaway.

## Synthetic phase limitation

The 30 synthetic P3a runs remain an unbounded A baseline. The generator does not provide FirstStage sampling observations, an MDPDE donor fit, or qualified donor provenance; B2/B4 and synthetic truth coverage under those bounds are therefore unavailable. See the [synthetic-phase report](../joint-global-width-bounds-synthetic-20261019/README.md). This is the main limitation on promotion evidence.

## Reproduction and validation

Use a fresh output path because the runner refuses to overwrite prior results.

```sh
cmake --build build --target joint_global_width_bound_experiment rhbm_tests --parallel 4
cmake --build build --target rhbm_gem --parallel 4
ctest --test-dir build -R 'rhbm_tests_joint_component_(numerical|contract|workflow)|joint_component_regression|joint_global_width_bound_experiment_test|joint_experiment_cleanup_test' --output-on-failure
python3 tests/experiments/run_joint_global_width_bound_paired.py \
  --phase real \
  --runner build/bin/joint_global_width_bound_experiment \
  --output /private/tmp/joint-global-width-bounds-real-reproduction \
  --model /Users/yslian/data/6z6u.cif \
  --simulation-map /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  --unsharpened-map /Users/yslian/data/emd_11103_unsharpened.map
```

The bounded-assembly regression and all six focused CTest entries passed. The production `rhbm_gem` target also built with the test-only research control compiled out. Final machine-readable comparisons are in [`paired-summary-final.json`](paired-summary-final.json); the earlier audit failure remains available for review at [`6z6u-unsharpened/pair-audit.json`](6z6u-unsharpened/pair-audit.json).
