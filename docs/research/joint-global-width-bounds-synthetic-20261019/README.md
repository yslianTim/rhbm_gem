# Global FirstStage Width Bounds: Synthetic Phase (P3a)

> **Archived workflow:** The paired A/B2/B4 executable and runner used for this study have been retired. The result tables and records below remain historical evidence; reproduction commands describe the former toolchain.

## Scope and outcome

The prescribed synthetic cases and seeds 1001–1005 were run with the existing Joint truth-controlled generator. The generator provides synthetic observations and scoring truth, but it does not provide FirstStage sampling observations, a FirstStage MDPDE fit, or a qualified FirstStage donor record. No contract-consistent FirstStage calibration could be produced from these inputs, so P95, B2, B4, and truth coverage under those bounds are unavailable. B2/B4 were not run.

For each case and seed, Group A used the existing unbounded Profile LM and the generator's original initialization. Synthetic truth was used only to score the completed fit. The 30 A runs completed without runner or audit errors. Qualification was 5/5 for well-conditioned, moderate, high-noise, and mixed-width; 2/5 for weak-identifiability; and 0/5 for A/C-correlation stress. Metrics for nonqualified runs are candidate metrics and are not formal Joint estimates.

| Synthetic case | A qualification | Mean full-population A RMSE | Mean full-population B RMSE (Å) | Mean full-population C RMSE | Mean prediction RMSE |
|---|---:|---:|---:|---:|---:|
| Well-conditioned | 5/5 | 0.00449704 | 0.000583062 | 0.000127657 | 0.000496954 |
| Moderate | 5/5 | 2.52155 | 0.208312 | 0.00995698 | 0.000496935 |
| Weak-identifiability | 2/5 | 31.1294 | 1.56187 | 0.0393754 | 9.95831e-6 |
| High-noise | 5/5 | 0.0451635 | 0.00581666 | 0.00128050 | 0.00496958 |
| Mixed-width | 5/5 | 2.89568 | 0.167897 | 0.0102334 | 0.000498853 |
| A/C-correlation stress | 0/5 | 1.91934e18 | 2349.39 | 4.34149e8 | 0.000497880 |

The RMSE columns average all five candidate outputs, including unqualified runs. They are descriptive only; the unqualified A/C-correlation candidates show why prediction fit alone cannot qualify the parameters. Per-seed status, metrics, total runtime, peak RSS, and canonical problem fingerprints are in [`paired-summary.json`](paired-summary.json).

## Reproduction

From the repository root, build the test-only runner and execute:

```sh
cmake --build build --target joint_global_width_bound_experiment --parallel 4
python3 tests/experiments/run_joint_global_width_bound_paired.py \
  --phase synthetic \
  --runner build/bin/joint_global_width_bound_experiment \
  --output /private/tmp/joint-global-width-bounds-synthetic-run \
  --model /Users/yslian/data/6z6u.cif \
  --simulation-map /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  --unsharpened-map /Users/yslian/data/emd_11103_unsharpened.map
```

The runner hashes the exact synthetic Joint problem for each seed. Case construction and the FirstStage limitation are implemented in `tests/experiments/joint_width_trust_experiment.cpp`; the paired orchestration and seed list are in `tests/experiments/run_joint_global_width_bound_paired.py`.
Use a new output path for each run; the runner refuses to overwrite prior results.

## Validation record

- Build: Release, AppleClang 21.0.0, OpenMP enabled.
- Build command: `cmake --build build --target joint_global_width_bound_experiment joint_width_trust_experiment --parallel 4`.
- Focused regression: `ctest --test-dir build -R 'rhbm_tests_joint_component_(numerical|contract|workflow)|joint_component_regression' --output-on-failure` — 4/4 passed.
- Paired-audit and experiment-cleanup contracts: `ctest --test-dir build -R 'joint_global_width_bound_experiment_test|joint_experiment_cleanup_test' --output-on-failure` — 2/2 passed.
- Python syntax check: `python3 -m py_compile tests/experiments/run_joint_global_width_bound_paired.py tests/integration/joint_global_width_bound_experiment_test.py` — passed.
- Synthetic A campaign: 30/30 completed; calibration records explicitly say unavailable, and no B2/B4 run was launched.

## Promotion assessment for P3a

This phase supplies an unbounded synthetic baseline and a reproducible reason that the requested bounded comparison is unavailable. It does not support or reject B2/B4. A synthetic FirstStage sampling path must be established before these workloads can contribute evidence for global-bound promotion.
