# Global P95 × 4 bounded LM: production candidate validation

> **Historical validation record:** This report predates the production promotion. The A/B2/B4 comparison test and runner referenced by its commands are retired; saved measurements remain historical evidence. See [the production promotion record](joint-bounded-active-set-production-promotion-20261020.md) for the current route.

## Decision

The Release, non-instrumented bounded candidate passed the checks below on the
6Z6U simulation, 6Z6U unsharpened map, and a truth-controlled synthetic fixture.
The public `EstimateJointComponents` default remains unchanged. The 37,406-target
bounded memory run and a bounded end-to-end CLI/downstream workflow were not
completed, so the production promotion gate remains open.

## Contract and source-level path

The candidate uses only qualified FirstStage widths from selected Joint target
atoms. It requires at least 20 qualified donors and at least 95% qualified-target
coverage, computes the deterministic Hyndman–Fan Type 7 P95, and uses the fixed
internal multiplier `k=4`. Invalid or unavailable calibration returns an
unavailable result; this path has no unbounded fallback. FirstStage estimates
are not modified; only the bounded Joint working initialization is projected.

The implementation responsibilities are:

1. `WidthCalibration` selects donors, performs the preflight, computes P95 and
   builds a feasible working initialization.
2. `JointComponentEstimator` passes the same global cap into FixedNeighbor
   orchestration and records calibration and projection evidence.
3. `Search` and `ActiveSetBoundedLM` solve the feasible damped subproblem and
   evaluate the trust ratio for the bounded step actually taken.
4. `Numerics`, component certification and assembly check A/C KKT, Width
   feasibility and projected Width KKT. Independent replay calculates its own
   gradients and checks prediction, residual, objective and bound consistency.
5. Persistence reuses existing Joint evidence fields; no schema or CLI option
   was added.

The bounded path is available in the Release library without
`RHBM_GEM_TEST_INSTRUMENTATION`. It is an internal candidate entry point at this
stage; the public workflow still calls the legacy unbounded path pending the
remaining validation gates.

## Release candidate results

All runs linked to the Release `rhbm_gem` target with instrumentation disabled.
OpenMP and BLAS inner threads were set to one for each run.

| Metric | 6Z6U simulation | 6Z6U unsharpened |
|---|---:|---:|
| Targets / observations | 1,418 / 180,478 | 1,418 / 189,377 |
| Qualified donors | 1,415 (99.788%) | 1,405 (99.083%) |
| P95 (Å) | 0.5689788338510722 | 0.8021089204681596 |
| Bmax = 4×P95 (Å) | 2.2759153354042887 | 3.2084356818726385 |
| Initial projection count / maximum (Å) | 0 / 0 | 1 / 1.606375 |
| Active bound atoms | 0 | 2 |
| Maximum log-width bound violation | 0 | 0 |
| Search sweeps / profile evaluations | 83 / 22,908 | 124 / 34,848 |
| Stop reason | `block-stationary` | `block-stationary` |
| Active-set solves / failures | 11,454 / 0 | 17,488 / 0 |
| Maximum active-set subproblem KKT | 1.96512e-19 | 1.62630e-19 |
| Projected Width KKT (threshold 1e-9) | 8.77651e-10 | 9.40273e-10 |
| Conditional / independent A/C KKT (threshold 1e-7) | 9.43457e-8 / 9.43457e-8 | 8.17416e-8 / 8.17416e-8 |
| Independent KKT consistency (threshold 1e-13) | 1.73472e-18 | 1.73472e-18 |
| Assembled consistency (threshold 2e-12) | 1.06581e-14 | 1.42109e-14 |
| Independent replay / RuntimeConvergence | Passed / Passed | Passed / Passed |
| Objective | 1.2330854370811315e-13 | 0.018417539368903946 |
| Prediction RMSE | 9.20650e-7 | See paired research record below |
| Peak RSS | 621,314,048 B | 1,403,518,976 B |
| Wall time | 304.70 s | 435.89 s |

The simulation's bound was inactive, truth coverage was 1,418/1,418, and its
calibration, serial-1032 width, objective, RMSEs and projected KKT match the
previous B4 research result within floating-point variation. In normalized
FirstStage truth units, A/B/C RMSE were 9.50818e-6 / 1.46475e-7 Å / 2.92104e-6.
The candidate returned serial-1032 width 0.500000003740957 Å.

The unsharpened map has no truth. The candidate returned serial-1032 at its
3.2084356818726385 Å cap; its objective, active atoms, KKT, sweeps and profile
evaluations match research B4. The research paired record reports prediction
RMSE 0.9260391 for B4; this is a residual description, not a parameter-accuracy
claim. Its A baseline was unqualified and is not treated as an estimate.
Against the research B4 record, the P95 and cap are exactly the same; the
serial-1032 width, two final active atoms, profile evaluations, sweeps,
projected KKT and replay status match. The objective differs by
1.43e-16 (research 0.018417539368903803; production
0.018417539368903946). The retained paired report does not store the full
per-atom Width vector, so vector-wise equality is not claimed.

The production candidate also retains point-level quality diagnostics. Both
maps had adequate local conditioning for all 1,418 targets. Parameter-stability
coverage was limited by the existing bounded perturbation policy: the simulation
assessed one of 138 blocks (6 solves, 24 evaluations), yielding 12 `Stable`
points and 1,406 `NotAssessed` points. On the unsharpened map, 139 of 140 blocks
were `NotAssessed`, one was `Unavailable`, and no perturbation solve ran; 12
points were consequently `Warning` and 1,406 `NotAssessed`. These statuses are
preserved as diagnostics; `RuntimeConvergence=Passed` does not imply every
target's parameter stability was assessed. This limited coverage remains part
of the promotion evidence and risk. The unavailable block reports
`baseline-local-kkt-failed`; the other 139 blocks were not probed because of the
representative-probe limit.

## Truth-controlled synthetic case

The fixture creates a 25-target map with known A/B/C, then obtains donor widths
through the real FirstStage sampling and MDPDE initialization. Calibration only
receives that FirstStage initialization. Results:

- Qualified donors: 25/25; P95 = 0.624761 Å; Bmax = 2.49904 Å.
- Truth coverage: 24/25 (96%); one 5 Å truth lies beyond the bound and its
  fitted Width reaches the bound. Maximum bound violation is zero.
- Active-set KKT = 1.51788e-18; projected Width KKT = 2.75557e-13;
  A/C KKT = 1.66533e-16; bound feasibility, replay and runtime convergence pass.
- A/B/C RMSE in the fixture's original map domain: 0.4 / 0.500191 Å / 0.0105113.
  The single broad atom accounts for the visible parameter bias.
- Prediction RMSE = 0.00217477 versus a truth-model residual RMSE of 0.00013511.
  The prediction error is 0.079684% of truth-signal RMS (2.72924). This records
  the bound-induced bias; it does not establish better parameter recovery.

The 6Z6U simulation is the bound-inactive control. The full synthetic campaign
was not rerun.

## Tests and remaining gates

The Release build used `cmake --build ... --parallel 4`. The focused CTest group
passed 7/7, including active-set LM, Width calibration, component numerical and
workflow tests, regressions, cleanup checks and this synthetic case. Failure
tests cover unavailable/invalid calibration, low donor populations, bounded
subproblem rank deficiency and infeasibility, projected KKT failure, replay
disagreement, recovery exhaustion and sweep budget exhaustion. A separate
regression now requires explicit `width-bound-feasibility` evidence; both its
pass and fail outcomes are tested.

Not verified in this candidate run:

- The 5,938,087-observation / 37,406-target bounded workload and its <=4 GB RSS
  target. Existing preparation-only memory figures do not establish this gate.
- The public CLI, result persistence through the full bounded workflow, or
  PostFitPeelingSamples/Local MDPDE/Group Fitting after a bounded Joint fit.
- A complete failure-injection campaign under the candidate Release path.

These gaps keep the default promotion gate closed.

## Reproduction commands

Build in the configured Release tree with no more than four compile jobs:

```sh
cmake -S . -B build
cmake --build build --target joint_width_repro_no_instrumentation joint_production_bounded_synthetic rhbm_tests --parallel 4
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 ctest --test-dir build -R 'rhbm_tests_joint_component_(numerical|contract|workflow)|joint_component_regression|joint_global_width_bound_experiment_test|joint_experiment_cleanup_test|joint_production_bounded_synthetic' --parallel 4 --output-on-failure
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 build/bin/joint_production_bounded_synthetic
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 build/bin/joint_width_repro_no_instrumentation /Users/yslian/data/6z6u.cif /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map /private/tmp/joint-production-bounded-simulation-quality-final
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 build/bin/joint_width_repro_no_instrumentation /Users/yslian/data/6z6u.cif /Users/yslian/data/emd_11103_unsharpened.map /private/tmp/joint-production-bounded-unsharpened-reasons-final
```

The two map/model inputs are local validation data. Machine-readable run
summaries and fingerprints were written under `/private/tmp`; they are not
committed to the repository.
