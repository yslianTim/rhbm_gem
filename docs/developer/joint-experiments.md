# Joint experiment and validation surface

The production estimator is FixedNeighbor with core size 12, forward serial
Gauss-Seidel, at most 200 complete sweeps, and at most one trusted accepted
bounded Active-set LM update per block visit for the Pure Exclusion policy.
Its global width ceiling is four times the qualified FirstStage target-width
P95. Global numerical convergence determines whether the Joint point is saved.
Local conditioning and bounded parameter stability provide separate
point-quality evidence. Pure Halo-Support Exclusion is the only current Joint
halo policy; FixedHalo is a retired historical method.
Profile LM is the bounded local profiled operation inside Active-set LM; the
unbounded FixedHalo route has been retired.

## Current reproducible tools

| Tool | Current question | Build option | Permanent owner |
| --- | --- | --- | --- |
| Joint benchmark | How do FixedNeighbor search, solve, workflow, and command paths use time and resources? | `RHBM_GEM_BUILD_BENCHMARKS` | `joint_benchmark_contract_test`, optional `joint_benchmark_smoke`, and Joint contract/workflow/numerical tests |
| Partial-selection measurement | What are the contributor-closure and halo measurements? | `RHBM_GEM_BUILD_BENCHMARKS` | PartialSelection tests |
| Bound-active Production synthetic | Does the calibrated Active-set LM satisfy known-truth bound, KKT, replay, and qualification checks? | `BUILD_TESTING` | `joint_production_bounded_synthetic` |
| Numerical oracles | Do retained local Profile, A/C, replay, and FullABC precision cases agree with independent references? | `RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS` | Local replay, precision, and width tests |

The benchmark workflow is `tests/integration/joint_benchmark.py`. It invokes
FixedNeighbor search and solve, the complete Joint/local/group workflow
benchmark, or the main command with explicit model and map inputs. It has no
custom core, reverse-order, or alternate solver option.

Current Python experiment drivers are `joint_benchmark.py` and
`joint_partial_selection.py`. Their C++ measurement programs live in
`tests/experiments/`.

## Validation tools

| Tool | Role | Build option | Permanent owner |
| --- | --- | --- | --- |
| Joint runtime harness | Runs public API fixtures, regression cases, and physical-input smoke checks | `BUILD_TESTING` | `joint_component_runtime_runner_test`, `joint_component_regression`, `joint_component_physical_smoke` |
| Production synthetic | Exercises the Global FirstStage P95 × 4 bound with truth-controlled active widths | `BUILD_TESTING` | `joint_production_bounded_synthetic` |
| Local precision controls | Compare local numerical operations with independent high-precision references | `RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS` | `joint_component_offline_tests` and local numerical tests |
| Joint cleanup contract | Keeps experiment and Python-driver directories within the current whitelist | `BUILD_TESTING` | `joint_experiment_cleanup_test` |
| Ordinary numerical and workflow tests | Verify convergence, local conditioning, stability, persistence, and downstream contracts | `BUILD_TESTING` | Joint contract, workflow, numerical, and persistence tests |

The Joint runtime harness is a validation tool, not a benchmark. Optional
precision controls are local numerical tests; they do not issue a production
rank certificate.

FixedHalo estimation tests and tools have been retired. The JSON/SQLite codec
still reads and round-trips historical FixedHalo results without exposing an
estimator path. The saved parameterization and solver evidence retain their
historical meaning.
The fixed-seed statistical study and its saved outcomes are historical; its old
fitting campaign used an earlier solver path and is not evidence for the
current Active-set LM. Only its deterministic noise generator contract remains
executable.

## Run current Production regressions

Configure and build the ordinary test surface:

```sh
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --parallel 4
```

Run the core Joint contract, workflow, numerical regressions, and the
bound-active synthetic:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 ctest --test-dir build -R 'rhbm_tests_joint_component_(contract|workflow|numerical)|joint_production_bounded_synthetic' --parallel 2 --output-on-failure
```

Enable independent high-precision and replay fixtures with
`-DRHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON`. Enable the supported benchmark
executables with `-DRHBM_GEM_BUILD_BENCHMARKS=ON`; both options require
`BUILD_TESTING=ON`. No retired research runner is needed for ordinary
Production regression.

## Historical decisions

| Topic | Conclusion | Current owner |
| --- | --- | --- |
| Fixed-B | Block updates and replay accounting were feasible, but Fixed-B was not promoted; FixedNeighbor superseded it. | FixedNeighbor tests |
| Core size | Repeated experiments selected core 12. | Production policy contract |
| Block order | Forward was adopted for serial Gauss-Seidel. | Production policy contract |
| Local update policy | OneAccepted is intrinsic to each block visit. | FixedNeighbor numerical tests |
| Sparse backend | EIGEN became the sole backend after bounded parity and resource comparisons. | Local solver tests and historical persistence tests |
| OperatorPcg and Schwarz | Investigated for large components but not adopted. | Historical research record |
| Global rank qualification | Removed from production, output, and resource planning. | Convergence, conditioning, and stability tests |
| Global width bound | Global FirstStage P95 × 4 selected for Production; A/B comparison runner retired. | Width calibration and bound-active Production regression |
| Halo policy comparisons | Pure Halo-Support Exclusion selected; Dynamic Halo and FixedHalo estimation retired. | Pure Exclusion tests and legacy persistence tests |

Global-rank and weak-direction research receipts remain historical. Their
thresholds, budgets, and rank outcomes do not describe current runtime behavior.
Detailed research conclusions and memory estimates are marked as retired in
[Joint component evidence](joint-component-evidence.md) and the
[Joint component memory census](joint-component-memory-census.md).
