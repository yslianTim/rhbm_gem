# Joint experiment and validation surface

The production estimator is FixedNeighbor with core size 12, forward serial
Gauss-Seidel, at most 200 complete sweeps, Profile LM, and at most one trusted
accepted update per block visit. Global numerical convergence determines
whether the Joint point is saved. Local conditioning and bounded parameter
stability provide separate point-quality evidence.

## Current reproducible tools

| Tool | Current question | Build option | Permanent owner |
| --- | --- | --- | --- |
| Joint benchmark | How do FixedNeighbor search, solve, workflow, and command paths use time and resources? | `RHBM_GEM_BUILD_BENCHMARKS` | `joint_benchmark_contract_test`, optional `joint_benchmark_smoke`, and Joint contract/workflow/numerical tests |
| Partial-selection measurement | What are the contributor-closure and halo measurements? | `RHBM_GEM_BUILD_BENCHMARKS` | PartialSelection tests |
| Statistical research | How do fixed-seed noise and position mismatch affect current outcomes? | `RHBM_GEM_BUILD_RESEARCH_TOOLS` | Statistical experiment contract and smoke tests |
| FullABC width study | What happens to local width profiles and derivative behavior in the saved study cases? | `RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS` | Local replay, precision, and width tests |

The benchmark workflow is `tests/integration/joint_benchmark.py`. It invokes
FixedNeighbor search and solve, the complete Joint/local/group workflow
benchmark, or the main command with explicit model and map inputs. It has no
custom core, reverse-order, or alternate solver option.

Current Python experiment drivers are `joint_benchmark.py`,
`joint_partial_selection.py`, and `joint_statistical_experiment.py`. Their C++
measurement programs live in `tests/experiments/`.

## Validation tools

| Tool | Role | Build option | Permanent owner |
| --- | --- | --- | --- |
| Joint runtime harness | Runs public API fixtures, regression cases, and physical-input smoke checks | `BUILD_TESTING` | `joint_component_runtime_runner_test`, `joint_component_regression`, `joint_component_physical_smoke` |
| Local precision controls | Compare local numerical operations with independent high-precision references | `RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS` | `joint_component_offline_tests` and local numerical tests |
| Joint cleanup contract | Keeps experiment and Python-driver directories within the current whitelist | `BUILD_TESTING` | `joint_experiment_cleanup_test` |
| Ordinary numerical and workflow tests | Verify convergence, local conditioning, stability, persistence, and downstream contracts | `BUILD_TESTING` | Joint contract, workflow, numerical, and persistence tests |

The Joint runtime harness is a validation tool, not a benchmark. Optional
precision controls are local numerical tests; they do not issue a production
rank certificate.

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

Global-rank and weak-direction research receipts remain historical. Their
thresholds, budgets, and rank outcomes do not describe current runtime behavior.
Detailed research conclusions and memory estimates are marked as retired in
[Joint component evidence](joint-component-evidence.md) and the
[Joint component memory census](joint-component-memory-census.md).
