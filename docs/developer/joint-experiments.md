# Joint experiment and validation surface

This inventory describes the current checkout. The production estimator is
FixedNeighbor with core 12, Forward serial Gauss-Seidel, at most 30 sweeps,
local Profile LM and at most one trusted accepted local update per block
visit.

## Current reproducible tools

| Tool | Current question | Build option | Permanent owner |
| --- | --- | --- | --- |
| Joint benchmark | How do current FixedNeighbor search, solve, workflow, postprocessing and command paths perform and use resources? | RHBM_GEM_BUILD_BENCHMARKS | joint_benchmark_contract_test, joint_benchmark_smoke, Joint contract/workflow/numerical tests |
| Partial-selection measurement | What are the current selection, contributor-closure and halo measurement results? | RHBM_GEM_BUILD_BENCHMARKS | PartialSelection_test.cpp and Joint workflow/numerical tests |
| Statistical research | How do fixed-seed noise and position mismatch affect the current estimator outcomes? | RHBM_GEM_BUILD_RESEARCH_TOOLS | joint_statistical_experiment_contract_test and joint_statistical_experiment_smoke |
| Offline numerical diagnosis | What does the current weak-halo and local numerical evidence show on immutable inputs? | RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS | joint_offline_diagnostic_smoke, joint_component_evidence_test, joint_component_offline_tests and joint_component_offline_regression |

The single benchmark workflow is tests/integration/joint_benchmark.py. It
calls the production FixedNeighbor benchmark for search and solve, the
postprocessing benchmark for workflow and postprocess, and the main command
for complete-command measurements. The benchmark policy has no custom core,
reverse order or alternative solver mode.

Current Python experiment drivers are joint_benchmark.py,
joint_partial_selection.py, joint_statistical_experiment.py and
joint_offline_diagnostic.py. Their C++ measurement programs live in
tests/experiments/. No other executable source in that directory is current.

The weak-halo diagnosis and noise/mismatch research retain machine-readable
results in docs/developer/figures/joint-validation/. Both answer current
diagnosis or scientific questions; other closed campaign results are not
retained as current-tree machine-readable evidence.

## Validation tools

| Tool | Role | Build option | Permanent owner |
| --- | --- | --- | --- |
| Joint runtime harness | Runs current public-API fixtures, regression cases and physical-input smoke checks | BUILD_TESTING | joint_component_runtime_runner_test, joint_component_regression and joint_component_physical_smoke |
| Joint offline audit | Produces current offline certification bundles for audit workflows | RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS | joint_component_offline_tests and joint_component_offline_regression |
| Joint cleanup contract | Keeps experiment and Python-driver directories within the current whitelist | BUILD_TESTING | joint_experiment_cleanup_test |
| Ordinary numerical and workflow tests | Permanently verify the frozen production contracts and current application paths | BUILD_TESTING | rhbm_tests_joint_component_contract, rhbm_tests_joint_component_workflow and rhbm_tests_joint_component_numerical |

The C++ runtime and audit harness sources live in tests/tools/. Python runtime,
audit and workflow smoke entry points are validation tools, not experiment
drivers. Files ending in _test.py and permanent C++ regression sources remain
ordinary tests; the experiment whitelist does not exclude them.

## Historical decisions

| Topic | Conclusion | Current owner |
| --- | --- | --- |
| Fixed-B | Block updates and replay accounting were feasible, but Fixed-B was not promoted; FixedNeighbor superseded it. | FixedNeighbor production tests |
| Core size | Repeated qualification selected core 12; larger tested cores did not justify their added local problem size. | Production policy contract |
| Block order | Forward was adopted as the serial Gauss-Seidel order. Reverse was a qualification diagnostic only. | Production policy contract |
| Local update policy | OneAccepted was selected and is intrinsic to the production block visit. | FixedNeighbor numerical tests |
| Sparse backend | EIGEN and SPQR passed bounded parity; SPQR was faster on completed 128/256 probes and EIGEN used less peak RSS. Both bounded 512 probes timed out. A later memory, dependency and maintenance policy selected EIGEN. | EIGEN numerical tests and historical persistence compatibility tests |
| OperatorPcg and Schwarz | Investigated for large connected components but not adopted. | FixedNeighbor with local Profile LM and EIGEN |
| Top-level LegacyCompact | The former global route was retired; the local FixedNeighbor primitive is named Profile LM. | FixedNeighbor production implementation |
| Qualification campaigns | Completed core, order, stationarity, prepared-block and local-route decisions are closed. Permanent tests own the frozen contracts. | Joint contract and numerical tests |

The detailed historical decisions are in
[joint-fixed-neighbor-experimental.md](joint-fixed-neighbor-experimental.md),
[joint-operator-search.md](joint-operator-search.md) and
[joint-component-evidence.md](joint-component-evidence.md). Historical
compatibility decoders preserve old result metadata without restoring retired
solver implementations or experiment modes.

Closed experiment implementations and machine-readable receipts are recoverable
from Git history and are intentionally not retained in the current tree. The
current tree keeps executable experiments only for current production
measurement, current feature measurement, current scientific research or
current offline diagnosis.
