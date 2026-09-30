# Joint capability and ownership guide

The estimator reports separate facts for command execution, numerical
qualification and scientific interpretation. A zero command exit code means the
result was saved; it does not mean the numerical checks passed. Target
convergence and full-parameter convergence remain distinct, especially when
halo parameters are weakly identifiable.

## Current Stage owners

| Stage | Question | Owner | Build category | Output meaning |
| --- | --- | --- | --- | --- |
| A | Why does this numerical problem behave this way? | [`joint_offline_diagnostic.py`](../../tests/integration/joint_offline_diagnostic.py) and `joint_offline_diagnostic` | OFFLINE | Per-start state, rank/spectrum, weak directions, precision and restart diagnosis; not a production acceptance gate |
| B | How does the estimator behave under controlled noise and position mismatch? | [`joint_statistical_experiment.py`](../../tests/integration/joint_statistical_experiment.py) and `joint_statistical_experiment` | RESEARCH | Fixed-seed scientific outcomes with process, numerical qualification and error summaries kept separate |
| C | Did the complete command path run, and what resources did it use? | [`joint_benchmark.py`](../../tests/integration/joint_benchmark.py) `--profile command` | BENCHMARK | Analysis, save/reload, JSON/CSV export and process/resource classification |

There is no combined Stage A/B/C runner. Stage A does not become a performance
benchmark when its cost is measured. Stage B is not a correctness test, and
Stage C resource completion does not replace persistence or numerical tests.

## Stage A: offline numerical diagnosis

The offline diagnostic keeps the actual returned state, target/halo roles,
rank and weak directions, precision comparison, and restart outcome. An
unavailable state remains unavailable; it is not converted to a generic
numerical failure. The [weak-halo diagnostic guide](joint-weak-halo-attribution.md)
is the current entry point and explains how to run and interpret the report.

## Stage B: statistical research

The statistical design, seeds, input definitions, metrics, failure denominator,
smoke and full-run commands belong to the [Stage B research guide](joint-noise-mismatch-validation.md).
Every condition remains in the receipt, including not-run and process failures;
the smoke is a tool check, not a statistical result.

## Stage C: complete-command benchmark

The [Joint benchmark guide](joint-benchmark.md) is the sole current entry point
for the command profile and its options. It keeps process/resource outcomes
separate from numerical qualification and scientific interpretation.

## Permanent correctness and evidence

Joint numerical correctness and persistence belong to C++ unit tests, offline
reference checks and the CLI workflow smoke. The benchmark measures those paths
without replacing their assertions. Refinement correctness has permanent owners
in `MDPDERegression_test.cpp` and `ProductionFitting_test.cpp`; its campaign
receipts do not own the production policy.

The historical weak-halo diagnosis, noise/mismatch results, resource campaign,
preflight, persistence comparisons, reports and provenance under
[`figures/joint-validation/`](figures/joint-validation/) remain intact. They are
historical evidence; active tools do not load those result files. Full Stage C
external-input measurements are unavailable unless the recorded model/map are
provided again. The limitations and measured results remain documented in the
[weak-halo diagnosis](joint-weak-halo-attribution.md),
[statistical experiment](joint-noise-mismatch-validation.md), and
[historical resource evidence](joint-component-evidence.md#6-historical-resource-envelope).
