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

The diagnostic replays the existing weak-halo fixture from its production start
and two fixed alternative starts. It retains the returned state, target/halo
parameter roles, rank and weak-direction evidence, independent precision
comparisons, local correction scans and restart outcome. An unavailable returned
state stays unavailable; collection completeness is reported separately.

```sh
cmake -S . -B build/offline -DBUILD_TESTING=ON \
  -DRHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON
cmake --build build/offline --target joint_offline_diagnostic
python3 tests/integration/joint_offline_diagnostic.py \
  --work-dir build/joint-diagnostic/run-01 \
  --executable build/offline/bin/joint_offline_diagnostic \
  --output build/joint-diagnostic/report/weak-halo.json
```

The small offline CTest smoke runs these three deterministic starts. The
historical Stage A report remains in
[`figures/joint-validation/weak-halo.json`](figures/joint-validation/weak-halo.json)
and [`weak-snapshot.json`](figures/joint-validation/weak-snapshot.json).

## Stage B: statistical research

The fixed design has two geometries, 20 PCG64 seeds per noisy condition, paired
initialization modes where specified, and the existing IID/correlated noise and
position-mismatch definitions. The smoke uses the noiseless control and first
seeded IID case and is marked `scope: smoke`; only a full run covers all 326
generated inputs and 450 fitted outcomes.

```sh
cmake -S . -B build/research -DBUILD_TESTING=ON \
  -DRHBM_GEM_BUILD_RESEARCH_TOOLS=ON
cmake --build build/research --target joint_statistical_experiment
python3 tests/integration/joint_statistical_experiment.py \
  --smoke --work-dir build/joint-statistical/smoke \
  --executable build/research/bin/joint_statistical_experiment
python3 tests/integration/joint_statistical_experiment.py \
  --work-dir build/joint-statistical/full \
  --executable build/research/bin/joint_statistical_experiment
```

Every condition remains in the receipt, including not-run and process failures.
Summary rows name attempted, completed, qualified, failed and unavailable counts;
error statistics retain their denominators and separate all-available from
converged-only estimates. The research smoke is a tool check, not a statistical
result.

## Stage C: complete-command benchmark

Use the unified command profile described in the [Joint benchmark guide](joint-benchmark.md):

```sh
python3 tests/integration/joint_benchmark.py \
  --profile command --case sample --build-dir build/debug \
  --model input.cif --map input.map --output build/joint-command.json
```

The profile records command completion, database persistence/reload, JSON/CSV
export, process error, timeout, RSS limit and not-run outcomes. Numerical
convergence and available values are reported independently. The built-in
`joint_benchmark_smoke` uses generated local data and checks this path without
running the historical resource campaign. Timeout and RSS limits can be supplied
to the benchmark command; old Stage C thresholds remain historical evidence.

## Permanent correctness and evidence

Joint numerical correctness and persistence belong to C++ unit tests, offline
reference checks and the CLI workflow smoke. The benchmark measures those paths
without replacing their assertions. Refinement correctness has permanent owners
in `MDPDEExperiment_test.cpp` and `ProductionFitting_test.cpp`; its campaign
receipts do not own the production policy.

The historical weak-halo diagnosis, noise/mismatch results, resource campaign,
preflight, persistence comparisons, reports and provenance under
[`figures/joint-validation/`](figures/joint-validation/) remain intact. They are
historical evidence; active tools do not load those result files. Full Stage C
external-input measurements are unavailable unless the recorded model/map are
provided again. The limitations and measured results remain documented in the
[weak-halo diagnosis](joint-weak-halo-attribution.md),
[statistical experiment](joint-noise-mismatch-validation.md), and
[historical command resource report](joint-command-resource-envelope.md).
