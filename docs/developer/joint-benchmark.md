# Joint benchmark

The current benchmark surface measures the qualified `FixedNeighbor` Joint
estimator. It is intentionally separate from the historical OperatorPcg and
Schwarz campaigns; those drivers are retired and their compact evidence is
kept under `docs/developer/figures/`.

## Current entry point

Use [`tests/integration/joint_benchmark.py`](/tests/integration/joint_benchmark.py)
with a benchmark build (`RHBM_GEM_BUILD_BENCHMARKS=ON`). The active profiles are:

| Profile | Scope | Driver |
| --- | --- | --- |
| `search` | FixedNeighbor search-only work, including local LegacyCompact profile work | `joint_fixed_neighbor_experiment --scaling-only` |
| `solve` | FixedNeighbor search, global replay, endpoint certification, and `RuntimeConvergence` | `joint_fixed_neighbor_experiment --case` |
| `workflow` | Joint workflow, postprocess, and persistence | `joint_postprocessing_benchmark` |
| `postprocess` | Joint postprocess, uncertainty, and persistence | `joint_postprocessing_benchmark` |
| `command` | CLI analysis, SQLite reload, and Joint export | `RHBM-GEM` |

For example:

```bash
python3 tests/integration/joint_benchmark.py \
  --profile search --case chain-8 \
  --build-dir build/qualification --output build/joint-search.json

python3 tests/integration/joint_benchmark.py \
  --profile solve --case chain-8 \
  --build-dir build/qualification --output build/joint-solve.json
```

The `search` profile is an attribution measurement and must not claim endpoint
qualification. Use `solve` when the result must include endpoint certification
and `RuntimeConvergence`. `workflow`, `postprocess`, and `command` are smoke
and persistence surfaces rather than substitutes for numerical qualification.

## FixedNeighbor production metadata

The default production policy is recorded in each report:

```text
search_method = FixedNeighbor
fixed_neighbor_core_atoms = 12
fixed_neighbor_block_order = forward
fixed_neighbor_maximum_sweeps = 30
fixed_neighbor_local_search = LegacyCompact
fixed_neighbor_policy = production
sparse_backend = EIGEN
```

FixedNeighbor visits structural blocks in Forward order. Each visit runs the
existing local LegacyCompact profile solver, accepts at most one trusted local
update, replays the candidate globally, and preserves the existing stationarity,
endpoint-certification, and `RuntimeConvergence` contracts. The benchmark does
not introduce a second numerical policy or an automatic route selector.

`--fixed-core-atoms` is available for search-only custom attribution cases.
The `solve` profile is restricted to the production core of 12 atoms. Custom
core sizes are not production defaults and do not change the estimator contract.

The wrapper records process wall time, sampled RSS, timeout/process status,
build metadata, source provenance, and the normalized numerical result. Use
`--repeat`, `--warmup`, `--timeout`, and the RSS options for repeatable small
measurements; generated output belongs in the build or output tree rather than
in the source tree.

The sparse backend is fixed to EIGEN in current builds. There is no backend
command-line option or CMake selector. SPQR appears only in the retained
historical qualification evidence and legacy provenance decoder.

## Current option boundary

The current benchmark deliberately has no options for an operator route,
preconditioner, Schwarz geometry, operator rank, or operator-factor ownership.
`joint_sparse_benchmark`, `joint_route_frontier.py`, and
`joint_schwarz_sweep.py` are retired. A request using their old route or option
names is rejected by the current parser or command catalog; it is not silently
mapped to FixedNeighbor.

The top-level `LegacyCompact` route is retired and is not a current benchmark
profile. The `LegacyCompact` numerical implementation itself remains active as
FixedNeighbor's local profile solver.

## Historical evidence

The following are closed investigations, not current benchmark instructions:

- global OperatorPcg/PCG search and preconditioner comparisons;
- Schwarz geometry and scaling campaigns;
- bounded-rank and operator-factor ownership experiments; and
- route-frontier comparisons between LegacyCompact, OperatorPcg, and
  FixedNeighbor.

Their compact summaries, manifests, and figures remain available for scientific
provenance. Raw repetitions and process telemetry are not required by the
current regression surface and are not regenerated as part of route cleanup.
Historical saved results remain readable through the legacy provenance decoder;
the old operator fields are compatibility data, not current solver controls.

## Regression contract

The parser and metadata contract are covered by
`tests/integration/joint_benchmark_test.py`. The registered smoke test runs the
small `chain-8` search and solve profiles. Numerical C++ tests cover ordinary,
multi-component, observable/nuisance, partial-selection, prepared-block,
stationarity, order, and persistence behavior separately from the benchmark
wrapper.
