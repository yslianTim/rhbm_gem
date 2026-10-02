# Joint benchmark

The current benchmark entry point is `tests/integration/joint_benchmark.py`. It
selects a focused profile, runs an existing C++ driver through the shared
process/resource support, and writes one versioned JSON result. It does not
reimplement estimator calculations or apply historical timing gates.

Build the optional drivers with:

```sh
cmake -S . -B build/debug -DBUILD_TESTING=ON -DRHBM_GEM_BUILD_BENCHMARKS=ON
cmake --build build/debug --target benchmarks_all rhbm_gem_cli -j4
```

Run a profile with an explicit case and output file:

```sh
python3 tests/integration/joint_benchmark.py \
  --profile fixed --case chain-8 --build-dir build/debug \
  --output build/joint-fixed.json
```

## Profiles

| Profile | Measures | Driver |
| --- | --- | --- |
| `prepare` | Joint problem construction and fixed-state basis preparation | `joint_sparse_benchmark` |
| `fixed` | Fixed-state operator preparation and one selected preconditioned step; the saved state is prepared in a separate, unmeasured process. Optional `--svd-mode legacy|values|auto` selects the compact-SVD path | `joint_sparse_benchmark` |
| `solve` | Joint search and returned-state assessment, with explicit Operator-PCG rank and Schwarz policies | `joint_sparse_benchmark` |
| `rank` | Bounded free-design rank evaluation; select `--rank-mode prototype` or `oracle` | `joint_sparse_benchmark` |
| `workflow` | Prepare, estimate, postprocess, and persist | `joint_postprocessing_benchmark` |
| `postprocess` | Uncertainty, peeling, summary, and persistence using the driver's fixed endpoint | `joint_postprocessing_benchmark` |
| `command` | CLI analysis, SQLite persistence/reload, and joint JSON/CSV export | `RHBM-GEM` CLI |

`prepare`, `fixed`, `solve`, and `rank` accept deterministic synthetic cases such
as `chain-8` and `cube-8`. Fixed and solve also accept frozen fixture cases such
as `baseline:first-stage-double` or `weak-1e-4:first-stage-double`. Workflow
and postprocess cases are `full`, `halo`, and `multi`.

The `fixed` and `solve` profiles accept these internal policy controls:

```text
--operator-rank auto|dense|spqr-bounds
--schwarz-core-atoms N
--schwarz-overlap-hops N
--schwarz-max-block-atoms N
--schwarz-storage-mib N
--schwarz-scratch-mib N
```

Defaults are `auto`, 128 core atoms, one overlap hop, 512 maximum block atoms,
512 MiB storage, and 256 MiB scratch. `auto` resolves to bounded rank on SPQR
and dense rank on Eigen. The benchmark writes the requested and resolved rank
backend, budget, search method, preconditioner, and all Schwarz limits into
`metadata.solver_policy`; the C++ result also records its effective policy.
These options belong to the internal benchmark drivers and do not add
production CLI switches.

For example, compare a small Schwarz core with core-only and one-hop overlap:

```sh
python3 tests/integration/joint_benchmark.py \
  --profile solve --case chain-8 --preconditioner schwarz \
  --operator-rank auto --schwarz-core-atoms 2 --schwarz-overlap-hops 0 \
  --schwarz-max-block-atoms 16 \
  --build-dir build/debug --output build/joint-schwarz-core-only.json
```

The command profile requires explicit model and map inputs:

```sh
python3 tests/integration/joint_benchmark.py \
  --profile command --case sample --build-dir build/debug \
  --model input.cif --map input.map --output build/joint-command.json
```

An unavailable CLI or missing input is recorded with `execution.status` set to
`unavailable`. It is not reported as a completed measurement.

## Result fields

Each result has `metadata`, `execution`, `problem`, `numerics`, `resources`, and
`result` sections. Metadata records the current commit, backend, build type,
compiler, benchmark source, case, and an input or fixture hash. `profile` and
`measurement_scope` identify the timer boundary. Profile-specific elapsed time
comes from the matching C++ phase; process wall time and peak RSS are also kept
for each run. Rank, objective, runtime convergence, and available active-face
atoms are retained as numerical results. Absolute local paths, process IDs, and temporary filenames are
not written to the result.

`execution.status` describes process completion, timeout, resource stop, or
work that did not run. `result.qualified` reports a solve's existing runtime
convergence assessment where one exists; it does not change process status or
create a CI performance gate. Repetitions retain individual outcomes, and a
failed or unavailable repetition prevents an aggregate elapsed/RSS value from
being reported as complete. Solve details retain the C++ search work counters,
rank status and resource use. Schwarz runs also include a partition summary
with block count, core and block size distributions, atom membership
distribution, graph workspace, storage, and scratch bound. Rank counters show
checks, time, entries, workspace, compact extraction count, and free-design
SVD count, so the bounded SPQR route can be checked directly.

The fixed-state compact-SVD, operator, rank, and action contracts remain owned
by permanent tests such as `Numerics_test`, `ProfileOperator_test`, and
`FreeDesignRank_test`. Search and bounded-work contracts remain in
`Search_test` and `OperatorSearch_test`. Independent replay and offline
certification remain in `joint_numerical_reference.py`, `joint_offline_support`,
and the opt-in offline tests. The benchmark measures these paths; it does not
replace their assertions.

Historical comparison results and archive provenance are summarized in the
[Joint evidence guide](joint-component-evidence.md). This guide covers current
benchmark profiles and measurement semantics.
