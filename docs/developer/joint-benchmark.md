# Joint benchmark

The current benchmark measures the production FixedNeighbor estimator and
current workflow costs. Its only workflow entry point is
[tests/integration/joint_benchmark.py](../../tests/integration/joint_benchmark.py).

## Profiles

| Profile | Measurement scope | Executable |
| --- | --- | --- |
| search | FixedNeighbor search and local Profile LM work; no endpoint qualification claim | joint_fixed_neighbor_experiment with --scaling-only |
| solve | Search, global replay, endpoint certification and RuntimeConvergence | joint_fixed_neighbor_experiment with --case |
| workflow | Joint fitting, target peeling, post-Joint local MDPDE, group inference and persistence | joint_postprocessing_benchmark |
| command | CLI analysis, save/reload and Joint export | RHBM-GEM |

Build the search, solve and workflow executables with
`RHBM_GEM_BUILD_BENCHMARKS=ON`, using no more than four jobs:

```sh
cmake -S . -B build/bench -DBUILD_TESTING=ON -DRHBM_GEM_BUILD_BENCHMARKS=ON
cmake --build build/bench --target benchmarks_all --parallel 4
```

The command profile uses the regular CLI build and requires `--model` and
`--map` inputs.

For small synthetic search and solve measurements:

    python3 tests/integration/joint_benchmark.py \
      --profile search --case chain-8 \
      --build-dir build/bench --output build/joint-search.json

    python3 tests/integration/joint_benchmark.py \
      --profile solve --case chain-8 \
      --build-dir build/bench --output build/joint-solve.json

These runs default to a 30-sweep diagnostic cap. Use `--maximum-sweeps 0` for
structural preparation only, `--maximum-sweeps 1` for one complete sweep, or a
larger cap up to 200 for a bounded multi-sweep measurement. These bounds apply
only to the experiment driver; they do not change the production estimator.
Aggregate telemetry is recorded by default. Add `--detailed` to retain
per-block and profile-trial records. The driver reports a resource census,
estimated preparation/global-state peaks, process RSS checkpoints before
preparation, after preparation, after the first sweep and after certification,
and the process peak RSS where the platform exposes it. Unsupported checkpoint
measurements are `null`.

For an experimental affected-row candidate replay comparison, add
`--incremental-replay-experiment` to a search or solve run. The experiment
measures the incremental objective/residual delta beside the unchanged full
component replay and compares acceptance decisions. It still uses full replay
as the authority and does not run a different estimator trajectory; results
are diagnostic evidence, not a production speedup claim.

Use workflow with `full`, `halo`, `multi`, `noisy`, `outlier`, `dense`,
`target-heavy`, or `connected` cases. The scaling cases reuse the same driver:
`full` has three targets, `target-heavy` has six inside the same spatial
envelope, `dense` keeps three targets and halves voxel spacing, `multi` has
three small components, and `connected` has one nine-target component. The
report records actual Joint observations, contributors and targets so these
workloads can be compared without assuming the row count. `halo` is partial
selection; `noisy` and `outlier` retain their deterministic perturbations.

## Measured postprocessing cases

One Release run used `build/joint-validation` on AppleClang with
`OMP_NUM_THREADS=1` and BLAS-family worker counts set to 1. Each case had one
measurement, a 600-second timeout and a 4 GiB RSS limit. All eight completed;
all had one observation-row map, one contributor map and one Joint-stage map.
Every selected target in these fixtures qualified for local evidence, and all
prepared samples were valid. There is no pre-change runtime baseline, so these
measurements establish current workload costs and mapping counts, not a measured
speedup.

| Case | Observations | Contributors / targets | Valid / prepared samples | Driver seconds | Postprocessing seconds | Process peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| full | 22,650 | 3 / 3 | 600 / 600 | 0.663 | 0.019 | 56.5 |
| halo | 15,967 | 5 / 1 | 200 / 200 | 0.437 | 0.008 | 50.8 |
| multi | 67,950 | 9 / 9 | 1,800 / 1,800 | 0.977 | 0.059 | 99.7 |
| noisy | 22,650 | 3 / 3 | 600 / 600 | 0.421 | 0.022 | 57.5 |
| outlier | 22,650 | 3 / 3 | 600 / 600 | 0.429 | 0.020 | 56.8 |
| dense | 181,122 | 3 / 3 | 600 / 600 | 1.932 | 0.082 | 263.0 |
| target-heavy | 27,185 | 6 / 6 | 1,200 / 1,200 | 0.654 | 0.035 | 78.1 |
| connected | 42,558 | 9 / 9 | 1,800 / 1,800 | 1.070 | 0.057 | 106.2 |

The `dense` run increased observations about 8.0× at three targets; postprocessing
was 0.082 s. `target-heavy` doubled targets while observations increased by
about 20%; postprocessing was 0.035 s. These are different synthetic shapes,
not controlled baseline comparisons. Process peak RSS includes fixture creation
and persistence; the 5 ms workflow RSS sample and payload estimates are reported
separately in the [memory census](joint-component-memory-census.md). No sample
fell outside the Joint observation domain in these generated cases; the
missing-coverage behavior is covered by the Joint postprocessing regression
fixture, not by a benchmark timing result.

Workflow JSON reports row, contributor and Joint-stage mapping build counts;
target/sample counts; context, peeling, charge-correction, local alpha/local
MDPDE, group alpha and WEB timings; and current-RSS checkpoints. The context
and workspace byte counts are estimates of container and sample payloads and
exclude allocator overhead. The 5 ms sampled workflow peak can miss short
spikes. Process peak RSS is the whole-process high-water mark, including fixture
construction and persistence, and does not attribute memory to a phase. RSS
availability is recorded separately; unsupported measurements are `null`.
Run the wrapper with explicit bounds for each case, for example:

```sh
python3 tests/integration/joint_benchmark.py \
  --profile workflow --case target-heavy \
  --build-dir build/bench --output build/joint-target-heavy.json \
  --timeout 600 --rss-limit 4294967296
```

These synthetic fixtures exercise workflow contracts but do not establish
large-map or full-atom-set resource requirements. A completed workflow report
preserves per-target Joint and local statuses; it does not imply every endpoint
qualified. Use command with a case name plus `--model` and `--map`. Reports and
temporary outputs belong in the build/output tree.

Production rank uses values-only projected-width SVD. The diagnostic weak
direction path is opt-in and has its own `ProjectedWidthDiagnostics` memory
estimate; the benchmark does not run it unless a separate diagnostic tool is
invoked.

## Fixed production policy

Every report records the frozen estimator policy:

| Setting | Value |
| --- | --- |
| Estimator | FixedNeighbor |
| Core | 12 atoms |
| Block order | Forward |
| Maximum sweeps | 200 |
| Local search | Profile LM |
| Provenance | `fixed_neighbor_local_search=profile-lm` |
| Accepted local updates | OneAccepted: at most one trusted update per block visit |
| Sparse backend | EIGEN |

Search and solve are the only FixedNeighbor measurement modes. The C++ driver
accepts --scaling-only for search measurement or --case for solve measurement.
It exposes no reverse-order mode, arbitrary core size or alternative policy.

Search reports search and structural-preparation work; solve also reports
primary endpoint, independent replay and rank-assessment time separately.
Reports state `reference_solve_calls=0` for production. Aggregate timings attribute local Profile LM work,
old-core and effective-response preparation, candidate replay and copying,
cache updates, sweep replay and global-state evaluation. Reports also include
block visits and solves, accepted updates, profile evaluations, full candidate
replays/copies and maximum block time. The Python wrapper records process
status, elapsed time, source/build metadata, endpoint evidence, estimated peak
bytes, observed process peak RSS, the diagnostic sweep cap and the unchanged
production solver policy. Use --repeat, --warmup, --timeout and
--rss-limit to control the measurement process.

Solve output also includes `assessment_stage_telemetry`: stage entry and
completion counts, dimensions, elapsed seconds, preflight status and byte
breakdowns, current RSS before/after, and process peak RSS. The two RSS fields
are labeled separately because process peak RSS is a process-wide high-water
mark, not a stage-local allocation. Missing platform readings are JSON `null`.
`production_reference_solve_calls` counts the retired global reference-solve
path; `reference_qr_factorization_calls_total` counts all calls to the
`ReferenceQR` primitive, including offline `DesignSpectrum()` oracle calls.
Production endpoint qualification no longer calls that oracle: its free-design
TiledQR/SVD is inside projected-width preparation. The large synthetic stage
results are documented in the [component memory census](joint-component-memory-census.md).

The wrapper modes describe different measurement scopes; they do not select
different estimator policies. There is no backend selector or alternate
solver route. EIGEN is the only current sparse backend.

## Permanent checks

joint_benchmark_contract_test checks the profile and metadata contract,
including rejection of arbitrary core sizes. joint_benchmark_smoke runs small
search and solve cases and checks the production policy. The Joint C++ contract,
workflow and numerical tests own numerical correctness and frozen thresholds;
a benchmark result does not replace those tests. Bounded-rank controls compare
the dense SVD oracle with the unchanged row, column-normalization and strict
threshold contract. Production structural evidence is a one-way rank-failure
check; matching and streaming estimates remain research diagnostics unless
their stronger proof requirements are met. The `joint:offline` SparseQR
research tests run Natural and COLAMD in isolated serial processes for
chain/cube sizes 128, 256, 512 and 1,024, with a memory preflight before any
factorization. Small SVD comparisons record pivot-rank disagreements; 1,024
cases above the planning budget report `ResourceLimit` and `NotRun`. These
tests characterize a research path and do not add a production solver option.

## Historical decisions

Completed qualification campaigns established core 12, Forward order and
OneAccepted. Fixed-B and global OperatorPcg/Schwarz were not promoted. The
bounded SPQR comparison passed numerical parity and favored SPQR on completed
128/256 wall-time probes, while EIGEN used less peak RSS; both 512 probes timed
out. A later memory, dependency and maintenance policy selected EIGEN.

Closed experiment implementations and machine-readable receipts are
recoverable from Git history and are intentionally not retained in the current
tree. The historical decisions are summarized in
[joint-fixed-neighbor-experimental.md](joint-fixed-neighbor-experimental.md),
[joint-operator-search.md](joint-operator-search.md) and
[joint-component-evidence.md](joint-component-evidence.md).
