# Joint benchmark

The current benchmark measures the production FixedNeighbor estimator and
current workflow costs. Its only workflow entry point is
[tests/integration/joint_benchmark.py](../../tests/integration/joint_benchmark.py).

## Profiles

| Profile | Measurement scope | Executable |
| --- | --- | --- |
| search | FixedNeighbor search and local Profile LM work; no endpoint qualification claim | joint_fixed_neighbor_experiment with --scaling-only |
| solve | Search, global replay, endpoint certification and RuntimeConvergence | joint_fixed_neighbor_experiment with --case |
| workflow | Complete in-memory Joint workflow and persistence | joint_postprocessing_benchmark |
| postprocess | Peeling, uncertainty, group processing and persistence | joint_postprocessing_benchmark |
| command | CLI analysis, save/reload and Joint export | RHBM-GEM |

Build the search, solve, workflow and postprocess executables with
RHBM_GEM_BUILD_BENCHMARKS=ON. The command profile uses the regular CLI build
and requires --model and --map inputs.

For small synthetic search and solve measurements:

    python3 tests/integration/joint_benchmark.py \
      --profile search --case chain-8 \
      --build-dir build/bench --output build/joint-search.json

    python3 tests/integration/joint_benchmark.py \
      --profile solve --case chain-8 \
      --build-dir build/bench --output build/joint-solve.json

Use workflow and postprocess with full, halo or multi cases. Use command with a
case name plus --model and --map. Reports and temporary outputs belong in the
build/output tree.

## Fixed production policy

Every report records the frozen estimator policy:

| Setting | Value |
| --- | --- |
| Estimator | FixedNeighbor |
| Core | 12 atoms |
| Block order | Forward |
| Maximum sweeps | 30 |
| Local search | Profile LM |
| Provenance | `fixed_neighbor_local_search=profile-lm` |
| Accepted local updates | OneAccepted: at most one trusted update per block visit |
| Sparse backend | EIGEN |

Search and solve are the only FixedNeighbor measurement modes. The C++ driver
accepts --scaling-only for search measurement or --case for solve measurement.
It exposes no reverse-order mode, arbitrary core size or alternative policy.

Search reports search_seconds; solve reports search_seconds and
assessment_seconds. Both can report peak RSS, sweeps, block solves, profile
evaluations, KKT and width-gradient metrics. The Python wrapper also records
process status, elapsed time, source/build metadata and the fixed solver policy.
Use --repeat, --warmup, --timeout and --rss-limit to control the measurement
process.

The wrapper modes describe different measurement scopes; they do not select
different estimator policies. There is no backend selector or alternate
solver route. EIGEN is the only current sparse backend.

## Permanent checks

joint_benchmark_contract_test checks the profile and metadata contract,
including rejection of arbitrary core sizes. joint_benchmark_smoke runs small
search and solve cases and checks the production policy. The Joint C++ contract,
workflow and numerical tests own numerical correctness and frozen thresholds;
a benchmark result does not replace those tests.

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
