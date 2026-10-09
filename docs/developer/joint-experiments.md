# Joint experiment lifecycle

This is the current-tree inventory for Joint and FixedNeighbor experiments.
Rows are experiment units rather than individual output files; the table is
kept deliberately regular so it can be consumed as a simple pipe-delimited
inventory. Status is one of `Active`, `Qualified`, `Closed`, or `Unknown`.

Current FixedNeighbor means the qualified production policy: core 12, Forward
serial Gauss-Seidel, 30 maximum sweeps, local `LegacyCompact`, and at most one
trusted accepted local update per block visit. The benchmark and route-frontier
drivers use that policy by default; only structural core-size overrides are
custom benchmark inputs.
Full/core128 and OneAccepted/core64 are historical evidence only; the latter
was superseded by the qualified core12 policy.

| Experiment | Status | Production relevance | Canonical evidence | Active driver | Active test | Raw artifacts |
| --- | --- | --- | --- | --- | --- | --- |
| Joint estimator runtime, persistence, and component regression | Active | Production correctness and provenance | `joint-component-runtime.md`; production C++ tests | `joint_component_runtime.py`; `joint_workflow_cli_smoke.py` | `joint_component_runtime_test.py`; runtime/physical smoke; persistence tests | No tracked execution telemetry |
| FixedNeighbor stationarity confirmation | Qualified | Locks the existing eta-only confirmation contract | `joint-fixed-neighbor-stationarity-r1/{README.md,analysis.json,summary.csv}` | `joint_fixed_neighbor_stationarity.py` | `joint_fixed_neighbor_stationarity_test.py` | Compact evidence; sweep telemetry is historical |
| FixedNeighbor historical local-update selection evidence | Qualified | Explains why the one-accepted update contract is intrinsic to the production route | `joint-fixed-neighbor-inexact-local-r1/`; `joint-fixed-neighbor-inexact-qualification-r2/` | none | `FixedNeighborBlockCoordinate_test.cpp` | Compact historical evidence; Full and TwoAccepted are alternatives from the completed selection study |
| FixedNeighbor outer-core qualification | Qualified | Justifies the adopted production core of 12 | `joint-fixed-neighbor-outer-core-qualification.md`; C3/C4 `README.md`, `analysis.json`, `summary.csv` | `joint_fixed_neighbor_outer_core_qualification.py` with `screen`, `endpoint`, and `frontier` phases | `joint_fixed_neighbor_outer_core_qualification_test.py` | r1-r5 current tree keeps compact evidence only; raw execution telemetry is historical |
| FixedNeighbor prepared-block behavior | Qualified | Preserves prepared mapping and symbolic-reuse contract under the current core12 policy | `joint-fixed-neighbor-experimental.md`; prepared-block structural contract | `joint_fixed_neighbor_prepared_block.py` | `joint_fixed_neighbor_prepared_block_test.py` | Small deterministic contract only; historical large-case output is compact evidence |
| FixedNeighbor order consistency | Qualified | Confirms production Forward order against the retained diagnostic | `joint-fixed-neighbor-order-r1/` | `joint_fixed_neighbor_order.py` | `joint_fixed_neighbor_order_test.py` | Compact comparison; raw pairs are not a runtime dependency |
| FixedNeighbor block acceptance / Fixed-B regression | Qualified | Protects replay-aware acceptance and block contracts | fixed-B campaign summaries; `joint-fixed-neighbor-experimental.md` | `joint_fixed_b_scaling.py`; C++ experiment | FixedB/FixedNeighbor component regression tests | Compact summaries; no production telemetry dependency |
| Global Joint benchmark and route frontier | Active | Current scalability evidence for explicit production routes; no auto-routing | `joint-benchmark.md`; route frontier summaries | `joint_benchmark.py`; `joint_route_frontier.py` | benchmark and route-frontier contract tests | Compact evidence is retained; superseded raw telemetry is pruned and new output belongs in the build tree |
| Superseded FixedNeighbor optimized frontier | Closed | Historical OneAccepted/core64 comparison; superseded by the qualified core12 production policy | `figures/joint-fixed-neighbor-optimized-frontier-r1/{README.md,analysis.json,summary.csv,campaign-manifest.json}` | none | none | Compact evidence retained; no active driver or raw dependency |
| Global Operator/SPQR/Schwarz scaling | Active | Current benchmark infrastructure, not a FixedNeighbor production policy | `joint-benchmark.md`; sparse backend docs | `joint_schwarz_sweep.py`; sparse benchmark tools | benchmark/scalability contracts | Canonical campaign summaries; outputs should be generated outside the source tree |
| Global factor ownership, bounded-rank, and projected-width/tail evidence | Qualified | Supports current Operator diagnostics without being a runtime gate | `joint-operator-search.md`; compact factor/rank/projected reports | analyzers are offline-only | no registered current test for the orphan analyzers | Compact evidence only; orphan analyzer tests have no current dependency |
| Historical FixedNeighbor endpoint/trajectory diagnostics | Closed | Negative evidence only; no production route or threshold change | `joint-fixed-neighbor-experimental.md` and compact diagnostic summaries | none | none after retirement | Raw endpoint decomposition, local certification, and trajectory output is removable |
| FixedNeighbor local Operator/Schwarz/hybrid branches | Closed | Explicitly not promoted; production remains LegacyCompact | negative evidence in `joint-fixed-neighbor-experimental.md` | none after retirement | none after retirement | Local hooks are removed from the current tree; historical implementation is recoverable from Git history |
| FixedNeighbor workspace residency, numeric reuse, and local attribution | Closed | No production wall-time or correctness benefit | `joint-fixed-neighbor-workspace-attribution.md` and compact summaries | none after retirement | none after retirement | Raw repetitions, wrappers, progress, and process logs are removable; no current policy switch remains |
| Historical FixedNeighbor CLI modes and campaign wrappers | Closed | Reproducible from Git history if needed; not a production surface | this inventory plus the closed-direction table in the canonical overview | none | none | Do not retain zero-caller mode telemetry or wrappers |
| Orphan factor-ownership/projected-tail/projected-width analyzer tests | Closed | No CMake registration, import, or current gate | corresponding compact figure reports | none | none | Test files are removable; canonical evidence remains |
| Generic experiment I/O, process monitoring, provenance, and manifest helpers | Active | Shared infrastructure for active benchmark/qualification drivers | `tests/integration/experiment_*.py` contracts | shared helpers | `experiment_support_test.py`; dependency tests | Generated output belongs in build/output directories |

## Retention policy

For a tracked campaign, retain the compact evidence that explains its result:

- `README.md`, when present;
- `analysis.json`, when present;
- `summary.csv`, when present; and
- `campaign-manifest.json`, when present.

These are defaults, not a requirement to invent a missing file. A small
representative result is retained only when the compact evidence cannot explain
the contract without it.

The following are not permanent source-tree artifacts by default:

- `runs.json` and other aggregate execution dumps;
- `individual-results/`;
- `*.progress.json`;
- `*-process/` and its `stdout.txt`/`stderr.txt`;
- per-run wrapper JSON; and
- duplicate raw result JSON.

Drivers should write new execution output below the build/output tree. Any
tracked exception must be named by the campaign README and justified by a
current analyzer or regression contract.

This cleanup only reduces the current checkout and prevents future source-tree
growth. It does not rewrite existing Git history. In particular, this work
does not use `git filter-repo`, BFG, force-push, reset to an old commit, or any
other history-rewriting operation.
