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
| `solve` | Joint search and returned-state assessment for LegacyCompact, OperatorPcg, or FixedNeighbor | `joint_sparse_benchmark` or `joint_fixed_neighbor_experiment` |
| `search` | LegacyCompact, OperatorPcg, or FixedNeighbor search only; no returned-state assessment or endpoint qualification | `joint_sparse_benchmark` or `joint_fixed_neighbor_experiment` |
| `rank` | Bounded free-design rank evaluation; select `--rank-mode prototype` or `oracle` | `joint_sparse_benchmark` |
| `workflow` | Prepare, estimate, postprocess, and persist | `joint_postprocessing_benchmark` |
| `postprocess` | Uncertainty, peeling, summary, and persistence using the driver's fixed endpoint | `joint_postprocessing_benchmark` |
| `command` | CLI analysis, SQLite persistence/reload, and joint JSON/CSV export | `RHBM-GEM` CLI |

`prepare`, `fixed`, `solve`, and `rank` accept deterministic synthetic cases such
as `chain-8` and `cube-8`. Fixed and solve also accept frozen fixture cases such
as `baseline:first-stage-double` or `weak-1e-4:first-stage-double`. Workflow
and postprocess cases are `full`, `halo`, and `multi`.

### Production FixedNeighbor route

FixedNeighbor is a supported production Joint search method, selected
explicitly through `FitOptions::joint_search_method` or
`potential_analysis --joint-search fixed-neighbor`. Its production policy is
core size 64, forward serial Gauss-Seidel order, at most 30 sweeps, and one
trusted accepted local `LegacyCompact` update per block visit. The benchmark
and historical experiment drivers still expose broader controls for
reproduction and attribution; those controls are diagnostic and are not
production defaults. Search-only measurements remain search-only and cannot
establish endpoint certification.

The current matched synthetic frontier shows the intended trade-off rather
than a performance gate: OperatorPcg is faster on the tested large cases,
while FixedNeighbor materially reduces peak RSS through bounded local block
factors. On cube-1024, the formal 600-second FixedNeighbor envelope timed out
near convergence; an uncapped diagnostic run completed at about 769 seconds
with peak RSS around 276 MiB, compared with about 221 seconds and 2681 MiB for
OperatorPcg. The timeout is a wall-time/resource trade-off, not a numerical
failure, and does not make FixedNeighbor the default. Search-only evidence
through larger sizes must not be described as full endpoint certification.

The repeated workspace-residency study is closed as Route C. It used the same
production policy on chain/cube × 256/512/1024 with one warmup and three
measured repetitions per case. All six numerical gates passed, while the
persistent-workspace wall-time gate was `not-material`: measured search changes
ranged from `-0.32%` to `+0.16%`, with higher persistent RSS in every case.
Exact numeric-factor reuse opportunities ranged from `1.17%` to `15.00%` and
were limited to initial-profile evaluation; they are not a production reuse
policy. No workspace-policy change, exact numeric-factor reuse, or candidate
replay optimization was made. The current large-case bottleneck is the local
derivative preparation/reduction path, which accounted for `83.7%` of
chain-1024 search and `81.5%` of cube-1024 search. See the [workspace
attribution](joint-fixed-neighbor-workspace-attribution.md) record.

The `fixed` and `solve` profiles accept these internal policy controls:

```text
--operator-rank auto|dense|spqr-bounds
--operator-rank-seconds S
--operator-rank-work-entries N
--operator-rank-workspace-mib N
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
Rank-budget defaults are 120 seconds, 100,000,000 work entries, and 256 MiB.
The `rank` profile records its actual `SpqrBounds` or dense-oracle path and
accepts the same three budget controls for the bounded prototype. These options
belong to the internal benchmark drivers and do not add production CLI
switches.

The `search` profile uses the same `SearchProfile(...)` path as the solve
profile for LegacyCompact and OperatorPcg and stops before
`AssessComponentSearch(...)`. The FixedNeighbor route uses the same synthetic
workload and stops before its final endpoint assessment. Each writes
`measurement_scope: "search-only"`; returned assessment and endpoint
qualification are not run. Search-only results do not establish runtime
convergence or endpoint qualification. For example:

```sh
python3 tests/integration/joint_benchmark.py \
  --profile search --case chain-512 --preconditioner schwarz \
  --operator-rank auto --schwarz-core-atoms 128 --schwarz-overlap-hops 1 \
  --schwarz-max-block-atoms 512 --build-dir build/joint-spqr \
  --output build/chain-512-search.json
```

### Prepared-block qualification

The prepared FixedNeighbor route has a dedicated matched search-only campaign:

```sh
cmake --build build/joint-spqr --target joint_fixed_neighbor_experiment -j4
python3 tests/integration/joint_fixed_neighbor_prepared_block.py \
  --build-dir build/joint-spqr \
  --output-dir build/joint-fixed-neighbor-prepared
```

The campaign fixes the production policy at 64-atom forward cores, a 30-sweep
budget, one accepted local `LegacyCompact` update and one Eigen thread. Its
chain/cube 256, 512 and 1024 cases are explicitly `fixed-neighbor-search-only`.
The gate checks the existing A/C KKT (`1e-10`), width-gradient (`1e-12`) and
complete-sweep eta confirmation (`1e-10`), then attributes preparation and
factor work with `prepared_block_count`, block/domain/mapping preparation
counters, symbolic factorizations, symbolic reuses and numeric factorizations.
The fresh-workspace symbolic count is an instrumentation control estimate;
`wall_time_gate` remains separate and is `not-run` unless a same-policy matched
wall-time control is supplied. These results must not be described as full
endpoint certification for the large search-only cases.

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
failed or unavailable repetition prevents `joint_benchmark.py` from reporting
its aggregate elapsed/RSS value as complete. Solve details retain the C++ search work counters,
rank status, certificate path, and resource use. Bounded rank diagnostics
include local-support coverage and its threshold/lower bound, the work stage,
charged entries, estimated total and remaining entries when the reconstruction
forecast is available, and the design, R-factor, reflector-nonzero and
reflector counts. Schwarz runs also include a partition summary
with block count, core and block size distributions, atom membership
distribution, graph workspace, storage, and scratch bound. Rank counters show
checks, time, entries, workspace, compact extraction count, and free-design
SVD count, so the bounded SPQR route can be checked directly.

When test instrumentation is enabled, solve results also include assessment
stage calls, completed calls, rows, columns, and seconds for primary and
reference evaluation, design spectrum, derivative preparation and reduction,
projected and normalized width spectra, correction/Jacobian spectrum, and the
total assessment. The benchmark snapshots the active and last completed stage
when each stage begins or ends. If a run times out, its individual result keeps
the last stage checkpoint and process resource observations; it remains a
timeout and is not treated as a completed assessment. This observer is
test/benchmark-only and is not part of the installed numerical API.

## Cost accounting

The sweep reports timing summaries across completed measurement repetitions.
Warmup runs are excluded. Timing and work-counter distributions retain their
minimum, median, and maximum; the canonical timing fields use the median. Peak
RSS uses the maximum. `measurements_complete` is true only when all requested
measurement repetitions completed, and scaling analysis uses complete rows.
For PCG, `measurements.pcg_iterations_per_solve.median` is the median of the
per-run medians. The pooled per-solve distribution remains available alongside
it. Counter variation across repetitions is kept in `measurements`.

The decomposed costs are defined per measurement run:

| Field | Definition |
| --- | --- |
| `topology_setup_seconds` | `partition_seconds` |
| `linearization_setup_seconds` | `operator_prepare_seconds + metric_seconds + local_seconds` |
| `damping_setup_seconds` | `factor_seconds` |
| `setup_seconds` | Topology + linearization + damping setup |
| `iterative_seconds` | `pcg_seconds` |
| `measured_seconds` | `search_seconds + assessment_seconds` |
| `wall_seconds` | Process wall time from the benchmark runner |

The setup sum includes `operator_prepare_seconds` once. Operator preparation
already covers design, fixed-factor, and rank work, so those subfields are
diagnostics and are not added again. `operator_seconds` sums operator
preparation and operator actions; Apply/Adjoint actions can also occur outside
PCG, so this sum is not the PCG cost. `operator_normals`,
`operator_applications`, and `operator_adjoints` count implementation calls.
`ApplyNormal` is one optimized operator action, not an implied Apply plus
Adjoint pair. The setup and PCG fractions divide their corresponding per-run
cost by `search_seconds`; they are null when that denominator is zero or
unavailable. These are descriptive measurements and have no hardware-specific
pass threshold.

## Six scalability axes

Do not use total solve wall time as a proxy for a single scaling property. The
current evidence separates six questions:

| Scalability axis | Measurement | Current evidence |
| --- | --- | --- |
| Rank certification | bounded rank work, certificate status, SPQR reconstruction | [rank frontier](figures/joint-rank-budget-frontier/campaign-manifest.json) and [local-support census](figures/joint-local-rank-witness-census/local-rank-witness-census.json) |
| Krylov iteration scaling | iterations per solve under a fixed Schwarz policy | [cube multi-block gate](figures/joint-search-scaling-r2/search-scaling-analysis.json) |
| Search operator throughput | operator actions, factorization, and search time | [factor ownership](figures/joint-factor-ownership-r1/factor-ownership-analysis.json) and [operator factor frontier](figures/joint-operator-factor-frontier-r1/factor-frontier-analysis.json) |
| Search factor memory | factor lifetime, owned-byte estimates, factor fill, and sampled process-tree RSS | [factor ownership](figures/joint-factor-ownership-r1/factor-ownership-analysis.json), [factor residency](figures/joint-factor-residency-r2/residency-analysis.json), and [cube memory analysis](figures/cube-search-memory-r1/cube-memory-analysis.json) |
| Assessment derivative reduction | row generation and projected-width QR cost | [projected-width comparison](figures/joint-projected-width-r1/projected-width-analysis.json) and [derivative scaling](figures/joint-derivative-scaling-r1/derivative-scaling-analysis.json) |
| End-to-end returned-state assessment | assessment completion, returned-state parity, timeout, and RSS | [compact acceptance](figures/joint-compact-acceptance-r2/acceptance-analysis.json) and [post-compact frontier](figures/joint-post-compact-frontier-r1/frontier-analysis.json) |

Search-only results do not establish assessment or endpoint qualification.
Assessment timeouts remain incomplete evidence even when search succeeded.
Production returned-state assessment now uses guarded compact stacked QR when
the singular-value decision is clear; the exact observation-scale Jacobian TSQR
remains the fallback at the one-ULP rank boundary. The 128/256 acceptance cases
passed full assessment and returned-state parity, and production-route 512
chain/cube assessments completed within the existing 600-second / 4-GiB
envelope. See the compact acceptance artifacts for per-case timings and RSS.

## Schwarz scaling sweep

`joint_schwarz_sweep.py` creates a configuration matrix and invokes the
single-run benchmark once per configuration. Its defaults are a small `chain`
sweep over 64 and 128 atoms, core sizes 32 and 64, overlap hops 0 and 1, and
one repetition. Larger workloads such as 512, 1024, 2048, and 4096 atoms must
be requested explicitly.

```sh
python3 tests/integration/joint_schwarz_sweep.py \
  --build-dir build/debug --output build/joint-schwarz-sweep.json \
  --csv build/joint-schwarz-sweep.csv \
  --topologies chain cube --atoms 512 1024 --cores 64 128 \
  --overlaps 0 1 2 --preconditioners identity diagonal schwarz \
  --operator-rank auto --repeat 1 --timeout 600
```

Schwarz configurations use the Cartesian product of core sizes and overlap
hops. Identity and Diagonal each produce one control per topology and atom
count, without multiplying across unused Schwarz settings. Each configuration
has its own versioned benchmark JSON under the aggregate output's sibling
`*_runs` directory. The aggregate JSON and optional CSV include configuration,
status, rows and free columns, rank backend and work, partition summaries,
PCG per-solve iteration counts, operator calls, timing distributions, search
and assessment time, wall time, peak RSS, objective, convergence, state
availability, and stop reason. `--warmup` defaults to zero and `--repeat`
defaults to one. An unavailable, timed out, or failed measurement keeps its
status and individual result; incomplete repetitions are not reported as a
complete campaign. Timing is descriptive and has no fixed performance pass
threshold.

## One-level scaling analysis

`joint_scaling_analysis.py` reads the sweep JSON and does not rerun the
estimator. It groups rows by topology, preconditioner, core size, overlap,
maximum block size, operator-rank mode, and sparse backend. Incomplete rows,
identity/diagonal controls, missing policies, and duplicate global sizes cannot
produce a one-level Schwarz verdict. The primary point value is
`measurements.pcg_iterations_per_solve.median`, the median of per-repetition
PCG solve medians. Total PCG iterations are supporting work only.
Every completed repetition must contain at least one PCG solve for the point to
be valid; missing per-solve telemetry is excluded from the valid-point count.
Campaign rows also carry `evidence_eligible` and
`evidence_exclusion_reasons`. The analyzer excludes rows explicitly marked
ineligible; unmarked inputs remain usable for general analysis.

At least four complete points with positive atom counts and positive per-solve
iteration medians are required for a slope verdict. `growth-observed` requires
both a log-log slope of at least 0.25 and an endpoint iteration growth ratio
of at least 1.5. Adjacent points within 15% of exact doubling also report their
iteration ratios; all adjacent pairs report interval log-log slopes. A
completed larger-size point with any `pcg-iteration-budget` repetition takes
the stronger `pcg-budget-limited` classification. That classification records
that coarse correction investigation is warranted; it leaves solver changes to
a separate design decision. `stable`
means these points do not meet the growth gate; it does not guarantee scaling
at larger sizes. `not-comparable` means the rows do not describe a single
one-level Schwarz policy.

The analysis separately reports PCG-solve, linearization, and damping-trial
scaling. Growth in these counts with stable per-solve iterations is diagnosed
as `nonlinear-work-growth`, not Krylov growth. If wall time or peak RSS grows
by at least 50% while PCG iterations remain stable, the diagnostic is
`cost-growth-with-stable-krylov`; investigate operator work, rank/setup,
memory bandwidth, or storage. These diagnostics do not alter the solver.

Example bounded tool smoke:

```sh
python3 tests/integration/joint_schwarz_sweep.py \
  --build-dir build/joint-spqr --output build/joint-scaling.json \
  --csv build/joint-scaling.csv --topologies chain --atoms 8 16 \
  --cores 4 --overlaps 1 --preconditioners schwarz --operator-rank auto \
  --warmup 1 --repeat 3 --timeout 120
python3 tests/integration/joint_scaling_analysis.py \
  --input build/joint-scaling.json --output build/joint-scaling-analysis.json
```

For a future campaign, keep each local Schwarz policy fixed while increasing
global size:

```sh
python3 tests/integration/joint_schwarz_sweep.py \
  --build-dir build/joint-spqr --output build/joint-scaling.json \
  --csv build/joint-scaling.csv --topologies chain cube \
  --atoms 256 512 1024 2048 4096 8192 --cores 128 --overlaps 1 \
  --preconditioners diagonal schwarz --operator-rank auto \
  --warmup 1 --repeat 3 --timeout 600
python3 tests/integration/joint_scaling_analysis.py \
  --input build/joint-scaling.json --output build/joint-scaling-analysis.json
```

Run identity controls at smaller sizes if desired. Interpret stable iterations
with rising wall time through operator/setup costs, and stable iterations with
rising solve, linearization, or damping counts as nonlinear search work.
Rising RSS with stable iterations is a memory/resource-scaling issue. Only a
repeated increase in per-solve PCG iterations under a fixed local policy is
the main evidence for a future coarse-space investigation.

### Archived initial rank-blocked campaign

The initial chain-128 diagonal attempt remains archived at
[`joint-schwarz-scaling`](figures/joint-schwarz-scaling/campaign-manifest.json).
Its three rank evaluations returned `unavailable` with `rank-work-budget` and
zero PCG solves. The refreshed campaign and its current evidence are recorded
below.

### Bounded-rank work frontier

The rank-only frontier is recorded in
[`campaign-manifest.json`](figures/joint-rank-budget-frontier/campaign-manifest.json),
with its normalized records in
[`rank-frontier.json`](figures/joint-rank-budget-frontier/rank-frontier.json)
and [`rank-frontier.csv`](figures/joint-rank-budget-frontier/rank-frontier.csv).
At 100M, chain-128 stopped in reconstruction after 99,995,408 charged entries;
the exact total estimate was 339,953,610, with 239,958,202 remaining. The first
tested complete point was 500M, which returned the existing rigorous
`rank-verified-full` result. Thus 100M is about 3.4 times below the measured
certificate work, rather than marginally below it.

| Case | Budget | Rank result | Charged entries | Estimated total | Rank time | Peak RSS |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| chain-128 | 100M | `rank-work-budget` | 99,995,408 | 339,953,610 | 3.58 s | 87.9 MB |
| chain-128 | 150M | `rank-work-budget` | 149,998,500 | 339,953,610 | 6.02 s | 87.9 MB |
| chain-128 | 250M | `rank-work-budget` | 249,995,604 | 339,953,610 | 12.29 s | 87.9 MB |
| chain-128 | 500M | `rank-verified-full` | 339,953,610 | 339,953,610 | 19.66 s | 88.0 MB |
| cube-128 | 500M | `rank-work-budget` | 499,995,688 | 1,758,700,534 | 15.22 s | 157.8 MB |
| chain-256 | 500M | `rank-work-budget` | 499,994,398 | 1,363,814,378 | 19.21 s | 141.2 MB |

The rank-only profile performed no PCG solves or assessment. Compact
extractions and free-design SVDs were zero in all six measurements. Doubling
chain size increases estimated SPQR reconstruction work by 4.01x; cube-128 estimates
5.17x the chain-128 work at the same atom count. The 500M samples for cube-128
and chain-256 stopped at the budget, but reported the complete reconstruction
forecast. The census results are recorded separately in
[`local-rank-witness-census.json`](figures/joint-local-rank-witness-census/local-rank-witness-census.json),
with per-case SPQR and oracle runs in the same directory.

### Local-support witness census

The diagnostic-only census covered every free column with disjoint exclusive
rows for chain and cube at 128 and 256 atoms. Each case had 128 or 256 groups
of two columns; the minimum certified local lower bound was 0.0465905, while
the threshold upper bounds ranged from 1.44e-11 to 8.27e-11. Independent dense
oracles returned full rank at every size. The four SPQR census runs had zero
compact extractions and zero free-design SVDs.

The bounded rank route now tries this sufficient-only local certificate after
structural checks and falls back to SPQR interval reconstruction when the
local test does not certify. These rank measurements do not establish PCG
scaling or support a two-level Schwarz conclusion.

### Current Operator-PCG search scaling

The original formal search-only campaign is recorded under
[`joint-search-scaling-r1`](figures/joint-search-scaling-r1/campaign-manifest.json).
The cube multi-block gate was then completed with a preselected 768-atom point;
the copied baseline points and official analyzer output are in
[`joint-search-scaling-r2`](figures/joint-search-scaling-r2/campaign-manifest.json).
It used SPQR with `auto` resolving to `SpqrBounds`, the rigorous
`LocalSupport` certificate, the production-equivalent rank budget of 120
seconds, 100,000,000 work entries, and 256 MiB. The fixed Schwarz policy used
128 core atoms, one overlap hop, a 512-atom maximum block, 512 MiB storage, and
256 MiB scratch. It requested chain and cube sizes 128, 256, 512, 1024, and
2048, one warmup and three measurements, one Eigen thread, a 600-second run
limit, and a 4 GiB RSS ceiling.

The `search` profile calls the same `SearchProfile(...)` implementation as a
full solve, then serializes its search result and telemetry and exits before
assessment. Its JSON says `measurement_scope: "search-only"`; assessment and
endpoint qualification are `not-run`, and runtime convergence is never
reported as passed. The permanent test checks zero assessments and zero
reference evaluations. Eighteen of 20 configurations were search-side
eligible. The original COLAMD cube-2048 diagonal and Schwarz runs hit the RSS
ceiling during search; their substage attribution and benchmark-only ordering
follow-up are recorded below. Search-side eligibility makes no
endpoint-correctness claim.

The Schwarz measurements were:

| Topology | Atoms | Blocks | PCG iters/solve | Rank s | Setup s | PCG s | Search s | Peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| chain | 128 | 1 | 3 | 0.393 | 1.355 | 0.841 | 3.376 | 146.5 MiB |
| chain | 256 | 2 | 4 | 0.792 | 2.782 | 2.163 | 7.254 | 257.5 MiB |
| chain | 512 | 4 | 4 | 1.591 | 5.757 | 4.381 | 14.917 | 476.0 MiB |
| chain | 1024 | 8 | 4 | 3.188 | 12.111 | 8.776 | 30.563 | 631.0 MiB |
| chain | 2048 | 17 | 4 | 6.401 | 26.612 | 17.616 | 62.961 | 1139.0 MiB |
| cube | 128 | 1 | 3 | 0.294 | 1.454 | 3.371 | 7.273 | 331.8 MiB |
| cube | 256 | 2 | 6 | 0.566 | 4.799 | 16.986 | 28.288 | 617.0 MiB |
| cube | 512 | 4 | 6 | 1.097 | 13.852 | 45.579 | 75.673 | 1269.4 MiB |
| cube | 768 | 7 | 6 | 1.630 | 23.925 | 101.152 | 160.107 | 3.71 GiB |
| cube | 1024 | 9 | 7 | 2.176 | 39.902 | 177.522 | 274.031 | 3.22 GiB |
| cube | 2048 | — | — | — | — | — | RSS limit | 4.03 GiB |

The full-range cube result includes the single-block 128-to-multi-block-256
transition, so the separate actual-block-count gate is the basis for a
coarse-space decision:

| Topology | Multi-block sizes | Actual block counts | Valid points | Slope | Growth ratio | Multi-block gate |
| --- | --- | --- | ---: | ---: | ---: | --- |
| chain | 256, 512, 1024, 2048 | 2, 4, 8, 17 | 4 | 0.0000 | 1.00 | stable |
| cube | 256, 512, 768, 1024 | 2, 4, 7, 9 | 4 | 0.0841 | 1.17 | stable |

The cube-768 point completed all three measurements at 3.71 GiB peak RSS.
The four-point cube multi-block gate is `stable`; the analyzer's unchanged
thresholds are four points, slope 0.25, and growth ratio 1.5. Together with the
stable chain gate, **current synthetic evidence does not warrant a two-level
Schwarz investigation.** The cube-2048 RSS failures concern factor memory and
do not reopen the Krylov decision. No coarse correction was implemented.

Rank work was 1.8%–11.7% of search time on completed rows. Each of five
nonlinear linearizations used one rank check per measurement; the rigorous
local-support certificate removed the rank-work blocker without changing the
rank budget. Rank-certificate continuation is not warranted. The full
diagonal/Schwarz result set, eligibility records, failures, and both gates are
in the machine-readable
[`campaign artifacts`](figures/joint-search-scaling-r1/).

### Cube search memory frontier and SPQR ordering

The follow-up campaign is recorded in
[`cube-search-memory-r1`](figures/cube-search-memory-r1/campaign-manifest.json).
All four available SuiteSparse orderings were measured on cube-512 and
cube-1024 with the diagonal preconditioner, SPQR, OperatorPcg, `auto` rank
(`SpqrBounds`), production-equivalent rank budgets, one Eigen thread, a
600-second timeout, and the original 4 GiB RSS limit. `DEFAULT` matched
COLAMD fill. `BEST` and METIS materially reduced fill and process RSS while
preserving rank status/certificate, accepted updates, evaluations, stop
reason, PCG iterations, objective, gradient, residual, and returned state
within floating-point tolerance.

| Ordering | Case | Numeric/fixed factor nnz | Exported factor storage | Symbolic s | Numeric s | Peak RSS | Search result |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| COLAMD | cube-512 | 143,684 / 130,980 | 413.9 MB | 0.069 | 1.694 | 1.631 GiB | complete; 6.97 PCG/solve |
| DEFAULT | cube-512 | 143,684 / 130,980 | 413.9 MB | 0.068 | 1.675 | 1.901 GiB | complete; 6.97 PCG/solve |
| BEST | cube-512 | 98,080 / 87,951 | 271.0 MB | 0.123 | 1.027 | 1.439 GiB | complete; 6.97 PCG/solve |
| METIS | cube-512 | 98,080 / 88,687 | 273.6 MB | 0.085 | 1.016 | 1.438 GiB | complete; 6.97 PCG/solve |
| COLAMD | cube-1024 | 519,728 / 481,584 | 1,472.4 MB | 0.119 | 7.984 | 3.860 GiB | complete; 7.8 PCG/solve |
| DEFAULT | cube-1024 | 519,728 / 481,584 | 1,472.4 MB | 0.123 | 8.048 | 3.854 GiB | complete; 7.8 PCG/solve |
| BEST | cube-1024 | 277,240 / 255,930 | 744.3 MB | 0.258 | 3.205 | 3.377 GiB | complete; 7.8 PCG/solve |
| METIS | cube-1024 | 277,240 / 257,141 | 750.4 MB | 0.150 | 3.201 | 3.203 GiB | complete; 7.8 PCG/solve |

The factor column is numeric-factor nnz / fixed-operator factor nnz. Exported
storage is an owned-factor proxy; it is not process RSS or SuiteSparse scratch.
Sampled process peak RSS is reported separately.

At cube-2048, no ordering completed search inside 4 GiB:

| Ordering / preconditioner | Status | Peak RSS | Active stage | Last completed stage |
| --- | --- | ---: | --- | --- |
| COLAMD / diagonal | RSS limit | 4.019 GiB | `spqr-fixed-factor` | `profile-evaluation` |
| COLAMD / Schwarz | RSS limit | 4.028 GiB | `spqr-fixed-factor` | `schwarz-partition` |
| DEFAULT / diagonal | RSS limit | 4.007 GiB | `spqr-fixed-factor` | `profile-evaluation` |
| BEST / Schwarz | RSS limit | 4.474 GiB | `spqr-fixed-factor` | `schwarz-partition` |
| METIS / diagonal | RSS limit | 4.460 GiB | `spqr-fixed-factor` | `profile-evaluation` |
| BEST / diagonal | RSS limit | 4.143 GiB | `rank-certificate` | `fixed-operator-factor` |

For the COLAMD baseline, diagonal and Schwarz failed at the same SPQR fixed
factor stage at nearly the same RSS. Schwarz partition completed before its
factorization failure; local Schwarz build and PCG did not start. This is a
fixed-factor fill/RSS blocker, not a Schwarz-local-memory or Krylov-iteration
failure. BEST/diagonal progressed through PCG and trial evaluation before
exceeding 4 GiB during a later rank certificate.

One diagnostic-only BEST/diagonal run used a 6 GiB cap. It completed at
5.672 GiB in 395.7 seconds with four accepted updates, five profile
evaluations, 30 PCG solves, and 7.8 mean iterations per solve. It is excluded
from the formal 4 GiB gate and provides no Schwarz block-count point. No
cube-2048 Schwarz point was obtained. The later preselected cube-768 campaign
closed the separate Krylov question with a stable four-point multi-block gate;
the 2048 failures are only memory/resource evidence. Do not infer that Schwarz
failed, and do not begin a two-level Schwarz investigation.

### Factor residency and operator representation

The [factor residency campaign](figures/joint-factor-residency-r2/residency-analysis.json)
records at most two concurrently resident global factors at the sampled peak:
an older profile/trial workspace and the next operator factor being built.
This confirms overlap during accepted-state transitions. The instrumentation's
owned-byte estimate is separate from the process-tree RSS measured by the
runner; at a resource stop it can describe stage-entry factors while another
factor construction is in progress.

The [operator factor frontier](figures/joint-operator-factor-frontier-r1/factor-frontier-analysis.json)
compares the dedicated exported `Fixed` control with a benchmark-only native
SuiteSparseQR representation. At cube-512/COLAMD, the native run matched the
current run's stop reason, accepted updates, rank certificate, all 30 PCG
iteration counts, and returned search state (maximum absolute state difference
2.22e-16). Search time changed from 75.58 s to 64.55 s, process RSS from
2.04 GiB to 1.89 GiB, and concurrent owned-byte estimate from 646 MB to
449 MB. At cube-1024/BEST, the corresponding figures were 153.53 s to 132.78 s,
2.94 GiB to 2.62 GiB, and 1,167 MB to 815 MB, with the same discrete search
trajectory and maximum state difference 2.22e-16. These are search-only
comparisons; they do not measure returned-state assessment or endpoint trust.

For cube-1024/COLAMD, one current-representation run exceeded the 4 GiB cap
while the native candidate finished only about 29 MiB below it; an earlier
current run completed at about 4.12 GB. That variability does not establish a
reliable process-RSS margin for the production ordering. All six cube-2048
representation/ordering retries stopped at the unchanged 4 GiB cap before any
PCG solve. The production OperatorPcg path now reuses the accepted profile
factor with copy-on-write workspace mutation. Native QR remains benchmark-only,
production ordering remains COLAMD, and the dedicated fallback retains the
exported `Fixed` representation. The accepted-factor measurements are in the
[factor ownership campaign](figures/joint-factor-ownership-r1/factor-ownership-analysis.json).

### Returned-state assessment attribution

Assessment has a separate measurement record in
[`joint-assessment-scaling-r1`](figures/joint-assessment-scaling-r1/campaign-manifest.json).
The stage times below are medians of measurement repetitions; derivative is
preparation plus reduction. The 512 cases ran one measurement with no warmup.
They completed rank and PCG search, then timed out during assessment. A timeout
is an incomplete assessment, not a pass.

| Case | Primary s | Reference s | Design spectrum s | Derivative s | Width s | Normalized s | Correction/Jacobian s | Total assessment s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| chain-128 | 0.114 | 1.392 | 0.938 | 18.127 | 0.215 | 0.194 | 0.250 | 21.340 |
| chain-256 | 0.236 | 6.684 | 4.134 | 104.090 | 0.998 | 0.866 | 1.199 | 118.334 |
| chain-512 | 0.464 | 37.168 | 20.199 | 42.140 prep; reduction timed out | — | — | — | timeout |
| cube-128 | 0.155 | 1.602 | 1.233 | 15.658 | 0.270 | 0.257 | 0.307 | 19.548 |
| cube-256 | 0.361 | 7.171 | 4.631 | 86.834 | 1.073 | 0.961 | 1.264 | 102.458 |
| cube-512 | 0.853 | 44.477 | 27.003 | 54.228 prep; reduction timed out | — | — | — | timeout |

Derivative reduction dominated every completed assessment, accounting for
66.5%–79.7% of total assessment time; it was also the active stage at both
512 timeouts. Total assessment grew 5.55x for chain and 5.24x for cube from
128 to 256. The 512 timeout rows retain the completed stages, active stage,
stage dimensions, 600-second process wall, and peak RSS (581.8 MiB for chain;
2.64 GiB for cube). Their stage checkpoints and per-run records are available
in the assessment campaign JSON, CSV, analysis, and `individual-results/`.

### Derivative-reduction scaling

The finer attribution campaign is recorded in
[`joint-derivative-scaling-r1`](figures/joint-derivative-scaling-r1/campaign-manifest.json).
It used returned-state assessment with Schwarz, `auto` rank (`SpqrBounds`),
one Eigen thread, a 600-second timeout, and 4 GiB RSS. The 128/256 cases use
three measurements; 512 uses one. Completed rows show median inclusive
reduction and median exclusive substages. The timeout values are cumulative
completed-tile work, not complete stage totals.

| Case | Reduction total | Rows | Projected QR | Jacobian QR | Norms | Status / active stage |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| chain-128 | 15.73 s | 1.06 s | 6.94 s | 7.11 s | 0.177 s | complete |
| chain-256 | 94.69 s | 4.29 s | 43.69 s | 44.38 s | 0.716 s | complete |
| chain-512 | — | 14.03 s* | 225.90 s* | 228.74 s* | 2.223 s* | timeout; `derivative-projected-qr` |
| cube-128 | 13.11 s | 0.945 s | 5.76 s | 5.90 s | 0.146 s | complete |
| cube-256 | 76.23 s | 3.72 s | 35.02 s | 35.60 s | 0.576 s | complete |
| cube-512 | — | 12.55 s* | 187.06 s* | 188.56 s* | 1.814 s* | timeout; `derivative-projected-qr` |

`*` denotes cumulative times observed before timeout. Chain-512 completed 22
QR appends per role before the next projected QR timed out; cube-512 completed
18 per role and was also interrupted on projected QR. Peak RSS was 0.545 GiB
for chain and 2.392 GiB for cube, so these were timeouts rather than memory
failures. On completed 128/256 cases, projected QR was 43.9%–46.1% and
Jacobian QR was 45.0%–46.9% of reduction time; rows were 4.5%–7.2% and norms
were 0.76%–1.12%. Inclusive reduction grew 6.02x for chain and 5.81x for cube
from 128 to 256. The exclusive sum is slightly below inclusive time because
outer-loop and checkpoint work is not assigned to a micro-stage.

### Compact Jacobian diagnostics

The identity and test-only QR prototype diagnostics are recorded in
[`joint-compact-jacobian-r1`](figures/joint-compact-jacobian-r1/compact-jacobian-analysis.json).
For `J = P - Z C`, the full identity is
`J^T J = P^T P - P^T Z C - C^T Z^T P + C^T Z^T Z C`. Across chain/cube at
8/32/128/256, maximum `||Z^T P||/(||Z|| ||P||)` was 1.40e-15; the full Gram
relative Frobenius error was 1.51e-15, absolute Frobenius error 1.35e-16, and
maximum elementwise error 7.63e-17. Omitting the measured cross terms changed
the Gram by at most 1.83e-15 relative on these fixtures.

The test-only QR stack `[R_P; -R_Z C]` matched rank on all eight lattice cases
and on cancellation/active-face fixtures. Across lattice cases, maximum
singular-value relative error was 3.05e-15, weak-direction projector
difference 2.11e-11, response-gradient difference 1.11e-17, and correction
difference 4.20e-15. For a near-exact active-face fit, the residual-normalized
`||Z^T r||/(||Z|| ||r||)` ratio reached 0.478 because `||r||` was near zero;
observation-scale normalization gave 6.60e-16.

Near-threshold rank fixtures cover clearly deficient/full-rank cases and the
floating-point boundary on both sides. The compact path accepts a result only
when its rank decision is outside the one-ULP boundary enclosure; ambiguous
cases fall back to the exact observation-scale Jacobian TSQR. Full Assessment,
AssessmentEvidence, convergence, and returned-state parity passed on the
chain/cube 8/32 endpoints plus cancellation, active-face, and rank-boundary
fixtures. Rank, objective, and KKT thresholds were unchanged.

The [compact assessment acceptance campaign](figures/joint-compact-acceptance-r2/acceptance-analysis.json)
reported 1.41x–1.57x median assessment speedups on chain/cube 128/256, with
peak RSS changes from -0.1% to +3.4%. Both production-route 512 cases were
rerun under the same 600-second / 4-GiB envelope and completed: chain-512 used
454.88 s assessment time and 627 MB peak RSS; cube-512 used 410.41 s and
2.56 GB. The guarded compact reduction is therefore the production assessment
path, with the current observation-scale TSQR retained as the exact fallback.

### Post-compact assessment frontier

The [preselected 640-atom frontier](figures/joint-post-compact-frontier-r1/frontier-analysis.json)
used the production compact-Jacobian assessment route, SPQR, OperatorPcg,
Schwarz, COLAMD, one Eigen thread, and the unchanged 600-second / 4-GiB
limits. Both cases completed search and entered returned-state assessment.

| Case | Status | Active stage | Last completed | Projected QR s | Compact Jacobian QR s | Wall s | Peak RSS |
| --- | --- | --- | --- | ---: | ---: | ---: | --- |
| chain-640 | timeout | `derivative-projected-qr` | `derivative-norms` | 376.54 | 0.00 | 600.03 | 0.62 GiB |
| cube-640 | timeout | `derivative-projected-qr` | `derivative-norms` | 266.81 | 0.00 | 600.03 | 3.63 GiB |

The first post-compact end-to-end blocker is projected-width QR time, not the
4-GiB RSS ceiling. The cube case has about 0.36 GiB of sampled RSS headroom.
The preselected 640 point timed out, so this campaign ran no 768-atom solve.

### Projected-width reduction comparison

The [projected-width campaign](figures/joint-projected-width-r1/projected-width-analysis.json)
compares production observation-tiled QR with a benchmark-only fixed-order
SPQR factorization of `[Z D]`, with every `Z` column before every `D` column.
The candidate checks the identity ordering and extracts the trailing compact
factor; it does not use a free column permutation or normal equations. The
current route remains the production path. Permanent tests check projected
Gram, spectra, ranks, norms, normalized spectra, weak-subspace projectors,
response, cancellation, near-collinearity, row permutation, scale, active-face,
and rank-boundary behavior. Ambiguous rank-boundary fixtures use the current
observation-tiled QR fallback.

| Case | Current assessment s | Candidate assessment s | Assessment speedup | Current projected QR s | Candidate reduction s | Current / candidate peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| chain-128 | 14.39 | 7.64 | 1.88x | 7.14 | 0.42 | 0.21 / 0.69 GiB |
| cube-128 | 13.79 | 8.58 | 1.61x | 5.83 | 0.62 | 0.52 / 1.03 GiB |
| chain-256 | 75.03 | 33.63 | 2.23x | 44.18 | 2.39 | 0.32 / 2.04 GiB |
| cube-256 | 68.32 | 35.89 | 1.90x | 35.52 | 3.11 | 1.06 / 2.24 GiB |

Full assessment and returned-state parity passed at 128 and 256. Candidate RSS
grew materially with size. At 512 it exceeded the same 4-GiB limit during
`derivative-projected-qr` for both chain and cube; production completed at
454.88 s / 410.41 s and 0.58 / 2.39 GiB sampled peak RSS. No 640 candidate was
run after its 512 RSS gate failed. The candidate remains benchmark-only; the
observation-tiled projected QR is still the production implementation.

The fixed-order `[Z D]` candidate is numerically valid and performance-promising
at 128/256, but its 512 chain and cube runs exceed the 4-GiB RSS limit. It is
memory-rejected and remains benchmark-only; this is a resource result, not a
mathematical failure.

### Projected Tail QR

The [projected-tail census](figures/joint-projected-tail-r1/projected-tail-summary.csv)
reuses the free-design factor `Z = Q [R; 0]`, applies its sparse `Q^T` transform
to `D`, and keeps the bottom block `T` from `Q^T D = [B; T]`. Since the
projected derivative is `P = Q_2 T`, the candidate factors only `T`; it restores
the tail QR column permutation before using the compact factor. Census counts
are from the actual sparse transformed representation, not a conceptual dense
matrix.

| Case | n | p | m | D nnz | Tail rows | Tail nnz | Tail density | Q transform s | Census peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| chain-128 | 58,300 | 256 | 128 | 65,920 | 58,044 | 3,998,050 | 53.81% | 0.041 | 0.345 GB |
| cube-128 | 47,792 | 256 | 128 | 65,920 | 47,536 | 4,610,601 | 75.77% | 0.074 | 0.471 GB |
| chain-256 | 116,540 | 512 | 256 | 131,840 | 116,028 | 15,421,154 | 51.92% | 0.167 | 0.823 GB |
| cube-256 | 93,431 | 512 | 256 | 131,840 | 92,919 | 16,330,509 | 68.65% | 0.354 | 0.871 GB |
| chain-512 | 233,020 | 1,024 | 512 | 263,680 | 231,996 | 53,360,618 | 44.92% | 0.861 | 2.840 GB |
| cube-512 | 183,040 | 1,024 | 512 | 263,680 | 182,016 | 58,348,670 | 62.61% | 2.096 | 2.392 GB |

The [assessment frontier](figures/joint-projected-tail-a4-r1/projected-tail-assessment-frontier.json)
compares the same search endpoint and assessment settings. The numerical
candidate uses QR and fill-reducing ordering after eliminating `Z`; it does not
form normal equations.

| Case | Current assessment s | Tail assessment s | Speedup | Current projected QR s | Tail reduction s | Current / tail peak RSS | Parity/status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| chain-128 | 14.54 | 8.10 | 1.79x | 7.18 | 0.86 | 0.189 / 0.527 GB | passed |
| cube-128 | 13.85 | 8.85 | 1.57x | 5.86 | 0.87 | 0.336 / 0.673 GB | passed |
| chain-256 | 75.56 | 36.78 | 2.05x | 44.38 | 5.46 | 0.281 / 1.580 GB | passed |
| cube-256 | 68.41 | 37.70 | 1.82x | 35.48 | 4.86 | 0.641 / 1.644 GB | passed |
| chain-512 | 439.89 | 175.26 | 2.51x | 295.38 | 30.70 | 0.439 / 4.104 GB | completed on repeat; narrow RSS headroom |
| cube-512 | 404.83 | RSS limit | — | 233.09 | symbolic stage did not complete | 1.461 / 4.635 GB | candidate stopped in `projected-tail-symbolic` |

Permanent fixtures pass projected Gram, singular spectrum, rank, column-norm,
normalized-width, weak-subspace, response, `P^T r`, compact-Jacobian, and
correction parity within their asserted tolerances. Full assessment/search
observable parity passes at 128/256 and for the completed chain-512 repeat.
At cube-512, tail extraction completed, but the tail QR symbolic stage crossed
4 GiB; the sparse Q transform was not the limiting stage. The candidate is
therefore not promoted. Production keeps observation-tiled projected QR as the
exact route and fallback.

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
