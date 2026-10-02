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
The `rank` profile accepts the same three budget controls for the bounded
prototype. These options belong to the internal benchmark drivers and do not
add production CLI switches.

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
rank status and resource use. Bounded rank diagnostics include the work stage,
charged entries, estimated total and remaining entries when the reconstruction
forecast is available, and the design, R-factor, reflector-nonzero and
reflector counts. Schwarz runs also include a partition summary
with block count, core and block size distributions, atom membership
distribution, graph workspace, storage, and scratch bound. Rank counters show
checks, time, entries, workspace, compact extraction count, and free-design
SVD count, so the bounded SPQR route can be checked directly.

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

## Current one-level scaling evidence

The current resource-bounded campaign is recorded under
[`joint-schwarz-scaling`](figures/joint-schwarz-scaling/campaign-manifest.json).
It requested chain and cube sizes 128, 256, 512, and 1024 with SPQR, rank mode
`auto`, one warmup, three measurements, a 600-second per-run timeout, and a
4 GiB process-tree RSS ceiling. The Schwarz policy was fixed at 128 core atoms,
one overlap hop, a 512-atom maximum block, 512 MiB storage, and 256 MiB scratch.

The first chain-128 diagonal control completed all three measurement
processes, but bounded rank was `unavailable` with `rank-work-budget` in every
repetition. Each process recorded zero PCG solves, so this row is ineligible for
iteration scaling. Its returned-state runtime convergence was `failed`; this
records an unqualified endpoint after search stopped at the rank budget, not a
PCG failure. The campaign stopped there. The remaining 15 configurations,
including every Schwarz configuration and the small identity controls, were
not run. No larger size was attempted and no rank budget was changed.

| Topology | Requested Schwarz sizes | Valid PCG points | Iteration slope | Endpoint ratio | Coarse gate |
| --- | --- | ---: | ---: | ---: | --- |
| Chain | 128, 256, 512, 1024 | 0 | — | — | `insufficient-evidence` |
| Cube | 128, 256, 512, 1024 | 0 | — | — | `insufficient-evidence` |

The one completed diagnostic row had median operator-linearization setup 3.787
s, search time 3.949 s, assessment time 21.911 s, and wall time 26.476 s. Peak
RSS was 178,749,440 bytes (about 170.5 MiB). It performed one linearization,
zero damping trials, zero accepted updates, and one profile evaluation per
measurement. These are single-configuration diagnostics, not scaling results;
PCG time and operator action counters were zero because no PCG solve ran, so
they are not solve-cost measurements. No iteration slope can be estimated.

The normalized rows and repetition statistics are in
[`scaling-summary.json`](figures/joint-schwarz-scaling/scaling-summary.json)
and [`scaling-summary.csv`](figures/joint-schwarz-scaling/scaling-summary.csv);
the gate output is in
[`scaling-analysis.json`](figures/joint-schwarz-scaling/scaling-analysis.json),
and the individual completed run retains full provenance. Current evidence is
resource-limited before a Krylov verdict can be established. It does not show
stable Krylov scaling or iteration growth and does not warrant a two-level
Schwarz investigation.

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
