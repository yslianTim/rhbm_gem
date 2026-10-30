# Joint bounded LM stability audit

Baseline: develop, 32609b2339736972c030a503778e868dc6647d67; clean tracked tree.
The ignored `output/diagnostics/7a6a_joint_components` evidence was available and read.

## Production responsibilities

| Responsibility | Implementation |
| --- | --- |
| FirstStage donor coverage, bounded initialization | `JointComponentEstimator.cpp`: `InitializeContributors`, `FitJointComponentsInitializedWithProgress`; `WidthCalibration.cpp` |
| Component partition and orchestration | `Components.cpp`: `OrchestrateFixedNeighborComponents`, `SolveFixedNeighborComponentView` |
| FixedNeighbor cores and effective local response | `FixedNeighborBlockCoordinate.cpp`: `SearchFixedNeighborComponent` |
| Constrained A/C profile, compact Width derivative | `Search.cpp`: `Profile::Get`, `linearize`; `Numerics.cpp`, `LinearSolve.cpp`, `TiledDerivative.cpp` |
| Bounded Width LM, damping, active bounds | `InstrumentedLM.hpp`: `minimizeOneStep`; `ActiveSetBoundedLM.cpp`: `SolveBoundedTrustRegionStep`, `SolveAtDamping` |
| Local candidate replay and acceptance | `Search.cpp`: `Profile::Trial`, `Accept`; initial A/C acceptance is separate from accepted Width updates |
| Global objective acceptance and safe state commit | `FixedNeighborBlockCoordinate.cpp`: candidate finite/feasibility/replay checks precede cache mutation; failed/rejected candidates retain accepted state |
| Recovery and sweep stationarity | same file: failure visits/retries, full sweep replay, original-coordinate A/C and projected Width KKT |
| Endpoint and independent verification | `CertifyFixedNeighborEndpoint`: streaming primary assessment, independent replay, trust and numerical qualification, memory preflight |
| Component qualification | `Components.cpp`: evidence/state creation; `JointComponentEstimator.cpp`: `JointComponentResult::RuntimeConvergence` |
| Assembly | `Components.cpp`: `AssembleComponents`: qualified component mapping, preserved coefficients, assembled replay/consistency |
| Downstream | `PotentialFittingWorkflow.cpp`: qualified snapshot gates Halo reconstruction, PostFitPeeling, Local and Group fitting |

## S0 baseline contracts

Controlled regression tests record (without changing production) the five-column
rank/bracket failure, near-unit Width radius collapse and coordinate-origin
discontinuity, and historical LM failure veto. Subsequent stages replace these
characterization assertions with the corrected contracts. Existing recovery,
rollback, nonfinite, endpoint and replay tests remain mandatory.

All build/test commands use at most four threads. Real-data runs use new temporary
SQLite files; diagnostic programs, databases and large output remain untracked.

## S1: increment coordinates and null sensitivity

The existing augmented QR solves in `z_i=S_i p_i`, with
`S_i=max(D_i, ||J_i||_2, sqrt(min_normal))`. The damping block is
`sqrt(lambda) D_i/S_i`; conversion back divides the solution by `S_i`.
Bounds, bound hits/releases, multipliers and KKT are checked directly in `p`.
No transformed bound tolerance is used. Scaling fills the existing augmented
workspace and adds only O(number of Widths) storage, with no dense fallback.

A direction is held at zero for this linearization only when its column norm is
below `64 epsilon` times the largest column norm AND its original gradient is at
most `64 epsilon` (normalized objective units). The full candidate gradient is
checked again, including held directions. This does not certify identifiability;
a later A/C or neighbor change recomputes the criterion. Small sensitivity with
nonstationary gradient is explicitly tested and remains movable.

The controlled five-Width case now takes one QR and zero damping brackets rather
than 81 failed QRs. Genuine rank deficiency still uses the existing bounded
positive-damping QR search. No SVD or additional production solver is needed.
S1 validation: Release build; 69 tests across ActiveSetBoundedLM, Search, Runtime,
and FixedNeighbor (including rollback/recovery) passed.

## S2: origin-independent trust radius

The initial radius is `factor * ||D||_2`, the weighted norm of a representative
unit log-width increment (production factor 0.1). It is independent of eta's
origin and limits the initial aggregate model displacement; bounds, replay,
weak-direction safeguards and actual objective reduction still control updates.
The trust-ratio regression explicitly selects a radius that reaches its active
bound rather than depending on the absolute starting coordinates.

Tiny steps use the original log-width increment infinity norm, not `||D p||`.
Production tiny-step success requires the original projected Width KKT plus valid,
feasible A/C. Nonstationary tiny steps return a numerical stop for block recovery.
A cosine-based early stop also cannot bypass production projected KKT.
The radius regression covers B=1, B=1.000000000626 and shifted reference origins,
finite evaluated candidates, objective improvement and final bound-aware KKT.
S2 validation: Release build and all 75 focused tests passed.

## S3: terminal qualification and recovery

`SearchResult::width_failed` separates LM subproblem/nonstationary tiny-step
failure from initial A/C profile acceptance. FixedNeighbor routes these failures
through the existing retry/rollback path even when the initial profile is trusted.
Failure logs retain initial profile trust, Width stop reason and solver history;
progress failed-block counts now include these block visits.

`accepted-state-integrity` is required component evidence: all accepted mutations
and sweep replays must remain intact, and no block recovery may remain pending.
Its value is the unresolved block count; its reason retains recovery attempts and
successes. Fatal stops cannot create qualified component state. Historical
`active-set-lm-subproblems` failures remain visible (including Failed status) but
are no longer themselves a convergence veto. Complete sweep stationarity, both
KKT gates, endpoint verification, independent replay, feasibility and terminal
integrity remain required. Missing integrity evidence rejects live assembly.

The existing extensible evidence serialization carries the new check without
changing JSON structure or introducing compatibility paths. The diagnostic test
facade now forwards actual terminal integrity instead of omitting it. Regression
coverage includes historical failure with qualified assembly, failed/missing
terminal evidence, endpoint/replay/KKT rejection, nonfinite state, successful and
exhausted recovery, and trusted initial A/C followed by an injected Width failure.

S3 validation: 77 focused tests and all four Joint CTest groups (contract,
workflow, numerical, persistence) passed. Current endpoint and replay veto
coverage is maintained in the FixedNeighbor and Numerics regression tests.
Production diagnostics also retain damping solve, bracket and QR totals in the
existing LM evidence reason.

## S4: final integration and regression evidence

Final Release build, four-thread ceiling throughout. Full CTest: **19/19 passed**
(1,062 GoogleTest executions plus the CLI/script checks). Final command:

```sh
OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 VECLIB_MAXIMUM_THREADS=4 \
  build/release/bin/RHBM-GEM potential_analysis --estimator joint-components \
  -a /Users/yslian/data/7a6a.cif -m /Users/yslian/data/emd_11668_unsharpened.map \
  --select-chain A -j 4 -d /tmp/rhbmgem-jce-fix/final.sqlite -k 7a6a-final
```

The final audit also found that absolute trust feasibility slack could permit an
oversized physical step when D and the radius are tiny. Feasibility slack is now
relative to the radius. Two regressions require a small weighted norm to retain a
meaningful physical step AND enforce the actual radius. The complete suite and
real case were rerun after this change.

| Metric | Baseline 32609b23 | Final |
| --- | ---: | ---: |
| FirstStage qualified donors | 1410/1418 | 1410/1418 |
| Initial Width coverage | passed | passed |
| FixedNeighbor blocks | 136 | 136 |
| Search sweeps | 200, budget exhausted | 111, block-stationary |
| A/C KKT | 5.73634534e-10 | 6.86763424e-8 |
| Independent projected Width KKT | 6.46334942e-6 | 9.57221247e-10 |
| Width upper bound (angstrom) | 3.0660309206224183 | unchanged |
| LM solves / failed subproblems | 27157 / 195 | 15104 / 0 |
| Recovery attempts / successes / unresolved | historical events obscured | 0 / 0 / 0 |
| Endpoint / Independent Replay | failed / failed | passed / passed |
| Available components | 0/1 | 1/1 |
| CLI result | 3 (application log) | 0 (captured subprocess return code) |
| Full-process wall seconds | 131.46 | 89.68 |
| CPU seconds (user + system) | 134.88 | 93.50 |
| Search seconds | 121.08 | 69.07 |
| Peak RSS bytes | unavailable | 1867923456 |

Baseline timing used macOS `/usr/bin/time -lp`: its final RSS query was denied by
the sandbox (`sysctl kern.clockrate`), so the wrapper returned 1 despite the
application's recorded CLI failure 3. Final timing used Python subprocess wait
and `getrusage(RUSAGE_CHILDREN)` (macOS byte units). Baseline timing overlapped
some initial test compilation; wall/CPU results are single-run comparisons, not
a controlled speedup claim. The earlier final-algorithm run also qualified in
111 sweeps (89.70 seconds, 1748434944-byte RSS). A baseline full-process RSS delta
cannot be established from these measurements; this comparison remains partial.

Same-compiler isolated solver benchmark (five samples, ten calls each):

| Subproblem | Baseline | Final |
| --- | ---: | ---: |
| Normal 4096x8 median microseconds | 12761.0 | 12838.6 |
| Normal QR / damping / bracket counts | 112 / 47 / 4 | 112 / 47 / 4 |
| Degenerate 5x5 median microseconds | 47.71 (failure) | 1.18 (valid) |
| Degenerate QR / damping / bracket counts | 81 / 81 / 80 | 1 / 1 / 0 |

The normal-case difference is about +0.6%; no meaningful normal-case slowdown was
observed. Final full-case diagnostics retain 31,859 damping solves, 6,507 bracket
solves and 31,921 QR factorizations. Baseline global QR/damping totals were not
serialized, so only the controlled subproblem counts are directly comparable.
No extra solver family or dense fallback was added; additional solver storage is
O(block Width count), not O(observations x Width count).

Original-endpoint local replay with the fixed production solver: LYS50 block 96
accepts a Width update (two profiles, zero failed subproblems, one QR); GLU163
block 55 now evaluates and accepts an update for both near-unit and exact-unit
initial widths. These are local regression evidence, not substitutes for the
full final run. Final original-coordinate gradients:

| Atom | Final Width | Absolute gradient |
| --- | ---: | ---: |
| LYS50 CB | 0.65880626 | 2.47e-13 |
| LYS50 CG | 0.51545363 | 1.09e-14 |
| LYS50 CD | 0.65307241 | 5.54e-15 |
| LYS50 CE | 0.65814955 | 9.34e-15 |
| LYS50 NZ (A=0) | 0.01396310 | 1.11e-75 |
| GLU163 OE2 (A=0) | 0.03382466 | 4.68e-24 |

The largest final projected Width KKT belongs to ASN12 CG. The raw gradient
infinity norm is 4.55916209e-8 because one Width is bound-active; it is not the
projected KKT used for qualification. All 1418 A/C/Width parameter states are
finite, A is nonnegative and Widths are positive and below the original bound.
Primary and independent projected KKT agree. Model/map hashes, normalization,
179949 fitting rows and 60636 excluded pure-Halo rows match the baseline. All
fitting rows have available predictions; excluded original rows remain excluded.

### Downstream output

All 1418 formal Joint target estimates passed. Halo reconstruction used 1400
eligible target donors and reconstructed all 734 halos. Peeling covered every
target (141800/141800 samples). The persisted run IDs link the Joint estimates,
peeling, Local results and Group summaries to the same run.

Local fitting completed: 1413 available, 5 unavailable (4 invalid-local-parameters,
1 local-mdpde-unqualified). Group fitting completed: 147 available inferences,
26 unavailable summaries (21 single-member, 5 singular-group-covariance). These
are downstream qualification limits, not missing execution; no claim is made
that every downstream fit succeeded. NZ/OE2 A=0 Widths remain unidentifiable;
passing optimization KKT does not confer statistical confidence.

### Coverage and remaining limits

Existing tests exercised constrained A/C, upper-active Width, weak Width,
rank deficiency, near singularity, high overlap, signed cancellation, noise-only,
objective acceptance/rejection, recovery/rollback, nonfinite objective, memory
fallback, endpoint verification and independent replay. Three seed widths are
covered by the existing initialization regression; one full real dataset was run.
Different maps and broad multistart behavior remain untested.

The first full-suite attempt imposed OMP_THREAD_LIMIT=1 and therefore failed
existing tests that explicitly require 2/4 threads. Rerunning all CTests serially
with a four-thread OpenMP ceiling passed every test, without changing those tests.

Selective block refinement was evaluated as unnecessary: the generic fixes
reached the original thresholds in 111/200 sweeps with zero unresolved failures.
It was not added. No diagnostic prototype, SQLite database or large output is
tracked. Raw evidence and compact validation are retained under
`/tmp/rhbmgem-jce-fix/`; the original ignored diagnostics and source datasets were
preserved. Baseline RSS and global factorization-count comparisons remain the
explicit performance evidence limits.
