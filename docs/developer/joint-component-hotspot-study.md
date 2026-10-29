# Joint Component Hotspot Study

## P0 baseline

The baseline was collected on `develop` at `e78f738d74b366363921118bae7ff10325648f5b`; the worktree was clean and `origin/develop` matched. The build is Release, AppleClang 21.0.0, Eigen 5.0.1, CMake 4.4.3, macOS ARM64 on an Apple M1 (T8103). Builds used `-j4`; benchmark and numerical-library thread limits were one. No dependency or Notion content changed.

`/Users/yslian/data/6z6u.cif` and `/Users/yslian/data/emd_11103_unsharpened.map` existed. The real run selected chains A and YA. The synthetic closed-loop benchmark and the new test-only attribution benchmark each used one warmup and five measured runs. Measurements below are medians with min–max ranges. Real runs were serial, each in a fresh database. The real output had 1,557 targets, 2,166 atoms including halo contributors, and 243,294 fitting observations.

The fine-grained attribution is behind `RHBM_GEM_TEST_INSTRUMENTATION` and is disabled by default. The test benchmark splits Eigen SparseQR `compute()` into its equivalent `analyzePattern()` plus `factorize()` calls to time them separately; `solve()` remains separate. This test-only measurement path retained the production output hashes, objective, KKT, qualification, and independent replay. Production continues to call `compute()`.

### Existing telemetry and timer scope

Production already updates `FixedNeighborWork`, `ProfileSearchWork`, `LinearSolveWork`, and `BoundedLmWork` when telemetry collection is enabled. Counts cover candidate replay paths, profile evaluations/derivatives, active-set work, SparseQR structure/factorization, refinement, and LM damping. P0 adds only test-build counters and timers for inner stages and pattern fingerprints.

Timers are nested. `candidate_replay` includes incidence preparation, prediction/residual, and objective delta; `incidence_preparation` includes ordering, lookup/count, prefix construction, and fill. Linear solve time includes row reduction and SparseQR; row reduction includes tiled row QR and sparse assembly; SparseQR inclusive time includes analysis, numeric factorization, and solve. Derivative preparation includes raw construction, tiled QR, rank SVD, and triangular solves; correction-matrix time itself includes its triangular solves. LM damping-search time includes active/free selection, augmented assembly, QR, solve, and KKT certification. These inclusive timers are not added together.

Optional attribution added about 0.36% median to the 300-target instrumented benchmark (1.16574 s on vs 1.16150 s off; ranges 1.16452–1.18539 and 1.16091–1.16440 s). Smaller sizes showed 0.21% lower, 0.03% higher, and 0.16% higher medians respectively, within run variation. This comparison isolates the new test-only timing/census work; it does not disable the existing coarse timers. The production executable has no new detailed timers.

### Synthetic baseline

| Targets | FixedNeighbor wall median (range) | Search median (range) | Profile search | Candidate replay | Candidate / wall | Attempts / affected | Sweeps; endpoint |
|---:|---:|---:|---:|---:|---:|---:|---|
| 12 | 39.202 ms (38.751–40.228) | 31.133 ms (30.720–31.889) | 25.273 ms | 2.130 ms | 5.43% | 4 / 4 | 4; qualified, independent replay passed |
| 48 | 186.981 ms (186.881–193.137) | 154.739 ms (154.595–159.787) | 124.902 ms | 11.330 ms | 6.06% | 20 / 20 | 5; qualified, independent replay passed |
| 240 | 923.514 ms (922.646–924.539) | 765.112 ms (764.721–766.568) | 614.637 ms | 57.658 ms | 6.24% | 105 / 105 | 5; qualified, independent replay passed |
| 300 | 1.15554 s (1.15473–1.15739) | 957.371 ms (956.570–958.538) | 768.671 ms | 72.239 ms | 6.25% | 130 / 130 | 5; qualified, independent replay passed |

Candidate copy, unchanged fast paths, and Full Replay memory fallbacks were effectively zero/zero/zero on these converged synthetic runs. Sweep-end Full Replay and replay audit were 25.64 ms and 0.253 ms at 300 targets; cache update was 0.976 ms. The candidate scratch peak was 233,512 bytes.

At 300 targets, the 130 affected replays touched 690,300 rows total (970–5,520 per candidate; 0.71–4.04% of 136,560 global rows) and 787,500 contributor entries (1,090–6,300 per candidate). Candidate copy was 0.036 ms; replay was 72.248 ms, split into 46.033 ms incidence preparation, 25.465 ms prediction/residual, and 0.696 ms objective delta. Preparation was 0.014 ms ordering, 22.393 ms row lookup/count, and 23.538 ms prefix/contributor construction. Thus structure preparation was 63.7% of candidate replay, with ordering negligible and the two support traversals about equally costly.

If all candidate replay vanished, Amdahl's bound is only 1.057–1.067× across the four synthetic sizes. If all incidence preparation vanished, the bound is 1.034–1.041× for complete FixedNeighbor wall time (3.3–4.0% maximum time reduction). Candidate replay is worth measuring, but incidence caching cannot meet the preferred 5% complete-workflow gate by itself.

### A/C SparseQR and refinement

For 300 targets there were 260 valid A/C solves, sparse structure builds, and SparseQR factorizations; 260 active-set iterations and zero releases. The linear solves took 289.263 ms total. Column normalization took 4.848 ms, selected sparse structure construction 6.741 ms, and ReduceSparseRows 233.292 ms. Within row reduction, row-major conversion took 4.181 ms, tile mapping 1.879 ms, tiled sparse-row QR 220.748 ms, and reduced sparse assembly 2.245 ms.

SparseQR analysis, numeric factorization, and solve took 2.915, 6.113, and 0.699 ms total. Together they were 3.37% of A/C linear-solve time and 0.84% of FixedNeighbor wall time; numeric factorization alone was 2.11% of linear-solve time. Initial and final projected KKT maxima were both `6.94e-17`.

The selected matrix pattern repeated at all 130 observed transitions. The actual reduced matrix passed to SparseQR repeated only 5/130 times; it changed 125/130 times. Values changed on all five same-pattern transitions. Selected dimensions were 970–5,520 rows, 4–24 columns, and 2,060–12,360 nonzeros; reduced dimensions were 4–144 rows and 10–782 nonzeros. Matching free faces or selected patterns therefore does not justify SparseQR factorization reuse.

Refinement was eligible 260/260 times and triggered 0 times; success/failure and refinement QR counts were both zero. This workload gives no reason to optimize or weaken conditional refinement.

### Profile derivative and bounded LM

The 300-target run made 130 derivative preparations and 130 reductions. It used 130 Tiled QR calls, 130 rank SVDs, and 390 triangular solves. Preparation took 235.939 ms and reduction 88.338 ms, together 41.6% of profile-search time and 27.8% of complete FixedNeighbor wall time. Tiled QR took 201.273 ms, raw construction 25.491 ms, rank SVD 10.497 ms, and triangular solves 0.445 ms. No cancellation fallback occurred in this search. A controlled cancellation-sensitive fixture triggered exactly one reference-order fallback on every measured run; its extra work had a 1.25 μs median (1.209–1.292 μs). The fallback is retained.

The 130 bounded-LM subproblems had no failures or releases. They used 1,754 damping solves/QR factorizations: 130 zero-damping solves, 546 bracket solves, and 1,078 bisection solves. The per-update histogram was 104 updates with one damping solve, 9 with 32–63, and 17 with 64–127. All damping search took 4.333 ms, dense QR 2.696 ms, and QR solve 0.580 ms: respectively 0.37%, 0.23%, and 0.05% of complete FixedNeighbor wall time. Even eliminating all LM work has a 1.004× workflow speedup ceiling.

The separate 4,096×8 bounded-LM benchmark took 12.467 ms per call (12.448–12.621), with 47 damping solves and 112 QR factorizations. A 30-target profile benchmark took 200.274 ms (200.230–208.512), including 99.924 ms derivative preparation, 30.257 ms derivative reduction, 68 damping solves, and 68 damping QR factorizations. This illustrates that derivative work remains material at smaller dimensions too, while bounded LM is not the end-to-end bottleneck.

### Real workload baseline

One warmup plus five measured runs produced identical endpoint values and trajectory counts. The measured rows were all 118 sweeps, 35,828 profile evaluations, and 14,680 accepted updates.

| Metric | Median (min–max) |
|---|---:|
| Complete CLI wall | 134.471 s (134.224–134.711) |
| Joint search | 108.909 s (108.869–109.357) |
| Peak RSS | 1,484,701,696 B (1,484,177,408–2,034,909,184) |
| Endpoint verification memory estimate | 52,661,461 B (same in all runs) |

All runs qualified and passed Independent Replay. Objective was `0.02195348206329014`, independent A/C KKT `4.39461e-9`, projected-width KKT `6.14044e-11`, and independent KKT consistency `6.94e-18`. The RSS range contains one high sample; preflight estimates were stable. The saved result reports 1,541/1,557 qualified width-calibration donors and 243,294 fitting rows. Fine-grained test-only stage counters are unavailable in the production CLI, so the real-workload attribution is limited to total search/stage costs and saved qualification evidence.

## P0 gate

The current priority is: (1) derivative preparation/reduction and its Tiled QR work; (2) ReduceSparseRows' tiled QR inside A/C solves; (3) Candidate Replay incidence preparation; then (4) SparseQR symbolic/numeric/solve time and bounded LM. Refinement was not triggered. The reduced SparseQR pattern changes on 96% of observed transitions, so symbolic reuse is not promising and numeric reuse is invalid because values change even on repeated patterns. Candidate replay is about 6.3% of complete synthetic FixedNeighbor time, but incidence preparation alone is capped near 4% workflow improvement. These gates drive the P1–P3 experiments recorded below.

## P1–P3 decisions and results

### P1 FixedNeighbor Candidate Replay

R1 confirmed that each affected replay rebuilt its sorted contributor incidence and performed two support traversals. At 300 targets this initially cost 46.033 ms of 72.248 ms candidate replay (63.7%); row lookup/count and contributor construction were each about 23 ms. Candidate replay was 6.25% of complete FixedNeighbor wall, so eliminating it entirely could save at most that amount.

R3 replaced the two per-support `lower_bound` lookups with direct indices in the existing `row_position` vector. After all block mappings are prepared, its informative-row entries are cleared and the vector is populated only for the last replayed block; the preceding block's entries are cleared before reuse. Atom ordering, contributor ordering, prediction/residual arithmetic, objective calculation, acceptance, cache commit, and independent replay remain unchanged. The vector already existed for block preparation and is accounted for in the global-state preflight estimate (`component_observations * sizeof(Eigen::Index)`); no allocation or preflight estimate was added. The largest candidate replay scratch remained 233,512 bytes.

In five independent 300-target benchmark processes, lower-bound versus direct-map medians were 1.165115 s versus 1.125677 s complete FixedNeighbor wall (−3.38%), 72.699 ms versus 32.780 ms candidate replay (−54.9%), and 87,851,008 B versus 87,801,856 B peak RSS. All five candidate-decision digests matched. The P1 production benchmark then compared one warmup plus five measured runs at each synthetic size:

| Targets | Candidate replay baseline → final | FixedNeighbor wall baseline → final | Wall change | Search baseline → final | Search change |
|---:|---:|---:|---:|---:|---:|
| 12 | 2.130 → 0.975 ms (−54.3%) | 39.202 → 37.812 ms | −3.55% | 31.133 → 29.759 ms | −4.41% |
| 48 | 11.330 → 5.001 ms (−55.9%) | 186.981 → 181.182 ms | −3.10% | 154.739 → 148.864 ms | −3.80% |
| 240 | 57.658 → 25.438 ms (−55.9%) | 923.514 → 893.422 ms | −3.26% | 765.112 → 735.608 ms | −3.86% |
| 300 | 72.239 → 32.056 ms (−55.6%) | 1.15554 → 1.12232 s | −2.87% | 957.371 → 921.170 ms | −3.78% |

All four sizes retained identical candidate decision digests, eta/beta hashes, objective, KKT, sweep counts, qualification, and Independent Replay. Preflight estimates and candidate scratch did not increase. Independent-process RSS on 300 targets was unchanged within noise. The smaller-than-5% synthetic workflow gains were stable across four sizes and came with a greater-than-50% reduction in the targeted stage, so this minimal reuse was promoted.

The real 6Z6U/EMD-11103 run with chains A and YA used one warmup and five fresh-database measured runs. It retained 1,557 targets, 118 sweeps, 35,828 profile evaluations, 14,680 accepted updates, objective `0.02195348206329014`, independent A/C KKT `4.39461e-9`, projected-width KKT `6.14044e-11`, endpoint qualification, and Independent Replay. CLI wall median was 114.499 s (114.422–118.010) versus 134.471 s (134.224–134.711) at P0, a 14.85% reduction; Joint search was 89.174 s (89.131–92.517) versus 108.909 s (108.869–109.357), an 18.12% reduction. Peak RSS was variable: final median 1,563,672,576 B (1,484,242,944–1,793,097,728) versus baseline median 1,484,701,696 B (1,484,177,408–2,034,909,184); the minima are essentially identical, the baseline had the larger maximum, and the preflight estimate stayed exactly 52,661,461 B. Since the change reuses an already-preflighted vector without adding storage, the RSS spread is treated as process-level measurement variation rather than an allocation increase.

R2 incidence caching was skipped after R3: direct mapping reduced all incidence preparation to 5.010 ms at 300 targets (0.45% of full FixedNeighbor wall), leaving less than a 0.5% workflow ceiling while a reusable cache would retain per-block structures. R4 candidate copying was also skipped: its measured ~0.04 ms was below 0.005% of complete wall. The lower-bound A/B selector and all test-only selector fields were removed; direct mapping is the sole production path.

### P2 Constrained A/C SparseQR and refinement

Q1/Q2 completed from the P0 stage breakdown. The 300-target case had 260 valid A/C solves and 260 sparse structure builds. Column normalization was 4.848 ms; selected sparse assembly 6.741 ms; row reduction 233.292 ms; SparseQR analysis/factorization/solve 2.915/6.113/0.699 ms. `ReduceSparseRows` dominated because tiled Householder work took 220.748 ms, while SparseQR itself was only 0.84% of complete FixedNeighbor wall. The free face and selected matrix pattern repeated across all 130 profile transitions, but the reduced pattern repeated only 5/130 times and changed 125/130 times. Numeric values changed on all five repeated-pattern transitions.

Q3 symbolic reuse was skipped: the actual factorization input rarely repeats, and a same-pattern fingerprint does not imply reusable numeric factors. The measured analysis time was only 2.915 ms total, so caching would add invalidation risk for a small maximum saving. Q5 refinement optimization was skipped because it was eligible on 260 solves but triggered 0 times; all initial KKT checks passed at `6.94e-17`. Refinement remains intact and strict KKT is unchanged.

Q4 tested reusing one `Eigen::HouseholderQR<MatrixXd>` object across `ReduceSparseRows` tiles, replacing per-tile construction while preserving the same `compute` arithmetic. On the 300-target instrumented benchmark, representative tiled-QR time changed from 222.881 ms to 215.794 ms (−3.2%); five-run complete search-wall medians changed from 1.125087 s to 1.120985 s (−0.36%). Numerical output and qualification passed, but the localized gain was below the 10% preference and the end-to-end change did not approach the 5% gate. The experiment was reverted. No P2 production optimization was promoted.

### P3 Profile Derivative and bounded LM

P3-D1 completed: 130 prepares and 130 reductions cost 235.939 and 88.338 ms at 300 targets. Prepare consisted of 25.491 ms raw derivative construction, 201.273 ms Tiled QR, 10.497 ms rank SVD, and 0.445 ms triangular solves; correction-matrix construction was 0.301 ms and its solve work is already included in the triangular-solve total. Reduction's derivative QR work is included in its 88.338 ms total. These counters are from the instrumented representative sample; the five-run stage and full-workflow medians/ranges are reported separately. The test counter labeled `workspace_resize_events` counts selected shape/setup events, including sparse `setFromTriplets` calls; it is not a process heap-allocation counter and is not used to claim a measured allocation count.

P3-D2 completed: the normal 300-target search triggered no cancellation fallback in 130 derivative preparations. A cancellation-sensitive controlled fixture triggered the retained reference-order fallback in all five measured runs; incremental extra time median was 1.25 μs (1.208–1.291 μs). The fallback and its numeric reference-order accumulation remain unchanged.

P3-D3 was skipped after inspection: `TiledQR` already reuses its tile scratch matrices within a derivative preparation, and `ReduceDerivative` reuses its local tile matrices. The 300-target component derivatives fit in one 8,192-row tile each, so cross-tile reuse did not apply to this synthetic workload. A profile-lifetime cache would require a broader API/lifetime change without evidence of repeated large-profile workspace creation. No persistent derivative workspace or memory-preflight change was introduced.

P3-L1/L2 completed: the 130 bounded-LM subproblems used 1,754 damping solves and the same number of QR factorizations (130 zero-damping, 546 bracket, 1,078 bisection). Damping search cost 4.333 ms, QR factorization 2.696 ms, and QR solve 0.580 ms, totaling 0.37%, 0.23%, and 0.05% of complete FixedNeighbor wall. The per-update QR histogram had 104 updates with one factorization, 9 with 32–63, and 17 with 64–127. No subproblem failed or released an active bound in this workload. Different damping values produce different augmented matrices; no evidence showed exact repeated numeric subproblems.

P3-L3 was skipped because the implementation already reuses bounded-LM augmented matrix, index vectors, and workspace. Assembly and active/free selection together were below 0.02% of workflow, with no measured repeated allocation cost. P3-L4 QR caching was skipped because no exact same free face, numeric matrix, damping, and threshold subproblem was identified; matching dimensions or active faces alone are insufficient. No P3 change passed the production gate, and Normal Equations remain absent.

## P4 Final validation and production audit

The 300-target comparison below uses the P0 instrumented Release baseline and the final production FixedNeighbor benchmark. Stage timings are nested and do not sum; SparseQR and refinement/derivative/LM substage times are from the instrumented representative sample. The full wall/search/candidate times are five-run medians. `stage_speedup` is baseline time divided by final time.

| Metric | Baseline | Final | Change / speedup |
|---|---:|---:|---:|
| FixedNeighbor complete wall | 1.15554 s | 1.118539 s | −3.20%; 1.033× |
| FixedNeighbor search | 957.371 ms | 920.006 ms | −3.90%; 1.041× |
| Candidate Replay | 72.239 ms | 31.995 ms | −55.7%; 2.258× |
| A/C SparseQR inclusive | 9.753 ms | 9.778 ms | +0.3%; 0.997× (unchanged solver) |
| Conditional refinement | 0 / 260 triggers | 0 / 260 triggers | Not triggered |
| Derivative prepare + reduce | 324.277 ms | 323.988 ms | −0.1%; 1.001× |
| Bounded-LM damping search | 4.316 ms | 4.328 ms | +0.3%; 0.997× |
| Profile evaluations | 266 | 266 | unchanged |
| SparseQR factorizations | 260 | 260 | unchanged |
| Damping QR factorizations | 1,754 | 1,754 | unchanged |
| Candidate scratch peak | 233,512 B | 233,512 B | unchanged |
| Global-state preflight peak | 36,716,466 B | 36,716,466 B | unchanged |
| Sweeps / accepted blocks | 5 / 130 | 5 / 130 | unchanged |
| Objective | `1.47904336e-10` | `1.47904336e-10` | unchanged |
| Global A/C KKT | `5.06726418e-10` | `5.06726418e-10` | unchanged |
| Projected width KKT | `1.87827531e-12` | `1.87827531e-12` | unchanged |
| Endpoint / Independent Replay | pass / pass | pass / pass | unchanged |

For synthetic 300 targets, five independent lower-bound versus direct-map benchmark processes measured Peak RSS medians of 87,851,008 B versus 87,801,856 B. The final 300-target closed-loop benchmark measured Peak RSS 88,227,840 B and the same preflight estimate. No persistent memory was added.

For the real A/YA workload, baseline versus final medians were CLI wall 134.471 s versus 114.499 s (−14.85%, workflow speedup 1.174×), and Joint search 108.909 s versus 89.174 s (−18.12%, stage speedup 1.221×). Baseline/final Peak RSS medians were 1,484,701,696 B and 1,563,672,576 B; run-to-run ranges overlap, with final minimum 1,484,242,944 B and final maximum below the baseline maximum. Endpoint-verification preflight estimate remained 52,661,461 B. All five final runs had the same 118 sweeps, 35,828 profile evaluations, 14,680 accepted updates, objective, independent A/C KKT, projected-width KKT, qualification, and Independent Replay result.

The final Release build completed with `cmake --build build/release -j4`. The serial CTest run passed 21/21 tests with numerical libraries capped at four threads. Coverage includes cancellation-sensitive derivative fallback; dense/reference derivative comparisons; SparseQR rank deficiency, active-set faces and conditional refinement; bounded trust-region damping, bracket and bisection; width-bound activation/release and projected KKT; incremental objective/acceptance and sweep replay; independent replay; component assembly; memory-budget fallback/rejection; non-finite inputs; and large-component memory planning. The 300-target production benchmark also passed endpoint qualification and Independent Replay.

Production audit: the only promoted implementation change is reuse of the existing row-position vector for affected-row direct lookup. FixedNeighbor Block-coordinate, constrained A/C SparseQR, nonnegative active-set, conditional strict-KKT refinement, and bounded active-set LM with width bounds and augmented QR remain the production algorithms. There is no second solver path, runtime experiment selector, or new persistent cache. Detailed attribution and cancellation hooks compile only in the test instrumentation target; their fields are consumed by the hotspot benchmark and regression tests. No third-party dependency or Notion content changed.

The lower-bound comparison switch, temporary per-tile QR reuse experiment, real-run databases/logs, and temporary runner scripts were removed after results were recorded. No unconsumed benchmark hook or experiment-only production branch remains.

Production remains on FixedNeighbor Block-coordinate, constrained A/C SparseQR with nonnegative active-set, and bounded active-set LM with width bounds and augmented QR.

## Follow-up: tiled QR performance study (2026-10-29)

This follow-up starts from `df32fbde2dc9816169c0b6292586b8021543afdf` on `develop`. It re-measures the current implementation; the historical 300-target QR figures above are not used as this study's baseline. The Release build uses the repository's configured AppleClang/Eigen toolchain. Benchmarks use one thread, one warmup, and five measured repetitions; tests/builds are capped at four threads. All counters added for this study are test-instrumentation-only.

### P0 workload census and baseline

`TiledQR::Append` already has an in-place Householder path and a reference-order path for cancellation-sensitive derivative accumulation. The A/C path column-normalizes the selected free-face sparse matrix, scans each 1,024-row tile to map its active columns, materializes that tile as a dense local block, factors it, transforms the matching RHS by Qᵀ, and retains the upper R rows. SparseQR then classifies rank and solves each active-set face. The profile preparation path assembles normalized free-design columns and raw derivative right-hand sides, streams them through an 8,192-row `TiledQR`, assesses free-design rank, and solves projection/correction systems. Derivative reduction builds the corrected profile Jacobian and streams it through a second `TiledQR`. The in-place Householder transform is already used where valid; this study does not repeat the prior Householder-object reuse experiment, which gained only 0.36% full search and was reverted.

The test-only attribution extension records exact structural and numeric row activity separately, dimensions and frequency of each factor shape, dense occupancy, materialization/factor/QᵀRHS/R extraction time, derivative stage time, and a conservative tile-temporary estimate. Tile records aggregate by shape, so instrumentation does not retain one record per observation or QR call. An on/off comparison (both profiling instrumentation and benchmark telemetry enabled together) measured 1.14064 s versus 1.11853 s search median at 300 targets, about 2.0% combined overhead; Peak RSS was 87,867,392 B versus 87,769,088 B. The on/off difference is the combined telemetry/timer/census cost, not a timer-only estimate.

| Synthetic targets | FixedNeighbor wall median (range) | Profile search median | A/C solves | Derivative prepares/reductions | Objective | Qualified / Independent Replay |
|---:|---:|---:|---:|---:|---:|---|
| 12 | 38.471 ms (38.466–38.515) | 25.656 ms | 8 | 4 / 4 | 1.54157728e-10 | pass / pass |
| 48 | 184.044 ms (183.766–184.358) | 127.805 ms | 40 | 20 / 20 | 1.51441738e-10 | pass / pass |
| 240 | 909.982 ms (909.465–910.885) | 631.328 ms | 210 | 105 / 105 | 1.48393658e-10 | pass / pass |
| 300 | 1.14064 s (1.13954–1.14300) | 790.713 ms | 260 | 130 / 130 | 1.47904336e-10 | pass / pass |

The 300-target instrumented run is the detailed census representative. A/C reduction handled 1,380,600 input rows over 1,500 tiles and 260 valid solves. Every design row was exactly nonzero; there were zero zero-design rows with nonzero RHS. The 3.09 million structural entries were all exact nonzeros. Dense occupancy averaged 9.43%, but row activity was 100%. Tile shapes were 1,200×(1,024×24), 240×(400×24), 40×(1,024×20), 10×(514×20), and 10×(970×4); the largest estimated tile temporary was 401,600 B. Thus this workload is sparse within rows but offers no exact-zero rows for row compaction.

| 300-target QR stage | Calls/tiles | Rows and factor columns | Materialization | Householder factor | QᵀRHS | R extraction | Stage total |
|---|---:|---|---:|---:|---:|---:|---:|
| A/C tiled row reduction | 1,500 tiles / 260 solves | 1,380,600 input rows; 4–24 active columns; 1 RHS | 3.87 ms | 206.73 ms | 13.47 ms | 0.89 ms | 234.92 ms |
| Derivative preparation | 130 QR calls/tiles | 690,300 rows; 970×4 (5), 4,610×20 (5), 5,520×24 (120); 12 RHS | 3.13 ms | 87.01 ms | 76.00 ms | <0.1 ms | 201.61 ms QR / 239.82 ms preparation |
| Corrected derivative reduction | 130 QR calls/tiles | 690,300 rows; 970×2 (5), 4,610×10 (5), 5,520×12 (120); 1 RHS | Included in Jacobian build | 24.78 ms | 4.63 ms | <0.1 ms | 94.08 ms reduction, including 35.31 ms corrected-Jacobian build |
| Projection-only width derivative | 26 tiles | 138,060 rows; 970×2 (1), 4,610×10 (1), 5,520×12 (24) | — | — | — | — | exact-nonzero in all 138,060 rows |

For preparation, free design contained 1,545,000 structural/exact entries and raw derivative 772,500; both had 690,300/690,300 exact-nonzero rows, with no row where both were zero. The corrected Jacobian had 8,189,000 structural and exact entries across 690,300 rows and no exact-zero row. Projection-only reduction had 1,637,800 structural/exact entries across 138,060 effective rows. The largest estimated temporary was 3,179,712 B for preparation, 1,148,256 B for corrected reduction, and 1,059,936 B for projection-only reduction. These are bounded tile workspaces; no observation-by-global-parameter dense matrix is retained. The synthetic benchmark's process RSS was about 88 MB.

Across the three A/C and derivative stages, measured median time was about 569 ms, roughly 49.9% of the 1.141 s FixedNeighbor wall. A/C Householder factorization is the single largest QR substage. In derivative preparation, factorization and QᵀRHS are comparable (87 and 76 ms); in derivative reduction, corrected-Jacobian construction (35 ms) is itself larger than the QR factorization (25 ms). The 300-target normal workload had no cancellation fallback; the separate cancellation fixture invoked the preserved reference-order fallback in all five measurements, with 1.42 μs median incremental cost. The independent dense reference fixtures covered 19 cases: narrow/wide matrices, 25/50/75/100% active rows, sparse/dense-like supports, a 2,053-row partial final tile, zero-design rows with nonzero RHS, active nonnegative A and signed C, full rank, rank deficiency, and near-rank deficiency. Full-rank solution relative error was at most about 1.8e-15; both deficient fixtures remained classified rank 23/24 and invalid. No tolerance changed.

The real files `/Users/yslian/data/6z6u.cif` and `/Users/yslian/data/emd_11103_unsharpened.map` existed and were run for chains A and YA with the current 1,024-row A/C tile. Five fresh-database measured CLI runs produced wall median 114.975 s (114.555–115.099), Peak RSS median 1,704,280,064 B (1,633,452,032–1,916,256,256), 118 sweeps, 35,828 profile evaluations, and 14,680 accepted updates. Every run completed endpoint certification, passed all component evidence, and passed Independent Replay. Objective was 0.02195348206329014, independent A/C KKT 4.39461e-9, independent width stationarity 6.14044e-11, and the memory preflight estimate was 133,985,764 B. Fine-grained census for this real CLI is unavailable: the CLI does not expose the test-only `ProfileSearchWork` breakdown. Therefore the QR stage timings and row-shape distribution above are attributed to the instrumented synthetic workflow; the real workflow is used for full end-to-end A/B and qualification, not for fine-grained stage attribution.

P0 answers: Householder factorization dominates A/C reduction; derivative preparation is the largest single complete derivative stage. The measured 300-target workflow has no exact-zero design rows, so A/C zero-row compaction has no production rows to remove. Preparation and corrected-Jacobian row compaction likewise have no exact-zero rows in this workload. Smaller local blocks have some cache/workspace opportunity, but tile candidates must be decided by full workflow and real-data results. Direct SparseQR is worth a controlled comparison because A/C Householder factor time is material and the selected matrix is sparse; it is evaluated in P3 below.

### P1 tile-size tuning

All sizes below use the 300-target workload, one warmup, five measured runs, one thread, and the same Release binary. Each candidate changed one tile-size parameter only. The persistent defaults remained A/C 1,024 and derivative preparation/reduction 8,192.

| A/C tile rows | Full FixedNeighbor median | `ReduceSparseRows` | Householder factor | SparseQR analysis/factor/solve | Estimated largest tile temporary |
|---:|---:|---:|---:|---:|---:|
| 128 | 1.17674 s | 206.97 ms | 163.45 ms | 8.08 / 56.68 / 4.51 ms | 50,368 B |
| 256 | 1.12548 s | 195.33 ms | 158.07 ms | 5.59 / 25.61 / 2.20 ms | 100,544 B |
| 512 | 1.10050 s | 186.68 ms | 156.30 ms | 3.86 / 12.63 / 1.16 ms | 200,896 B |
| 1,024 (baseline) | 1.14486 s | 237.32 ms | 209.22 ms | 2.97 / 6.14 / 0.65 ms | 401,600 B |
| 2,048 | 1.14209 s | 242.42 ms | 215.82 ms | 2.01 / 2.83 / 0.36 ms | 803,200 B |

Tile 512 increased A/C tile count from 1,500 to 2,760. It cut Householder time 25.3% and row-reduction time 21.3%; SparseQR time increased because the tiled reduction retained more rows. Full synthetic search improved 3.88%. At 12, 48, and 240 targets the same 512-row choice improved full FixedNeighbor wall by 3.38%, 3.73%, and 3.50%, respectively. Objectives, rank, KKT, sweeps, qualification, and Independent Replay passed; floating-point hashes changed as expected when row grouping changed. This synthetic-only benefit did not hold on the real workflow: five fresh-database 6Z6U/EMD-11103 A/YA runs at 512 rows took 115.544 s median (115.142–119.833), versus the 1,024-row median of 114.975 s (114.555–115.099), about 0.5% slower. The candidate used 122 sweeps and 14,772 accepted updates versus 118 and 14,680 at baseline. Its objective was `0.021953482063290126` versus `0.02195348206329014`; independent A/C KKT was `2.87957e-9` versus `4.39461e-9`, width stationarity `3.33512e-11` versus `6.14044e-11`; all evidence and Independent Replay passed. The preflight estimate remained 133,985,764 B. RSS ranges varied substantially between runs and do not support attributing their lower candidate median to a sub-megabyte tile-workspace change. Since the real full workflow did not accelerate, tile 512 was rejected and the 1,024-row production default restored.

| Derivative tile rows | Preparation median | Full FixedNeighbor median | Corrected-Jacobian reduction median | Estimated largest prepare / reduce tile temporary |
|---:|---:|---:|---:|---:|
| 256 | 236.91 ms | 1.13722 s | 98.77 ms | 154,560 / 54,592 B |
| 512 | 225.65 ms | 1.12058 s | 94.78 ms | 302,016 / 107,840 B |
| 1,024 | 241.69 ms | 1.14453 s | 94.81 ms | 596,928 / 214,336 B |
| 2,048 | 282.33 ms | 1.20132 s | 99.76 ms | 1,186,752 / 427,328 B |
| 4,096 | 260.55 ms | 1.16932 s | 98.41 ms | 2,359,488 / 852,064 B |
| 8,192 (baseline) | 241.59 ms | 1.14486 s | 93.67 ms | 3,179,712 / 1,148,256 B |

Preparation tile 512 reduced that stage by 6.6% and full synthetic wall by 2.1% in this one target size, below the preferred stage gate and without real-data evidence. Smaller reduction tiles added 5.4–6.5% reduction time; 512/1,024 had no meaningful end-to-end change, while 2,048 and 4,096 were slower. Candidate measured RSS ranged up to 94,355,456 B for the reduction-512 experiment versus 88,047,616 B in its baseline sample; because this exceeded the 5% preference and lacked end-to-end benefit, no derivative tile change was promoted. All tile-size runs passed qualification and Independent Replay; no tolerance was changed. No CLI or persistent benchmark tile-size selector remains.

### P2 zero-row compaction

Skipped without a Production prototype. In the measured production-like 300-target A/C workload, all 1,380,600 rows were exact-nonzero design rows, so compacting exact-zero rows would remove no QR work. Free-design, raw-derivative, corrected-Jacobian, and projection-only derivative rows were also all exact-nonzero in the P0 census. The 19 controlled A/C fixtures demonstrate the separate case where zero-design rows carry nonzero RHS: at 25/50/75/100% active-row ratios, the full RHS is retained and full-rank solutions agree with an independent dense QR reference to at most about 1.8e-15. The identity `||Ax-y||² = ||A_I x-y_I||² + ||y_Z||²` permits omitting those rows from the solve, but their constant residual contribution must remain in the full objective/certificate. Since the real CLI does not expose row census and no measured production opportunity was established, no row-index cache or compaction path was retained.

### P3 direct Eigen SparseQR comparison

A temporary test-only direct path passed the same normalized selected sparse matrix, free columns, and RHS directly to the existing Eigen SparseQR rank/solve code. It was removed after the A/B. The input memory estimate was bounded by the observed 5,520×24 matrix with at most 12,360 stored values (about 0.15 MB for values and indices, plus outer pointers and object overhead); this experiment did not factor any larger component. Peak RSS for the 300-target run was 88,424,448 B, below 1% above the same-build baseline.

| Targets | Tiled A/C wall | Direct SparseQR wall | Change | Qualification / replay |
|---:|---:|---:|---:|---|
| 12 | 38.471 ms | 38.398 ms | −0.19% | pass / pass |
| 48 | 184.044 ms | 184.513 ms | +0.25% | pass / pass |
| 240 | 909.982 ms | 911.207 ms | +0.13% | pass / pass |
| 300 | 1.14064 s | 1.14113 s | +0.04% | pass / pass |

For 300 targets, direct mode reduced row-reduction from 234.92 ms to about 1.01 ms, but SparseQR grew from 2.90/6.10/0.70 ms analysis/factor/solve to 21.97/192.49/19.38 ms. The direct factorization's reported R had at most 214 stored entries; the input had up to 12,360. Eigen does not expose the full internal Q reflector storage as a comparable public fill metric, so factor/solve time and process RSS were used alongside R and input nnz. Thus the tiled Householder stage disappeared but numeric SparseQR became much more expensive; the complete-workflow difference stayed within measurement noise and varied direction across target sizes. The direct path also passed the 19 controlled cases: full-rank dense-reference error stayed below 1.8e-15, active A/signed C remained feasible, and rank-deficient and near-rank-deficient inputs were still rejected at rank 23/24. This is not a Production replacement: there is no stable end-to-end gain, and the maintained compressed-QR path keeps a smaller SparseQR factorization. The temporary direct selector and implementation were fully removed.

### P4 derivative row compaction and P5 architecture study

P4 compaction was skipped because the P0 derivative census found zero rows where free design and raw derivative were both exact zero, zero exact-zero corrected-Jacobian rows, and zero exact-zero projection-only rows. No “small value” threshold was introduced. Residual and cancellation fallback contracts remain intact.

P5 received the required mathematical review; no prototype was justified. The first QR works on free design and raw derivative RHS to assess free-design rank and compute projection coefficients. The nonzero-residual correction then changes the profile Jacobian. The second QR factors that corrected Jacobian into the reduced LM system and transforms its residual RHS. Those factors serve different matrices and the second factorization cannot be removed by reusing only the first R or QᵀRHS result. A one-pass alternative would need to preserve a larger orthogonal representation and still apply the correction while retaining rank and width-derivative semantics; no simpler exact streamed representation was established. The two derivative stages total 333.90 ms (29.3% of the 300-target instrumented FixedNeighbor wall), but the current workflow uses one tile per preparation/reduction, so a persistent workspace or broader QR API had no demonstrated benefit. The existing cancellation detection, reference-order accumulation, and QR fallback remain unchanged.

### P6 final validation and audit

No production optimization passed the real full-workflow gate. Consequently, baseline and final Production metrics are identical by source: A/C rows remain 1,024; derivative preparation/reduction remain 8,192; the current single constrained A/C SparseQR and bounded active-set LM path remain. The P0 five-run baseline set is also the final five-run set because rejected experiments were removed and no Production code changed after it. A post-cleanup real CLI run independently completed in 114.755 s with the same 118 sweeps, 35,828 profile evaluations, 14,680 updates, objective, KKT and preflight; every evidence item passed, including Independent Replay.

| Metric | Baseline | Final Production | Evidence |
|---|---:|---:|---|
| Synthetic 300-target FixedNeighbor wall | 1.11853 s median (1.11737–1.11968) | same | Closed-loop benchmark, five measured runs |
| Synthetic Profile Search | 790.713 ms median (instrumented) | same | Hotspot benchmark, five measured runs |
| A/C row reduction / Householder / QᵀRHS | 234.92 / 206.73 / 13.47 ms | same | Instrumented five-run benchmark |
| A/C SparseQR analysis / factor / solve | 2.90 / 6.10 / 0.70 ms | same | Instrumented representative breakdown |
| Derivative preparation / QR / factor / QᵀRHS | 239.82 / 201.61 / 87.01 / 76.00 ms | same | Instrumented representative breakdown |
| Derivative reduction / factor / QᵀRHS | 94.08 / 24.78 / 4.63 ms | same | Instrumented representative breakdown |
| QR work counts | 1,500 A/C tiles; 260 SparseQR; 130 prepare; 130 reduce; 26 projected-width tiles | same | Instrumented representative breakdown |
| Synthetic process Peak RSS / memory preflight | 88,276,992 B / 36,716,466 B | same | Closed-loop five-run maximum / estimate, preflight passed |
| Sweeps / accepted updates | 5 / 130 | same | closed-loop benchmark |
| Objective / A/C KKT / width KKT | `1.47904336e-10` / `5.06726418e-10` / `1.87827531e-12` | same | qualification and Independent Replay passed |
| Real A/YA CLI wall | 114.975 s median (114.555–115.099) | same five-run baseline; final sanity 114.755 s | five baseline runs plus one post-cleanup run |
| Real Peak RSS / memory preflight | 1,704,280,064 B median / 133,985,764 B | same preflight; final sanity 1,850,720,256 B RSS | process RSS varied; all five baseline and final sanity passed |

`cmake --build build/release -j4` succeeded. `ctest --test-dir build/release -j4 --output-on-failure` passed 21/21 tests. The build emitted three existing missing-field-initializer warnings in `tests/data/DataObjectPersistence_test.cpp`; that file was not changed. The final `rhbm_joint_optimization_benchmark` reported 12.4676 ms bounded-LM median, 200.034 ms profile-search median, and 19.1829/5.8166 ms derivative prepare/reduce medians, all with replay passed. The 300-target fixed-neighbor closed-loop benchmark passed endpoint qualification, Independent Replay, and memory preflight.

Audit: there is no QR optimization in Production from this follow-up, no new runtime solver selector, no direct SparseQR shadow path, no tile-size selector, and no new persistent cache. The maintained test-only census is consumed by `rhbm_joint_hotspot_attribution_benchmark`; the direct SparseQR study code is absent. The existing Production solver remains a single understandable pipeline. The real CLI does not expose fine-grained QR stage timing or row census, which limits real-data attribution; it does expose full workflow result evidence and preflight. No Notion document was read or changed, no third-party dependency was added, and no branch or worktree was created.

| Stage | Status | Decision |
|---|---|---|
| P0 Workload Census | completed | Baseline and exact row/shape/time census committed |
| P1 Tile Size Tuning | completed, no promotion | 512-row A/C helped synthetic but slowed real; retain 1,024/8,192 |
| P2 Sparse Row Compaction | skipped | No exact-zero production-like rows to remove |
| P3 Direct SparseQR | completed, no promotion | Numeric factor/solve cost erased row-reduction savings |
| P4 Derivative Row Compaction | skipped | No exact-zero derivative rows in observed workload |
| P5 Derivative Architecture | research-only | Two QR stages have distinct mathematical roles; no simple exact merge |
| P6 Validation/Cleanup | completed | Release, 21 tests, main benchmarks, and real acceptance passed |

The best next direction for A/C remains reducing dense Householder factorization work without raising SparseQR cost; current exact-zero compaction and tile-size candidates do not meet the real workflow gate. For Profile Derivative, the preparation QR's Householder and QᵀRHS are both substantial, while derivative reduction also spends more time constructing corrected Jacobians than factoring them. Tile-size changes affect cache/workspace footprint and QR count materially but did not produce real-data acceleration. Direct SparseQR is not worthwhile for the measured workload. A new QR kernel is not justified yet; this round should close with the current numerical path and the retained test-only census.

Hotspot attribution instrumentation and its dedicated benchmark were removed after study completion. The numerical regression coverage is maintained by the Joint Component GoogleTest suite. Historical performance data remain preserved in this report.
