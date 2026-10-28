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

Pending gated experiments. Production remains on FixedNeighbor Block-coordinate, constrained A/C SparseQR with nonnegative active-set, and bounded active-set LM with width bounds and augmented QR.
