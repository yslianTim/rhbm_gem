# Joint Incremental Prediction Replay Study

This record compares candidate replay methods with the production Full Replay.
The production implementation remains unchanged and uses Full Replay.

## N0: baseline and workload census

Baseline `develop` HEAD: `ad49ee2c2984d17657628e639542a9ccc18b381b`.
The tree was clean, and a fresh fetch confirmed `origin/develop` at the same
commit. Release configuration used AppleClang 21, C++20, Eigen 5.0.1, system
dependencies, OpenMP enabled, and a four-worker build. CTest ran serially.

The original test configuration passed all 19 CTest groups in 22.73 s. After
adding a test-only candidate observer, the Joint numerical group passed again;
the observer-on and observer-off final eta, beta, prediction, residual, and
objective matched bitwise.

The frozen `chain`, 3-target FixedNeighbor diagnostic is a small synthetic
fixture, not a production cryo-EM workload. It uses 1,425 informative rows,
1,545 support entries, 3 targets, one block, and 1,425 affected-row entries.
Contributor counts per informative row were 1/1/2 (min/median/max), averaging
1.084. Its three-sweep search accepted 3 candidates, rejected 0, and skipped
0 unchanged candidate replays. The search stopped at the configured sweep
budget, so endpoint certification and Independent Replay were not attempted.

| Baseline measure | Result |
| --- | ---: |
| Candidate Full Replay time | 0.1566 ms |
| Local search time | 2.255 ms |
| Sweep-end Full Replay time | 0.1543 ms |
| Global state evaluation time | 0.3047 ms |
| FixedNeighbor search time | 3.111 ms for the recorded run; 3.289 ms median across five runs |
| Candidate replay share | 5.04% |
| Amdahl upper bound if candidate replay were free | 1.053× |
| Estimated peak bytes | 377,037 |
| Benchmark process peak RSS | 49,332,224 bytes for the entire benchmark executable |

Initial eta was `[-0.5978369, -0.5978369, -0.5978369]`. Final eta was
`[-0.693144, -0.693154, -0.693147]`; final beta was
`[2.00002, 0.199998, 1.99993, 0.200015, 1.99999, 0.200005]`. Final objective
was `2.01991e-10`; global A/C KKT was `3.33067e-16` and projected Width KKT was
`6.70651e-10` at the last sweep. These are diagnostic outputs, not a qualified
endpoint.

Streaming trace digests (including each candidate's accepted-before state,
parameters, prediction, residual, objective, decision, and replay flag):

- Final state: `9f490134fa1515af4eb77a39884ce778078f90a61b094c3f4d3c523bf17a12cc`
- Candidate sequence: `7edb83c9df15a984c087e1f13998c04d3ea76068e8ece82dfa3aed29bc26752e`

This fixture makes Full Replay a small fraction of search time. It only
establishes the instrumentation baseline; larger replay-only fixtures are
needed to evaluate the affected-row methods.

## N1: deterministic replay comparison

The test-only microbenchmark uses the same accepted-before parameter and
prediction state for every candidate method. It applies one fixed candidate
per structural block and does not let any candidate result change a later
candidate. Timings include candidate evaluation; D1/E also include rebuilding
the block-local contributor incidence for every candidate. Fixture/sequence
construction and full-reference verification are reported separately, rather
than folded into one method's timing. The candidate sequence is synthetic and
is not a solver trajectory.

| Targets | Rows | Blocks | A Full Replay | B1 | B2 | C | D1 | E | D1 incidence prep |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 1,425 | 1 | 0.0520 ms | 0.1488 ms | 0.1482 ms | 0.1306 ms | 0.1619 ms | 0.1438 ms | 0.0784 ms |
| 30 | 13,710 | 3 | 0.5397 ms | 0.6913 ms | 0.6941 ms | 0.4694 ms | 0.7581 ms | 0.5330 ms | 0.3127 ms |
| 300 | 136,560 | 26 | 5.2585 ms | 2.7318 ms | 2.7099 ms | 0.5531 ms | 2.7880 ms | 0.6242 ms | 0.3680 ms |

Times are five-sample medians per candidate. D1 replay-only evaluation at 300
targets was 2.3637 ms per candidate; together with incidence preparation it
was 1.89× faster than Full Replay. D1 was slower at 3 and 30 targets. The
largest fixture therefore approaches but does not reach the suggested 2×
candidate replay threshold. C was 9.51× faster and E 8.42× faster at 300
targets, but their incremental objectives differed from Full Scan for all 26
candidates (maximum absolute difference `4.33681e-17`, maximum 50 ULP).

At 300 targets, B1 and B2 both had Full Scan objectives bitwise equal to A,
but their prediction and residual arrays were not bitwise equal. B1 differed
on 3,928 of 3,550,560 row comparisons, with maximum prediction error
`2.77556e-17`; B2 differed on 4,492, with maximum error `5.55112e-17`.
Residual maximum error was `5.55978e-14` relative and 256 ULP on near-zero
residuals. B2 produced more prediction mismatches and twice B1's maximum
absolute prediction error, confirming the two update expressions are not
interchangeable.

D1 reproduced prediction, residual, and Full Scan objective bitwise on all
3,550,560 informative-row comparisons and 26 candidate objectives at 300
targets. The test also verifies that each affected row's incidence contains
every contributor in `layout.full_atoms` order, including original support
traversal order. Estimated temporary bytes were 321,824 for D1, 88,512 for
B1/B2/C, and 2,184,960 for A. The 300-target process peak RSS was 72,056,832
bytes, including the synthetic input, reference verification, and every
strategy. Fixture plus candidate-sequence construction took 0.0804 s, and
full-reference verification took 0.7491 s.

N1 establishes D1 as the only candidate with bitwise-identical prediction,
residual, and Full Scan objective. It is a candidate for shadow and closed-loop
validation, not yet a production optimization: total FixedNeighbor search
speedup, acceptance agreement, long-run behavior, and larger-component
performance remain unmeasured.

## N2: Full-Replay-controlled shadow search

The test-only shadow ran on the deterministic 30-target chain with five
sweeps. Production Full Replay made all 15 candidate decisions; all 15 were
accepted. The shadow candidates never changed the production state or later
local searches. The final state, sweep objectives, sweep replay errors, and stop
reason matched the observer-off baseline bitwise. Existing tests separately
cover unchanged candidates and block recovery/retry behavior.

| Shadow method | Prediction mismatches | Residual mismatches | Objective mismatches | Acceptance disagreements | Maximum absolute prediction / residual / objective error |
| --- | ---: | ---: | ---: | ---: | ---: |
| B1 | 9,163 | 9,163 | 10 | 0 | `1.11022e-16` / `1.11022e-16` / `1.85288e-22` |
| B2 | 9,238 | 9,238 | 10 | 0 | `1.11022e-16` / `1.11022e-16` / `7.94093e-23` |
| C | 9,163 | 9,163 | 15 | 0 | `1.11022e-16` / `1.11022e-16` / `3.33067e-15` |
| D1 | 0 | 0 | 0 | 0 | 0 / 0 / 0 |
| E | 0 | 0 | 15 | 0 | 0 / 0 / `3.33067e-15` |

B1 and B2 maximum prediction differences were 3 and 4 ULP, respectively. Their
maximum residual difference was 536,870,912 ULP because some residuals are
near zero; maximum relative residual error was `9.35544e-8`. Their first
prediction and residual mismatch was informative row 1,507; the first Full Scan
objective mismatch occurred at candidate 5. C and E's largest objective
difference was 71,661,654,325 ULP near a zero objective, with maximum relative
objective error `1.22105e-5`; C first differed at candidate 0, E at candidate 0.
These are 15 accepted updates without a shadow refresh. D1's minimum objective acceptance margin was
`1.0955168297673593e-16`, or `0.00010955168294349627` of its replay enclosure.
There was no acceptance disagreement in this fixture, but it had no rejected
candidate and its narrowest margin means this small run alone cannot establish
that B/C/E are safe around the decision boundary. N3 must check whether refresh
can contain their accumulated drift before deciding whether they proceed; D1
proceeds to closed-loop testing.

## N3: long-run drift, refresh, and rollback

The repeatable synthetic stress used 30 targets, 13,710 rows, 15,450 support
entries, and the first overlapping 12-atom block (5,520 affected rows). It made
503 attempts on that same block, with 474 scripted accepted updates and every
17th attempt rejected (29 rollbacks). The sequence includes alternating signed
charge coefficients, repeated small updates, width swings, updates within
`1e-12` of the upper width bound, and near-converged coefficient perturbations.
The 11-update refresh interval represents one synthetic sweep. This is a
deterministic arithmetic stress, not a FixedNeighbor LM trajectory: the update
and reject schedule is scripted, and the reported acceptance check is the
objective-only nonincrease/replay-enclosure proxy. All 29 rejects preserved the
accepted shadow state bitwise for each strategy. Every full refresh restored
prediction, residual, and objective exactly (zero post-refresh mismatches).

The timings below cover candidate calculation, accepted-state updates, and
scheduled Full Replay refresh; independent reference replay and error scanning
are outside the timed interval. This was one run per case, so the times are
diagnostic rather than a stable performance claim. An 11-update sweep refresh
is shown last.

| Strategy | Refresh interval | ms/attempt | Max prediction error | Max residual error | Max objective error | Max prediction/residual ULP | Max objective ULP |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| A | Full Replay each candidate | 0.533 | 0 | 0 | 0 | 0 / 0 | 0 |
| B1 | 1 | 1.393 | `1.39e-16` | `2.22e-16` | `5.55e-17` | 32 / 4 | 1 |
| B1 | 10 | 0.858 | `1.67e-16` | `2.22e-16` | `5.55e-17` | 32 / 10 | 1 |
| B1 | 50 | 0.814 | `5.00e-16` | `5.00e-16` | `1.11e-16` | 32 / 32 | 2 |
| B1 | 100 | 0.808 | `8.33e-16` | `8.33e-16` | `1.11e-16` | 54 / 60 | 2 |
| B1 | 11 (sweep) | 0.850 | `1.39e-16` | `2.22e-16` | `5.55e-17` | 32 / 10 | 1 |
| B2 | 1 | 1.321 | `5.55e-17` | `5.55e-17` | `5.55e-17` | 2 / 2 | 1 |
| B2 | 10 | 0.842 | `1.11e-16` | `1.11e-16` | `5.55e-17` | 4 / 6 | 1 |
| B2 | 50 | 0.819 | `2.22e-16` | `2.22e-16` | `5.55e-17` | 9 / 16 | 1 |
| B2 | 100 | 0.794 | `2.78e-16` | `2.78e-16` | `5.55e-17` | 25 / 20 | 1 |
| B2 | 11 (sweep) | 0.838 | `1.11e-16` | `1.11e-16` | `5.55e-17` | 4 / 6 | 1 |
| C | 1 | 1.083 | `1.39e-16` | `2.22e-16` | `3.11e-15` | 32 / 4 | 56 |
| C | 10 | 0.615 | `1.67e-16` | `2.22e-16` | `3.22e-15` | 32 / 10 | 58 |
| C | 50 | 0.578 | `5.00e-16` | `5.00e-16` | `2.61e-15` | 32 / 32 | 47 |
| C | 100 | 0.594 | `8.33e-16` | `8.33e-16` | `2.78e-15` | 54 / 60 | 50 |
| C | 11 (sweep) | 0.616 | `1.39e-16` | `2.22e-16` | `2.89e-15` | 32 / 10 | 52 |
| D1 | 1 / 10 / 50 / 100 / 11 | 1.411 / 0.936 / 0.880 / 0.883 / 0.922 | 0 | 0 | 0 | 0 / 0 | 0 |
| E | 1 | 1.187 | 0 | 0 | `3.11e-15` | 0 / 0 | 56 |
| E | 10 | 0.716 | 0 | 0 | `3.22e-15` | 0 / 0 | 58 |
| E | 50 | 0.665 | 0 | 0 | `2.61e-15` | 0 / 0 | 47 |
| E | 100 | 0.670 | 0 | 0 | `2.83e-15` | 0 / 0 | 51 |
| E | 11 (sweep) | 0.701 | 0 | 0 | `2.89e-15` | 0 / 0 | 52 |

B1 and B2 prediction/residual drift increased as refresh intervals grew. For
example, B1's maximum prediction error rose from `1.39e-16` with refresh every
update to `8.33e-16` at 100 updates; B2 rose from `5.55e-17` to `2.78e-16`.
The incremental objective in C/E differed on 489–494 of 503 attempts even with
refresh every update: refresh clears carried state error, but the next
incremental objective still has its own rounding difference. C's maximum was
`3.22e-15` (58 ULP); E had bitwise-exact predictions/residuals but the same
objective issue. No objective-only acceptance disagreements occurred in this
scripted sequence; that does not establish trajectory safety near real LM
acceptance boundaries. The measured average was 0.533 ms/attempt for A; even
the fastest refresh setting per strategy remained slower (B1 0.808, B2 0.794,
C 0.578, D1 0.880, E 0.665 ms/attempt). Refresh did not make these methods a
net win at this component size.

Estimated peak bytes were 438,720 for A/B1/B2/C and 540,208 for D1/E. Candidate
scratch alone was 88,512 bytes for B1/B2/C, 320,848 for D1/E, and 219,360 for
A. Process peak RSS for the entire stress executable was 26,656,768 bytes.
Existing tests cover unchanged candidates, recovery retry/exhaustion, and
resource-limited endpoint certification. The stress runner does not inject a
thrown exception or exercise a real solver trajectory.

## N4: closed-loop D1 search

A test-only candidate replay callback let D1 control candidate objective,
consistency, acceptance, accepted prediction/residual, and later block response.
The production build has no callback field or strategy branch. On the
three-target chain with two-atom cores and practical convergence criteria, the
Full Replay and D1 searches both completed five sweeps, ten block solves, ten
accepted candidates, and 20 local profile trials. Both stopped with
`block-stationary`, converged, and passed endpoint certification and
Independent Replay.

Level A passed. Every candidate's parameter snapshot, accepted-before state
digest, candidate prediction/residual digest, objective values, local profile
trial counts, acceptance flags, local objectives, width/gradient updates, LM
ratios, and block record matched bitwise. Sweep objectives, cache-replay
errors, global A/C KKT, projected Width KKT, sweep states, and final
eta/beta/prediction/residual/objective also matched bitwise.
Final objective was `2.0199079729737389e-10`; global A/C KKT was
`4.3341136235497402e-9`, and projected Width KKT was `1.6235413013987454e-10`.
Independent Replay passed for both. This is a single-component, two-block
fixture; the full regression suite remains part of final acceptance. The test
is for correctness only and makes no speed claim.

## N5: larger synthetic performance and memory study

No cryo-EM `.mrc`/`.star` or large observation data files are present in the
project workspace, so N5 uses the reproducible synthetic `chain` workload. The
largest closed-loop case has 1,000 targets and 455,060 observations; the
5,938,087-observation workload was not available. Each closed-loop size had one
warm-up and five measured repeats. Timings are medians; ranges are the five
measured runs. The diagnostic runs five search sweeps and disables endpoint
certification, so these are FixedNeighbor search measurements, not complete
algorithm runtimes.

| Targets / rows / blocks | Block prep A → D1 | A search median (range) | D1 search median (range) | A→D1 search change | A→D1 candidate replay | D1 added scratch estimate | Process peak RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 / 1,425 / 1 | 0.125 → 0.117 ms | 4.953 ms (4.744–5.094) | 5.606 ms (5.320–9.124) | 13.2% slower | 0.267 → 0.852 ms | 81,768 B | 23,904,256 B |
| 30 / 13,710 / 3 | 1.220 → 1.242 ms | 100.144 ms (99.400–101.688) | 104.198 ms (102.393–154.713) | 4.0% slower | 7.948 → 11.824 ms | 321,824 B | 40,419,328 B |
| 300 / 136,560 / 26 | 11.749 → 11.939 ms | 1.690 s (1.636–1.704) | 1.364 s (1.333–1.427) | 19.3% faster (1.24×) | 694.612 → 373.164 ms | 321,824 B | 97,730,560 B |
| 1,000 / 455,060 / 84 | 39.721 → 39.540 ms | 10.863 s (10.502–10.961) | 6.704 s (6.606–7.080) | 38.3% faster (1.62×) | 7.514 → 3.393 s | 321,824 B | 244,400,128 B |

The 1,000-target run accepted 420 candidates in both methods. All five
Full-Replay and D1 runs had bitwise-identical block/sweep traces and final
state. A's candidate replay occupied 69.1% of search time; replacing it with
D1 left local search and sweep/global work essentially unchanged. The estimated
baseline peak was 121,371,048 B; D1 raised it to 121,692,872 B. The process
RSS is for the benchmark process, including both methods and repeated runs.
At 3 and 30 targets D1 was slower; at 3 targets the absolute median increase
was 0.65 ms. There were no rejected candidates in these closed-loop searches,
so N5 did not produce a standalone rollback timing. N3 separately verified
bitwise rollback for 29 scripted rejected updates. The search diagnostic does
measure sweep-end replay; at 1,000 targets its median was 88.9 ms for A and
88.2 ms for D1. It does not run endpoint certification or Independent Replay;
N4's converged fixture passed those checks.

The 1,000-target replay-only study measured all replay candidates five times
with a fixed accepted state and 84 fixed block candidates. Full objective scan
is included except for C/E, whose incremental objectives are experimental.
D1's incidence preparation is included in its total candidate time.

| Method | Preparation | Candidate time | Temporary / index memory | Prediction, residual, objective |
| --- | ---: | ---: | ---: | --- |
| A | 0 | 17.525 ms | 7,280,960 B | Reference |
| B1 | 0 | 7.992 ms | 88,512 B | 12,960 prediction and residual mismatches; full-scan objective exact |
| B2 | 0 | 7.962 ms | 88,512 B | 14,776 prediction and residual mismatches; full-scan objective exact |
| C | 0 | 0.587 ms | 88,512 B | Same prediction/residual differences as B1; incremental objective differs on all 84 candidates |
| D1 | 0.391 ms/candidate incidence build | 8.009 ms total (7.553 ms replay) | 321,824 B | Bitwise exact across 38,225,040 row comparisons and 84 objectives |
| D2 | 4.264 ms global index build once | 7.655 ms (7.706 ms/candidate amortized over 84) | 96,320 B scratch + 11,880,488 B retained CSR | Bitwise exact across 38,225,040 row comparisons and 84 objectives |
| E | 0.386 ms/candidate incidence build | 0.658 ms total (0.213 ms update) | 321,824 B | Prediction/residual exact; incremental objective differs on all 84 candidates |

B1/B2 prediction maximum absolute errors were `2.78e-17`/`5.55e-17`; residual
errors reached 256 ULP near zero. C/E's incremental objective maximum absolute
error was `3.64e-17` (42 ULP). They are not production candidates. The replay
benchmark process peak RSS was 197,754,880 B across verification and all
methods; it is not a per-method RSS measurement. D2's index-build peak estimate
was 19,616,508 B, including its temporary row counts/cursors. Across 84
candidates, D2 repaid its 4.264 ms build after 13 candidates (0.144 of a
84-candidate sweep) relative to D1, but saved only 0.303 ms per candidate after
amortization and retained 11.88 MB. That incremental speed is not worth the
global index memory and lifetime; D2 is rejected.

D1 meets the approximate 2× candidate-replay and 15% FixedNeighbor search
thresholds at 1,000 targets, and meets the total-search threshold at 300
targets. At 3 and 30 targets its regression is material in relative terms,
though the 3-target cost is under 1 ms. This supports the bounded, on-demand
D1 implementation with a Full-Replay memory fallback; it does not establish
performance on the unavailable 5.9-million-observation data. N6 adds and
preflights that production path, then reruns the suite and closed-loop checks.

## Roadmap status

| Stage | Status | Result |
| --- | --- | --- |
| N0 | Complete | Baseline telemetry and deterministic state/candidate digests recorded. |
| N1 | Complete | A/B1/B2/C/D1/E measured at 3, 30, and 300 synthetic targets; D1 passed bitwise candidate replay checks. |
| N2 | Complete | Five-method shadow; Full Replay retained every decision; D1 remained bitwise exact, while B/C/E accumulated drift. |
| N3 | Complete | 503-attempt refresh/drift stress at 1/10/50/100/sweep; B methods drifted more with longer intervals, C/E retained objective error, D1 stayed exact; all scripted rejections rolled back bitwise. |
| N4 | Complete | D1 controlled candidate results in closed loop and reached Level A on the converged two-block fixture; Independent Replay passed. |
| N5 | Complete | Five-run closed-loop comparison through 1,000 synthetic targets; D1 reached 1.62× search speedup with bitwise-identical traces; D2 was exact but rejected for a small gain and 11.88 MB retained index. |
| N6 | Pending | |
