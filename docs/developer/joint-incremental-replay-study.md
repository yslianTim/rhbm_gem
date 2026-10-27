# Joint Incremental Prediction Replay Study

This record preserves the historical candidate replay study and records the
later E1 Production promotion. At the N0 baseline, Production used Full
Replay; the historical N6 checkpoint adopted D1 with a Full Replay memory
fallback. The current Production decision is recorded in M0–M4 below.

> Historical note: the N0–P4 sections below record earlier study checkpoints.
> Their statement that Production remained D1 describes that earlier state and
> is superseded by the M0–M4 promotion record at the end of this document.

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
though the 3-target cost is under 1 ms. N6 implements and preflights the
bounded, on-demand D1 production path, then reruns full search and correctness
checks. Neither stage establishes performance on the unavailable
5.9-million-observation data.

## N6: production decision and final acceptance

Production now builds a per-candidate compressed incidence for the current
block's affected rows from core plus neighboring atoms. Contributors are
ordered by `layout.full_atoms`, then by their existing per-atom support order.
The replay uses the same kernel, multiply/add expression, and informative-row
objective sum order as `ReplayInto()`. It retains only affected prediction and
residual values; accepted candidates update those rows, while rejected
candidates leave accepted state untouched. Full Replay remains available when
the D1 scratch estimate exceeds the memory budget, and continues to verify
each sweep and endpoint.

The resource preflight now counts the maximum contributor entries and local
atoms for any block and includes incidence, offsets, widths, and affected-row
vectors in the Global State Evaluation stage. The production path uses the
existing 3.2 GB planning budget. A test forces the D1 plan over budget while
leaving Full Replay within budget; search falls back before creating the D1
incidence and returns a bitwise-identical result. On the 1,000-target fixture,
the observed D1 scratch was 233,512 B and the preflight increase over the
previous Full Replay estimate was 277,848 B.

The production-path closed-loop run repeated each size five times after a
warm-up. The A side uses the test-only Full Replay callback; D1 uses the
production path. Every run compares all block and sweep decisions and final
state bitwise. Endpoint certification is disabled for timing, so this remains
a FixedNeighbor search benchmark rather than total algorithm runtime.

| Targets | A search median (range) | D1 search median (range) | Search change | A → D1 candidate replay | D1 scratch / preflight peak | Process peak RSS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 4.740 ms (4.736–4.838) | 5.140 ms (5.134–5.478) | 8.4% slower | 0.260 → 0.660 ms | 58,976 / 447,445 B | 23,461,888 B |
| 30 | 95.447 ms (95.061–97.998) | 98.697 ms (98.372–99.877) | 3.4% slower | 7.803 → 10.492 ms | 233,512 / 3,930,461 B | 35,143,680 B |
| 300 | 1.632 s (1.616–1.638) | 1.314 s (1.300–1.328) | 19.5% faster (1.24×) | 679.116 → 357.100 ms | 233,512 / 36,716,466 B | 94,076,928 B |
| 1,000 | 10.669 s (10.525–14.231) | 6.506 s (6.495–6.561) | 39.0% faster (1.64×) | 7.420 → 3.309 s | 233,512 / 121,648,896 B | 248,119,296 B |

The largest run accepted 420 of 420 candidates for both methods and retained
bitwise-identical search traces. The single A timing outlier at 1,000 targets
is reflected in its range; the median comparison remains 1.64×. At small
sizes, the candidate replay itself is slower, but the total search penalty is
0.40 ms at 3 targets and 3.25 ms at 30 targets. At 300 and 1,000 targets the
search wins exceed the 15% engineering threshold. The 1,000-target candidate
replay exceeds 2× speedup. This is enough to adopt D1 without a size-based
heuristic: the memory preflight controls whether it runs, and Full Replay is
used when D1's scratch would exceed budget.

The converged two-block comparison in N4 ran production D1 against a Full
Replay reference for ten accepted candidates and twenty local profile trials.
It remained Level A: parameters, predictions, residuals, objectives, candidate
and sweep traces, KKT values, stop reason, endpoint certification, and
Independent Replay qualification matched. Final objective was
`2.0199079729737389e-10`, global A/C KKT `4.3341136235497402e-9`, and projected
Width KKT `1.6235413013987454e-10`. The full 19-group CTest suite passed
serially in 21.96 s after a four-worker build. The new memory-fallback case and
all existing component, endpoint, persistence, workflow, and CLI tests passed.

The final choice is D1 with memory-bounded Full Replay fallback. D2, B1/B2, C,
and E remain benchmark-only; D2's index was exact but its 3.8% amortized replay
gain did not justify 11.88 MB at 1,000 targets, while B/C/E failed bitwise
prediction or objective equality. The real 5.9-million-observation case and
complete end-to-end runtime were unavailable; these are the remaining limits
on the performance claim. No FirstStage, Local MDPDE, Group Fitting, KKT,
solver acceptance, bound, convergence, or Independent Replay rules changed.

## Incremental Objective study

### R0: numerical contract and current Production baseline

The study started on `develop` at `6f41026cde46153d537a8f41d7ab749fcd14f4c0`.
The working tree was clean, `develop` and the locally recorded
`origin/develop` ref were equal, and there were no local commits ahead or
behind. No fetch was performed in this run, so this records the local tracking
ref rather than a fresh server-side query. Production contained D1 and its
preflighted Full Replay fallback; incremental objective was available only in
test experiments (`ReplayObjectiveStrategy::Incremental`). There was no
Production objective cache, error bound, conditional scan, or global CSR path.

The numerical levels for this study are:

- **Level A — Bitwise identity:** prediction, residual, objective, and the
  compared search decisions are identical to Full Replay.
- **Level B — Decision-safe numerical equivalence:** a conservative bound
  proves the existing decision for every value in the bound; otherwise Full
  Scan is used. Sweep verification, endpoint certification, and Independent
  Replay still pass.
- **Level C — Empirical agreement:** measured results are close without a
  decision-safety proof. This is research evidence only and is not eligible
  for Production.

#### Acceptance and replay contract

For each visited block, `local_objective_before` is half the sum of squared
accepted residuals on its affected rows, divided by `scale²`. The local solver
returns only a valid, available, feasible, KKT-qualified, trusted endpoint.
`local_objective_after` is its qualified local certificate objective divided
by `scale²`. If the local objective increases, the candidate is left unchanged;
on a recovery visit the same condition exhausts that block recovery. This
local comparison has no additive roundoff allowance in the outer block code.

After local qualification, D1 computes the candidate predictions for all
affected rows and scans every informative row in the production order to form
the candidate objective. It then computes
`global_delta = candidate_objective - accepted_objective` and
`consistency_error = abs(global_delta + local_before - local_after)`. The
consistency limit is `1e-12 + 2e-12 * max(abs(candidate_objective),
abs(accepted_objective))`; the tolerance therefore depends on the absolute
objective values. A consistency failure stops the search with
`block-objective-replay-failed`. If `global_delta > 0`, the block is accepted
only when the local objective did not increase and the global increase is at
most that same enclosure; otherwise it remains unchanged, or terminates a
recovery retry. A candidate whose core parameters are exactly unchanged skips
replay. Accepted prediction, residual, objective, and parameter state are
updated only after these checks, so the candidate replay itself cannot mutate
the accepted cache before acceptance.

At each completed sweep, a Full Replay reconstructs predictions, residuals,
and objective. Prediction cache error must be no greater than
`2e-12 + 2e-13 * max(1, max(abs(full_prediction)))`; objective error must be no
greater than the objective enclosure above, referenced to the replay
objective. A failure stops with `block-cache-replay-failed`. This check
compares the cache but does not overwrite it. Global A/C and projected Width
KKT are evaluated from a separately streamed state evaluation. Stationarity
then invokes endpoint certification; a successful production endpoint still
requires the existing primary and independent replay trust checks. The strict
and practical KKT, LM, width-bound, convergence, and Independent Replay
criteria are unchanged by this study.

The existing `BlockObjectiveReplayEnclosure` is an engineering comparison
tolerance, not a derived floating-point error bound. It cannot certify an
incremental objective by itself. The Full Scan fallback is required whenever
a candidate decision cannot be proved inside an eventual conservative error
interval. Nonfinite or overflowed intermediate values cannot be treated as a
fast acceptance.

#### Baseline measurements

The Release build in `build/replay-study` used system dependencies and
`BUILD_TESTING=ON`; the build had no work to do. Serial CTest passed all
19/19 groups in 22.24 s. No thresholds or Production code changed.

The 1,000-target replay benchmark used 455,060 informative rows, 515,000
support entries, 84 candidates, and 38,225,040 prediction/residual comparisons
per method. It reported 8.898 s for independent Full Replay verification.
D1 matched all 84 objectives and all prediction/residual values bitwise. The
incremental E objective differed on all 84 candidates: maximum absolute error
`3.64292e-17`, relative error `7.18793e-15`, and 42 ULP. Its candidate median
was 0.648 ms including incidence preparation (0.209 ms for replay/update);
D1 was 7.878 ms including incidence preparation (7.439 ms replay). The
process RSS was 197,787,648 B for the benchmark containing every method and
reference verification, not per-method RSS.

The 1,000-target closed-loop timing used five measured runs after warm-up,
five sweeps per run, and no endpoint certification. Full Replay took a median
10.393 s of search and D1 6.401 s (1.624x); candidate replay was 7.248 s and
3.283 s, respectively. Both accepted 420/420 candidates and the benchmark
verified identical traces. D1 scratch was 233,512 B and the preflight peak
estimate was 121,648,896 B. Process peak RSS was 242,663,424 B across both
methods and all repetitions. This is a FixedNeighbor synthetic benchmark,
not an end-to-end algorithm result or a 5.9-million-row measurement.

For a converged two-block fixture, the existing N4 study recorded ten accepted
candidates, twenty local profile trials, bitwise-identical Full Replay and D1
traces, endpoint certification and Independent Replay passing, objective
`2.0199079729737389e-10`, global A/C KKT `4.33411e-9`, projected Width KKT
`1.62354e-10`, and a narrowest observed D1 objective acceptance margin of
`1.09552e-16` (about `1.09552e-4` of the replay enclosure). This is one small
fixture and supplies no general acceptance-bound proof.

### R1: incremental objective strategy comparison

R1 extends only test support and benchmark code. The 30-target `chain`
workload used 13,710 observations, 15,450 support entries, a 12-atom core,
5,520 affected rows, and 503 scripted attempts (474 accepted and 29 rejected).
Prediction used exact affected-row replay. Every candidate was checked against
Full Replay outside the timed interval. The stress includes signed charge
coefficients, alternating small updates, width swings and upper-bound
proximity. It is deterministic arithmetic stress, not an LM trajectory.

The objective alternatives were: current direct sum of squared differences;
scale each difference before summing; Kahan or Neumaier accumulation of the
unscaled or scaled differences; recompute changed fixed-size tiles of 256,
512, 1,024, or 2,048 rows and sum tile values; and update a fixed binary tree
above the same tiles. The tiled caches only commit after a scripted accepted
candidate. Rejected candidates discard their sparse cache updates. A focused
test confirms that discarded changes leave subsequent cache evaluations
unchanged and that accepted updates continue to match a fresh Full Scan within
the observed reduction-order difference.

At an 11-attempt refresh interval, direct incremental objective took
`0.675 ms/attempt`; scale-per-term took `0.675 ms`; Kahan ranged
`0.691 ms`; Neumaier ranged `0.677 ms`. They differed from Full Scan on
489–490 of 503 attempts, with maximum absolute error `2.89e-15` for direct and
`2.94e-15` for the compensated variants (52 and 53 ULP respectively). Refresh
every accepted attempt also left 494/503 objective mismatches, with maximum
`3.11e-15` and 56 ULP. More frequent refresh removes carried cache error but
does not remove the next candidate's incremental arithmetic difference.
Compensation did not reliably reduce the objective error and added a small
measured cost.

Across tile sizes, cache evaluation took about `0.075 s` for 503 attempts
(`0.150 ms/attempt`), compared with `0.675 ms/attempt` for the entire direct
incremental candidate path. Total tile-method time was `0.820–0.822 ms/attempt`
at the 11-attempt interval; the hierarchical tree was `0.820–0.825 ms/attempt`.
The direct tile cache used 56–432 B at 2,048–256 rows; tree cache used 128–1,024
B. Sparse update scratch was 112–864 B for tiles and 224–3,456 B for trees.
The strategies differed on 481–491 of 503 objectives, with maximum absolute
error `2.39e-15–2.55e-15` (43–46 ULP). The tree did not pay for itself on this
fixture: tile recomputation plus propagation was slower than direct
incremental accumulation, and regrouping the reductions still differed from
the production Full Scan.

No scripted objective-only acceptance disagreement occurred for these
strategies, and all 29 scripted rejects preserved the accepted state. This is
empirical evidence only; it does not prove decision safety near an LM
acceptance boundary. The observed signed-delta cancellation ratio was 2.81 on
this fixture, but it does not bound other workloads. Underflow/overflow and
nonfinite behavior were not covered by this stress. For normal finite
arithmetic, a first-order analysis can bound each residual square and a
positive sequential sum with a `gamma_n` term, but the incremental difference
also subtracts nearly equal squares and carries prior objective rounding.
Kahan and Neumaier only address part of the signed-term summation. The current
scalar cache does not carry the residual-scale data or accumulated error state
needed to certify those terms for every candidate. No empirical tolerance is
promoted to a certified bound.

Direct incremental objective is the fastest simple strategy and is carried
into the incidence comparison as an experimental candidate. None of the R1
strategies is eligible for Production before an acceptance error interval
that covers both the candidate and accepted objectives is implemented and
validated. R1 changed no Production path or numeric threshold. The full
serial CTest suite passed 19/19 in 22.02 s after this test-only change.

### R2: E1 on-demand incidence versus E2 global CSR

R2 reuses the direct incremental objective from R1 and keeps exact
affected-row prediction. E1 builds the existing local incidence for every
candidate. E2 builds one global row-to-contributor CSR, then filters each
affected row to the block's core and neighbor atoms. The E2 CSR stores a
32-bit position in `layout.full_atoms` and a 32-bit position in that atom's
existing support list. It reads the squared distance from `domain.atoms` and
does not copy squared-distance values. Construction checks the row-count,
contributor-count, offset, support-index, atom-index, and byte arithmetic for
overflow. Contributor fill order follows `layout.full_atoms` and original
support order. E2 remains test-only; this study did not add a Production
selector or change the `StructuralPartition` path.

The benchmark used five timed repeats over 84 candidates at 1,000 targets and
26 candidates at 300 targets. E2 prediction and residual matched Full Replay
bitwise on all 38,225,040 row comparisons at 1,000 targets and 3,550,560 at
300 targets. E1 and E2 share the same incremental objective calculation, so
both differed from Full Scan on every candidate: maximum absolute error was
`3.64292e-17` (42 ULP) at 1,000 targets and `4.33681e-17` (50 ULP) at 300.

| Targets | Actual support entries | E1 incidence prep / candidate | E1 replay / candidate | E1 total / candidate | E2 CSR build once | E2 replay / candidate | E2 amortized / candidate | Break-even candidates | E2 CSR resident | E2 build scratch | E2 candidate scratch |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 300 | 154,500 | 0.365 ms | 0.201 ms | 0.623 ms | 0.724 ms | 0.214 ms | 0.242 ms | 2 | 2,328,488 B | 2,321,520 B | 91,020 B |
| 1,000 | 515,000 | 0.380 ms | 0.209 ms | 0.647 ms | 2.916 ms | 0.225 ms | 0.260 ms | 7 | 7,760,488 B | 7,736,020 B | 97,320 B |

The amortized E2 time includes its CSR build divided across the measured block
candidates. It was 2.58x faster than E1 at 300 targets and 2.49x at 1,000.
E2's candidate replay itself was about 6–8% slower than E1's replay because it
filters global row entries, but avoids E1's 0.365–0.380 ms per-candidate local
incidence build. The CSR break-even was below one 26- or 84-candidate sweep.
The incremental objective cache adds no per-row resident data; the accepted
scalar objective is already part of the state.

At 1,000 targets the E2 index had 3,640,488 B of row offsets and 4,120,000 B
of contributor entries. At 300 targets those were 1,092,488 B and 1,236,000 B.
The reported E2 build scratch consists of row counts, a fill cursor, and the
informative-row bitmap. E2's build peak estimate, including the accepted
prediction/residual state, was 6,834,968 B at 300 targets and 22,777,468 B at
1,000 targets. The full benchmark process peak RSS was 78,086,144 B and
205,602,816 B, respectively; those process peaks include the synthetic input,
both CSR variants, reference checks, and all methods, not a per-method peak.
E1's candidate scratch was 321,824 B; E2's was 91,020–97,320 B. Accepted state
is not changed during candidate evaluation, so rollback retains no extra full
copy in either path.

The 5,938,087-observation real workload and its measured support-entry count
were unavailable. No fixed average contributor count was used to extrapolate
its CSR memory. E2 advances to R3 as an experimental candidate because its
replay-time savings were repeatable at both measured synthetic sizes and its
scratch remains bounded by actual index sizes. Promotion still depends on the
R3 decision-safety gate; these replay-only timings do not establish a
FixedNeighbor or end-to-end speedup.

### R3: acceptance safety gate outcome

R3 did not qualify an incremental objective for Production. No conservative
bound was implemented that covers all of: the prior accepted objective's
accumulated error, changed residual squaring and subtraction, signed delta
accumulation, scaling, the candidate objective, the accepted and candidate
Full Scan values, the global delta, local/global consistency, and the
objective-magnitude-dependent replay enclosure. The `gamma_k` model can bound
a positive sequential sum under finite-arithmetic assumptions, but using a
generic `gamma_k` term as an acceptance tolerance would not bound the
incremental update or the difference between both Full Scan decisions. The
cache would also need to carry and refresh its error interval without
overwriting evidence of an already accepted bad update.

The measured error is small in absolute terms but not uniformly small
relative to observed decisions: direct E differed from Full Scan on 489–494
of 503 attempts in the 30-target stress, reached `3.22e-15` at a longer refresh
interval, and E/E2 differed on every candidate at 300 and 1,000 targets. The
earlier D1 closed-loop fixture recorded a minimum acceptance margin of
`1.09552e-16`. Therefore the zero-disagreement scripted checks cannot establish
that a candidate whose estimate differs by `1e-15` remains on the same side of
every production acceptance condition.

The R3 gate's safe disposition is consequently Full Scan for every objective
decision: Fast Accept `0`, Fast Reject `0`, Full Scan fallback `100%`. These
are the required disposition counts, not measurements from a partially
implemented conditional gate. No error-bound violations are reported because
no candidate error interval was emitted; the count is not zero. An unknown or
nonfinite interval also remains on the Full Scan path. This has no speed
benefit over D1 and does not justify adding an objective cache to Production.

Production remains D1 with its existing preflighted Full Replay memory
fallback, sweep-end Full Replay check, endpoint certification, Independent
Replay, local/global consistency, and convergence rules. No threshold,
acceptance branch, recovery behavior, or accepted-state cache was changed.
R4 is diagnostic only: it may measure E1/E2 on a real D1 trajectory, but it
cannot let an unqualified incremental estimate control later search state.
R5's Production outcome is D unless a future study supplies a conservative
decision interval and validates it over the actual solver path.

### R4: closed-loop diagnostic on D1 trajectories

The benchmark added a test-only observer that compares E1/E2 and C estimates
at each candidate generated by a D1-controlled search. None of those estimates
updates the accepted state or controls a later candidate. It compares the
incremental objective proxy with the candidate's D1 Full Scan and records C's
delta-prediction, residual, and objective differences. The benchmark ran five
timed repeats and five sweeps at each size. These are search-only runs; endpoint
certification and Independent Replay were disabled.

| Targets | Full Replay search median | D1 search median | Full Replay candidate replay | D1 candidate replay | Candidate attempts / accepted | E1/E2 objective mismatches; max absolute / relative / ULP | C prediction mismatches; max absolute / relative / ULP | C residual mismatches; max absolute / relative / ULP |
| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |
| 3 | 4.749 ms | 5.146 ms | 0.260 ms | 0.659 ms | 5 / 5 | 5/5; `7.59e-16` / `5.60e-11` / 493,194 | 147; `5.55e-17` / `1.59e-16` / 1 | 147; `5.55e-17` / `5.75e-10` / 4,194,304 |
| 30 | 94.740 ms | 97.620 ms | 7.752 ms | 10.394 ms | 15 / 15 | 15/15; `3.33e-15` / `5.12e-13` / 4,360 | 2,464; `5.55e-17` / `3.18e-16` / 2 | 2,464; `5.55e-17` / `1.83e-8` / 134,217,728 |
| 300 | 1.608 s | 1.293 s | 671.469 ms | 353.133 ms | 130 / 130 | 126/130; `2.94e-15` / `3.35e-13` / 1,811 | 24,785; `5.55e-17` / `3.18e-16` / 2 | 24,785; `5.55e-17` / `1.31e-6` / 8,589,934,592 |
| 1,000 | 10.412 s | 6.434 s | 7.249 s | 3.290 s | 420 / 420 | 414/420; `1.07e-14` / `5.14e-13` / 3,123 | 82,797; `5.55e-17` / `3.18e-16` / 2 | 82,797; `5.55e-17` / `8.02e-7` / 4,294,967,296 |

E2 uses the same exact prediction and objective estimate as E1 here, so its
objective error and decision-proxy counts are the same. Across these traces,
the E1/E2 acceptance proxy disagreed with D1 zero times, and the C proxy did
also. Every candidate in these four runs was accepted, however; this did not
exercise rejection or near-boundary decisions. Zero proxy disagreements are
therefore not evidence of a safe gate. The policy fallback remains Full Scan
for every candidate (100% by policy; no conditional gate was implemented).

C prediction differences were at most two ULP, while its residual ULP and
relative errors grew near zero residuals; those relative and ULP measures are
ill-conditioned around zero. The absolute prediction/residual difference was
at most `5.55e-17`. The largest E1/E2 objective difference was `1.07e-14`.
These measured values are Level C evidence and do not establish error bounds.

The D1 search was 38.2% faster than Full Replay at 1,000 targets, 19.6%
faster at 300 targets, and 3.0% and 8.4% slower at 30 and 3 targets. These
five-repeat FixedNeighbor measurements do not include endpoint qualification.
At 3 targets, the diagnostic D1 run reported a stationarity stop; that is not
a qualified endpoint because certification and Independent Replay were off.
The other D1 runs stopped at their five-sweep budgets. The benchmark did not
let E1, E2, or C control the trajectory, so their final parameter, width, KKT,
and convergence differences are not measurable here. No E1/E2/C convergence
or Independent Replay pass is claimed. R4 passes as a diagnostic study, but
no incremental objective qualifies for R5 Production promotion.

### R5: performance decision and Production outcome

All closed-loop sizes used five timed repetitions and five sweeps. For the
1,000-target Full-Replay-versus-D1 comparison, the search medians were
10.412323 s and 6.433959 s, respectively: D1 was 38.2% faster. Candidate
replay accounted for 7.249199 s of Full Replay search and 3.290498 s of D1
search. D1's candidate replay was 51.1% of its total search time, giving a
2.05x Amdahl upper bound only if that whole candidate-replay cost vanished.
This is a bound on D1 search, not a measured incremental-objective speedup.
The benchmark reported a 121,648,896 B memory preflight peak and 233,512 B
affected-row candidate scratch for D1. Process peak RSS was 242,302,976 B for
the benchmark containing both methods and all five repetitions.

The replay-only 1,000-target comparison used five timed repeats and 84
candidates. These per-candidate values are not a FixedNeighbor total-search
comparison:

| Method | Candidate replay median | Preparation / amortized index cost | Incremental memory observed | Numerical result | Production |
| --- | ---: | ---: | ---: | --- | --- |
| D1 | 7.877 ms (including 0.380 ms incidence preparation) | On-demand incidence per candidate | 321,824 B candidate scratch; no global index | Prediction, residual, and objective bitwise match Full Replay | Yes; only production candidate path |
| E1 | 0.647 ms (0.209 ms replay plus 0.380 ms incidence preparation) | On-demand incidence per candidate | Same 321,824 B candidate scratch; no retained objective cache | Prediction/residual exact; objective differed on 84/84 candidates, max `3.64e-17` | No; no decision-safety bound |
| E2 | 0.260 ms amortized, including CSR build divided over 84 candidates (0.225 ms replay) | 2.916 ms one-time CSR build | 7,760,488 B resident CSR; 7,736,020 B build temporary; 22,777,468 B estimated build peak including accepted state | Prediction/residual exact; objective differed on 84/84, max `3.64e-17` | No; no decision-safety bound |
| C | 0.562 ms in the replay-only comparison | No global index | 88,512 B experiment candidate scratch | Prediction drift; objective differed on 84/84, max `3.64e-17` | No; prediction drift also needs control |

E2 was 2.49x faster than E1 in amortized replay-only candidate time and broke
even after seven candidates in this 84-candidate synthetic fixture. That
supports its replay-only index hypothesis, but not a Production memory
decision: neither E1 nor E2 passed R3, and the 5,938,087-observation input and
its actual support count were unavailable. Process RSS in the E1/E2 replay
benchmark was 205,602,816 B across the fixture, reference checks, methods, and
repetitions; it is not a per-method resident-memory measurement.

The conditional Full Scan fallback rate is 100% by the safe R3 policy because
no validated error interval exists. Fast Accept and Fast Reject counts are
both zero by policy; no conditional-gate execution counts or error-bound
violations exist to report. There is no accepted incremental objective cache
in Production and thus no Production objective drift. Experimental direct
incremental objective showed measured differences up to `3.22e-15` in the
503-attempt study and up to `1.07e-14` on the 1,000-target D1 shadow. No
conservative bound was emitted, so bound violations are N/A rather than zero.

The qualified D1 two-block fixture passed convergence, endpoint certification,
and Independent Replay; global A/C KKT was `4.33411e-9`, projected Width KKT
was `1.62354e-10`, and objective was `2.0199079729737389e-10`. Its complete
Full Replay and D1 state traces, including final eta/beta/prediction/residual
and objective, matched bitwise. The N5 1,000-target five-sweep comparison also
had bitwise-identical block/sweep traces and final state. The R4 shadow runs
did not let E1/E2/C change parameters; final A/B/C/Width deltas for those
methods are therefore not applicable, not zero. Their KKT, convergence
qualification, endpoint certification, and Independent Replay remain
unmeasured. At 1,000 targets the non-qualifying five-sweep D1 diagnostic ended
with global A/C KKT `2.77528e-10` and projected Width KKT `5.63438e-13`, but
stopped at its sweep budget without endpoint certification or Independent
Replay. No claim that the incremental methods preserve final estimates to two
decimal places is supported by these shadow comparisons.

R5 selects Outcome D. Production remains D1 with its existing Full Replay
memory fallback, sweep-end replay, endpoint certification, Independent
Replay, and unchanged acceptance and convergence criteria. No incremental
objective Production branch was added, so there was no unqualified branch to
remove. The direct, compensated, tiled, tree, E1/E2, and C mechanisms remain
test-only research comparisons. They do not meet the 15% Production promotion
criterion because no safe conditional gate was implemented and no valid
incremental-objective FixedNeighbor total was measured. No real large-data
benchmark or complete end-to-end algorithm timing was available.

Direct incremental accumulation was the fastest simple R1 objective update;
Kahan and Neumaier did not consistently reduce error, and tile/tree cache
recompute was slower on the measured fixture. E2 was fastest in replay-only
incidence comparison, while E1 had lower retained memory. Conditional Full
Scan is the only proven-safe objective decision policy from this study, but
at 100% fallback it delivers no speedup over D1.

R5 ran after the final code change from R4. Serial CTest remained 19/19 green
in 20.66 s. Git changes for this stage are documentation only.

## Roadmap status

| Stage | Status | Result |
| --- | --- | --- |
| N0 | Complete | Baseline telemetry and deterministic state/candidate digests recorded. |
| N1 | Complete | A/B1/B2/C/D1/E measured at 3, 30, and 300 synthetic targets; D1 passed bitwise candidate replay checks. |
| N2 | Complete | Five-method shadow; Full Replay retained every decision; D1 remained bitwise exact, while B/C/E accumulated drift. |
| N3 | Complete | 503-attempt refresh/drift stress at 1/10/50/100/sweep; B methods drifted more with longer intervals, C/E retained objective error, D1 stayed exact; all scripted rejections rolled back bitwise. |
| N4 | Complete | D1 controlled candidate results in closed loop and reached Level A on the converged two-block fixture; Independent Replay passed. |
| N5 | Complete | Five-run closed-loop comparison through 1,000 synthetic targets; D1 reached 1.62× search speedup with bitwise-identical traces; D2 was exact but rejected for a small gain and 11.88 MB retained index. |
| N6 | Complete | Adopted preflighted production D1 after Level A closed-loop/Independent Replay checks; 19/19 CTest groups passed; Full Replay remains the memory fallback and final reference. |

| Incremental stage | Status | Result | Commit |
| --- | --- | --- | --- |
| R0 | Complete | Acceptance and numerical-level audit; D1 baseline measured, Production unchanged. | `205c0815` |
| R1 | Complete | Compared direct, scaled, compensated, tiled, and tree objectives; no reliable error bound. | `7584c5f8` |
| R2 | Complete | E2 beat E1 in replay-only tests but added 7.76 MB CSR at 1,000 targets; no Production integration. | `2deba209` |
| R3 | Complete, candidate unqualified | No conservative decision bound; safe disposition is D1 Full Scan on every candidate. | `2624c585` |
| R4 | Complete, diagnostic only | E1/E2/C shadowed D1 traces; errors observed and no safe gate or candidate-controlled convergence established. | `46a4c288` |
| R5 | Complete | Outcome D; retain D1, report replay-only gains and leave incremental methods test-only. | This final study commit |

## Closed-loop Incremental Objective reconsideration

### E0: refreshed D1 qualification baseline

This follow-up starts from `develop` HEAD `0ae5f00a1bacbd910d8338da6609ee8724f14d83`,
with a clean tree and `origin/develop` at the same locally recorded commit. The
Release Ninja build in `build/replay-study` uses `/usr/bin/c++`, system
dependencies, OpenMP enabled, and `BUILD_TESTING=ON`. The baseline build
completed with `--parallel 4`; serial CTest passed all 19/19 groups in 22.47 s.
OpenMP and common BLAS thread limits were set to four. No repository source
changes were part of this baseline.

The D1 qualification run used the deterministic synthetic `chain` workload,
initial log-width `log(0.55)` for every atom, zero initial A/C coefficients,
width upper bound 2.2, 12-atom cores, the existing `PracticalJointConvergenceCriteria`,
and a 30-sweep maximum. Every component stopped at `block-stationary` before
the budget. Endpoint Certification, Independent Replay, runtime convergence,
and parameter-stability assessment passed. These are one-run qualification
records; they are not repeated performance measurements. The known synthetic
truth is A=2, C=0.2, Width=0.5 for each target.

| Targets | Observations / support entries | Blocks | Sweeps | Candidate accept / reject | Full Objective | Global A/C KKT | Projected Width KKT | Endpoint / Independent Replay | Preflight peak |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| 3 | 1,425 / 1,545 | 1 | 4 | 4 / 0 | `2.0199078132865162e-10` | `3.33e-16` | `3.89e-15` | pass / pass | 447,445 B |
| 30 | 13,710 / 15,450 | 3 | 5 | 15 / 0 | `1.517076030684615e-10` | `1.60e-9` | `1.88e-11` | pass / pass | 3,930,461 B |
| 300 | 136,560 / 154,500 | 26 | 5 | 130 / 0 | `1.4790433590012254e-10` | `5.07e-10` | `1.88e-12` | pass / pass | 36,716,466 B |
| 1,000 | 455,060 / 515,000 | 84 | 5 | 420 / 0 | `1.4790435904528611e-10` | `2.78e-10` | `5.63e-13` | pass / pass | 121,648,896 B |

The final D1 parameter medians and the single-run qualification timing were:

| Targets | A median | C median | Width median | Search / wall time | Qualification-process peak RSS |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 1.9999904 | 0.20000529 | 0.49999991 | 6.88 / 8.13 ms | 23.6 MB |
| 30 | 2.0000014 | 0.20000169 | 0.50000007 | 111.69 / 133.11 ms | 35.7 MB |
| 300 | 2.0000001 | 0.20000037 | 0.50000001 | 1.294 / 1.519 s | 101.3 MB |
| 1,000 | 2.0000005 | 0.19999991 | 0.50000003 | 6.434 / 7.187 s | 275.7 MB |

The RSS column belongs to the one-run qualification process, which ran all
three methods and retained candidate diagnostic references; it is not a
per-method Production RSS measurement. These times are qualification-run
records, not repeated performance results.

All widths were feasible, the existing global convergence and parameter-
stability criteria passed, and local-conditioning assessment covered every
block. The 1,000-target A/C/Width medians were 2.0000005 / 0.19999991 /
0.50000003; P95 absolute errors against the known truth were `6.15e-5`,
`1.48e-5`, and `3.01e-6`. Its A/C/Width RMSEs were `4.36e-5`, `1.05e-5`, and
`2.13e-6`. The 300-target RMSEs were `4.36e-5`, `1.05e-5`, and `2.14e-6`; the
30-target values were `4.45e-5`, `1.05e-5`, and `2.21e-6`; the 3-target values
were `4.06e-5`, `9.16e-6`, and `2.14e-6`. This is a deterministic uniform
truth fixture, not real cryo-EM evidence or a low-signal precision study.

The qualification process peak RSS was 23.6 MB, 35.7 MB, 101.3 MB, and 275.7
MB at the four sizes, respectively. These values cover a diagnostic process
that ran all three methods and retained per-candidate reference diagnostics;
they are not per-method Production RSS values. The 5,938,087-observation real
input remains unavailable.

The prior N4 fixture remains the explicit Full Replay versus D1 converged
comparison: 3 targets, two blocks, five sweeps, ten accepted candidates,
bitwise-identical final state, Endpoint Certification, and Independent Replay
passed. The broader baseline above supplies a qualified D1 endpoint at each
requested synthetic size. The old N5 timing runs remain search-only and are
not reused as qualification evidence.

### E1: controlled closed-loop search

The test-instrumented build now has three actual search controls: D1 Full Scan,
E1-Direct, and E1-Guarded. E1 computes the exact affected-row prediction and
residual patch, then adds the affected squared-residual difference to its own
currently accepted objective cache. That value controls candidate objective,
global delta, local/global consistency, acceptance, and the next block's state.
Parameters and prediction/residual rows are committed only after those checks
pass. A rejected or failed block leaves the last accepted state available to
the existing retry path. The candidate observer is a separate diagnostic mode;
it copies full vectors and runs a full objective reference, and its cost is not
included in Performance Mode. E1-Direct, E1-Guarded, refresh intervals, and
their counters compile only with `RHBM_GEM_TEST_INSTRUMENTATION`; the ordinary
Production build still selects D1 and has no E1 runtime selector.

Each qualified method below stopped at `block-stationary`, before its sweep
budget, with unchanged convergence criteria. Both E1 methods passed Endpoint
Certification, Independent Replay, parameter-stability assessment, and runtime
convergence. E1-Direct's cached final objective is shown separately from the
independently recomputed objective.

| Targets | Candidate accept / reject (D1, Direct, Guarded) | Sweeps | Final cached objective D1 / Direct / Guarded | Independently recomputed objective | Global A/C KKT | Projected Width KKT | Endpoint / Independent Replay | Max |ΔA| / |ΔC| / |ΔWidth| vs D1 |
| ---: | --- | ---: | --- | ---: | ---: | ---: | --- | --- |
| 3 | 4/0, 4/0, 4/0 | 4 | `2.0199078133e-10` / `2.0199153986e-10` / `2.0199078133e-10` | `2.0199078132865188e-10` | `3.33e-16` | `3.89e-15` | pass / pass | 0 / 0 / 0 |
| 30 | 15/0, 15/0, 15/0 | 5 | `1.5170760307e-10` / `1.5170575065e-10` / `1.5170760307e-10` | `1.5170760306846197e-10` | `1.60e-9` | `1.88e-11` | pass / pass | 0 / 0 / 0 |
| 300 | 130/0, 130/0, 130/0 | 5 | `1.4790433590e-10` / `1.4791233916e-10` / `1.4790433590e-10` | `1.4790433590012404e-10` | `5.07e-10` | `1.88e-12` | pass / pass | 0 / 0 / 0 |
| 1,000 | 420/0, 420/0, 420/0 | 5 | `1.4790435905e-10` / `1.4789356718e-10` / `1.4790435905e-10` | `1.4790435904529177e-10` | `2.78e-10` | `5.63e-13` | pass / pass | 0 / 0 / 0 |

The ordinary-chain E1 search therefore did not change candidate decisions,
sweep count, independently recomputed final objective, or final A/C/Width on these four fixtures. The
independent references compared each E1 candidate with the full objective at
that same candidate and its actual accepted-before state; they did not compare
different D1 and E1 trajectories. All 1,000-target blocks received local
conditioning assessment. The separate width-upper-bound test and initial
widths 0.35, 0.55, and 0.9 also passed for D1 and E1, with the same final
parameters.

The final parameter medians (and therefore D1/E1 parameter deltas) were:

| Targets | A median | C median | Width median | D1/E1 max difference distribution, for each parameter |
| ---: | ---: | ---: | ---: | --- |
| 3 | 1.9999904 | 0.20000529 | 0.49999991 | mean / median / P95 / P99 / max = 0 |
| 30 | 2.0000014 | 0.20000169 | 0.50000007 | mean / median / P95 / P99 / max = 0 |
| 300 | 2.0000001 | 0.20000037 | 0.50000001 | mean / median / P95 / P99 / max = 0 |
| 1,000 | 2.0000005 | 0.19999991 | 0.50000003 | mean / median / P95 / P99 / max = 0 |

The per-candidate audit found rounding drift but no reference acceptance
disagreement. ULP counts are large because the objective and some decision
deltas are close to zero; the absolute differences give the useful scale.

| Targets | Direct max absolute / relative error | Direct max ULP | Guarded Full Scan fallbacks | Guarded max relative raw-estimate error after refreshes | Min decision margin | Reference decision disagreements, Direct / Guarded |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | `7.59e-16` / `3.76e-6` | 29,344,345,814 | 1/4 (25.0%) | `5.60e-11` | `4.12e-18` | 0 / 0 |
| 30 | `3.33e-15` / `1.22e-5` | 71,661,907,826 | 3/15 (20.0%) | `6.12e-12` | `1.10e-16` | 0 / 0 |
| 300 | `1.29e-14` / `5.41e-5` | 309,611,272,926 | 28/130 (21.5%) | `2.43e-11` | `1.00e-17` | 0 / 0 |
| 1,000 | `2.01e-14` / `7.30e-5` | 417,490,110,813 | 86/420 (20.5%) | `3.27e-11` | `6.40e-18` | 0 / 0 |

Guarded's raw incremental estimate still has rounding error on fast candidates;
fallback candidates use a full objective scan. The guard band is empirical and
does not certify an error bound. A smaller 3-target band of 0.01 fell back on
1/10 candidates, the 0.1 band on 2/10, and the 1.0 band on 10/10. All three
settings converged and had zero reference decision disagreements on that
fixture. The 0.1 policy is the measured Guarded configuration above.

No outer candidate was rejected in these converged chain runs; the 300-target
qualification captured 260/260 accepted local profile trials as well. The
forced-recovery E1 test still passed after an injected untrusted profile
failure. A numeric candidate-rejection fixture was not established, so rejection
decision coverage remains incomplete. The signed-cancellation/near-zero
amplitude fixture produced 42/42 accepted replay candidates and no reference
decision disagreement, but both D1 and E1 reached the 500-sweep budget without
qualification: global A/C KKT was `2.22e-16`, while projected Width KKT was
`1.22801358e-9`, just above the unchanged `1e-9` requirement. Neither method
received Endpoint Certification there. Its E1 max candidate/cache error was
`7.87e-16`; this was a shared convergence limit, not an E1-only state drift.

### E2: cache reconciliation and refresh intervals

For each periodic refresh, the experiment retained the exact affected-row
state, recomputed the full replay, compared prediction/residual and objective,
and refreshed only the objective after those checks passed. The full replay
audit at every sweep end remains in place and precedes its optional objective
refresh. The 300-target direct run had 130 accepted updates and qualified under
all intervals:

| Accepted updates between refreshes | Refresh count | Max cache error at refresh | Max candidate objective error, absolute / relative / ULP | Refresh time | Decision disagreements / refresh failures / positive true deltas |
| ---: | ---: | ---: | ---: | ---: | --- |
| 10 | 13 | `1.02e-14` | `1.02e-14` / `2.06e-12` / 14,484 | 72.2 ms | 0 / 0 / 0 |
| 50 | 2 | `8.00e-15` | `1.29e-14` / `1.23e-9` / 7,398,953 | 11.1 ms | 0 / 0 / 0 |
| 100 | 1 | `8.00e-15` | `1.29e-14` / `5.14e-5` / 309,611,272,926 | 5.56 ms | 0 / 0 / 0 |
| Every sweep | 5 | `8.01e-15` | candidate-reference scan not attached | Included in mandatory sweep replay | 0 / 0 / 0 |

The periodic intervals each ended with Endpoint Certification and Independent
Replay. Every-sweep reconciliation saw at most 26 accepted updates since the
preceding refresh. The every-candidate reference is D1 Full Scan. No true
objective increase outside the existing acceptance enclosure was observed;
therefore no cache refresh masked such an event. No rollback was triggered, and
the two checkpoint representations (full accepted state versus parameters with
state reconstructed by replay) were not compared. The refresh interval study
is complete; E2 checkpoint/rollback qualification is partial. No extra
checkpoint buffer was allocated, and the measured cache/scratch paths stayed
within the existing preflight budget.

### E3/E4: guard policy and estimate qualification

The E1 errors did not change an acceptance decision in the audited cases,
including minimum reference margins from `4.1e-18` to `1.0e-16`. This is
evidence for these synthetic fixtures, not a general safety guarantee. Guarded
reduced observed raw-objective error after its fallbacks and sweep refreshes,
but its 20–25% fallback rate consumed some of Direct's gain.

For the 3, 30, 300, and 1,000-target qualified chain fixtures, maximum absolute
parameter differences between D1 and either E1 method were zero for A, C, and
Width. Their mean / median / P95 / P99 / maximum absolute errors against known
truth are below, followed by A/C/Width RMSE and maximum relative errors.

| Targets | A error distribution | C error distribution | Width error distribution | A / C / Width RMSE | Max relative error A / C / Width |
| ---: | --- | --- | --- | --- | --- |
| 3 | `3.31e-5 / 2.43e-5 / 6.54e-5 / 6.54e-5 / 6.54e-5` | `7.37e-6 / 5.29e-6 / 1.48e-5 / 1.48e-5 / 1.48e-5` | `1.69e-6 / 1.67e-6 / 3.31e-6 / 3.31e-6 / 3.31e-6` | `4.06e-5 / 9.16e-6 / 2.14e-6` | `3.27e-5 / 7.41e-5 / 6.62e-6` |
| 30 | `4.01e-5 / 4.53e-5 / 6.30e-5 / 6.71e-5 / 6.71e-5` | `9.38e-6 / 1.09e-5 / 1.47e-5 / 1.48e-5 / 1.48e-5` | `2.00e-6 / 2.21e-6 / 3.14e-6 / 3.68e-6 / 3.68e-6` | `4.45e-5 / 1.05e-5 / 2.21e-6` | `3.35e-5 / 7.40e-5 / 7.37e-6` |
| 300 | `3.93e-5 / 4.35e-5 / 6.15e-5 / 6.16e-5 / 6.30e-5` | `9.42e-6 / 1.04e-5 / 1.47e-5 / 1.48e-5 / 1.48e-5` | `1.93e-6 / 2.13e-6 / 3.01e-6 / 3.02e-6 / 3.14e-6` | `4.36e-5 / 1.05e-5 / 2.14e-6` | `3.15e-5 / 7.40e-5 / 6.29e-6` |
| 1,000 | `3.92e-5 / 4.37e-5 / 6.15e-5 / 6.16e-5 / 6.45e-5` | `9.41e-6 / 1.05e-5 / 1.48e-5 / 1.48e-5 / 1.50e-5` | `1.92e-6 / 2.14e-6 / 3.01e-6 / 3.02e-6 / 3.21e-6` | `4.36e-5 / 1.05e-5 / 2.13e-6` | `3.22e-5 / 7.51e-5 / 6.43e-6` |

At 1,000 targets, A/C/Width maximum absolute truth errors were `6.45e-5`,
`1.50e-5`, and `3.21e-6`, all below 0.01 in their stored parameter units.
This verifies the two-decimal scale on the uniform synthetic truth fixture;
it does not establish it for low-signal or weakly identifiable real data. The
FixedNeighbor state has two fitted coefficients per target (A and C) plus
Width; it has no separate B coordinate, so a ΔB statistic is not applicable
to this implementation.

The active-width-upper-bound test passed with the existing projected KKT
threshold. The near-zero amplitude/signed-charge fixture remained unqualified
for both methods as described above. E1 qualification for rank-deficient or
near-singular local blocks, weakly identifiable widths, multiple outer
rejections, and a broad set of noisy/overlapping initializations remains
incomplete. Consequently E4 is partial even though the standard synthetic
fixtures pass.

### E5: five-run performance and memory

Performance Mode used five measured runs per method and one warm-up, with the
same synthetic input, initial parameters, and 4-thread cap. It attached no
per-candidate diagnostic Full Scan. It ran through Endpoint Certification,
Independent Replay, local-conditioning assessment, and parameter-stability
assessment. Search medians exclude the final checks; wall medians include
them. The range is the five search measurements.

| Targets | Candidates | D1 search median (range) | E1-Direct search median (range) | E1-Guarded search median (range) | Search improvement Direct / Guarded | Wall median D1 / Direct / Guarded | Guard fallbacks | Peak RSS D1 / Direct / Guarded |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 4 | 4.383 ms (4.363–4.498) | 4.283 ms (4.264–4.301) | 4.288 ms (4.284–4.313) | 2.3% / 2.2% | 5.502 / 5.392 / 5.392 ms | 1/4 | 23.64 / 23.66 / 23.66 MB |
| 30 | 15 | 97.616 ms (97.450–98.103) | 94.117 ms (93.679–96.205) | 94.726 ms (94.502–96.997) | 3.6% / 3.0% | 118.712 / 115.119 / 115.888 ms | 3/15 | 35.75 / 35.80 / 35.75 MB |
| 300 | 130 | 1.291 s (1.285–1.293) | 1.013 s (1.004–1.017) | 1.066 s (1.064–1.206) | 21.5% / 17.5% | 1.514 / 1.236 / 1.289 s | 28/130 | 93.47 / 94.58 / 94.40 MB |
| 1,000 | 420 | 6.413 s (6.408–6.431) | 3.372 s (3.372–3.381) | 3.990 s (3.986–4.193) | 47.4% / 37.8% | 7.165 / 4.125 / 4.743 s | 86/420 | 247.61 / 246.27 / 248.09 MB |

Endpoint Certification and Independent Replay medians at 1,000 targets were
71 ms and 37 ms, respectively; at 300 targets they were 21 ms and 11 ms.
Candidate evaluation fell from 3.284 s on D1 to 0.243 s on either E1 method at
1,000 targets. Guard fallback consumed 0.621 s there and 0.060 s at 300
targets. Sweep-end full replay was already part of every method's search; the
objective scalar refresh added no second replay. At 1,000 targets D1 candidate
evaluation accounts for 51.2% of search time, so the Amdahl limit is 2.05x;
Direct measured 1.90x and Guarded 1.61x on search time.

All three methods used the same memory preflight: 447,445 B at 3 targets,
3,930,461 B at 30, 36,716,466 B at 300, and 121,648,896 B at 1,000. Affected
row scratch was 58,976 B at 3 targets and 233,512 B for the larger sizes. E1
kept on-demand incidence and one objective scalar; it added no CSR, full-state
checkpoint, or rollback buffer. The 3.2 GB budget was not approached on these
synthetic sizes. The approximately 37,406-target / 5,938,087-observation input
was unavailable, and no end-to-end full algorithm timing is claimed.

### E6: Production decision

The controlled search establishes that E1 can run to qualified endpoints on
the standard 3–1,000-target synthetic fixtures and can deliver greater than
15% search improvement at 300 and 1,000 targets. It does not complete the
requested hard-case qualification: the low-signal width boundary fixture is
unqualified for both methods, outer rejected-candidate reference coverage is
absent, rank-deficient/weak-width E1 coverage is incomplete, and objective
checkpoint rollback was not exercised. The Guarded band remains empirical.
The real cryo-EM component was also unavailable.

Production remains D1 affected-row exact recompute plus Full Objective Scan.
No Production acceptance, KKT, width-bound, convergence, endpoint, independent
replay, or memory threshold changed. Sweep-end replay, Endpoint Certification,
Independent Replay, and the 3.2 GB preflight remain in the formal path. No E1
runtime selector or third-party dependency was added. The test-only E1,
Guarded, refresh, and qualification benchmarks remain as research evidence.
The final CMake build without test instrumentation passed, and serial CTest
passed 19/19.

### E1 follow-up roadmap

| Stage | Status | Result | Commit |
| --- | --- | --- | --- |
| E0 | Complete | Refreshed D1 baseline qualified on 3, 30, 300, and 1,000 synthetic targets. | `2ffe6b48` |
| E1 | Complete for tested fixtures | Incremental objective controlled closed-loop decisions; all audited candidate decisions matched their reference; full CTest passed. | `92acb92c` |
| E2 | Partial | 10/50/100-update and every-sweep reconciliation measured; checkpoint representation and rollback were not exercised. | `92acb92c` |
| E3 | Complete for tested fixtures | Guard bands 0.01/0.1/1.0 compared; use remains test-only and empirical. | `92acb92c` |
| E4 | Partial | Standard, active-bound, multi-initial, and known-truth qualifications passed; hard weak-signal/rank/rejection coverage remains. | `92acb92c` |
| E5 | Partial | Five-run certified performance covered 3/30/300/1,000 synthetic targets; real large data and end-to-end runtime unavailable. | `92acb92c` |
| E6 | Complete | Retain D1 Production; repository report updated; no Notion edits or push. | This report commit |

## E1 promotion follow-up

### P0 — Historical Full Replay baseline

Commit ancestry confirms that `6787a777da83dd820e72899a9bf4f53b91e6ad93`
(`test(joint): benchmark large replay workloads`) is the direct parent of
`6f41026cde46153d537a8f41d7ab749fcd14f4c0` (`perf(joint): use preflighted
affected-row replay`). The latter is the D1 production introduction. The F0
commit's production candidate path calls `ReplayInto` for a full prediction,
residual, and objective replay; its test-only candidate observer does not
replace that path. Thus F0 is the last committed Production Full Replay
baseline before D1, rather than the current-code F1 reference.

The source diff from F0 to the D1-introduction commit adds exact affected-row
candidate prediction, memory census/planning, and a memory-bounded fallback to
the existing full replay. It does not change the local A/C solver, bounded
Width LM, candidate consistency/acceptance conditions, Global KKT tests,
convergence criteria, Endpoint Certification, or Independent Replay. For an
affected row, contributors are evaluated in the existing global atom order;
the objective still scans all informative rows. Unaffected rows reuse the
accepted residual, and the existing sweep-end Full Replay checks that cache.
The memory fallback calls the original full replay directly. This is a replay
implementation and resource-planning change, with possible last-bit floating
order effects covered by replay consistency checks; no other intended
optimization or convergence semantic change was introduced by D1.

F0 was exported with `git archive` into `/private/tmp`, then built without
switching `develop` or creating a branch/worktree. The only archive-local
change added a diagnostic test; no historical solver or acceptance code was
modified. Release settings used C++20, CMake 4.4.3, Ninja 1.13.2, AppleClang
21.0.0, Boost 1.92.0, GTest 1.18.0, system Eigen, and OpenMP 5.1, with
`BUILD_TESTING=ON`, `RHBM_GEM_ENABLE_UMAP=OFF`, and a four-worker build. The
targeted archived F0 test ran successfully.

The signed-cancellation input is the existing `chain` 3-target synthetic
workload after the exact observation mutation in
`IncrementalObjectiveHandlesSignedChargeCancellationAndLowAmplitude`. It has
1,425 rows and 1,545 support entries. Its canonical SHA-256 is
`0b5c3798029192eb18c137134db3214e9647fa87868433366bd744210cce34eb`. The
generator source is byte-identical between F0 and current HEAD; the current
test suite now pins this digest so future runs can detect fixture drift.

All runs used atom order 1/2/3, initial width 0.55, upper bound 2.2, core size
2, objective scale `5.678970475776449`, `PracticalJointConvergenceCriteria`,
and the listed sweep budget. No KKT or stopping threshold was changed.

| F0 sweep budget | Outcome | Final objective | Global A/C KKT | Projected Width KKT | Endpoint certification | Last 100 sweeps |
| ---: | --- | ---: | ---: | ---: | --- | --- |
| 200 | Unqualified, `block-sweep-budget` | `4.0589644398391108e-10` | `2.2204460492503131e-16` | `1.2280135797126945e-9` | Not reached | Objective and parameters stationary; 0 width movement |
| 500 | Unqualified, `block-sweep-budget` | `4.0589644398391108e-10` | `2.2204460492503131e-16` | `1.2280135797126945e-9` | Not reached | Objective and parameters stationary; 0 width movement |
| 1,000 | Unqualified, `block-sweep-budget` | `4.0589644398391108e-10` | `2.2204460492503131e-16` | `1.2280135797126945e-9` | Not reached | Objective and parameters stationary; 0 width movement |

The practical projected Width KKT limit is `1e-9`; F0 remains above it by
`2.280135797126945e-10`. Prediction and objective replay differences were
zero on each sweep. Each budget produced the same final parameters:

```text
eta   = [-0.69313238734324734, -0.69315371810798543, -0.69314737995827003]
width = [0.50000739666305904, 0.49999673123666483, 0.49999990030084757]
beta  = [1.8623863311018945e-09, -0.20000155320669047,
         1.9999352041743306, 0.20001477549941502,
         1.9999905326991614, -0.19999470127506408]
```

The first 24 sweeps accepted 42 changed block updates; all later block attempts
were unchanged. There were no candidate rejections, recovery retries, active
width bounds, nonfinite values, or cache drift. The last 100 sweeps in every
run had exactly the same objective, A/C KKT, Width KKT, and parameters. This is
direct evidence that the signed-cancellation/near-zero Width failure existed
in historical Production Full Replay before D1; its CTest success must not be
read as a qualified endpoint.

### P1 — Signed-cancellation failure attribution

The same frozen input digest and scale were run against four implementations:
archived F0, current F1 using `FullReplayForTesting` for candidate replay, D1
using exact affected-row prediction plus full objective scan, and E1-Direct
using exact affected-row prediction plus the incremental objective. F1 is a
current-code reference measurement and is not substituted for the archived F0
run. Each implementation ran with sweep budgets 200, 500, and 1,000.

| Budget | F0 outcome | F1 outcome | D1 outcome | E1 outcome | Full Replay objective in all modes | E1 cached objective at stop | Global A/C KKT | Projected Width KKT |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| 200 | Unqualified | Unqualified | Unqualified | Unqualified | `4.0589644398391108e-10` | `4.0589565781233565e-10` | `2.2204460492503131e-16` | `1.2280135797126945e-9` |
| 500 | Unqualified | Unqualified | Unqualified | Unqualified | `4.0589644398391108e-10` | `4.0589565781233565e-10` | `2.2204460492503131e-16` | `1.2280135797126945e-9` |
| 1,000 | Unqualified | Unqualified | Unqualified | Unqualified | `4.0589644398391108e-10` | `4.0589565781233565e-10` | `2.2204460492503131e-16` | `1.2280135797126945e-9` |

Every run stopped at `block-sweep-budget`; convergence and Endpoint
Certification were not reached, and Independent Replay was not attempted.
For each budget, F1, D1, and E1 had the same final parameters as F0. The maximum
per-parameter E1-to-D1 and F1-to-D1 differences in A/C and Width were zero in
this build. Their final truth RMSEs were also identical: A `3.7807092865e-5`,
C `9.1068530959e-6`, and Width `4.6692389841e-6`. Full objective and KKT were
constant throughout the final 100 sweeps. All methods had zero parameter
movement and zero active Width bounds during those sweeps.

E1's final cached objective differed from Full Replay by
`7.8617157543253364e-16`, while the largest candidate incremental-objective
difference was `7.8691393678997912e-16`. Both are within the existing
`BlockObjectiveReplayEnclosure` (approximately `1e-12` here). Prediction
replay error was zero. E1 made zero candidate-decision disagreements against
the Full Replay objective reference; it had 42 accepted changed candidates,
zero candidate rejections, and zero recovery retries at each budget. The other
block attempts ended unchanged with `block-objective-increase` before a global
candidate replay: 358, 958, and 1,958 block records at budgets 200, 500, and
1,000 respectively. These are local no-update outcomes, not rejected global
candidate transactions. LM local stop reasons, candidate decisions, and each
sweep's full objective, KKT values, parameter changes, replay errors, and
active-bound count were written to:

- F0: `/private/tmp/rhbm-gem-f0-signed-trace.csv`
- F1/D1/E1: `/private/tmp/rhbm-p1-signed-cancellation.csv`

The current diagnostic trace contains all 5,100 sweeps, 10,200 block records,
30,843 profile-trial records, and 378 candidate events. All three modes passed
the 436,948-byte memory preflight with no Full Replay memory fallback. No
nonfinite state or resource failure occurred.

This is a **Legacy Baseline Limitation — Non-blocking for E1 Promotion**. It
fails the unchanged projected Width KKT limit in historical F0, current F1,
D1, and E1, with identical final objective, KKT, and parameters. E1 shows no
additional failure in this case. The fixture is unqualified; the passing test
records and asserts that status and does not certify it as a successful solve.

### P2 — Numerical non-regression qualification

The archived F0 qualification harness and current benchmark used the same
versioned synthetic input stream. A new regression test pins the four standard
fixture SHA-256 values, and the hard-case comparison test pins the additional
fixture hashes below. F0 was the archived 6787a777 solver; F1 was a current-code
Full Replay candidate reference; D1 was the current exact affected-row
prediction plus full-objective scan; E1 was current exact affected-row
prediction plus incremental objective. All runs used the unchanged practical
convergence criteria, bounds, initial states, core sizes, and sweep budgets.

| Targets | SHA-256 | F0 outcome / sweeps | D1 and E1 outcome / sweeps | F0 objective; D1/E1 independent objective | Global A/C KKT; Width KKT | Truth RMSE A / C / Width |
| ---: | --- | --- | --- | --- | --- | --- |
| 3 | `4b3381463b3e7d97099342a9de86df01bf435cf2d64b8ea96af9ac0ee2eab176` | Qualified / 4 | Qualified / 4 | `2.0199078132865162e-10`; `2.0199078132865188e-10` | `3.33e-16`; `3.89e-15` | `4.0637447e-5 / 9.1600727e-6 / 2.1408647e-6` |
| 30 | `13cb46a3575b92368d8393f168cf45cc7db9615e0e2cce06b7f97326a8580c57` | Qualified / 5 | Qualified / 5 | `1.517076030684615e-10`; `1.5170760306846197e-10` | `1.6023889e-9`; `1.8806623e-11` | `4.4464113e-5 / 1.0451471e-5 / 2.2118189e-6` |
| 300 | `bfb80c5d252cee06f7e7f57906c39123adc007d37bccbe34fe546390268d8544` | Qualified / 5 | Qualified / 5 | `1.4790433590012254e-10`; `1.4790433590012404e-10` | `5.0672642e-10`; `1.8782753e-12` | `4.3627912e-5 / 1.0462239e-5 / 2.1372762e-6` |
| 1,000 | `3d81218fbf38f8c661551631fbe4038694c67254b26283f51d10d49299ee33a3` | Qualified / 5 | Qualified / 5 | `1.4790435904528611e-10`; `1.4790435904529177e-10` | `2.7752783e-10`; `5.6343818e-13` | `4.3588368e-5 / 1.0456485e-5 / 2.1346573e-6` |

Every run reached Endpoint Certification and Independent Replay. The
independent objectives differ from F0 by at most `5.66e-24`; D1 and E1 final
A/C/Width parameter differences were all exactly zero in these runs. E1's
maximum candidate objective error against a full-residual calculation was
`7.59e-16`, `3.33e-15`, `1.29e-14`, and `2.01e-14`, respectively. There were
zero candidate acceptance disagreements. All tested candidates were accepted
(4/15/130/420); no candidate rejection was exercised by these standard inputs.
E1's pre-promotion objective cache was not reconciled on every sweep in this
benchmark; its largest final sweep replay difference was `1.08e-14` at 1,000
targets. The audited P4 benchmark below includes the proposed per-sweep
reconciliation after each successful Full Replay audit. It remains test-only
because P3 did not pass its rejection-state-safety gate.

The separately hashed active upper-width-bound input
(`4b3381463b3e7d97099342a9de86df01bf435cf2d64b8ea96af9ac0ee2eab176`)
also qualified on F0, F1, D1, and E1. Each had one active upper bound, objective
`2.02540027678003e-10`, global A/C KKT `1.1910380459667635e-9`, projected Width
KKT `4.4614645311469303e-11`, Endpoint Certification, and Independent Replay.
The bound was feasible and the KKT used the existing projected-bound
criterion. The three initial-width runs `.35`, `.55`, and `.9` qualified on
F0, F1, D1, and E1 in 7, 5, and 23 sweeps. Their F0 objectives were
`2.0199080158137056e-10`, `2.0199079729737389e-10`, and
`2.0199079727070732e-10`; current F1/D1/E1 objectives matched to rounding,
with zero E1-to-D1 parameter differences.

| Hard fixture | SHA-256 prefix | F0 | F1 / D1 / E1 | Objective, KKT and E1 audit | Assessment |
| --- | --- | --- | --- | --- | --- |
| High-overlap cube, 8 targets | `2a1894dc5dc5f7789fbc802477e77e08f520a6c25f0c7eac6b42a4caf3f717b2` | Qualified | Qualified / Qualified / Qualified | Objective `1.1487848300043995e-10`; A/C KKT `6.9202588e-10`; Width KKT `1.4783064e-11`; 10 E1 candidates, zero disagreements, maximum error `2.30e-15` | Pass |
| Rank-deficient constant-support | `fd39aa13434e60a2193641f97ff9bdd05669fb6bc1b57a4d52f454bdacaa7151` | `block-recovery-exhausted` | Same / Same / Same | Objective `0.49999999999999944`; A/C KKT `0.4828122672`; Width KKT `0`; no parameter movement; local rank failure preceded E1 candidate evaluation | Legacy unqualified; E1 added no failure |
| Near-singular overlapping support | `fc9e6c0ec7f746d9de6221d1b680334bd870d7797012a666c32cdcdb467454ef` | `block-recovery-exhausted` | Same / Same / Same | Objective `0.24948812932689915`; A/C KKT `0.6897680220`; Width KKT `0`; no parameter movement; local KKT failure preceded E1 candidate evaluation | Legacy unqualified; E1 added no failure |
| Noise-only low signal | `0d2ecf601cb21d6ce3d4e23f1aa7f37b78f17fe7ec58a19c49b397374ad99d2f` | `block-sweep-budget` | Same / Same / Same | Objective `3.5558894899614887e-8`; A/C KKT `4.4687101e-13`; Width KKT `1.5473044e-12`; 2 active width bounds; final changes A/C `2.33e-13`, log Width `1.7644e-6`; 184 E1 accepts, zero disagreements, maximum candidate error `1.39e-22` | Legacy unqualified; E1 added no failure |
| Weak-width direction | `aa1d50f8e6dbcdaa22e401be857e849573ca0c6457d714583029003f23f9ef9c` | Qualified | Qualified / Qualified / Qualified | Objective `1.1716558769804921e-32`; A/C KKT `1.1102230e-16`; Width KKT `5.5511151e-17`; 2 E1 candidates, zero disagreements, maximum error `6.14e-32` | Pass |

The rank-deficient and near-singular fixtures are common solver limitations,
not incremental-objective decision failures: all four modes retain the same
finite accepted state and full-replay objective, and both failures occur
before any outer candidate replay. The noise-only fixture has a passing
endpoint assessment and Independent Replay but does not meet search
convergence by its 100-sweep budget; those passing checks do not change its
unqualified classification. F0, F1, D1, and E1 return exactly the same
objective, KKT, parameters, and active-bound count for it. No convergence
limit was relaxed.

The P2 targeted tests passed for the digest-pinned standard cases, active
upper bound, three initial widths, high overlap, rank deficiency, near
singularity, weak Width, low signal, and the P1 signed-cancellation baseline.
Local conditioning tests independently classify exact singularity and
ill-conditioning. The outer-transaction rejection and recovery cases remain
for P3 because none of these qualification fixtures naturally rejects a
global candidate. These results support numerical non-regression; they do not
by themselves qualify Candidate Rejection/Recovery state safety.

### P3 — Candidate state, recovery, and memory safety

The accepted-state update remains after both the local/global objective
consistency check and the existing objective acceptance check. The candidate
prediction, residual, and objective are prepared separately; a rejection
returns or continues before the accepted cache and parameters are updated. A
successful sweep still performs a full replay audit before any test-configured
objective refresh. The P4 production-path reconciliation has not been applied.

The D1/E1 retry test now records the accepted state after the first sweep,
checks the retry's `objective_before` and candidate's `accepted_before` against
independent Full Replay, and compares the recovered final parameters and
objective. Both modes recovered after one retry with equal final parameters.
The exhausted-recovery test verifies that both modes retain the last complete
accepted parameters and a Full Replay-consistent objective. The existing
unchanged-candidate test still exercises the no-replay fast path.

The E1 memory tests exercise the incremental strategy under a budget that
requires the memory-bounded Full Replay candidate fallback. The fallback
records only Full Replay candidates, and each accepted-before and candidate
objective matches an independent Full Replay. A one-byte preflight budget
returns `resource-limit` before an accepted state or candidate is created.
Nonfinite observations are rejected before accepted-state creation. Endpoint
certification's resource-limited assessment behavior remains covered by the
existing certification test. A separate internal-boundary test pairs finite
`1e200` observations with the ordinary fixture scale, producing an infinite
initial objective; D1 and E1 both stop at `block-cache-replay-failed`, perform
no candidate replay or acceptance, and retain the initial parameters. This is
an adversarial inconsistent-context check, not a valid constructed
`JointProblem` workload.

No natural outer candidate rejection was found in 512 deterministic
observation/scale/initial-width perturbations and 196 signed-cancellation
factor/width perturbations. The pinned standard and hard fixtures also
produced zero candidate decision disagreements and zero rejected changed
candidates. No test injects a false rejection. This is missing evidence for
the requested objective-increase and local/global-consistency rejection
transactions; the recovery retry above is a local profile recovery and is not
counted as an outer candidate rejection. The near-boundary candidate evidence
currently consists of the qualified E1 candidate audits and their full
objective comparison, not a rejected boundary case. Consequently P3 is
partially qualified and is not a Promotion pass. No KKT, width, or convergence
threshold was changed.

### P4 — Promotion decision and audited performance

P3 did not complete its required natural outer-rejection evidence. Production
therefore remains D1 (exact affected-row prediction with a full objective
scan); no E1 production candidate update or runtime selector was promoted.
This is a Partial Qualification decision, not a numerical regression finding.

For a full-cost E1 measurement, the closed-loop benchmark now enables the
incremental objective together with the existing full Replay at every sweep
and objective-cache reconciliation after the replay audit passes. The
reconciliation timer includes cache prediction/objective comparison and the
objective refresh. Each mode ran one warm-up plus five measured Release runs
serially, using practical criteria and the original bounds. All benchmark runs
converged, passed Endpoint Certification and Independent Replay, and had zero
candidate rejection or incremental-vs-Full-Replay decision disagreement.

| Targets | D1 search median (s) | E1 audited search median (s) | Search change | D1 complete wall (s) | E1 complete wall (s) | Wall change | D1/E1 process peak RSS (MiB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 0.004358 | 0.004443 | 2.0% slower | 0.005467 | 0.005596 | 2.4% slower | 22.57 / 22.57 |
| 30 | 0.097897 | 0.094575 | 3.4% faster | 0.118995 | 0.115968 | 2.5% faster | 34.05 / 34.16 |
| 300 | 1.293647 | 1.008356 | 22.1% faster | 1.519822 | 1.232017 | 18.9% faster | 90.01 / 89.99 |
| 1,000 | 6.434934 | 3.381034 | 47.5% faster | 7.187297 | 4.133619 | 42.5% faster | 244.86 / 245.05 |

The 300- and 1,000-target search speedups exceed the requested 15% guide. The
small fixture's wall time regressed by 0.13 ms; the 30-target wall time
improved by 3.0 ms. E1's per-run median objective reconciliation costs were
0.000422, 0.005066, 0.050565, and 0.169161 seconds for 3, 30, 300, and 1,000
targets. Endpoint Certification medians were 0.000219, 0.002256, 0.020984,
and 0.070278 seconds; Independent Replay medians were 0.000113, 0.001236,
0.010803, and 0.036284 seconds. E1 objective-cache replay errors were at most
`7.96e-24` and all memory preflights passed. The measured RSS is process-wide
maximum resident set size across the warm-up and five repetitions.

For each size, the D1/E1 independent objectives and final A/C/Width estimates
were equal at the precision recorded in P2. Their known-truth RMSEs were also
identical to the P2 values. These synthetic workloads do not represent the
5,938,087-observation production dataset; that full workload remains
unmeasured. The E1 timings include closed-loop search, every-sweep Full Replay
audit and reconciliation, Endpoint Certification, and Independent Replay.

The Release build completed and the full serial CTest suite passed 19/19. The
production library was built without test-only instrumentation; the audited E1
measurements use the test benchmark policy and do not change the production
default. No third-party dependency, optimization threshold, or convergence
criterion changed. The repository is not pushed.

### M0 — Promotion baseline at `94ac5f3`

The current `develop` source was measured with the Release, test-instrumented
closed-loop benchmark before changing the production path. Each D1/E1 mode ran
one warm-up and five timed runs, serially. Both modes completed endpoint
certification and independent replay on all four standard chain fixtures;
each stopped with `block-stationary`. D1 used a full objective scan for each
candidate. E1 used exact affected-row prediction, incremental objective
updates, and a full replay audit and objective refresh after every sweep.

| Targets | D1 / E1 independent objective | Mean A / C / width (D1 and E1) | Global A/C KKT | Projected width KKT | Sweeps | D1 / E1 search median (s) | D1 / E1 wall median (s) | Preflight peak (bytes) | D1 / E1 process peak RSS (bytes) |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | `2.0199078132865188e-10` / same | 1.9999831 / 0.20000604 / 0.49999942 | `3.330669074e-16` | `3.885780586e-15` | 4 | 0.004412251 / 0.004281458 | 0.005630500 / 0.005387041 | 447,445 | 23,724,032 / 23,756,800 |
| 30 | `1.5170760306846197e-10` / same | 1.9999989 / 0.20000029 / 0.49999995 | `1.602388888e-9` | `1.880662293e-11` | 5 | 0.097959709 / 0.095761459 | 0.119049083 / 0.117569250 | 3,930,461 | 35,749,888 / 35,651,584 |
| 300 | `1.4790433590012404e-10` / same | 1.9999998 / 0.20000006 / 0.49999999 | `5.067264247e-10` | `1.878275313e-12` | 5 | 1.290707000 / 1.009375957 | 1.515713083 / 1.235971167 | 36,716,466 | 100,040,704 / 94,289,920 |
| 1,000 | `1.4790435904529177e-10` / same | 2.0000000 / 0.20000002 / 0.50000000 | `2.775278306e-10` | `5.634381850e-13` | 5 | 6.433934333 / 3.376282041 | 7.186736791 / 4.130182332 | 121,648,896 | 256,966,656 / 256,966,656 |

The final A/C/width distributions matched between D1 and E1 for these fixtures;
the table records their means. The previous E1 `objective_reconciliation`
values (0.000408, 0.005168, 0.050663, and 0.169116 seconds) included the
overlapping global-state evaluation (0.000405, 0.005039, 0.050317, and
0.168025 seconds). They are not objective-refresh costs. M0 now times the
full-replay consistency audit separately and starts the reconciliation timer
only after the audit passes and after global-state evaluation has completed.

The acceptance inequalities also constrain the positive-delta branch: if the
local objective is nonincreasing and the local/global discrepancy is enclosed
by `epsilon`, then `global_delta <= local_after - local_before + epsilon <=
epsilon`. Thus its separate rejection condition is unreachable for finite
values satisfying those checks; it remains in place to preserve the explicit
acceptance contract and protect floating-point boundary behavior. The
consistency failure retains its existing terminal `block-objective-replay-failed`
semantics. The audit also found that `WithinBlockObjectiveReplay(Inf, Inf)` is
true under IEEE comparisons, so a nonfinite candidate objective cannot rely
on the enclosure check alone. M1 will explicitly qualify candidate parameters,
prediction, residual, and objective for finiteness before commit.

A post-fix 1,000-target timing check separated the E1 median sweep replay
(0.085563 s), replay audit (0.000881 s), global-state evaluation (0.168027 s),
and objective-only reconciliation (0.000000125 s). This is a timer-boundary
check, not a replacement for the final Production Release benchmark.

### M1 — E1-Direct Production promotion

`develop` now uses E1-Direct in the ordinary Production candidate path. The
exact affected-row replay still recomputes prediction and residual from all
contributors; the candidate objective uses the accepted cache plus the
affected residual-square delta. D1's affected-row full objective remains only
in the instrumented test reference. If the existing scratch preflight rejects
the affected-row workspace, the candidate uses Full Replay and its complete
objective instead.

Candidate parameters and replay values are prepared and checked before any
accepted-state write. Accepted prediction and residual values update only on
affected rows; parameters swap into place after all acceptance checks. A
nonfinite candidate parameter, prediction, residual, or objective now follows
the terminal `block-objective-replay-failed` path without updating the
accepted state. The replay enclosure helper also rejects nonfinite values,
including the previously problematic `Inf <= Inf` comparison.

Every successful sweep performs a Full Replay prediction, residual, and
objective audit. Only a passing audit reaches objective-scalar reconciliation;
the scalar is set from the already computed replay objective. Production work
telemetry reports candidate replay counts, sweep replay and audit times,
global-state evaluation time, reconciliation count/time, and maximum replay
errors. Endpoint Certification and Independent Replay remain in the existing
completion path.

Validation at this stage: the Production Release CLI and instrumented test
targets both compiled; all 39 `JointFixedNeighborBlockCoordinateTest` cases
passed. The group includes D1/E1 reference comparisons, qualified standard and
hard fixtures, width bounds, recovery, unchanged candidates, nonfinite initial
state, and memory fallback. The signed-cancellation case retains its shared
F0/D1 legacy width-KKT limitation; E1 introduces no new failure there. No
numerical acceptance, KKT, width, or convergence threshold changed.

### M2 — Candidate transaction and failure-path tests

A test-only candidate fault callback now injects faults after local candidate
calculation and before acceptance. It is compiled only with
`RHBM_GEM_TEST_INSTRUMENTATION`; it does not select a Production strategy or
alter local solvers or acceptance limits. The added cases cover a local
objective increase, objective-consistency failure, NaN and infinite objective,
nonfinite candidate parameters, and nonfinite affected prediction/residual.
For failures after a prior accepted update, the tests snapshot the exact
accepted-before state for the failing candidate and verify parameters,
prediction, residual, and objective remain bitwise unchanged. The local
objective-increase test captures state before the sweep's separately required
Objective Reconciliation, since that later audit is allowed to refresh the
objective scalar.

The consistency and nonfinite faults run against both D1 Full Scan and E1
Incremental Direct. They terminate with the existing
`block-objective-replay-failed` reason and do not masquerade as an ordinary
candidate rejection. The positive-delta test checks the enclosure boundary
and a value immediately outside it; it supports the acceptance-control-flow
analysis, not a claim about every possible floating-point input. These
fault-injection tests are Transaction Safety and branch-coverage evidence,
not Numerical Non-regression evidence and not proof that natural inputs can
never produce different decisions.

After these additions, all 42 `JointFixedNeighborBlockCoordinateTest` cases
passed. The fault-injection cases passed with their expected state-preservation
and failure semantics. The earlier P3 “missing evidence” classification for
natural candidate rejections is superseded for branch-coverage purposes; no
natural-rejection fixture is claimed or required.

### M3 — Retired replay experiment cleanup

Production FixedNeighbor now has one candidate objective path: exact
affected-row prediction and residual recomputation followed by the direct
objective cache delta. The only D1 full-objective candidate calculation is a
test-only reference flag. It is absent from the Production build. Scratch
preflight still selects Full Replay when affected-row replay cannot fit; a
fallback candidate uses the full prediction, residual, and objective result.

Removed the IncrementalGuarded path, empirical guard-band and refresh-interval
controls, cache-trust state, candidate full-scan fallback counters, and
periodic duplicate Full Replay checks. Sweep replay audit and Objective
Reconciliation remain mandatory: the already computed Full Replay is checked
before the cached scalar is refreshed. The 3.2 GB planning budget and all
acceptance, KKT, width, and convergence thresholds are unchanged.

Removed `JointReplayExperiment.{hpp,cpp}`, which implemented Delta Prediction
B1/B2/C, external CSR/incidence replay variants, compensated or scaled
incremental sums, and tiled/tree objective caches. Removed the replay-only and
drift benchmark executables and sources. Replaced the strategy-comparison
benchmark with a Release benchmark linked to the uninstrumented `rhbm_gem`
Production library; it performs one warm-up and five measured closed-loop
runs and reports candidate, sweep replay, audit, reconciliation, global
evaluation, endpoint, Independent Replay, memory preflight, and process RSS.
The separate Joint Optimization benchmark remains because it measures other
solver work and did not depend on the retired replay framework.

The broad experiment shadow tests and guard/refresh-interval matrices were
replaced by a focused Production E1 vs D1 Full Objective Reference test. It
checks every E1 affected-row prediction/residual against independent Full
Replay, checks the incremental objective and acceptance decision, verifies
the final D1/E1 state on the standard fixture, and confirms every successful
sweep reconciles after audit. Active-width, multiple initial-width,
high-overlap, weak-width, rank-deficient, near-singular, low-signal,
signed-cancellation, recovery, unchanged, fault-injection, and memory-fallback
coverage remains. The focused FixedNeighbor group now has 36 tests and all
36 passed after cleanup. This is numerical non-regression and branch-safety
evidence for the covered fixtures; it does not claim a universal floating
point error bound.

### M4 — Production Release qualification and measured performance

The `build/release` tree is an explicit Release build with
`BUILD_TESTING=OFF` and `RHBM_GEM_TEST_INSTRUMENTATION` absent from the
`FixedNeighborBlockCoordinate.cpp` compile command. The benchmark executable
uses the same ordinary `rhbm_gem` Production library without test
instrumentation. Its source calls `SearchFixedNeighborComponent` directly,
requires a converged and qualified endpoint plus passing Independent Replay,
and reports one warm-up followed by five measured runs. Objective values are
reported on the normalized scale (`independent replay objective / scale²`).

For a fair D1 reference, a Release copy of revision `b6a69d85` was built
outside the repository with only the candidate Objective calculation changed
back to the Full Objective Scan. It used the same current benchmark driver,
full sweep audit, reconciliation, endpoint certification, and Independent
Replay. Both executables were uninstrumented, used the same synthetic
workloads, and had the same thread limits. This comparison measures the
candidate Objective calculation difference; it is not the historical
test-only strategy selector.

| Targets | D1 search median [range] | E1 search median [range] | Search change | D1 wall median [range] | E1 wall median [range] | Wall change | D1 / E1 peak RSS | Memory preflight |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 4.365 [4.358, 4.394] ms | 4.283 [4.277, 4.305] ms | 1.9% faster | 5.259 [5.252, 5.338] ms | 5.176 [5.167, 5.244] ms | 1.6% faster | 24.23 / 24.20 MB | 447,445 B, pass |
| 30 | 97.526 [97.505, 97.606] ms | 94.213 [94.201, 94.471] ms | 3.4% faster | 116.382 [116.335, 116.480] ms | 113.071 [113.047, 113.374] ms | 2.8% faster | 35.42 / 35.52 MB | 3,930,461 B, pass |
| 300 | 1.287 [1.287, 1.295] s | 1.007 [1.007, 1.010] s | 21.8% faster | 1.490 [1.489, 1.498] s | 1.210 [1.209, 1.214] s | 18.8% faster | 88.31 / 88.56 MB | 36,716,466 B, pass |
| 1,000 | 6.418 [6.416, 6.419] s | 3.369 [3.365, 3.373] s | 47.5% faster | 7.099 [7.096, 7.105] s | 4.050 [4.047, 4.058] s | 42.9% faster | 233.83 / 233.65 MB | 121,648,896 B, pass |

The following are E1 median component times in seconds. Sweep Full Replay,
replay audit, global-state evaluation, and Objective Reconciliation have
separate non-overlapping timers. Endpoint certification excludes its separately
reported Independent Replay.

| Targets | Candidate evaluation | Sweep Full Replay | Replay audit | Global state | Objective reconciliation | Endpoint certification | Independent Replay |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0.000454 | 0.000206 | 0.000002 | 0.000405 | 0.000000083 | 0.000103 | 0.000108 |
| 30 | 0.007026 | 0.002561 | 0.000025 | 0.005037 | 0.000000042 | 0.001021 | 0.001071 |
| 300 | 0.072423 | 0.025655 | 0.000253 | 0.050424 | 0.000000042 | 0.010188 | 0.010716 |
| 1,000 | 0.241705 | 0.085553 | 0.000885 | 0.168241 | 0.000000124 | 0.033997 | 0.035934 |

Every measured D1 and E1 endpoint was qualified, passed Independent Replay,
and stopped with `block-stationary`. The 3/30/300/1,000-target normalized
objectives were respectively `2.01990781e-10`, `1.51707603e-10`,
`1.47904336e-10`, and `1.47904359e-10`; the final means and KKT values matched
between D1 and E1 to the reported precision. E1 used 4/5/5/5 sweeps and
4/15/130/420 affected candidate evaluations. No benchmark candidate used the
memory fallback; the fallback and over-budget preflight are covered by the
FixedNeighbor tests. D1/E1 preflight estimates were identical. RSS varies by
process allocation and remains comparable at all four sizes.

The wider numerical suite retains the original D1/E1 comparisons: standard,
active-width, high-overlap, weak-width, and multiple initial-width fixtures
remain qualified without a meaningful estimate or objective regression.
Rank-deficient and near-singular fixtures retain the same recovery-exhausted
baseline result; the low-signal fixture retains the same sweep-budget stop.
Signed cancellation retains the known F0/D1 width-KKT limitation, with no new
E1 failure. The fault-injection cases are transaction and branch-coverage
evidence, not a claim that all possible natural floating-point inputs have
identical decisions.

A test-only sweep replay perturbation now exercises the cache audit failure
before reconciliation. The injected mismatch produces
`block-cache-replay-failed`, leaves the reconciliation count at zero, leaves
`objective_cache_refreshed` false, and preserves eta, beta, prediction,
residual, and objective bitwise at the sweep boundary. This specifically tests
the check-before-refresh ordering; it does not simulate a natural E1 numerical
drift.

The repository contains no 37,406-target / 5,938,087-observation production
input, so runtime and peak RSS at that actual size remain unmeasured. The
1,000-target result is synthetic and is not extrapolated to that workload.
Release build, all 19 serial CTest groups, the focused 37-case FixedNeighbor
suite, transaction fault tests, endpoint certification, Independent Replay,
and memory-preflight coverage passed. M4 did not change numerical thresholds,
add dependencies, update Notion, or push commits.

### Final Production decision

E1-Direct is the Production candidate Objective path. Normal candidates use
the accepted Objective cache plus exact affected-row residual-square deltas;
the D1 Full Objective Scan exists only as an instrumented test reference.
Scratch-memory preflight still selects the full trusted replay fallback before
affected-row allocations. Each successful sweep Full Replays prediction,
residual, and objective; the cache audit must pass before the objective scalar
is reconciled. Endpoint Certification and Independent Replay remain mandatory.
The measured Release build has a small 3-target improvement, a modest
30-target improvement, and stable 21.8% / 47.5% search-time reductions at 300
and 1,000 targets. There is no E1-specific numerical or transaction-safety
blocker in the covered fixtures. No universal floating-point error bound is
claimed.
