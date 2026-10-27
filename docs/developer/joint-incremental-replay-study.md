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

## Roadmap status

| Stage | Status | Result |
| --- | --- | --- |
| N0 | Complete | Baseline telemetry and deterministic state/candidate digests recorded. |
| N1 | Complete | A/B1/B2/C/D1/E measured at 3, 30, and 300 synthetic targets; D1 passed bitwise candidate replay checks. |
| N2 | Complete | Five-method shadow; Full Replay retained every decision; D1 remained bitwise exact, while B/C/E accumulated drift. |
| N3 | Complete | 503-attempt refresh/drift stress at 1/10/50/100/sweep; B methods drifted more with longer intervals, C/E retained objective error, D1 stayed exact; all scripted rejections rolled back bitwise. |
| N4 | Pending | |
| N5 | Pending | |
| N6 | Pending | |
