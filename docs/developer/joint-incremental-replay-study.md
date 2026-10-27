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

## Roadmap status

| Stage | Status | Result |
| --- | --- | --- |
| N0 | Complete | Baseline telemetry and deterministic state/candidate digests recorded. |
| N1 | Complete | A/B1/B2/C/D1/E measured at 3, 30, and 300 synthetic targets; D1 passed bitwise candidate replay checks. |
| N2 | Complete | Five-method shadow; Full Replay retained every decision; D1 remained bitwise exact, while B/C/E accumulated drift. |
| N3 | Pending | |
| N4 | Pending | |
| N5 | Pending | |
| N6 | Pending | |
