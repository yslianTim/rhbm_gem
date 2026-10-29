# Joint Component Projected Width KKT Promotion Study

## Decision

**Partial Qualification.** Production remains `practical-v2`: global A/C projected KKT `1e-7`, bound-aware projected Width KKT `1e-9`. The T4 candidate (`1e-7` / `1e-5`) passed numerical qualification on the tested successful synthetic cases and stayed within the stated two-decimal truth error. It saved one sweep only on the 30-target workload. The result does not justify promotion: the 30-target Parameter Stability classification changed from `Stable` to `Unstable`, stability coverage did not include weak-width or active-width blocks, and the required continuation, multi-start, scale/noise, and real-map evidence is incomplete.

This is an offline study record. It does not change the production convergence policy, policy identifier, persistence schema, or runtime controls.

## Scope and reproducibility

- Initial branch: `develop`; initial HEAD and `origin/develop`: `6c6f6a0e27cf5f3cd2d08b93377e57b13f96f9a6`; initial working tree clean.
- Production criteria were verified as A/C `1e-7`, Width `1e-9`; provenance remains `practical-v2`.
- Release build: AppleClang 21.0.0, arm64, `Release` (`-O3 -DNDEBUG`), CMake/Ninja. Build parallelism was at most 4.
- Each benchmark policy used the same deterministic `chain` synthetic input, initial values, partition, upper bound, solver budgets, and one-thread numerical-library settings (`OMP_NUM_THREADS=1`, BLAS/MKL/VECLIB threads 1), with one warmup and five measured runs. Medians and min–max ranges are reported. The CPU model was not re-read for this run.
- The initial Release CTest run passed 19/19 when run with `OMP_NUM_THREADS=2` and CTest parallelism 2. `OMP_NUM_THREADS=1` hit an existing Eigen thread-count expectation in an unrelated test; the final verification used the passing setting.
- Threshold injection, sweep trace, and CSV export were temporary test-only research changes. They were removed after measurement. No third-party dependency, Notion document, Production instrumentation, branch, or worktree was added.
- Synthetic benchmark results are not real cryo-EM validation. No usable real map fixture was present in the repository; real-map workload is **Unavailable**.

## Production dependency map

| Path | Role | Tolerance impact |
|---|---|---|
| `PracticalJointConvergenceCriteria` / `StrictJointConvergenceCriteria` | Defines practical `{1e-7, 1e-9}` and strict `{1e-10, 1e-12}` KKT pairs | Source of production policy values |
| `CalculateWidthBoundKkt` and `ClassifyFixedNeighborStationarity` | Computes bound-aware Width stationarity and combines it with A/C KKT, feasibility, completed sweeps, and block-failure vetoes | Determines whether a search endpoint may be attempted |
| `CertifyFixedNeighborEndpoint` → `VerifyJointEndpoint` | Independently qualifies the endpoint | Certification must receive the same policy as search |
| Independent Replay and assembled evidence (`AssessmentEvidence`) | Replays predictions/objective/residual and reports global qualification | Must use the candidate criteria consistently; all benchmark candidate runs did |
| `StabilityProbePriority`, `StabilitySensitivityScore` | Scores local blocks relative to the practical KKT scales | Can change representative-block order/selection if policy scales change |
| `ParameterStability.cpp` | Qualifies baseline and perturbed profiles and classifies sensitivity | Its internal practical-width gate also affects whether a probe finishes |
| `Components.cpp` / `JointResultJson.cpp` | Writes and validates candidate policy provenance | `practical-v2` remains correct because Production stays at `1e-9` |

The study preserved initialization, A/C solver, LM, partition/order, candidate acceptance, objective monotonicity, rollback/recovery, width calibration/bound, replay, halo reconstruction, peeling, and group fitting.

## W0 — Baseline and bottleneck

All baseline endpoints passed Endpoint Certification and Independent Replay, had one certification attempt, and had no failed blocks. Peak RSS is the benchmark process high-water mark.

| Targets | Wall median (range) | Search median | Sweeps / accepted blocks | Final A/C KKT | Final projected Width KKT | Objective | Peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 3.697 ms (3.663–3.812) | 2.835 ms | 3 / 3 | `2.78e-16` | `6.71e-10` | `2.01991e-10` | 24.0 MB |
| 30 | 105.888 ms (105.471–107.545) | 87.183 ms | 5 / 15 | `1.60e-9` | `1.88e-11` | `1.51708e-10` | 36.8 MB |
| 300 | 962.872 ms (956.380–965.638) | 762.076 ms | 4 / 104 | `2.99e-8` | `1.11e-10` | `1.47923e-10` | 89.4 MB |
| 1,000 | 3.202 s (3.185–3.366) | 2.523 s | 4 / 336 | `1.64e-8` | `3.33e-11` | `1.47923e-10` | 237.0 MB |

### Sweep evidence

| Targets | Sweep | A/C projected KKT | Projected Width KKT | Interpretation |
|---:|---:|---:|---:|---|
| 3 | 1 / 2 / 3 | `1.67e-16` / `1.67e-16` / `2.78e-16` | `3.62e-3` / `2.27e-5` / `6.71e-10` | T4 still needs sweep 3; no sweep saved |
| 30 | 1 / 2 / 3 / 4 / 5 | `8.40e-3` / `2.53e-4` / `5.61e-6` / `9.46e-8` / `1.60e-9` | `5.80e-4` / `3.38e-6` / `7.18e-8` / `1.112e-9` / `1.88e-11` | At sweep 4 A/C is eligible but T0 Width is 11.2% above `1e-9`; this is the only observed B case and T4 saves one eligible sweep |
| 300 | 1 / 2 / 3 / 4 | `2.66e-3` / `8.00e-5` / `1.77e-6` / `2.99e-8` | `6.20e-5` / `3.78e-7` / `7.17e-9` / `1.11e-10` | At sweep 3 Width meets T4 but A/C does not; A/C is the bottleneck |
| 1,000 | 1 / 2 / 3 / 4 | `1.46e-3` / `4.38e-5` / `9.72e-7` / `1.64e-8` | `1.86e-5` / `1.13e-7` / `2.15e-9` / `3.33e-11` | At sweep 3 Width meets T4 but A/C does not; A/C is the bottleneck |

At 30 targets the sweep-4 confirmation was available, four full sweeps had completed, no block failed, and T4 could proceed to Certification. It did so with one attempt, which passed; Independent Replay also passed. T0 then performed one more sweep. There was no certification retry increase. The observed tail sweep reduced Width KKT from `1.112e-9` to `1.881e-11`, but its parameter and prediction changes remained very small (see W2).

## W1 — Five-policy controlled benchmark

Every listed endpoint was qualified under its own policy and passed Independent Replay; each used one Certification attempt and zero failed blocks. The baseline-to-candidate speedup is the ratio of wall medians. RSS values are rounded to MB.

| Targets | Width threshold | Sweeps | Wall median (range) | Speedup | Peak RSS |
|---:|---:|---:|---:|---:|---:|
| 3 | `1e-9` | 3 | 3.697 ms (3.663–3.812) | 1.000× | 24.0 MB |
| 3 | `1e-8` | 3 | 3.872 ms (3.676–3.892) | 0.955× | 24.1 MB |
| 3 | `1e-7` | 3 | 3.892 ms (3.788–6.035) | 0.950× | 24.1 MB |
| 3 | `1e-6` | 3 | 3.840 ms (3.650–3.917) | 0.963× | 24.1 MB |
| 3 | `1e-5` | 3 | 4.012 ms (3.711–4.287) | 0.922× | 24.2 MB |
| 30 | `1e-9` | 5 | 105.888 ms (105.471–107.545) | 1.000× | 36.8 MB |
| 30 | `1e-8` | 4 | 90.765 ms (90.491–91.296) | 1.167× | 36.9 MB |
| 30 | `1e-7` | 4 | 88.103 ms (87.995–89.205) | 1.202× | 36.9 MB |
| 30 | `1e-6` | 4 | 88.661 ms (88.185–89.380) | 1.194× | 38.1 MB |
| 30 | `1e-5` | 4 | 88.807 ms (88.452–91.918) | 1.192× | 36.9 MB |
| 300 | `1e-9` | 4 | 962.872 ms (956.380–965.638) | 1.000× | 89.4 MB |
| 300 | `1e-8` | 4 | 959.923 ms (944.718–963.116) | 1.003× | 88.9 MB |
| 300 | `1e-7` | 4 | 960.854 ms (944.373–968.138) | 1.002× | 88.8 MB |
| 300 | `1e-6` | 4 | 947.651 ms (944.297–970.771) | 1.016× | 88.1 MB |
| 300 | `1e-5` | 4 | 955.958 ms (945.062–965.085) | 1.007× | 88.1 MB |
| 1,000 | `1e-9` | 4 | 3.202 s (3.185–3.366) | 1.000× | 237.0 MB |
| 1,000 | `1e-8` | 4 | 3.205 s (3.195–3.318) | 0.999× | 238.8 MB |
| 1,000 | `1e-7` | 4 | 3.201 s (3.199–3.240) | 1.000× | 234.2 MB |
| 1,000 | `1e-6` | 4 | 3.229 s (3.223–3.314) | 0.992× | 243.3 MB |
| 1,000 | `1e-5` | 4 | 3.253 s (3.222–3.373) | 0.984× | 236.6 MB |

Only 30 targets saved a sweep. T4 wall time there fell 16.1% (1.192× speedup); search time fell 19.5% (87.183 ms to 70.169 ms). At 3, 300, and 1,000 targets there was no sweep reduction and no repeatable wall-time improvement. The small median movements at 300 and 1,000 overlap the measured ranges. Memory stayed in the same range; there is no evidence of a material increase.

The measured benchmark includes search and endpoint certification. Parameter Stability is separately timed below and must be considered for end-to-end interpretation; the 30-target diagnostic added about 5 ms in the full-policy run, consuming part of the ~17 ms search-plus-certification saving. These diagnostic timing samples were single runs and are not a paired five-run end-to-end result.

## W2 — Accuracy and reference audit

The deterministic synthetic ground truth is A=2, C=0.2, B=0.5, with low-amplitude nuisance noise. The application error threshold was 0.005. Errors below are max absolute errors across targets. For 3, 300, and 1,000 targets, T0 and T4 parameter vectors matched exactly in the saved CSVs. At 30 targets, T1–T4 reached the same sweep-4 endpoint; the T4 endpoint is compared to the T0 sweep-5 reference.

| Targets | Policy | Max A error | Max C error | Max B error |
|---:|---|---:|---:|---:|
| 3 | T0 / T4 | `6.5353e-5` | `1.4822e-5` | `3.3070e-6` |
| 30 | T0 | `6.7098e-5` | `1.4793e-5` | `3.6832e-6` |
| 30 | T1–T4 | `6.7098e-5` | `1.4787e-5` | `3.6832e-6` |
| 300 | T0 / T4 | `6.5510e-5` | `1.5113e-5` | `3.2833e-6` |
| 1,000 | T0 / T4 | `6.5510e-5` | `1.5113e-5` | `3.2833e-6` |

At 30 targets, the T0 vs T4 truth-error RMS / P95 / P99 values were:

| Parameter | T0 RMS / P95 / P99 | T4 RMS / P95 / P99 | Max candidate-reference difference | RMS candidate-reference difference |
|---|---:|---:|---:|---:|
| A | `4.4464e-5 / 6.2970e-5 / 6.7098e-5` | `4.4506e-5 / 6.2970e-5 / 6.7098e-5` | `3.7998e-6` | `9.6350e-7` |
| C | `1.0451e-5 / 1.4714e-5 / 1.4793e-5` | `1.0454e-5 / 1.4714e-5 / 1.4787e-5` | `3.1544e-7` | `7.9974e-8` |
| B | `2.2118e-6 / 3.1438e-6 / 3.6832e-6` | `2.2150e-6 / 3.1438e-6 / 3.6832e-6` | `2.6068e-7` | `6.6106e-8` |

At 30 targets, objective changed from `1.51707603e-10` to `1.51721796e-10` (absolute `1.4193e-14`, relative about `9.4e-5`). Across 13,710 prediction rows, max / RMS / P95 / P99 absolute prediction change was `3.2855e-7 / 3.6282e-8 / 1.2596e-7 / 1.5385e-7`. With the same observations, residual changes have the same magnitudes. No large target outlier was found in the measured vectors. The remaining truth-error margin is large relative to 0.005, but this synthetic result does not establish identifiability in difficult workloads.

**Shadow continuation: Not Run.** An earlier probe reused candidate widths as initial values but did not carry over the candidate A/C state; that is not a valid continuation from the candidate endpoint and is excluded as evidence. The available endpoint path in this study did not expose a full-state continuation input. No continuation result is claimed.

Local Conditioning at the selected 30-target blocks was `WellConditioned`; the three selected blocks were full rank. This does not cover weak-width geometry. Weak-width and noise-only hard cases completed numerical qualification and replay but reported Parameter Stability as `NotAssessed`; they are not accuracy evidence.

## W1-A / W1-B — Parameter Stability impact

The comparison used the same three representative block IDs at each tested size. Weak-width coverage and active-width coverage were both zero because the `chain` fixtures contained no such selected blocks.

| Targets / policy | Selected blocks | Local conditioning | Overall stability | Conditioning time | Stability time | Peak memory |
|---|---|---|---|---:|---:|---:|
| 30 T0 | `block-1`, `block-2`, `block-3` | WellConditioned | Stable | 14.60 ms | 245.13 ms | 23,717,495 B |
| 30 T4 solver-only, T0 selection scale | same | WellConditioned | Unavailable | 14.39 ms | 37.57 ms | 23,717,495 B |
| 30 T4 full policy, T4 selection/assessment scale | same | WellConditioned | Unstable | 14.89 ms | 250.29 ms | 23,717,495 B |
| 300 T0 | `block-1`, `block-7`, `block-21` | WellConditioned | Unstable | 153.69 ms | 208.03 ms | 53,378,250 B |
| 300 T4 solver-only, T0 selection scale | same | WellConditioned | Unstable | 154.38 ms | 208.33 ms | 53,378,250 B |
| 300 T4 full policy, T4 selection/assessment scale | same | WellConditioned | Unstable | 182.99 ms | 223.37 ms | 53,378,250 B |

On 30 targets, solver-only T4 could not complete its probes under the unchanged practical-width assessment gate, so `Unavailable` is not a reduced diagnostic cost. Passing T4 criteria through the full-policy assessment completed the probes but changed the result from T0 `Stable` to T4 `Unstable`. On 300 targets, IDs and overall `Unstable` classification did not change. These small samples do not show a block-selection change, but they do show that the final stability result can change and that an unchanged downstream gate can produce `Unavailable`. Do not count solver-only stability time as a valid end-to-end saving.

The evidence is not sufficient to require a permanent decoupling change: selection sets matched in both tested fixtures, and neither fixture measured weak-width or active-width coverage. If a future promotion is pursued, the selection scale and the assessment's qualification scale must be specified and tested separately; this study does not validate a new full-policy stability contract.

## W3 — Robustness matrix

| Workload | Result under T4 | Classification / limitation |
|---|---|---|
| 3 / 30 / 300 / 1,000 target chain | Certification and Independent Replay passed; no failed blocks | Qualified; only 30-target saved one sweep |
| High-overlap cube, 8 targets | 4 sweeps, one certification, replay passed; A/C KKT `6.07e-8`, projected Width KKT `1.30e-9`; max truth errors A/C/B `6.97e-5 / 1.37e-5 / 4.46e-6` | Qualified for this fixture; not a broad overlap or multi-start study |
| Active-width upper bound | Bound-aware projected KKT and feasibility checks passed under T4 | Qualified for the existing upper-bound fixture; ordinary gradient was not substituted for projected KKT |
| Signed cancellation | T4 KKT-only qualification test passed with replay | Qualified for the tested fixture |
| Noise-only | 2 sweeps, Certification and replay passed; A/C KKT `3.53e-9`, Width KKT `8.55e-13` | Qualified but Parameter Stability `NotAssessed` |
| Weak-width direction | 2 sweeps, Certification and replay passed; truth errors effectively zero for this seeded fixture | Qualified but Parameter Stability `NotAssessed`; not evidence for general weak-width accuracy |
| Rank-deficient constant support | T4 did not qualify; recovery exhausted, A/C KKT `0.482812`, no certification/replay pass | Correctly unqualified |
| Near-singular overlapping support | T4 did not qualify; recovery exhausted, A/C KKT `0.689768`, no certification/replay pass | Correctly unqualified |
| Near-boundary separate fixture | Not Run | Unavailable |
| Multi-start, changed noise levels, true-width scales, A/C amplitudes | Not Run | Required evidence gap |
| Real cryo-EM map | Not Run | No repository fixture available |

Certification failure, retry/recovery, independent replay veto, rollback, resource-limited outcomes, and full end-to-end diagnosis were not exhaustively enumerated across the new threshold matrix. Existing safeguards were not modified, and the rank-deficient/near-singular checks above remained vetoed, but the full promotion gate remains incomplete.

## W4 — Gate assessment and cleanup

| Gate | Assessment |
|---|---|
| G1 Numerical validity | **Passed on exercised successful fixtures, incomplete overall.** Candidate search and endpoint used the same threshold; all four chain sizes passed Certification and replay, while rank-deficient and near-singular cases remained unqualified. Full retry/resource/recovery matrix was not completed. |
| G2 Estimation accuracy | **Passed on measured identifiable chain fixtures, incomplete overall.** Truth errors were below `0.005`, and the 30-target reference differences were small. Full-state shadow continuation, multi-start and scale/noise studies are absent; stability changed from Stable to Unstable in one fixture. |
| G3 Efficiency and robustness | **Not established.** Only 30 targets saved one sweep; 3 targets did not improve, and 300/1,000 targets were A/C-limited. No significant memory increase was measured, but the broad-workload benefit and full-policy stability impact are insufficient for promotion. |

The only allowed conclusion supported by the collected evidence is **Partial Qualification**. Production remains Width KKT `1e-9`; policy ID remains `practical-v2`.

Cleanup restored every temporary experiment change in the benchmark, block-coordinate implementation/header, stability implementation/header, and block-coordinate test. Temporary threshold switches, per-sweep tracing, CSV exporters, and stability test hooks are absent. No candidate Production branch or orphaned helper remains. Existing production KKT, projected upper-bound, certification, replay, stability, conditioning, recovery, and persistence tests remain in the source tree. There is one Production convergence policy and no legacy `1e-5` policy in the runtime.

## Final verification

After cleanup, the full Release build passed with `--parallel 4`. Full CTest passed **19/19** in 18.85 seconds with `OMP_NUM_THREADS=2`, BLAS/MKL/VECLIB threads 1, and CTest parallelism 2. This includes the KKT, projected upper-bound, endpoint/replay, stability, conditioning, recovery, JSON persistence, and integration test groups.

The restored Production benchmark used one warmup and five measured runs, with single-thread numerical libraries:

| Targets | Search median (range) | Wall median (range) | Sweeps / accepted blocks | Certification / replay | Peak RSS |
|---:|---:|---:|---:|---|---:|
| 30 | 88.932 ms (86.741–89.844) | 108.255 ms (105.245–108.994) | 5 / 15 | 1 / passed | 35.4 MB |
| 300 | 758.839 ms (746.926–1,253.718) | 958.731 ms (946.558–1,468.253) | 4 / 104 | 1 / passed | 90.7 MB |
| 1,000 | 2.509 s (2.500–2.562) | 3.191 s (3.179–3.236) | 4 / 336 | 1 / passed | 239.0 MB |

The 300-target run had one high timing outlier; its median and final qualification remained consistent with the earlier baseline. The five-policy closed-loop matrix above used temporary test-only injection and the final cleanup removes that injection path. Peak RSS comparisons are process high-water marks and show no consistent candidate increase. Remaining limitations are the missing full-state continuation, multi-start/noise/scale matrix, dedicated near-boundary case, and real cryo-EM fixture.
