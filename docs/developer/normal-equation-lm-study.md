# Bounded LM Normal-Equation Study

## Decision

**Do not promote the normal-equation solver. Keep augmented QR as the production solver.**

The normal-equation route passed the tested KKT, endpoint, and replay checks on the requested real data. It did not meet the performance requirement: over five independent full runs, its LM subproblem was 6.72 times slower at the median, while complete Joint wall time was 0.57% slower. A tall, narrow microbenchmark was faster with normal equations, but that workload did not represent the actual component solver calls.

This report records the temporary solver implementation and tests. The experiment-only implementation, selector, instrumentation, and benchmark changes were removed after the decision. Production code remains on the original QR path.

## Setup

- Repository branch: `develop`; starting HEAD: `3d156ccfc4551b8986d9da11bc2c007784553711`.
- Real input: `/Users/yslian/data/6z6u.cif`, map `/Users/yslian/data/emd_11103_unsharpened.map`, selected chains `A` and `YA`.
- Real workflow selected 1,557 targets, 2,166 contributors, 243,294 observations, and one component.
- Hardware: Apple M1 (T8103), macOS ARM64; Apple Clang 21; Release build.
- CMake: `BUILD_TESTING=ON`, Python bindings off, experimental feature on, UMAP off, system dependencies. OpenMP and numerical-library thread limits were all set to 4. Build commands used at most `-j4`; benchmark workflows used one job.
- No third-party dependency was added. No Notion content was changed.

The N0 audit confirmed the current convention: the residual step solves `J·δ ≈ -r`; damping adds `λD²`; upper-bound multiplier signs and trust-region acceptance remain those of the existing solver. The experiment did not change the damping search, trust radius, acceptance ratio, A/C profile solver, or endpoint checks.

## Roadmap results

| Stage | Result | Evidence |
|---|---|---|
| N0 | Complete | Existing contract, workflow, and numerical CTest groups passed before the experiment. QR timings and numerical baselines were collected. |
| N1 | Complete for evaluation | A test-only normal-equation route formed `JᵀJ` and `Jᵀr`, used a scaled free-face LLT solve, and kept the QR implementation available as fallback. It was not retained. |
| N2 | Pass | Analytic damped KKT, active-bound activation/release, all-active and tiny-residual cases, rank-deficient damping fallback, ill-conditioning guard, scaled variables, invalid scaling, and tiny trust radius passed. The full three-group CTest invocation passed (3/3). |
| N3 | Fail for promotion | Three five-run synthetic comparisons and five independent QR/hybrid runs on the requested real input showed slower FixedNeighbor and LM timings; complete Joint time had no measurable improvement. Correctness checks passed. |
| N4 | Complete for evaluation | The conditioning-gated LLT route passed all 53,584 real-data attempts with no QR fallbacks. Singular and ill-conditioned unit cases fell back to QR. A microbenchmark still incurred one guard fallback per solve. |
| N5 | Keep QR | Performance target was not met. All experimental solver, test selector, instrumentation, and benchmark code was removed; this report is retained. |

## N0 QR baseline

The baseline benchmark used the existing FixedNeighbor workload with a warmup and five measured runs. Values below are the same test-instrumented build used for the A/B comparison.

| Targets | Search median | FixedNeighbor wall median | LM subproblem median | Sweep / replay |
|---:|---:|---:|---:|---|
| 12 | 32.211 ms | 33.743 ms | 0.195 ms | 4 sweeps; qualified; independent replay passed |
| 48 | 161.723 ms | 167.829 ms | 0.783 ms | 5 sweeps; qualified; independent replay passed |
| 240 | 799.143 ms | 829.159 ms | 3.890 ms | 5 sweeps; qualified; independent replay passed |

The 4096×8 bounded-step microbenchmark took 15.609 ms per QR solve; QR factorization accounted for 11.558 ms. Derivative preparation and reduction on a 13,710-row, 60-column workload took 19.965 ms and 5.822 ms respectively. These numbers are useful for establishing an upper-bound case for the proposed solver, not as a representative Joint-workflow baseline.

## N2 numerical checks

The temporary hybrid implementation normalized the free system by `D`, certified its scaled residual, checked a singular-value rank threshold and an estimated roundoff bound, and then recomputed the original-J-space KKT residual. Any failed check retried the same damping candidate with augmented QR. Damping search and active-set activation/release code were shared with QR.

The unit cases passed the original KKT tolerances. In particular, the rank-deficient zero-damping case proceeded to positive damping, while an ill-conditioned face used the QR fallback and matched the direct QR step. No tolerance was relaxed.

## N3 synthetic A/B

Each solver/magnitude entry below is a median over five measured runs after one warmup. Both modes used the same workload and endpoint qualification.

| Targets | Metric | QR | Normal-equation hybrid | Change |
|---:|---|---:|---:|---:|
| 12 | LM subproblem | 0.195 ms | 1.696 ms | 8.69× slower |
|  | FixedNeighbor wall | 33.743 ms | 35.478 ms | +5.14% |
| 48 | LM subproblem | 0.783 ms | 7.028 ms | 8.98× slower |
|  | FixedNeighbor wall | 167.829 ms | 174.765 ms | +4.13% |
| 240 | LM subproblem | 3.890 ms | 33.956 ms | 8.73× slower |
|  | FixedNeighbor wall | 829.159 ms | 869.115 ms | +4.82% |

All six synthetic endpoints were qualified and passed independent replay. Per-size objective, global A/C KKT, projected-width KKT, and sweep counts matched to the reported precision:

| Targets | Objective | Global A/C KKT | Projected-width KKT | Sweeps |
|---:|---:|---:|---:|---:|
| 12 | `1.54157728e-10` | QR `1.94e-16`; hybrid `1.67e-16` | `9.99e-16` | 4 |
| 48 | `1.51441738e-10` | `1.2669e-9` | `1.1746e-11` | 5 |
| 240 | `1.48393658e-10` | `5.6655e-10` | `2.3480e-12` | 5 |

The small workflow differences are adverse in every size. The normal route had no QR fallbacks in these closed-loop runs. Peak RSS differences were small: +32 KiB, 0, and +192 KiB for 12, 48, and 240 targets respectively. Memory preflight passed for every run.

The separate 4096×8 microbenchmark favored hybrid: 2.429 ms versus 15.609 ms, or 6.43× faster. It averaged 59 damping candidates per solve, accepted 58 normal-equation candidates, and fell back once due to the conditioning guard. This narrow, tall matrix is not representative of the production component's repeated active-set solves.

## N3 requested real-data A/B

Five independent workflow invocations per mode were run in alternating QR/hybrid order. Each invocation opened the same CIF/map, selected chains A and YA, used one worker, saved to a fresh temporary SQLite database, and removed that database/output directory on completion.

| Metric | QR median (range) | Hybrid median (range) | Change |
|---|---:|---:|---:|
| Complete Joint command wall | 144.244 s (143.711–145.885) | 145.066 s (142.492–145.624) | +0.57% |
| LM subproblem | 0.154832 s (0.150342–0.160937) | 1.041209 s (0.987727–1.051925) | +572.48% / 6.72× slower |
| Peak RSS median | 1,627,242,496 B (1,624,522,752–1,963,573,248) | 1,624,604,672 B (1,560,461,312–1,775,386,624) | −0.16% by median |

The RSS ranges are wide and overlap; the small median difference is not evidence of a memory reduction. No stable peak-RSS increase was observed in these runs, but the memory measurement alone does not change the performance decision.

The real-data endpoint was consistent across modes:

- 1,557 targets, 2,166 contributors, 243,294 observations, one component, 117 sweeps, and 35,525 profile evaluations.
- Both modes were qualified; A/C and target runtime KKT statuses passed; independent KKT consistency and independent replay passed.
- Objective was `0.0219534820633` in both; absolute difference was `2.78e-17`.
- Maximum absolute difference across A/C, B, and log-B parameters was `7.60e-6`; width-gradient maximum absolute difference was `1.45e-11`.
- QR accepted 14,557 updates; hybrid accepted 14,613. Both stopped qualified after the same 117 sweeps.
- 151 local-conditioning blocks were assessed in each run. Trajectory stability was observed; perturbation stability was unavailable.
- The saved Joint result does not retain prediction vectors, so prediction RMS difference was unavailable. Independent replay passed for both modes.
- Hybrid attempted and accepted 53,584 normal-equation candidates with zero QR fallbacks. QR performed 53,894 QR factorizations.

The measured normal factorization and solve account for only a small part of the 0.887 s median LM-subproblem increase. The experimental timing did not isolate each singular-value guard call; the remaining time is within the full bounded-solve path. This overhead is a likely explanation for the mismatch between the tall microbenchmark and the real component result, but that attribution was not separately measured.

## N4 guard and fallback observations

The temporary guard rejected rank-deficient/ill-conditioned scaled systems and non-finite or uncertified solves. It preserved the existing QR fallback at the same damping value. The microbenchmark incurred one `normal-equation-conditioning-guard` fallback per solve; the full real workflow incurred none. The singular and ill-conditioned unit tests verified both fallback reasons and QR behavior. No LDLT route was added because the LLT plus QR experiment already failed the production performance gate.

## N5 decision and cleanup

The solver did not satisfy the N5 speed objective. The full real-data Joint workflow did not improve by 5%; it was 0.57% slower at the median, within run-to-run variation. The LM subproblem was substantially slower on the real input and all three representative FixedNeighbor sizes. A tall microbenchmark speedup does not offset those production-shape results.

Therefore, the production solver remains augmented QR. The experimental C++ branch, selector, test-only instrumentation, experimental tests, and benchmark modifications were removed. Only this research report is retained. No solver dependency or runtime flag was added.
