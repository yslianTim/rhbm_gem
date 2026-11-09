# Joint LM simplification and validation

## P0 baseline

The initial `develop` HEAD was `22e6de6d10279b2f811c97fe68b9915960bd0c95`, equal to the local `origin/develop` tracking ref. The working tree was clean. The checkout used for this study was the existing `develop` checkout.

### Data flow and numerical boundaries

Let `n` be the number of observation rows in a profile, `m` the number of Width variables in the block, and `p ≤ 2m` the number of free A/C coefficients.

| Stage | Data and work |
| --- | --- |
| `SearchProfile` / `Profile::Get` | Evaluates the sparse `n × 2m` basis and derivative, solves constrained A/C and caches an `n`-entry residual. Each trace row stores an `Endpoint`, with coefficients and Width gradients, without the full residual. |
| `BoundedLmDriver::Initialize` / trial evaluation | Baseline allocates two `n`-entry normalized residual buffers. `Evaluate` writes one buffer and the driver computes its stable norm. Accepted trials swap the buffers. |
| `Profile::Linearize` | Reuses the cached evaluation at the accepted coordinates. `PrepareDerivative` builds sparse free-design and Width derivative data and the `p × m` corrections. `ReduceDerivative` feeds tiles of Jacobian and normalized residual into Householder QR and retains up to `m × m` and `m` response entries. The discarded orthogonal residual tail still contributes to actual reduction. |
| `BoundedLmDriver::Advance` | Pivoted QR of the compact model supplies a scaled gradient and the model-step norm, shared across retries. Each retry requests a feasible step and compares full normalized residual norms for actual reduction. |
| `SolveBoundedTrustRegionStep` | Creates one workspace per call and reuses it across damping candidates. For a free face of size `f`, pivoted QR factors a compact matrix with `f` columns and up to `m + f` rows. The KKT check uses a residual and gradient at the feasible step. Baseline recomputes that residual for its tolerance scale. |
| `Profile::AssessTrial` / `Accept` | Runs Replay on proposed updates and may check large steps for weak Width identifiability. The accepted Width stationarity and objective survive rejected evaluations even when the candidate replaces `cached`. |
| `SearchResult` / BCD | BCD distinguishes Width failure from a trusted initial A/C profile, reconciles objective changes independently, and retries failed blocks. Complete-sweep stationarity and Endpoint Certification qualify the accepted component. |

Keep the independent Replay and certification calculations, the feasible damped QR solve and its rank checks, and the full-residual actual reduction. The production Profile driver uses `ftol=1e-14`, `xtol=gtol=1e-12` and a proposal threshold of `1e-4`. Practical BCD uses `1e-7` for global A/C KKT and `1e-9` for projected Width KKT. Profile's accepted Width check uses that same Width threshold.

### Reproducible inputs

`JointComponentSearchTest.ProfileSearchPreservesAcceptedTrajectory` records the fixed three-target chain, initial widths `.46, .49, .52`, upper bound `.53` and three-update budget. The baseline accepts all four evaluations including the initial point, makes three derivative evaluations and bounded solves, and stops at `accepted-update-budget`. Normalized objective is `2.019907813286478e-10`, Width KKT is `3.4194869158454821e-14`, and Independent Replay and endpoint qualification pass.

| Update | Radius | Damping | Actual decrease | Predicted decrease | Ratio |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | .029053145019303075 | 0 | .99992667830883186 | .98251500091694766 | 1.0177215384758853 |
| 2 | .058106290038606151 | 0 | .99362367168962051 | .99361807323160711 | 1.0000056344164465 |
| 3 | .00050343055423290634 | 0 | 9.6486963085951061e-7 | 9.6486491817342245e-7 | 1.0000048842962361 |

The second trace uses the same chain with widths `.41, .89, 1.3`, upper bound `2.2`, and the default evaluation and update budgets. It has eight evaluations, seven derivatives and bounded solves, six accepted updates, 112 damping QRs and one final rejected trial. The final accepted objective is `2.0199078132865648e-10`, Width KKT is `7.0610184366159956e-14`, and Independent Replay and endpoint qualification pass. Its terminal reason is `native-lm-stop`, with `width_failed=false`. The two traces retained their status semantics from P0 to final, with `-1` becoming `running` and `2` becoming `step-tolerance` in diagnostic JSON.

Use the `JointComponentSearchTest`, `ActiveSetBoundedLMTest`, `JointBcdSearchTest` and `JointComponentRuntimeTest` suites for the fixed retry, tiny-step, boundary, rank, replay and recovery fixtures. `TransientWidthFailureRecoversAndAssemblesQualifiedComponent` requires one recovered failure, qualified certification, and a converged assembled result.

For timing, use `SyntheticJointWorkload("chain", count)` and time a production-library `SearchBcdComponent` call, excluding input construction. Use counts 3, 30, 300 and 1,000, initial `eta=log(.55)`, upper bound `log(2.2)`, 12 core atoms, 12 maximum sweeps, practical convergence and default profile budget. BCD sets one accepted Width update per visit. The benchmark's fixed cap applies only to this workload, which qualified within four sweeps at every size. For standalone Profile timing, use the same chain with 12 targets and an eight-update budget, and time five calls per process. Use one OpenMP thread, a fresh process for each sample and macOS `getrusage(RUSAGE_SELF).ru_maxrss` in bytes.

Baseline Release configuration and validation used `cmake --preset release -DBUILD_TESTING=ON`, `cmake --build build/release --parallel 4` and `OMP_NUM_THREADS=4 ctest --test-dir build/release --output-on-failure --parallel 4`. All 21 CTest entries passed.

## P1 and P2 changes

P1 confines the `ActiveSetStep` copy and reset to the test build in `ProfileLmSearch.cpp`. The production work counters and maximum subproblem KKT still use the original step. `ActiveSetBoundedLM.cpp` takes the ordinary `.norm()` of the residual already formed for the gradient at the same feasible step. A separate 36,000-input Eigen check included rectangular matrices, exact cancellation and scales from `1e-140` to `1e140`. Both ways of forming the residual norm produced identical double values in that check.

P2 makes `BoundedLmDriver.hpp` store the accepted residual norm and keep the candidate norm local to each retry. `ProfileLmSearch.cpp` computes `(cached.residual/context.scale).stableNorm()` without materializing a normalized vector. Division still precedes the norm per coefficient, which matters when the unnormalized norm would overflow or lose subnormal information. The Profile supplies the accepted compact linear model through `Get(eta)` and updates accepted Width stationarity only in `Accept`. Actual reduction continues to use all observation rows.

The change removes two Eigen vector allocations with `16n` bytes of coefficient storage per nonempty driver initialization. The Profile's full evaluation remains available for derivatives and Replay. The normalization comparison found zero bit differences in 60,410 sampled expressions, spanning 1–16,379 entries and scales from `1e-290` to `1e290`. A separate 25-case check disabled Eigen heap allocation during the expression. `Search_test.cpp` adds checks for initial and trial evaluation failures, Profile budget exhaustion and materialized norm parity around Eigen's 4,096-entry internal block boundary, including subnormal and near-overflow inputs.

For P1 and P2, all 129 fixed Profile trace records matched P0 exactly at 17-digit precision, excluding timing and RSS. The comparison included accepted eta, candidate step, radius, damping, reductions, ratio, accepted and rejected decisions, objective, Width and active-set KKT, replay, each damping candidate's rank and working set, and endpoint qualification. The work counts and stop reasons also matched. The existing compact-vs-full model, retry, boundary, tiny-step, and recovery tests passed at both stages.

## P3 Driver QR decision — Reject

The candidate used `Jc.transpose() * (rc / full_residual_norm)` for the scaled gradient and `Jc * step` for the model-step norm. It removed one compact pivoted QR per successful derivative evaluation. Active-set damping still used augmented QR. The five Joint CTest entries passed. In the fixed Profile traces, the largest absolute prediction difference was `4.44e-16` and the largest ratio difference was `1.71e-13`. The two traces kept their accepted states and endpoint qualification. All four BCD workloads kept the same counters, objectives, KKT values and certification.

A grid of 2,400 compact models covered 3 and 12 columns, scales of `1e-140`, `1`, and `1e140`, and four conditioning regimes. The table gives the maximum QR vs direct difference in each regime.

| Compact model | Scaled gradient absolute | Scaled gradient relative | Model-step norm relative | Prediction absolute | Ratio absolute |
| --- | ---: | ---: | ---: | ---: | ---: |
| Normal | `5.55e-16` | `1.17e-15` | `4.38e-16` | `2.60e-18` | `2.84e-14` |
| Near rank deficiency | `6.66e-16` | `1.55e-15` | `4.38e-16` | `3.47e-18` | `5.68e-14` |
| Rank deficient | `5.55e-16` | `1.19e-15` | `5.53e-16` | `2.17e-18` | `1.82e-12` |
| Tiny gradient | `9.86e-32` | `8.78e-16` | `5.14e-16` | `8.67e-18` | `1.14e-13` |

Two controlled `Advance()` inputs changed decisions with the same supplied bounded step and full residual norms.

| Boundary | Pivoted QR | Direct products |
| --- | --- | --- |
| Proposal threshold `1e-4` | Model-step norm `1.4901161193847654e-6`, prediction `2.2204460492503123e-12`, ratio `0.00010000000000000003`, accepted | Model-step norm `1.4901161193847658e-6`, prediction `2.2204460492503139e-12`, ratio `0.000099999999999999964`, rejected |
| Gradient precision with nearly dependent columns | Scaled gradient `2.2204460492503136e-16`, returned `Running` | Scaled gradient `2.2204460492503131e-16`, returned `GradientPrecisionLimit` |

The proposal fixture used the same actual reduction of `2.2204460492503131e-16`, zero damping and pre-update radius `3.48361350192021e-6`. Both returned `ReductionTolerance`, but only the QR version advanced eta. The gradient fixture held reduction and step tolerances at zero to expose the precision test and accepted the same first step in both versions. These are deliberately threshold-adjacent driver fixtures with a supplied step. They demonstrate sensitivity to arithmetic order without asserting a frequency for production data. The local speedup did not establish a BCD speedup. The final driver retains the pivoted QR and no alternate path.

The proposal fixture used this compact model, a full initial norm of `1`, and a trial norm of `nextafter(1., 0.)`.

```text
Jc = [-.85807777455308432  -.82167451797840085  -.86171356973434066]
     [0                    .69604280936766361  -.71740127394732167]
     [0                   0                    -.67707825218251239]
p  = [2.3271962539772264e-6, -1.4326159676191485e-6, -1.8335562010032557e-6]
```

## P4 workspace decision — Reject

The candidate gave each Profile one active-set workspace while retaining the standalone solver entry point. Only the scratch buffers carried over. The context, rank checks, damping search and working-set initialization still ran per call. The five Joint CTest entries passed. Separate sequences on two concurrent threads agreed with fresh workspaces through size changes, a rank-deficient system recovered by damping, an exhausted damping bracket, invalid bounds, release from an upper bound, and later valid solves. The two Profile traces and four BCD workloads also matched P2.

Each allocation sample covers 1,000 groups of seven solves. Resize times are medians and ranges of three runs. They count actual resizing and reserving of the workspace buffers, excluding factorization and returned-step allocations. Peak is the workspace payload, excluding allocator metadata. Resize time includes allocator and timing overhead.

| Compact columns / radius | Allocations fresh → reused | Requested bytes fresh → reused | Peak bytes per workspace | Resize time fresh → reused |
| --- | ---: | ---: | ---: | ---: |
| 3 / 10 | 84,000 → 12,000 | 2,240,000 → 320,000 | 320 | 1.95 [1.91, 1.96] → 0.28 [0.28, 0.28] ms |
| 3 / .2 | 98,000 → 38,000 | 3,584,000 → 2,240,000 | 416 | 2.55 [2.50, 2.59] → 1.41 [1.36, 1.48] ms |
| 12 / 10 | 84,000 → 12,000 | 14,840,000 → 2,120,000 | 2,120 | 2.01 [1.98, 2.11] → 0.29 [0.27, 0.30] ms |
| 12 / .2 | 98,000 → 38,000 | 32,312,000 → 27,080,000 | 3,368 | 2.59 [2.53, 2.59] → 1.47 [1.43, 1.51] ms |

The damped shapes reallocate the augmented buffer when switching from the undamped matrix, even with reuse. Existing scratch reuse within a bounded solve already covers the damping bracket and bisection. The measured BCD workloads have one bounded solve per block visit, so they have no cross-solve workspace allocation to amortize. Reuse did not produce a stable benefit in the standalone Profile workload or BCD, and would keep the scratch buffers live during trial evaluation and Replay. The final code keeps workspace ownership in `SolveBoundedTrustRegionStep`.

## P4 diagnostic contract — Promote

The only readers of the legacy `lm_status` integer were `tests/support/JointRuntimeJson.hpp` and the trajectory test. The JSON builder is reached through the test-support Profile and component helpers. Searches of tracked production, integration and persistence code found no consumers of the integer. The Profile trace now carries the internal `BoundedLmStatus` in test builds and the test JSON emits names such as `running` in diagnostic schema version 3. The trajectory test checks both the enum and its serialization.

`Numerics.hpp` and `ProfileLmSearch.cpp` keep the accepted eta and scaling diagonal copies, radius, damping, proposed acceptance and detailed working-set trace in the instrumented build. Production still keeps the step, reductions and ratio for BCD failure diagnostics, alongside endpoint evidence and trust. `BoundedLmWork` continues to carry solve counts and KKT summaries into BCD. `ActiveSetStep` continues to report the feasible step, damping, working set, and solver status. The public persisted Joint result continues to use its existing contract.

## Numerical results and work counts

The production BCD runs produced the following exact P0–P4 results, including both experimental builds. The terminal reason was `block-stationary` in all four cases. Each had accepted-state integrity, independent Replay and qualified Endpoint Certification, with no recovery required. The separate transient Width failure regression exercised one failure and one successful retry, then verified assembled and persisted qualification.

| Targets | Sweeps | Normalized objective | Global A/C KKT | Projected Width KKT | Accepted local updates | Profile evaluations | Bounded solves | Driver QRs retained / P3 | Damping QRs |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 3 | `2.0199078132865863e-10` | `3.8858e-16` | `4.0945e-13` | 3 | 6 | 3 | 3 / 0 | 3 |
| 30 | 4 | `1.5171102971714166e-10` | `4.4909e-8` | `5.2731e-10` | 12 | 24 | 12 | 12 / 0 | 12 |
| 300 | 4 | `1.4790860082213613e-10` | `1.4191e-8` | `5.2646e-11` | 104 | 208 | 104 | 104 / 0 | 162 |
| 1,000 | 4 | `1.4790861411606352e-10` | `7.7723e-9` | `1.5792e-11` | 336 | 672 | 336 | 336 / 0 | 394 |

The standalone 12-target Profile recorded 7 evaluations, 5 derivatives and 4 accepted updates per call, with 7 bounded solves and 72 damping QRs at every stage. It executed 5 Driver QRs per call with the retained model and zero in P3. Driver QR counts follow the code path, once per successful linearization. In these BCD runs every block visit made one successful linearization and accepted one update. The 129 fixed trace records were identical again in the final build.

## Performance and memory

Three fresh-process measurements per cell used the Release production library on arm64 macOS with AppleClang 21.0.0 and one OpenMP thread. The benchmark interleaved all six binaries and varied their order between repetitions. Cells show the median [minimum, maximum]. The separate Driver microbenchmark measured 40,000 individual `Advance()` calls on a 12-column dense model with a supplied feasible step and excluded initialization. The Profile timings include the actual bounded solver and Replay. BCD timings include complete sweeps and endpoint verification.

### Wall time per call

| Workload | P0 | P1 | P2 | P3 direct | P4 workspace | P4 final |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BCD, 3 targets (ms) | 6.48 [6.31, 7.82] | 6.51 [6.48, 6.51] | 7.07 [6.45, 7.45] | 6.61 [6.12, 7.08] | 7.07 [6.52, 7.61] | 6.71 [6.22, 7.29] |
| BCD, 30 targets (ms) | 92.64 [90.64, 99.57] | 93.07 [90.87, 93.66] | 93.41 [90.59, 93.58] | 91.92 [91.69, 93.13] | 94.86 [94.17, 95.30] | 95.35 [92.36, 95.39] |
| BCD, 300 targets (ms) | 992.09 [965.77, 1023.88] | 972.87 [951.31, 978.27] | 981.61 [969.23, 998.70] | 1017.20 [977.38, 1020.21] | 975.53 [963.74, 1015.01] | 965.57 [956.92, 1038.23] |
| BCD, 1,000 targets (ms) | 3281.39 [3220.55, 3371.72] | 3298.38 [3294.10, 3407.21] | 3265.22 [3264.67, 3286.11] | 3264.50 [3191.13, 3336.73] | 3293.66 [3256.58, 3625.73] | 3226.97 [3225.55, 3278.75] |
| Profile, 12 targets (ms) | 27.28 [26.66, 27.39] | 27.80 [26.93, 28.78] | 27.27 [27.15, 28.55] | 27.40 [26.82, 29.01] | 27.47 [26.50, 28.10] | 27.88 [26.64, 28.36] |
| Driver, 12 columns (µs) | 2.074 [2.037, 2.147] | 2.061 [2.046, 2.152] | 2.146 [2.141, 2.166] | 0.506 [0.504, 0.519] | 2.128 [2.116, 2.208] | 2.143 [2.097, 2.145] |

### Peak process RSS (MiB)

| Workload | P0 | P1 | P2 | P3 direct | P4 workspace | P4 final |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| BCD, 3 targets | 22.8 [22.8, 22.8] | 22.7 [22.7, 22.8] | 22.7 [22.7, 22.8] | 22.6 [22.6, 22.7] | 22.7 [22.7, 22.8] | 22.6 [22.6, 22.7] |
| BCD, 30 targets | 33.2 [33.0, 33.5] | 33.0 [33.0, 33.1] | 33.0 [32.9, 33.1] | 32.8 [32.8, 32.8] | 32.9 [32.8, 33.1] | 32.8 [32.8, 32.9] |
| BCD, 300 targets | 81.2 [81.1, 86.3] | 81.3 [81.1, 82.3] | 81.2 [81.2, 82.0] | 81.6 [81.5, 83.7] | 81.2 [81.1, 84.2] | 81.1 [81.0, 84.2] |
| BCD, 1,000 targets | 210.2 [205.3, 215.6] | 216.4 [216.0, 219.3] | 215.5 [206.2, 215.7] | 215.8 [215.3, 215.9] | 206.4 [203.5, 216.0] | 217.0 [215.6, 217.2] |
| Profile, 12 targets | 28.1 [28.0, 28.1] | 28.1 [28.1, 28.8] | 28.2 [28.2, 28.3] | 28.1 [28.1, 28.2] | 28.1 [28.1, 28.2] | 28.1 [28.1, 28.2] |
| Driver, 12 columns | 20.1 [20.0, 20.1] | 20.0 [20.0, 20.1] | 20.1 [20.0, 20.1] | 20.1 [20.0, 20.1] | 20.0 [20.0, 20.1] | 20.0 [20.0, 20.0] |

The P3 candidate reduced isolated Driver time about 4.2× relative to P2. Its BCD and standalone Profile times were within the variation across runs. The Driver fixture's norm-only adapter recomputes the compact response during linearization, so the P0-to-P2 Driver microbenchmark also includes that fixture difference. The production Profile and BCD measurements provide the comparable end-to-end timings. The final build showed no repeatable BCD speedup or slowdown in these samples.

Peak RSS varies with process allocation layout. At 1,000 targets, the initial and interleaved runs gave opposite directions for the median RSS difference. Pooling those six runs per build gives 213.0 [205.3, 221.3] MiB for P0 and 215.8 [202.7, 217.2] MiB for the final build. The ranges overlap and do not establish a persistent peak-RSS change. Total process allocation counts and allocator overhead were unavailable. The specific buffer reductions were two full residual vectors in P2 and two test-trace vectors per populated trial in the final P4 build. P3 skipped the compact QR allocation and used a gradient vector, while the measured P4 workspace allocation data above covers only the bounded solver's scratch.

## Build and verification commands

Both build configurations compile `rhbm_gem` and `rhbm_gem_test`. CTest retains its one-thread OpenMP settings for the Joint suites. The filtered Release command covered P1, P2 and both experimental builds. The final Release run passed 21 CTest entries with zero failures or skips, including 1,117 GoogleTest cases. The focused Debug run passed three CTest entries with zero failures or skips, including 248 GoogleTest cases.

```sh
export OMP_NUM_THREADS=4
cmake --preset release -DBUILD_TESTING=ON
cmake --build build/release --parallel 4
ctest --test-dir build/release --output-on-failure --parallel 4 -L 'joint:'
ctest --test-dir build/release --output-on-failure --parallel 4

cmake --preset debug -DBUILD_TESTING=ON
cmake --build build/debug --parallel 4 --target rhbm_gem rhbm_tests
ctest --test-dir build/debug --output-on-failure --parallel 4 \
  -R 'rhbm_tests_joint_component_(numerical|contract|workflow)$'
```

## Git and cleanup

All work was performed on the existing `develop` checkout. The final code HEAD is `685f9f7820b3e754cf7adf91b9c42f526097c94d`. The report commit is the final checkout HEAD and can be resolved with `git log -1 --format=%H -- docs/developer/joint-lm-simplification-validation.md`.

| Commit | Purpose |
| --- | --- |
| `bae97830` | `docs(joint): establish LM data flow and validation baseline` |
| `05e9f480` | `refactor(joint): eliminate redundant bounded LM diagnostics and KKT work` |
| `9a39c1c6` | `refactor(joint): use residual norms in bounded LM driver` |
| `685f9f78` | `refactor(joint): prune retired LM diagnostic contracts` |
| This report commit | `docs(joint): record LM simplification validation and decisions` |

The final working tree is clean and `develop` is five commits ahead of the local `origin/develop` tracking ref recorded at the start. No push was made. No branch or worktree was created, no Notion content was modified, and no dependency was added. Experimental solver variants, benchmark binaries, and temporary instrumentation were removed. The retained path preserves the formal numerical policies and their qualification thresholds.
