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

The second trace uses the same chain with widths `.41, .89, 1.3`, upper bound `2.2`, and the default evaluation and update budgets. It has eight evaluations, seven derivatives and bounded solves, six accepted updates, 112 damping QRs and one final rejected trial. The final accepted objective is `2.0199078132865648e-10`, Width KKT is `7.0610184366159956e-14`, and Independent Replay and endpoint qualification pass. Its terminal reason is `native-lm-stop`, with `width_failed=false`.

Use the `JointComponentSearchTest`, `ActiveSetBoundedLMTest`, `JointBcdSearchTest` and `JointComponentRuntimeTest` suites for the fixed retry, tiny-step, boundary, rank, replay and recovery fixtures. `TransientWidthFailureRecoversAndAssemblesQualifiedComponent` requires one recovered failure, qualified certification, and a converged assembled result.

For timing, use `SyntheticJointWorkload("chain", count)` and time a production-library `SearchBcdComponent` call, excluding input construction. Use counts 3, 30, 300 and 1,000, initial `eta=log(.55)`, upper bound `log(2.2)`, 12 core atoms, 12 maximum sweeps, practical convergence and default profile budget. BCD sets one accepted Width update per visit. The benchmark's fixed cap applies only to this workload, which qualified within four sweeps at every size. For standalone Profile timing, use the same chain with 12 targets and an eight-update budget, and time five calls per process. Use one OpenMP thread, a fresh process for each sample and macOS `getrusage(RUSAGE_SELF).ru_maxrss` in bytes.

Baseline Release configuration and validation used `cmake --preset release -DBUILD_TESTING=ON`, `cmake --build build/release --parallel 4` and `OMP_NUM_THREADS=4 ctest --test-dir build/release --output-on-failure --parallel 4`. All 21 CTest entries passed.
