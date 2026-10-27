# Joint component evidence

This file is the compact historical record for research decisions. The
production behavior is described in [Joint component runtime](joint-component-runtime.md).
Values below apply only to the stated datasets and revisions; they are not
general accuracy or performance guarantees.

## Current production contract

The internal `EstimateJointComponents()` entrypoint and the normal
`joint-components` workflow use
Target-only Joint fitting with Pure Halo-Support Exclusion at 2.5 Å. Halo
support rows are excluded from fitting; Halo parameters do not enter the Joint
objective. The solver uses disjoint FixedNeighbor cores of at most 12 Targets,
forward serial Gauss–Seidel sweeps, at most 200 complete sweeps, and at most
one trusted accepted update per block visit.

Target widths are constrained by one global bound,
`Bmax = 4 × Q95(qualified FirstStage Target widths)`. Calibration requires at
least 20 qualified donors and 95% Target coverage. The local A/C solve retains
nonnegative A and signed C constraints. Bounded Active-set LM, monotone global
objective acceptance, complete-sweep confirmation, projected A/C and Width
KKT, endpoint verification, and Independent Replay determine whether a
component is qualified. Unqualified estimates are not published downstream.

After qualification, Halo A/B/C values are reconstructed from the nearest
four qualified Target estimates for Original-domain Peeling. Conditioning and
stability are read-only local diagnostics; they do not change the accepted
Target estimates or independently deny a qualified endpoint. Current Joint
results use JSON schema v19, target-only Pure Exclusion, geometry-local kNN-4
Halo reconstruction, and no legacy decoding.

## Search policy and block solver

The fixed-width conditional solver was numerically feasible but was not
promoted; the nonlinear FixedNeighbor search was selected. Repeated
correctness-checked chain and cube experiments selected core size 12 over the
historical core-64 control. Median search time improved by 33.0%–64.3% in the
four measured cases, with every compared run passing its search and replay
checks:

| Case | Core 12 | Core 64 | Measured reduction |
| --- | ---: | ---: | ---: |
| Chain, 512 atoms | 61.18 s | 171.26 s | 64.3% |
| Cube, 512 atoms | 151.02 s | 322.81 s | 53.2% |
| Chain, 1,024 atoms | 167.73 s | 346.28 s | 51.6% |
| Cube, 1,024 atoms | 490.85 s | 732.64 s | 33.0% |

Forward order is the production serial Gauss–Seidel order. In one qualified
6Z6U comparison, forward and reverse runs used 79 and 68 sweeps and produced
objectives `0.006497513583164643` and `0.006497513583172617`; this is a
single-dataset order check, not a general equivalence result. OneAccepted was
selected over Full and TwoAccepted local update policies. Fixed-B, reverse
order, core 64, and the alternative update policies are not production
options. Source: [commit `de211df5`](https://github.com/yslianTim/rhbm_gem/commit/de211df5a57647fccc0c8abee2262f82ccdb0c53).

## Linear solver and qualification choices

The production constrained A/C solve uses Eigen SparseQR with normalized free
columns, a nonnegative-amplitude active set, bounded local QR refinement, and
the existing strict local KKT limit of `1e-10`. The Block 21 study used a full
rank `930 × 24` design with a free-face condition indicator near `6.62e10`.
SparseQR, dense QR, SVD, active-face enumeration, and direct free-face QR all
missed the strict gate; the best recorded scaled KKT was `4.65661e-10`. A
50-digit reference was accurate in its own arithmetic, but rounding its
coefficients did not give a double candidate that also passed the independent
functional-kernel certificate. No multiprecision production fallback or
parameter rotation was adopted. The strict solve failure remains an explicit
qualification failure, not a relaxed threshold. Source: [commit `aebf67c6`](https://github.com/yslianTim/rhbm_gem/commit/aebf67c64a5d80071daef698aa887a36d04d7c87).

The former global-rank certification path was removed because its dense and
projected-width workspaces were not suitable for the largest planned shape:
the retired estimates included `67.2–78.4 GB` for the earlier BDC-SVD shape,
`189.5 GB` for dense design spectrum, and `39.1 GB` for a combined
projected-width stage. These were planner estimates for retired algorithms,
not measured allocations or current production requirements. Production keeps
the local rank and conditioning checks needed by each block and makes no
global identifiability claim. Sources: [commit `d9fcbd9e`](https://github.com/yslianTim/rhbm_gem/commit/d9fcbd9ec219ddab04d383cf24ce72b0950a49a7) and [commit `de211df5`](https://github.com/yslianTim/rhbm_gem/commit/de211df5a57647fccc0c8abee2262f82ccdb0c53).

The global width cap and bounded solver were promoted without relaxing the
local A/C solve, endpoint replay, or practical global qualification limits. A
post-promotion 6Z6U chain-A simulation run used 1,415 qualified donors among
1,418 Targets, `Q95 = 0.5689788338510722 Å`, and `Bmax = 2.2759153354042887 Å`.
It converged in 83 sweeps and 22,908 profile evaluations; Independent Replay
and runtime qualification passed. Candidate objective was `1.2330854370810358e-13`,
prediction RMSE was about `9.21e-7`, A/C KKT was `9.43e-8`, and projected Width
KKT was `8.78e-10`. The two recorded Release runs took 298–305 s and measured
586–617 MB peak RSS. They used 1,418 Targets and 180,478 observations; they do
not close the separate 37,406-Target full-workflow resource gate. Sources:
[commit `de211df5`](https://github.com/yslianTim/rhbm_gem/commit/de211df5a57647fccc0c8abee2262f82ccdb0c53) and [commit `975ed2c7`](https://github.com/yslianTim/rhbm_gem/commit/975ed2c72bfb35df695a387770f6ef6d796bba44).

## Halo policies and FullABC promotion decision

Pure Halo-Support Exclusion is the production policy: it removes each retained
Halo's recorded support rows, including the exact 2.5 Å boundary, from the
Target fitting domain. It does not require Halo FirstStage parameters and does
not infer Halo estimates from the excluded rows. The original observation
domain and mappings remain available for reconstruction and Peeling.

An earlier Dynamic Halo-mean study did not justify adopting a changing Halo
mean. On its recorded 6Z6U case, the dynamic E2 result had an observation
objective 58.4% above the fixed E0 result, with Halo, Target, and full-field
prediction RMSEs 5.73×, 5.04×, and 6.49× E0. Forward/reverse E2 runs both
qualified with exact replay, but the Target parameter changes did not improve
the field fit. This is one historical model and does not establish that the
policies are interchangeable. Separately, a 450-fit noise/mismatch pilot from
the retired Fixed-Halo/Profile-LM path had 65 qualified and 385 unqualified
fits; it is deliberately unbalanced and is not a population success rate or a
validation of current Pure Exclusion Active-set LM. Source: `develop` commit
`c84eeb42c0abaac314fce20309114313ada5af40`.

FullABC and Recovery were evaluated as opt-in research paths and were not
promoted; their solver paths have since been removed from production source.
Two different P4 reports used different samples and gates; their results must
not be pooled:

| Report | Stated sample and result | Decision and limitation |
| --- | --- | --- |
| P4 promotion gate, [`da63a718`](https://github.com/yslianTim/rhbm_gem/commit/da63a718459e7ca2191484679862f15126acbab6) | P2 had 5 local scenes with 0/5 FullABC-qualified results. A boundary Recovery endpoint replayed consistently but had `A ≈ 8.35e122`, `B ≈ 0.0934 Å`. In a separate three-short-bond seed comparison, component-local Recovery Halo-signal RMSE was `0.00652` versus `0.07363` for kNN-4; boundary Recovery Peeling RMSE was about `1.91e57`. | No promotion: replay and objective agreement did not exclude a pathological endpoint. Five local scenes had no Group members; the long Chain-A run was interrupted before a peak-RSS result, and full save/reload integration was not measured. |
| P4 production feasibility, [`969f3371`](https://github.com/yslianTim/rhbm_gem/commit/969f3371e0a43af9ca136c119ee939a46635721b) | A distinct P3 synthetic physical set had 10/10 qualified Target endpoints and 8/10 qualified FullABC endpoints; the exact-3 and exact-4 observation free-Halo cases failed qualification. In the 8 qualified pairs, Halo-signal and Peeling RMSE improved, while Local/Group estimates improved in only 4/8. One high-noise A case had RMSE `1.413` and maximum absolute error `3.150`. | No promotion: small synthetic improvement did not establish downstream reliability. The complete 5,938,087-row Halo FullABC workflow, its actual Halo count, process RSS, persistence, and UMAP/Two Stage integration were not measured. |

The production Halo reconstruction remains geometry-local kNN-4 Target Median.
It uses only qualified Target estimates and has deterministic Atom-ID tie
breaking. FullABC failures, Recovery outputs, and the reports' small-sample
signal gains do not change this contract.

## Local conditioning, stability, and large-component limits

Local Conditioning assesses the free A/C face, conditional projected-width
directions, and the full profile derivative within each FixedNeighbor block.
Parameter Stability probes at most eight representative blocks, with bounded
local solves and `±1e-3` log-width seeds. A full-rank matrix can still be
ill-conditioned; the current roundoff-amplification warning is
`κ × ε_double >= 1e-8`. These diagnostics are conditional local evidence, not
covariance, global identifiability, or endpoint qualification. They do not
modify accepted estimates or turn an unrun probe into `Stable`. Sources:
[commit `32cb6362`](https://github.com/yslianTim/rhbm_gem/commit/32cb6362a3d988d7fa77f715515351a9071d58dd), [commit `ddfabe47`](https://github.com/yslianTim/rhbm_gem/commit/ddfabe470f37886086c800c4e657b4b6395aea6a), and [commit `08151497`](https://github.com/yslianTim/rhbm_gem/commit/08151497f68f6c49b250c97f94f4dbe562dc092f).

The deterministic large-shape checks use 5,938,087 rows, 37,406 contributors,
74,811 support memberships, and 3,118 blocks. One preparation/zero-state
verification run measured 622,460,928 bytes peak process RSS, with preparation
and endpoint estimated peaks of 972,989,100 and 625,051,683 bytes. A separate
preflight measured 1,721,008,128 bytes peak RSS. These runs did not execute a
complete search, local conditioning and stability over all blocks, or
downstream work. The planner's 3.2 GB estimate budget and 4.0 GB hard limit
are not measured RSS guarantees. The full physical-data resource gate remains
open. Sources: [commit `d9fcbd9e`](https://github.com/yslianTim/rhbm_gem/commit/d9fcbd9ec219ddab04d383cf24ce72b0950a49a7) and [commit `de211df5`](https://github.com/yslianTim/rhbm_gem/commit/de211df5a57647fccc0c8abee2262f82ccdb0c53).

## Peeling and downstream workflow

Post-Joint Local charge correction was fixed to use Original-domain support
and row lookup. In the recorded 6Z6U B/C/D comparisons, Joint A/B/C, Halo
fields, and Peeling were unchanged; Local and Group results improved, but
case D still had Local A/B RMSE `0.709793 / 0.021414`. The evidence supports
the domain-correctness fix, not complete scientific recovery or a change to
the Local MDPDE/Group Fitting contracts. An independent small synthetic
workflow also covers Pure Exclusion → kNN-4 Halo reconstruction → Original
Peeling → own-charge correction → Local MDPDE → Group Fitting. Sources:
[fix `a9684293`](https://github.com/yslianTim/rhbm_gem/commit/a968429320eebaf79e1e960b576eea83f7dc8f8d), [validation `7cc9a35e`](https://github.com/yslianTim/rhbm_gem/commit/7cc9a35e8ae2e4ad497cebf092c4cdaedd0f1852), and [workflow validation `975ed2c7`](https://github.com/yslianTim/rhbm_gem/commit/975ed2c72bfb35df695a387770f6ef6d796bba44).

The persisted result retains the Original/Fitting domain identities, exclusion
mappings, width calibration and global bound, accepted solver provenance,
qualification evidence, and Halo donor provenance. Historical JSON/SQLite
decoding remains part of the compatibility contract. The Joint objective does
not estimate a global parameter covariance; downstream Local MDPDE and Group
Fitting retain their own qualification rules.

## Bounded Active-set LM and Profile Derivative optimization baseline

Frozen baseline: `8275deea7ea21c4891ef2a882239273005313511` on `develop`.
The local Release build used AppleClang 21.0.0.21000334, C++20, Eigen 5.0.1,
Ninja, macOS Darwin 25.6.0 arm64 (Apple M1), OpenMP 5.1, Boost 1.92.0, and
SQLite 3.54.0. CMake used `BUILD_TESTING=ON`,
`RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE=ON`, `RHBM_GEM_ENABLE_UMAP=OFF`,
`RHBM_GEM_DEP_PROVIDER=SYSTEM`, and `CMAKE_BUILD_TYPE=Release`. Build and
runtime thread limits were four; CTest ran serially. The benchmark was built
with test-only instrumentation and run with five warmup calls followed by five
samples of five calls each. Its full numeric traces are hashed from the raw
IEEE-754 bit patterns; timing is diagnostic-build data and should only be
compared with the same configuration.

The Bounded LM fixture is a deterministic `4096 × 8` dense subproblem with a
nonzero damping bracket and bisection. The Profile Derivative fixture is a
30-Target chain with 13,710 rows and 60 derivative columns, reduced in two
tiles at the existing 8,192-row limit. The separate synthetic Joint timing is
a 3-Target, 3-sweep diagnostic run; it is a local regression fixture, not a
large-workload performance claim.

| Baseline measure | Result |
| --- | ---: |
| Bounded LM total, median (5 samples) | 17.059 ms/call (16.851–17.173) |
| Damping solves / bracket / bisection | 59 / 4 / 54 per call |
| Dense QR factorizations | 136 per call |
| Dense QR factorization time | 11.620 ms/call |
| Explicit temporary allocation sites | 1,581 per call |
| Estimated temporary workspace peak | 328,584 bytes, excluding QR internals |
| Derivative preparation, median (5 samples) | 20.233 ms/call (19.807–20.507) |
| Derivative reduction, median (5 samples) | 5.977 ms/call (5.860–6.007) |
| Tiled QR time during preparation | 16.986 ms/call |
| Triplet vector capacity growths | 32 per preparation |
| Estimated maximum tiled workspace | 5,898,240 bytes |
| Synthetic Joint search, median (5 samples) | 3.317 ms/call (3.263–3.386) |
| Benchmark process peak RSS | 44,728,320 bytes |

Baseline Bounded LM trace SHA-256:
`290e1a65097e34b1a1198064bc49c2c357a83271c709ce6b156ccbb1eec798ab`.
Profile preparation/reduction bit-pattern SHA-256:
`8c2f90707e1f5b97058dc6708c58b136737878f3795d5c731f2150555255c9d7`.
The joint contract, workflow, and numerical CTest groups all passed (3/3,
5.69 s). Deterministic Bounded LM tests pin result bits and record the damping
candidate sequence, per-iteration active masks, steps, gradients, KKT
tolerances, free-face ranks, and trust-region feasibility decisions; a release
fixture verifies the active-face transition bitwise. A three-row fixture
confirms cancellation detection and reference-order recomputation in both full
and projected-only derivative preparation.
