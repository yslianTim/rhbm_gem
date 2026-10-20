# FixedNeighbor Joint: current contract and evidence

This guide records the current production route and the historical studies that
selected its search policy. Retired global-rank certification code and its
campaign tools are not part of the current estimator.

Width-trust and global-bound comparison runners, Dynamic Halo implementations,
and completed Halo strategy campaigns have also been retired. Their numerical
decisions remain in `docs/research/`; current benchmarks use the Production
FixedNeighbor path, while optional numerical audits retain independent A/C,
Profile, width, and replay references. Fixed-Halo compatibility remains a
supported Profile LM path with dedicated tests.

## Current production contract

`FixedNeighbor` is the sole Joint estimator selected by
`--estimator joint-components`. Selected non-hydrogen atoms fix observation
rows. An unselected halo enters the contributor model only when its support
contains at least three rows and three exactly distinct squared distances.
All retained atoms use FullABC. Rows and observations remain fixed after halo
selection; excluded halo signal remains in the residual and can shift target
estimates. See the [halo eligibility guide](joint-observable-halo.md).

| Policy | Current value |
| --- | --- |
| Mutable core | 12 atoms, disjoint blocks |
| Block order | Forward serial Gauss-Seidel |
| Maximum sweeps | 200 complete sweeps |
| Local search | Bounded Active-set LM (`fixed_neighbor_local_search=active-set-lm`) |
| Accepted local update | At most one trusted update per block visit |
| Acceptance | Existing trust checks and roundoff-safe monotone global objective |
| A/C stationarity | Global projected KKT `<= 1e-7` |
| Width stationarity | Global projected Width KKT `<= 1e-9` |
| Sweep confirmation | Prior complete sweep and log-width change `<= 1e-6` |
| Endpoint | Primary streaming evaluation and independent streaming replay |

A confirmed stationary endpoint becomes the component state only after
numerical endpoint verification. Replay checks feasibility, prediction,
residual, objective, A/C KKT, and projected Width KKT consistency at the
accepted parameters; it does not solve a second global A/C problem. Recoverable numerical
verification failures can be retried at a bounded interval. Invalid primary
state and endpoint-verification resource limits stop the search. Maximum
sweeps means budget exhaustion and does not imply convergence. If the last
sweep has not been checked, one final bounded check still requires that sweep
to meet stationarity.

After convergence, each fixed-neighbor block receives local conditioning
evidence at the final accepted state. It reports free A/C, conditional
projected-width, and full profile-derivative conditioning. Bounded parameter
perturbation probes run only for up to eight sensitive blocks, with explicit
profile-evaluation, update, workspace, and CPU budgets. These diagnostics do
not change accepted parameters, objective, or convergence status. Their
warnings are point-quality evidence, not global identifiability claims.

The objective, constrained A/C solve, global objective acceptance, and
independent endpoint replay retain their existing contracts. The full
derivative keeps the nonzero-residual coefficient correction. Pure Exclusion
Production bounds every target width by `4 * Q95` of qualified FirstStage
target widths; calibration needs at least 20 donors and 95% target coverage.
Projected Width KKT certifies an active upper boundary. Calibration failure
persists an unavailable result and does not fall back to unbounded search.
Fixed-Halo compatibility retains the unbounded Profile LM path.

For accepted trials, the local state cache is rebased from the full candidate
replay already used for acceptance. Sparse predictions, transpose products,
derivative reductions, and independent scalar replay use compensated double
arithmetic. The 50-digit reference remains confined to offline precision
controls.

## Memory planning and large components

`ResourceWork` preserves preflights for construction, structural construction,
FixedNeighbor preparation, global streaming evaluation, endpoint verification,
reference solves, local conditioning, and parameter stability. The planning
budget is 3,200,000,000 estimated peak bytes below a 4,000,000,000-byte hard
limit. These are allocation estimates, not observed process RSS. Current local
conditioning and perturbation workspace follow block dimensions; the global
rank factors and projection coefficients from the retired pipeline are gone.

A deterministic shape with 5,938,087 rows, 37,406 contributors, 74,811 support
memberships, and 3,118 blocks previously reported 2,892,736 prepared-mapping
bytes and an 876,119,666-byte streaming-evaluation estimated peak. A direct
no-sweep test process recorded 474,546,176 bytes peak RSS. The run disabled
sweeps and endpoint verification; it does not measure a full physical-data fit.
The exact preparation estimate depends on actual affected rows and support
memberships. The more recent 37,406-target bounded preflight is recorded in the
[bounded resource validation](../research/joint-bounded-giant-component-resource-validation-20261020.md);
it also did not run a full fit, so the 4 GB full-workflow gate remains open. See
the [memory census](joint-component-memory-census.md) for the measurement scope
and retired rank estimates.

The benchmark driver reports structural counts, estimates, available RSS
checkpoints, and process peak RSS. Preparation-only, one-sweep, and bounded
multi-sweep modes are measurement harness settings; they do not represent a
physical full-search result. Per-block timings and profile trials require the
`--detailed` option.

## Historical KKT and reference agreement evidence

This section records the retired coefficient-agreement contract. Production
uses the versioned reference-free contract above. KKT is a first-order
optimality check and is not numerically equivalent to coefficient agreement.

For fixed widths and an unchanged active set, let `X_F` contain the free A/C
columns and let `g_F = X_Fᵀr` be the least-squares gradient on that face. In the
unconstrained full-rank case, `g_F = (X_FᵀX_F)(β_F - β_F*)`, giving

```text
||β_F - β_F*||₂ ≤ ||g_F||₂ / σ_min(X_F)².
```

With nonnegative A constraints, this bound applies only on a stable free face;
active coordinates also require a constrained error bound using projected KKT
and the active multipliers. The runtime projected KKT value is a normalized
residual, while reference agreement is a scaled coordinate-wise infinity
difference. They are not interchangeable, and numerical rank alone does not
bound their relationship without conditioning and scaling information.

The 168-atom diagnosis provides this historical continuation evidence at the
same input and initial widths. The 150-sweep row is a diagnostic run that
continued beyond the 107-sweep candidate; it is not an offline regular
certificate or a replay of the new production contract.

| Sweeps | Projected A/C KKT | Primary width gradient | Reference width gradient | Scaled A/C difference | Endpoint trust |
| ---: | ---: | ---: | ---: | ---: | --- |
| 30 | 3.6848e-5 | 7.3085e-7 | 1.4159e-6 | 1.0196e-3 | Failed: `kkt-failed` |
| 107 | 4.5219e-11 | 8.3490e-13 | 1.6499e-12 | 1.2573e-9 | Failed: `reference-disagreement` |
| 150, continued diagnostic | 8.3766e-14 | 9.1791e-15 | 1.2882e-14 | 2.8179e-12 | Passed: `trusted` |

Under the legacy contract, at sweep 107, the KKT and primary width-gradient checks passed, but the
coefficient difference remained above `1e-10`, and the reference width gradient
remained above `1e-12`. Continuing the ordinary accepted search reduced both
disagreements enough to pass the existing endpoint checks. This case therefore
shows that KKT did not imply the legacy reference agreement check, and that
the observed disagreement was recoverable through more sweeps. It provides no
evidence of an irreducible floating-point failure.

The saved diagnosis does not contain the active-face `σ_min(X_F)` or condition
number, so it cannot establish how much conditioning explains the sweep-107
coefficient difference. Because projected KKT and scaled coefficient metrics
use different normalizations, these recorded values cannot certify coefficient
accuracy under the current reference-free contract.

## Historical decisions

| Topic | Conclusion |
| --- | --- |
| Fixed-B | Conditional block updates and replay accounting were numerically feasible, but Fixed-B was not promoted as the estimator. Nonlinear FixedNeighbor superseded it. |
| Core size | Repeated correctness-checked qualification selected core 12. Against the historical core-64 control, measured median search improved by 33.0%–64.3% across chain/cube at 512 and 1024 atoms. Larger cores did not justify their added local problem size. |
| Block order | Forward was adopted as the production serial Gauss-Seidel order. Reverse was a qualification diagnostic; the largest fully checked parameter differences stayed below existing interpretation bounds. |
| Stationarity | Eta-only confirmation was selected with the historical `1e-10` threshold and requires a previous complete sweep. The bounded production route uses practical-v1 projected A/C KKT `1e-7`, projected Width KKT `1e-9`, and log-width change `1e-6`; the stricter values remain in legacy numerical controls. |
| Prepared blocks | Structural block, domain and mapping preparation is checked by a deterministic test. Preparation adds no convergence criterion and does not reuse numeric factors across sweeps. |
| Local update | OneAccepted was selected over Full and TwoAccepted local updates; the at-most-one trusted accepted update is now intrinsic to production. |
| Optimized/frontier variants | The historical OneAccepted/core-64 and alternative-route frontiers were superseded by the qualified core-12 production policy. |

The repeated core-size measurements that explain the adopted policy are:

| Case | Core 12 median search | Historical control 64 | Improvement |
| --- | ---: | ---: | ---: |
| chain-512 | 61.18 s | 171.26 s | 64.3% |
| cube-512 | 151.02 s | 322.81 s | 53.2% |
| chain-1024 | 167.73 s | 346.28 s | 51.6% |
| cube-1024 | 490.85 s | 732.64 s | 33.0% |

Every measured run passed search/replay correctness. The 64-atom values are
historical controls, not the current production default.

The table is a compact historical comparison only; core 64 is not a current
default or benchmark option.

## Closed directions

These results are retained as compact negative evidence only. They do not
define a production option or a current CTest runtime contract.

| Direction | Decision | Reason retained |
| --- | --- | --- |
| Local `OperatorPcg` | Closed | The pure local route failed the existing local gate: `assessment_local = 1.458662526e-10 > 1e-10`, despite passing global KKT and endpoint certification in the diagnostic case. |
| Local Schwarz geometry/tuning | Closed | No pure Schwarz candidate passed the unchanged endpoint/`RuntimeConvergence` gate on both topologies; the best attribution candidate was not promoted. |
| Operator-to-Legacy hybrid polish | Closed | The selected route failed the 512 promotion gate: cube-512 eta confirmation was `1.547e-10 > 1e-10`; there is no production hybrid method or fallback. |
| Workspace residency | Closed | Persistent workspace had no material wall-time benefit and higher RSS in the measured cases; the current route keeps fresh numeric factorization. |
| Exact numeric-factor reuse | Closed | Exact reuse opportunities were only `1.17%`–`15.00%` and occurred at initial-profile evaluation; no factor-cache policy was promoted. |
| Candidate replay localization | Closed | Replay and accepted-state semantics remain the existing production contract; no localized replay optimization was adopted. |
| Derivative tile-size tuning | Closed | No tested tile reached the 10% improvement gate; production QR/scratch policy is unchanged. |
| Historical core 64/128/256 campaigns | Superseded | Their compact results are controls for the core12 decision; they are not current defaults or separate permanent drivers. |
| Endpoint decomposition/local certification/attribution | Diagnostic only | These answered completed research questions and had no production invariant after qualification; their campaign drivers and runtime tests were retired. |

The negative evidence is intentionally numeric where it affects a decision.
The retired FixedNeighbor-local OperatorPcg, Schwarz and hybrid hooks, as well
as the global OperatorPcg/Schwarz route, are removed from the current tree.
The EIGEN sparse backend remains current because FixedNeighbor's conditional
profile solve uses it. There is no hidden FixedNeighbor
policy switch. Closed experiment implementations and machine-readable receipts
are recoverable from Git history and are intentionally not retained in the
current tree.
