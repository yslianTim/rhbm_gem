# FixedNeighbor Joint: current contract and evidence

This is the concise overview for the production FixedNeighbor route. Historical
execution telemetry and retired campaign implementations remain recoverable
from Git history; they are not part of the current production surface.

## Current production contract

`FixedNeighbor` is the sole Joint estimator. It is selected by
`--estimator joint-components`; there is no separate global search-route
option or automatic routing layer.

Before FixedNeighbor builds components, selected non-hydrogen atoms fix the
observation rows. Unselected halo candidates enter the model only when their
support has at least three rows and at least three exactly distinct squared
distances. All retained atoms use FullABC. The rows and observations stay fixed
after pruning, and components are rebuilt from retained supports. Excluded
signal remains in the residual, so this contributor rule changes the additive
model and can shift target estimates. Selection eligibility and the full
numerical qualification below are separate gates. See the
[halo eligibility guide](joint-observable-halo.md) and the
[target-recovery study](../research/joint-halo-pruning-validation.md).

| Policy | Current value |
| --- | --- |
| Outer core | 12 atoms, disjoint mutable blocks |
| Block order | Forward serial Gauss-Seidel |
| Maximum sweeps | 200 |
| Local search | Profile LM (`fixed_neighbor_local_search=profile-lm`) |
| Local update contract | At most one trusted accepted Profile LM update per block visit |
| Acceptance | Existing replay and roundoff-safe monotone acceptance |
| Candidate stationarity | A/C KKT `≤ 1e-10` and width-gradient infinity norm `≤ 1e-12` |
| Confirmation | A previous complete sweep exists and eta-change infinity norm is `≤ 1e-10` |
| Endpoint | Primary streaming evaluation and independent streaming replay must pass |
| Certification retry | Valid failed endpoints continue; retry every 10 confirmed sweeps while all stationarity requirements hold; any failed requirement resets eligibility |
| Exhaustion | One bounded assessment certifies the last accepted state if its final sweep has not already been certified; a trusted state may remain available while search is incomplete |

The unchanged thresholds are A/C projected KKT `1e-10`, width stationarity
`1e-12`, and eta confirmation `1e-10`. A confirmed stationary endpoint is
checked at the accepted parameters by a primary streaming evaluation and an
independent scalar replay. The replay verifies prediction, residual, objective,
KKT and width-gradient consistency; it does not solve a second global A/C
problem. Conditional A/C optimality, endpoint verification and full numerical
qualification are separate evidence. Full qualification also requires design
rank, projected-width rank and numerical identifiability. A threshold-sensitive
or unavailable rank cannot pass. If rank resources are unavailable, the
endpoint may be verified while qualification and `RuntimeConvergence` remain
unavailable. Valid failed endpoints continue normal
updates and retry certification at a bounded interval. Budget exhaustion does
not increase the sweep count or mark search complete: one final assessment runs
on the accepted state only if its final sweep has no certification. Runtime
convergence remains independent of search completion. Prepared structural
block mappings and workspace reuse do not
reuse numeric factors, accepted states, or replay results across sweeps. The
route preserves the objective, local A/C solve, width parameterization, serial
update semantics and OneAccepted behavior. A sparse unit-weight local A/C solve
that fails its existing KKT gate can use the bounded free-face QR residual
refinement described in the [Joint runtime contract](joint-component-runtime.md);
it keeps the same objective, rank threshold and KKT threshold. The solver
provenance remains v5 and FixedNeighbor policy v4, separate from schemas 3–8
and their stored evidence. See the
[Joint runtime contract](joint-component-runtime.md) for the implementation
details.

## Memory planning and large components

FixedNeighbor preflights its actual component rows, support memberships,
affected rows, neighbor lists and largest block before preparing block views.
The accepted global sweep evaluation uses support iteration and vector
reductions; it avoids a full component design and width-derivative matrix. This
does not make the local profile solve or rank assessment matrix-free. Production
endpoint verification uses O(N+P) streaming passes and does not construct a
global design matrix.

The old per-block row map would use 148,119,642,128 bytes (148.120 GB, about
137.9 GiB) for 5,938,087 rows and 3,118 blocks, before any other data. The
current prepared mapping stores two indices per affected row and a
`LocalSupportIndex` per retained support membership. Its exact reported bytes
depend on `sum(block.affected_rows)` and the actual support memberships and
capacities; observations and atom counts alone cannot determine that number.
The deterministic large chain fixture has 74,811 support memberships and
3,118 blocks; it reports 2,892,736 prepared-mapping bytes. Its global-evaluation
estimated peak is 876,119,666 bytes, and a direct no-sweep test process recorded
474,546,176 bytes peak RSS (0.475 GB, 0.442 GiB). The fixture disables sweeps
and endpoint certification and does not materialize a global dense matrix.

The planner rejects estimated peaks above 3,200,000,000 bytes (3.200 GB,
2.980 GiB), below the nominal 4,000,000,000-byte limit. This is an allocation
estimate with safety margin, not observed peak RSS. It does not establish that a
full physical search fits within 4 GB: real support count and neighborhood shape
are required, and rank and target assessments remain unavailable when their
own stage estimates exceed budget.

The benchmark driver reports the structural census, estimated preparation and
global-state peaks, RSS before and after preparation, RSS after the first sweep
and after certification when those stages run, and process peak RSS. Its
preparation-only (`--maximum-sweeps 0`), one-sweep and bounded multi-sweep
modes are diagnostic harness settings. Aggregate counters are the default;
per-block timings and local profile trials require `--detailed`. No physical
map/model input is present in this workspace, so driver measurements are
synthetic and cannot establish physical full-search feasibility.

The 5,938,087-row, 37,406-atom single-component census continues through
streaming endpoint verification and the combined rank preflight. The synthetic
input has 74,811 support memberships, 3,118 blocks and 37,406 targets. The
all-active free-design estimate is 91,400,872,907 bytes; mixed and all-interior
faces estimate 167,840,033,907 and 266,666,536,155 bytes. All exceed the
3.2 GB planning budget and factorization is `NotRun`. The earlier
189,484,131,163-byte Dense DesignSpectrum and 39,065,926,939-byte projected-
width estimates describe the retired two-stage planner. The current combined
stage accounts for the free-design QR/SVD, projection coefficients and
projected-width QR/SVD. Derived Target Identifiability is estimated at
482,239,947 bytes and has no dense factor workspace. A former global covariance
preflight estimated 418,399,065,443 bytes and rejected before factorization;
that path has been removed. Current local conditional covariance uses an
atom-local two-column design, as documented in the
[memory census](joint-component-memory-census.md). In this shape, structural
rank deficiency leaves no trusted
rank evidence, so target results remain unavailable even though the
evidence-memory preflight passes. See the
[component memory census](joint-component-memory-census.md) for the allocation
breakdown and RSS scope.

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
| Stationarity | Eta-only confirmation was selected with the frozen `1e-10` threshold and requires a previous complete sweep. Current tests retain the KKT `1e-10` and width-gradient `1e-12` gates. |
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
