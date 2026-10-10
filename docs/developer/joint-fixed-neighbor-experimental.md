# FixedNeighbor Joint: current contract and evidence

This is the concise overview for the production FixedNeighbor route. Historical
execution telemetry and retired campaign implementations remain recoverable
from Git history; they are not part of the current production surface.

## Current production contract

`FixedNeighbor` is the sole Joint estimator. It is selected by
`--estimator joint-components`; there is no separate global search-route
option or automatic routing layer.

| Policy | Current value |
| --- | --- |
| Outer core | 12 atoms, disjoint mutable blocks |
| Block order | Forward serial Gauss-Seidel |
| Maximum sweeps | 30 |
| Local search | `LegacyCompact` |
| Local update contract | At most one trusted accepted `LegacyCompact` update per block visit |
| Acceptance | Existing replay and roundoff-safe monotone acceptance |
| Convergence | Existing A/C KKT, width-gradient and eta-confirmation checks |
| Endpoint | Existing trust, certification and `RuntimeConvergence` semantics |

The numerical thresholds are unchanged: A/C KKT `1e-10`, width-gradient
`1e-12`, and eta confirmation `1e-10`. Prepared structural block mappings and
workspace reuse do not reuse numeric factors, accepted states, or replay
results across sweeps. The route therefore preserves the objective, A/C solve,
width parameterization, serial update semantics, replay, endpoint
certification, persistence and `RuntimeConvergence` contract. See the
[Joint runtime contract](joint-component-runtime.md) for the implementation
details.

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
The EIGEN sparse backend remains current because FixedNeighbor and its local
LegacyCompact profile search still use it. There is no hidden FixedNeighbor
policy switch. Closed experiment implementations and machine-readable receipts
are recoverable from Git history and are intentionally not retained in the
current tree.
