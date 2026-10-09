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

## Current qualification evidence

| Evidence | Decision | Canonical record |
| --- | --- | --- |
| Stationarity confirmation | Eta-only confirmation is safe under the existing `1e-10` change threshold; confirmed endpoints pass the existing gates through the qualified cases | [stationarity summary](figures/joint-fixed-neighbor-stationarity-r1/README.md) |
| Historical local-update selection | Full and TwoAccepted were evaluated as historical alternatives; OneAccepted was selected and is now intrinsic to FixedNeighbor. All four 256/512 full-endpoint cases passed, while the 1024 comparison remains search-only | [local-work qualification](figures/joint-fixed-neighbor-inexact-qualification-r2/README.md) |
| Outer-core C3 endpoint gate | Cores 12, 16 and historical control 64 passed endpoint, trust, `RuntimeConvergence`, replay and parity checks on chain/cube 256/512 | [outer-core qualification](joint-fixed-neighbor-outer-core-qualification.md) |
| Outer-core C4 repeated frontier | Core 12 passed the correctness gate and exceeded the 10% search-improvement requirement on every measured topology/size | [C4 summary](figures/joint-fixed-neighbor-core-size-r5/README.md) |
| Prepared-block behavior | Structural preparation and mapping counts are checked separately from numerical solving; no new convergence criterion was introduced | [prepared-block qualification](joint-benchmark.md) |
| Forward/reverse order | All five pairs confirmed convergence; maximum eta difference was `3.75e-11` and maximum scaled A/C difference was `2.09e-12`; full endpoint/runtime parity passed through 512 | [order summary](figures/joint-fixed-neighbor-order-r1/README.md) |
| Fixed-B replay floor | The roundoff-aware acceptance contract and its regression evidence remain unchanged | [Fixed-B requalification](figures/joint-fixed-b-requal-r1/analysis.json) |

The C4 measurements that explain the adopted core are:

| Case | Core 12 median search | Historical control 64 | Improvement |
| --- | ---: | ---: | ---: |
| chain-512 | 61.18 s | 171.26 s | 64.3% |
| cube-512 | 151.02 s | 322.81 s | 53.2% |
| chain-1024 | 167.73 s | 346.28 s | 51.6% |
| cube-1024 | 490.85 s | 732.64 s | 33.0% |

Every measured run passed search/replay correctness. The 64-atom values are
historical controls, not the current production default.

The former OneAccepted/core64 optimized frontier is superseded by the
qualified core12 production policy. Its compact comparison is retained as
historical evidence; it is not an active route or benchmark default.

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
policy switch. The old campaign wrappers, per-run JSON, progress files and
process logs are not required to interpret the current route. The current repository keeps canonical
README/analysis/summary/manifest evidence and the focused regression contracts;
raw execution data is historical and recoverable from Git history when needed.
