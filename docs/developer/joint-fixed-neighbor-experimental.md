# Fixed-neighbor joint search: experimental evidence

FixedNeighbor remains an internal experimental search route. It is not part of the production main flow and does not change the default `LegacyCompact` method.

## Fixed-B qualification

The cube-256 fixed-B plateau was classified as a global roundoff lock. The local conditional solves continued to return valid certificates and nonzero updates, but the global replay could increase by about `9.91e-16`, within the existing objective replay enclosure of about `1e-12`. The existing strict global acceptance test rejected those updates, leaving the sweep KKT at `7.0516e-9`.

The selected repair makes block acceptance roundoff-aware only when the local conditional objective does not increase, the local/global delta identity is enclosed by the replay error, and any positive global replay delta is within the existing replay enclosure. The global KKT target remains `1e-10`.

After that repair, fixed-B chain/cube controls passed at 256 atoms in both orders. Forward scaling remained stable-looking through 1024 atoms: chain required `5, 5, 5` sweeps at 256/512/1024; cube required `9, 11, 11`. The cube trend has only three multi-block sizes, so the formal four-point growth gate was not evaluated. The 1024 cube run reached about 1.48 GB peak RSS.

Detailed evidence is in:

- [fixed-B floor attribution](figures/joint-fixed-b-floor-n1/analysis.json)
- [fixed-B coupling diagnostics](figures/joint-fixed-b-coupling-n1/coupling.json)
- [fixed-B requalification](figures/joint-fixed-b-requal-r1/analysis.json)
- [fixed-B scaling](figures/joint-fixed-b-scaling-r1/analysis.json)

## Nonlinear endpoint experiment

The nonlinear experiment used disjoint 128-atom cores, serial Gauss-Seidel updates, local `LegacyCompact` width search, and the existing global replay and endpoint certification. Global `LegacyCompact` and `OperatorPcg` returned runtime-converged references for all six cases.

| Case | Sweeps | FixedNeighbor objective | Global A/C KKT | Width gradient infinity norm | Scaled A/C difference vs LegacyCompact | Endpoint | Runtime |
|---|---:|---:|---:|---:|---:|---|---|
| chain-32 | 1 | `1.5342583e-10` | `1.80e-16` | `2.11e-17` | `2.18e-16` | passed | passed |
| cube-32 | 1 | `9.8793029e-11` | `1.25e-16` | `2.19e-17` | `1.51e-14` | passed | passed |
| chain-128 | 1 | `1.4868379e-10` | `9.71e-17` | `7.46e-18` | `1.21e-16` | passed | passed |
| cube-128 | 1 | `8.7674368e-11` | `1.11e-16` | `1.56e-16` | `1.47e-16` | passed | passed |
| chain-256 | 6 | `1.4847246e-10` | `3.99e-12` | `1.60e-14` | `7.62e-12` | passed | passed |
| cube-256 | 11 | `8.2223286e-11` | `5.66e-11` | `2.19e-13` | `1.38e-10` | **failed** | **failed** |

All six FixedNeighbor objectives were within the experiment's objective parity tolerance. Cube-256 also met the global A/C KKT target, but the existing endpoint checks failed: inner coefficient difference was `1.05e-9` against `1e-10`, width stationarity was `1.15e-12` against `1e-12`, and local correction was `1.24e-9` against `1e-10`. Its endpoint trust result was `reference-disagreement`. Objective parity and KKT alone do not certify this endpoint.

An exploratory continuation past the initial stationary point reduced cube-256 KKT to `1.97e-13` by sweep 14, then produced the same KKT at sweep 15. Endpoint certification had still not passed through sweep 14. This did not identify a supported repair, so no sweep-budget or certification semantics were changed.

The F2 endpoint gate therefore remains failed. Fixed-neighbor large-component nonlinear scaling, order sensitivity, and route promotion were skipped. Coordinated/shared-parameter blocks remain deferred: the observed cube-256 failure is an endpoint-assessment mismatch, not evidence of sweep growth, converged order dependence, or boundary oscillation.
