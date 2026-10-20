# P4 Assessment: Global FirstStage-Calibrated Width Bounds

> **Historical decision record:** This assessment predates the 2026-10-20 production promotion. The A/B2/B4 comparison executable and paired runner are retired; its results below remain evidence for the decision history. See [the production promotion record](joint-bounded-active-set-production-promotion-20261020.md) for current policy.

## Decision

**Insufficient Evidence for production promotion.** Keep the production numerical default unchanged. The research-only B2/B4 bounded solver passed its active-bound KKT and replay checks on one unsharpened 6Z6U map, but the prescribed synthetic campaign could not produce contract-consistent FirstStage donors, and the truth-bearing simulation did not activate either bound.

This is not evidence that the bounded method is numerically invalid. It is evidence that its accuracy and generality have not yet been established well enough to change production behavior.

## Evidence by criterion

| Criterion | Finding | Consequence |
|---|---|---|
| Calibration provenance | P95 uses only qualified FirstStage target widths. B2 and B4 equal 2×P95 and 4×P95. Truth leakage checks pass on both real maps. | Calibration contract is implemented and audited for these inputs. |
| Bounded solver and endpoint evidence | On the unsharpened map, B2 and B4 hit the cap, passed independent replay, had zero feasibility/sign/complementarity violations and zero failed active-set subproblems, and passed projected-width stationarity below 1e-9. | Supports numerical correctness for this bounded case. |
| Truth-bearing simulation | A/B2/B4 all qualified and achieved 100% truth coverage, but neither cap became active and all three returned the same serial-1032 width to the reported precision. | Compatibility evidence; no measurable cap effect. |
| Unsharpened real map | A stopped with recovery exhausted, width 6,055.801 Å, and failed replay. B2/B4 limited that width to 1.604218/3.208436 Å and qualified numerically. The map has no truth, and residual RMSE/objectives were nearly unchanged. | Shows runaway prevention, not parameter accuracy or general benefit. |
| Synthetic validation | The 30-seed P3a campaign supplied A-only runs. The generator lacks FirstStage sampling observations, MDPDE donor fits, and qualified donor provenance, so B2/B4 were correctly not run. | The planned multi-seed truth-controlled bounded comparison is still missing. |
| Production scope | The upper-bound control and paired runner are test-only. The production library builds with the test instrumentation disabled; no production default or public CLI changed. | The experiment can remain isolated while evidence is incomplete. |

## Interpretation

On the real simulation map, B2/B4 were inactive, so the near-identical parameter RMSEs cannot distinguish bounded from unbounded fitting. On the unsharpened map, the caps prevented extreme widths and the bounded endpoints passed the numerical gates after the assembled assessment was given the same research bound as the component solver. The fits still have no truth-based accuracy score, and their observed-data objective and residual RMSE are close to the unqualified A candidate. The bounded behavior is promising as a stability intervention, but these data do not establish a production accuracy improvement.

The initial unsharpened pair audit failed only because its global assembled stationarity check used the unbounded replay context. The check was not relaxed: the test-only research bound was propagated into assembled evidence, a regression was added, and A/B2/B4 were rerun. Both the initial failure record and the corrected passing audit are preserved in the [P3b report](joint-global-width-bounds-real-20261019/README.md).

## Evidence needed before reconsidering promotion

The next decision should wait for a truth-controlled synthetic path that creates FirstStage samples and qualified MDPDE donors under the same calibration contract. It should run paired A/B2/B4 over multiple seeds and include cases where the bound is active, then compare truth-based parameter error and prediction fit while retaining the existing replay, stationarity, and feasibility gates. The current synthetic generator limitation and P3a results are documented in the [synthetic report](joint-global-width-bounds-synthetic-20261019/README.md).

## Supporting records

- [P3b real-map report and paired summaries](joint-global-width-bounds-real-20261019/README.md)
- [Final real-map machine-readable summary](joint-global-width-bounds-real-20261019/paired-summary-final.json)
- [P3a synthetic report](joint-global-width-bounds-synthetic-20261019/README.md)
