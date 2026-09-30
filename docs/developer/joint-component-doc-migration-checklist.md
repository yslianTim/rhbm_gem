# Joint document migration checklist

This checklist records ownership migration for possible PR-H3 cleanup.
No document or evidence artifact is deleted by PR-H1 or PR-H2. A historical
document is a deletion candidate only when its unique current contract and
historical evidence have owners and its remaining unique information is none.

| Historical document | Unique current contract moved to | Historical evidence moved to | Remaining unique info? |
|---|---|---|---|
| [joint-abc-profile-experiment.md](joint-abc-profile-experiment.md) | [Runtime guide](joint-component-runtime.md) and [component contract](joint_abc_components_contract.md) | [Canonical evidence](joint-component-evidence.md), fixture catalog | none |
| [joint-abc-coverage-experiment.md](joint-abc-coverage-experiment.md) | [Runtime guide](joint-component-runtime.md) and [component contract](joint_abc_components_contract.md) | [Canonical evidence](joint-component-evidence.md), fixture catalog | none |
| [joint-abc-certification-experiment.md](joint-abc-certification-experiment.md) | [Certification contract](joint_abc_certification_contract.md) | [Canonical evidence](joint-component-evidence.md), transition receipts | none |
| [joint-abc-components-experiment.md](joint-abc-components-experiment.md) | [Component contract](joint_abc_components_contract.md) | [Canonical evidence](joint-component-evidence.md), transition receipts | none |
| [joint-bounded-search-rank.md](joint-bounded-search-rank.md) | [Sparse backend](joint-component-sparse-backend.md), [operator search](joint-operator-search.md) | [Canonical evidence](joint-component-evidence.md), [validation report](joint-bounded-search-rank-validation.md) | none |
| [joint-initial-rank-diagnostic.md](joint-initial-rank-diagnostic.md) | [Observable-halo guide](joint-observable-halo.md) | [Canonical evidence](joint-component-evidence.md) | none |
| [joint-component-v1-acceptance.md](joint-component-v1-acceptance.md) | [Runtime guide](joint-component-runtime.md), formal contracts | [Canonical evidence](joint-component-evidence.md), v1 receipts | yes — exact lane configuration and acceptance record |
| [joint-component-partial-selection-acceptance.md](joint-component-partial-selection-acceptance.md) | [Runtime guide](joint-component-runtime.md) | [Canonical evidence](joint-component-evidence.md), acceptance receipt | yes — complete control and resource-run details |
| [joint-component-compact-svd-acceptance.md](joint-component-compact-svd-acceptance.md) | [Compact SVD guide](joint-component-compact-svd.md) | [Canonical evidence](joint-component-evidence.md), compact receipts | yes — exact historical case table and artifact checks |
| [joint-component-sparse-acceptance.md](joint-component-sparse-acceptance.md) | [Sparse backend](joint-component-sparse-backend.md) | [Canonical evidence](joint-component-evidence.md), sparse verification files | yes — case-level measurements and verification record |
| [joint-component-reference-acceptance.md](joint-component-reference-acceptance.md) | [Compact SVD guide](joint-component-compact-svd.md) | [Canonical evidence](joint-component-evidence.md), archive inventory | yes — per-case timings and archive status |
| [joint-fixed-actions-validation.md](joint-fixed-actions-validation.md) | [Fixed-actions guide](joint-fixed-actions.md) | [Canonical evidence](joint-component-evidence.md), comparison receipts | yes — detailed A/B/C run matrix |
| [joint-operator-search-validation.md](joint-operator-search-validation.md) | [Operator-search guide](joint-operator-search.md) | [Canonical evidence](joint-component-evidence.md), campaign manifest | yes — run-level deadline and comparison record |
| [joint-bounded-search-rank-validation.md](joint-bounded-search-rank-validation.md) | [Sparse backend](joint-component-sparse-backend.md) | [Canonical evidence](joint-component-evidence.md) | yes — detailed test and prototype outcomes |
| [joint-component-tiled-backend.md](joint-component-tiled-backend.md) | [Profile-operator guide](joint-profile-operator.md) | [Canonical evidence](joint-component-evidence.md), measurement JSON | yes — case-level measurements and parity record |
| [joint-observable-acceptance.md](joint-observable-acceptance.md) | [Observable-halo guide](joint-observable-halo.md) | [Canonical evidence](joint-component-evidence.md) | yes — detailed incomplete 6Z6U and fixed-case results |
| [joint-analysis-consolidation.md](joint-analysis-consolidation.md) | [Stage integration](joint-stage-integration.md), [data I/O architecture](architecture/dataobject-io-architecture.md) | [Canonical evidence](joint-component-evidence.md), [consolidation benchmark record](joint-analysis-consolidation-benchmark.json) | yes — migration sequence and detailed acceptance record |
| [joint-command-resource-envelope.md](joint-command-resource-envelope.md) | [Benchmark guide](joint-benchmark.md) | [Canonical evidence](joint-component-evidence.md), resource manifest | yes — exact input hashes and phase records |

The following remain current owners or active research guides and are not
PR-H3 deletion candidates: the runtime, benchmark, sparse backend, compact SVD,
fixed-actions, profile-operator, operator-search, observable-halo,
target-estimability, and stage-integration guides; both formal mathematical
contracts; and the weak-halo and noise/mismatch guides.

Before a PR-H3 deletion, update any inbound links to a retained owner and remove
only files whose row is marked none. Machine-readable receipts, fixtures, and
inventories require their own row-level evidence review; this checklist does
not authorize their deletion.
