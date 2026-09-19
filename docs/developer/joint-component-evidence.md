# Joint component evidence and retired experiments

The maintained entry point is the [runtime guide](joint-component-runtime.md).
Historical comparison runners were retired after the transition acceptance;
the numerical and observation contracts remain regression responsibilities.

## Retained conclusions and controls

| Research group | Conclusion and retained protection |
| --- | --- |
| Observation matching and estimated-neighbor sweep | Observation/sampling geometry and frozen support matter. Keep sampler, strict-boundary, contributor and post-sum quantization contracts. |
| Matched joint A/C and Experiments A/B/C | These comparison workflows are not prerequisites for joint LS. Retain their provenance and observation contracts. The composite MDPDE negative results motivate future robustness research; they do not disprove heterogeneous MDPDE. |
| Fixed-B and initial joint-ABC profile | The constrained linear subproblem remains part of variable projection. Keep independent double/float32 recovery and full profile Jacobian checks at nonzero residual. |
| Coverage | Nine immutable datasets survive in the fixture catalog. Weak, active-face, zero and duplicate cases retain their original limitations rather than becoming success cases. |
| Certification/stabilization | Guarded reference/replay prevents acceptance of untrusted trials. Keep near-0.02/narrower in both precisions; derivatives and boundary certification have separate offline tools. |
| Exact components, local certification and extraction | Preserve structural partition, same-state and actual-state assembly, shared parent normalization/rank context, local evidence and failure isolation. Local certification never overwrites historical global qualification. |

Selected historical counterexamples are retained here so that retirement does
not erase the limits of the comparisons:

| Frozen comparison | Result and limitation |
| --- | --- |
| Observation-matched, truth neighbors | Compact prediction versus sampler RMSE was `5.17001e-16`; float32 residual was `2.91318e-8`. This isolates observation consistency, not unknown-neighbor convergence. |
| Estimated-neighbor sweep, eight checkpoints | Aggregate A/B/C errors improved by approximately 19–21 / 25–26 / 4.7–4.9 times, but C worsened for 70–76 of 168 atoms. Aggregate improvement is not per-atom recovery. |
| Matched joint A/C, eight checkpoints and two starts | All 16 joint runs exhausted 2,000 evaluations without qualification, despite full design rank 336. The frozen-neighbor controls qualified 1,335/1,344 blocks. Rank and lower C error do not establish joint stationarity. |
| Experiment A, baseline-best28 | At alpha=0, A/C RMSE was `0.011471954 / 0.00063732864`; alpha=0.1 reached `0.0091138914 / 0.00049314988` with qualification. Alpha=0.5 and 1 exhausted their budgets. |
| Experiment B, baseline-best28 | Alpha=0 qualified at A/C RMSE `0.0091707311 / 0.00032313699`; the positive-alpha runs exhausted 100 evaluations without qualification. |
| Experiment C, two starts | Both exhausted 100 evaluations without qualification; there was no independent repeat. One endpoint had A/C RMSE `0.017916446 / 0.00011271024` and effective sample fraction about `0.8983%`. Without a common-alpha control, this does not isolate the cause or disprove heterogeneous MDPDE. |

These values describe the frozen historical reports at the final baseline
commit, not newly measured runtime benchmarks. The inventory below identifies
the corresponding original report paths and hashes.

The unchanged baseline recorded 128 paired branches, 64 required global regular
pairs and 192 component-local scopes, of which 128 were regular. The 96 regular
globally restricted scopes use a separate certification scope. These are numerical
qualifications for the recorded fixtures and policy, not cross-platform or
statistical generalization claims.

The [fixture catalog](../../tests/fixtures/joint_component/catalog.json) binds
inputs and compact expected states to original archive members and hashes.
The [certification contract](joint_abc_certification_contract.md) and
[component contract](joint_abc_components_contract.md) retain the mathematical
thresholds and historical scope definitions. Current public objectives are
normalized by the parent observation scale squared; archived experiment
objectives remain raw half-RSS.

## Transition acceptance

The [acceptance receipt](figures/joint-component-evidence/transition-acceptance.json)
records exact scientific comparisons after excluding timing, peak RSS and process
completion manifests:

- All 1,248 migrated search/state/context records match the pre-extraction oracle.
- All 2,096 historical global/restricted audit records match their frozen oracle.
- All 3,712 baseline audit records, including separate local scopes, match; local
  results never replace global/restricted evidence.

The complete fresh matrix used the frozen baseline executable. The migrated
search driver reproduced its states and contexts exactly, and the migrated
offline tools separately reran the four representative local certificates.
The receipt distinguishes these runs; it does not claim that the migrated
offline driver reran the entire audit matrix.

Per-case results remain reviewable in the compact [global pair table](figures/joint-component-evidence/global-pairs.csv),
[global/restricted failure matrix](figures/joint-component-evidence/global-and-restricted-scopes.csv)
and [local scope table](figures/joint-component-evidence/local-scopes.csv).
Public objective scaling is checked separately against nonzero residuals and
constant rows; historical raw half-RSS values were not rewritten.

The [post-retirement validation receipt](figures/joint-component-evidence/validation.json)
records 22 passing default CTests and 3 passing extended/offline CTests in a
clean source copy without Git metadata, historical archives or external original
inputs. It also records repository guards, testing-disabled installation,
consumer execution and retained `solve`/`refine` entry-point checks. Validation
was performed on the recorded macOS/compiler configuration, not across platforms.

## Historical retrieval

Complete reports, source snapshots, matrices and traces remain in Git at:

- Final runtime baseline: `ba2449faf1e1b163f1b4209c3ef386ff6516de43`.
- Pre-extraction component reference: `6f30510c06def630226a3dab25e9bcc579d72687`.

Use an isolated checkout of the full commit to reproduce a historical runner.
To retrieve one archived artifact without restoring the old working tree:

```sh
git show ba2449faf1e1b163f1b4209c3ef386ff6516de43:docs/developer/figures/joint-abc-components/search-records-a.tar.gz > /tmp/search-records-a.tar.gz
```

The [retired-artifact inventory](figures/joint-component-evidence/retired-artifacts.json) records the exact paths, sizes and SHA-256
hashes. Verify the retrieved bytes against that inventory. No Git history is
rewritten. Deleting the matrices reduces the current checkout; it does not
remove their objects from repository history.

The older certification report referred to `scientific-records.tar.gz` and
`frozen-diagnostics/frozen-endpoint-records.tar.gz`, which are absent from the
final baseline tree. They are not required inputs or promised recoverable
artifacts. The retained pre-extraction reference is the actual component
archive identified in the catalog.
