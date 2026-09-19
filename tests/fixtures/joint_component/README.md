# Joint component runtime fixture catalog

`catalog.json` identifies nine immutable observation/support snapshots, their
source commits, archive members and SHA-256 hashes. Each package contains the
existing `dataset.json`, `snapshot.json`, `voxels.csv`, `contributors.csv` and a
compact `cases.json`. The latter retains the eight recorded initial-width
vectors, selected endpoint/search assertions and historical local qualification
checks. It does not regenerate an oracle using the implementation under test.

Search references come from the pre-extraction exact-component archive at
`6f30510c06def630226a3dab25e9bcc579d72687`. Component-local certificate references
come from `ba2449faf1e1b163f1b4209c3ef386ff6516de43`; they supplement, never replace,
historical global or globally restricted qualifications.

The Python runtime runner verifies the package and member hashes before
extracting into its build/work directory. No external CIF/MRC or historical
experiment archive is required. Numeric states use scaled `1e-10` comparison;
accepted/rejected decisions, stopping reasons, ranks, counts and availability
must match exactly. Independent scalar forward/KKT/gradient replay supplies an
additional endpoint check.

Default cases are baseline first-stage in both precisions, near-0.02/narrower in
both precisions, and weak-1e-4, active-a, zero-signal and duplicate first-stage
double. Extended tests add heterogeneous-168/first-stage-float32, weak-1e-2 and
near-0.10 first-stage double. `--dataset NAME --all-starts` selects the original
four starts in both precisions when initialization changes require that check.

The selected historical **local** certificate expectations are:

| Dataset / representative start | Expected evidence or limitation |
| --- | --- |
| heterogeneous-168 / first-stage float32 | Regular in the frozen local scope; extended lane only. |
| baseline / first-stage double | Regular; independent recovery and quantization controls also apply. |
| weak-1e-2 / first-stage double | Regular for this recorded case; no blanket weak-signal claim. |
| weak-1e-4 / first-stage double | `resolution-unverified` derivative; not regular even though numerical rank passes. |
| active-a / first-stage double | Local correction fails and derivative is `transition-unverified`; not regular. |
| zero-signal / first-stage double | Width identification and derivative fail despite zero residual. |
| near-0.10 / first-stage double | Regular in the frozen local scope. |
| near-0.02 / narrower double and float32 | Guarded must reject untrusted trials and retain the trusted state; selected local scopes are regular. |
| duplicate / first-stage double | No usable state after rank failure; endpoint evidence remains unavailable. |

These frozen expectations do not make the runtime return an offline certificate;
its offline evidence remains `NotRun`. `cases.json` retains exact per-case checks
and failure/availability lists for the other stored starts and precisions.

See the [runtime guide](../../../docs/developer/joint-component-runtime.md) and
[evidence index](../../../docs/developer/joint-component-evidence.md).

`simulation-contract.json` preserves the original coverage input hashes and
heterogeneous v1 width contract for `simulation_contract.resolve`; it is not
a dependency on the retired coverage registry or runner.
