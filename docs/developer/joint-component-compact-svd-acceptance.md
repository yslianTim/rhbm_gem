# Joint compact SVD acceptance

The historical compact campaign has been verified and passes its numerical and
performance gates. Its baseline is `ce58c89747d4091f91967e7f4bd9c7890d51c203`,
production fingerprint
`9718531067e72cd3ecf180d12ae8033b8e39360b7fe1d49d24ae847de4240fa0`.
The measured candidate matches `e32919f3`, fingerprint
`f1ac45c36945e8e89e588af5c305158a7723b268a21c5c7ab0b22a3ae676a755`.
Git source proofs, historical harness hashes, raw samples, inputs, matrices,
audits, replay and completed SQLite/JSON/CSV exports were checked. Reaggregation
also rechecks the exported objective directly, since it is already normalized.

| Historical SPQR fixed-state metric | Auto / legacy median time | Required maximum |
|---|---:|---:|
| Single 128 reference + derivative | 0.268861 | 1.10 |
| Heterogeneous 168 reference + derivative | 0.730425 | 1.10 |
| Single 512 reference SVD | 0.009068 | 0.70 |
| Single 512 free-design SVD | 0.005008 | 0.70 |
| Single 512 reference + derivative | 0.184551 | 0.70 |

All three modes have three fixed-state samples for both sparse backends, and
all numerical audit/replay comparisons pass. Historical Single 128 commands
pass numerical, convergence and persistence/export checks. The historical
Single 512 candidate completed three times (median 289.238 seconds), with
runtime convergence and exports. Its `ce58c897` command baseline stopped at
600 seconds; no completed-command speedup is assigned to that censored run.

These historical results are distinct from the latest endpoint-reference
simplification. Its independently rebuilt `e32919f3` baseline and `c9e0f8c9`
candidate have fresh complete-command measurements in the
[endpoint-reference acceptance](joint-component-reference-acceptance.md).
That record also supplies the latest fixed-state numerical controls.

The [implementation guide](joint-component-compact-svd.md) describes unchanged
rank policies and reproduction. The [evidence index](figures/joint-reference-acceptance/README.md)
contains the immutable historical campaign together with the new evidence and
relocatable reaggregation instructions. The earlier
[verification record](figures/joint-compact-acceptance/verification.json) retains
the original build/test logs; the new acceptance record contains the latest
regression results.
