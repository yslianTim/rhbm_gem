# Joint local A/C KKT fixture

`block-82.json` captures the initial local profile evaluation that failed during
the 6Z6U simulation workflow. It contains only the affected local support
geometry, effective response, initial widths, atom/row identities, and numerical
policy; it does not include the CIF or map.

The captured block has 881 affected rows, an 881 × 24 local design with 3,368
nonzero entries, numerical rank 24, and free-face rank 22. The workflow reports
the failure at sweep 2, block 82 (one-based indexing).

The fixture preserves the pre-fix coefficients and diagnostics. To verify the
current solver against it without rebuilding the full Joint problem, enable the
offline audit target and run `joint_local_ac_kkt_recovery` and
`joint_local_ac_kkt_accuracy`. The fixture records doubles with round-trip
precision.

`../joint_local_replay/block-82-trial.json` captures the rejected Local Profile
LM trial from the same block, including its exact trial A/C/log-B state and the
primary and replay evaluations. It was captured on AppleClang 21.0.0, arm64,
with IEEE-754 binary64 and fast math disabled. The fixture contains only the
local support geometry and effective response. With Joint offline audits
enabled, `joint_local_replay_fixture` reloads it and recomputes both paths
without rerunning 6Z6U. Cross-platform floating-point differences are reported
by the runner and are not treated as a bitwise contract.
