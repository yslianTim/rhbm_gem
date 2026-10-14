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
