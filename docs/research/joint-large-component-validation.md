# Joint large-component validation scope

## 6Z6U end-to-end audit

After the P0–P5 commits, configure with
`RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON`, build with at most four workers,
and run the trajectory and halo-width comparison against the supplied CIF,
map, and adjacent `.simulation.json` metadata. For the full accepted-update
and LM-trial trace, use `joint_width_trajectory MODEL MAP OUTPUT all 10000`.
Use a fresh output directory for each experiment.

The instrumented `joint_width_trajectory` executable runs the actual
FixedNeighbor component workflow on the supplied map and model. With the
`all` selector, it records up to 10,000 profile events for all 2,044
contributors. Its summary reports component and target-state availability,
formal runtime convergence, independent replay status, candidate objective,
objective by completed sweep, maximum accepted width and log-width step,
local A/C conditioning, projected A/C KKT, global width gradient, LM
rejections, recovery attempts and successes, elapsed wall and CPU time, and
process peak RSS. It lists every observed accepted width at least 10,000 Å
with its target/halo role, plus the 20 accepted atoms with the largest A/C
cancellation ratios and their local A/C condition indicators. The
profile-event and accepted-update capture completion flags qualify captured
maxima; a capped trace must be treated as a lower bound.

This instrumented runtime includes detailed per-profile diagnostics. The
reported wall time therefore includes diagnostic overhead and is not a
production-runtime benchmark. Component solve preflight is not exposed by the
public result and is reported as unavailable; the tool reports the separate
problem-construction preflight estimate. SparseQR factorization counts are
not exposed by FixedNeighbor telemetry; production high-precision fallback
calls are zero because no such fallback is enabled. The independent endpoint
replay and the saved component state remain the only basis for reporting
formal target estimates.

## Large synthetic controls

The repository contains separate controls for preparation, streaming
verification, and target mapping:

| Registered test | Coverage | Limit |
| --- | --- | --- |
| `joint_component_mapping_memory_scalability` | `ProductionScaleSingleComponentPreparationStaysBelowBudget`: 5,938,087 observation rows, 37,406 contributors, 3,118 blocks; checks memory preflight and bounded prepared mappings, then separately certifies a supplied accepted state | Zero search sweeps; no complete fit |
| `rhbm_tests_joint_component_numerical` | Includes `EndpointReplayUsesStreamingVerificationForLargeComponents` | Streaming endpoint verification on a synthetic 37,406-atom shape; not a 6Z6U fit |
| `rhbm_tests_joint_component_workflow` | `LargeComponentMappingWorkIsLinearInPointsAndMemberships` | Target evidence mapping complexity; does not run FixedNeighbor fitting |

The preparation test's row count is a synthetic dimension, not a physically
allocated 5.9-million-row map. These tests validate memory planning, mapping,
and endpoint paths separately; they do not establish that a full fit at this
scale completes within the resource budget. No global dense design spectrum,
projection matrix, rank QR/SVD, covariance, SPQR, or Operator-PCG path is part
of these controls.

Final verification uses a Release build with offline audits enabled, no more
than four build or test workers, the actual registered CTest names, and a
single full CTest pass after the P0–P5 commits. Full 6Z6U fitting and the
fixed-width model comparisons are separate runtime experiments; neither is
replaced by fixture or synthetic tests.
