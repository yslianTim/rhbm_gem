# Pure Exclusion Baseline Before FixedHalo Retirement

This is the pre-refactor baseline for removing the FixedHalo execution path.
It is anchored to `develop` at `cddf3f990331d80bcb87e9e04dca6ce5cc43fe53`
(2026-10-20). The local branch and the GitHub `develop` ref both resolved to
that commit before edits began.

## Deterministic geometry baseline

`PureHaloSupportExclusionTest.ExcludesExactBoundarySupportAndPreservesOriginalDomain`
uses a 1-D map whose values are 1 through 11, target atom `1` at x=0, and
non-target atom `2` at x=5. The map spacing is 0.5 Å, so the non-target has
one support incidence exactly at 2.5 Å (`r²=6.25 Å²`). That boundary row is
excluded.

| Field | Baseline |
| --- | --- |
| Selected target IDs | `["1"]` |
| Original contributors | `["1", "2"]` |
| Halo candidate / retained / zero-support counts | `1 / 1 / 0` |
| Original row IDs | `"0"` through `"10"` |
| Fitting row IDs | `"0"` through `"9"` |
| Excluded row identity | `["10"]` |
| Original / fitting / excluded rows | `11 / 10 / 1` |
| Fitting-to-original contributors | `[0]` |
| Fitting-to-original rows | `[0,1,2,3,4,5,6,7,8,9]` |
| Target support rows in fitting domain | `[0,1,2,3,4,5,6,7,8,9]` |
| Observation scale | `sqrt(506) = 22.494443758403985` |

The target's original support is rows 0–10, with squared distances
`[6.25,4,2.25,1,0.25,0,0.25,1,2.25,4,6.25]`. The halo support is row 10
only. The test checks exact row identities, provenance counts, mappings,
boundary inclusion, and that the scale is computed from all original
observations.

## Production and downstream baseline owners

The existing deterministic tests are the frozen contracts used for parity
checks after geometry decoupling:

| Contract | Test owner |
| --- | --- |
| Boundary inclusion, exact exclusion, original/fitting rows and mappings | `PureHaloSupportExclusionTest.ExcludesExactBoundarySupportAndPreservesOriginalDomain`; `CapturesOriginalAndFittingIdentitiesSeparately` |
| Multiple Halo/Target support overlap, zero-support Halo, empty-domain and target-coverage rejection, deterministic order | Remaining construction cases in `PureHaloSupportExclusion_test.cpp` and `HaloObservationExclusion_test.cpp` |
| Original-domain scale and constant-row objective behavior | `JointComponentRuntimeTest.PublicObjectivesUseParentScaleAndIncludeConstantRows` |
| Calibration, bounded search, replay, KKT and numerical behavior | `rhbm_tests_joint_component_numerical`; `joint_production_bounded_synthetic` |
| Halo Median Reconstruction, peeling, local fitting and Group Fitting | `PureHaloSupportExclusionTest.RunsFullPureWorkflowWithGroupAndGlobalHaloSources`; `RunsPureWorkflowWithoutHalos`; `PeelingUsesOriginalRowsAndReconstructedHalo` |
| Persistence and reload | `PureHaloSupportExclusionTest.PersistsVersionedDualDomainAndMedianProvenance`; `joint_persistence_contract_test` |
| Production API and CLI outcomes | `PureHaloSupportExclusionTest.ProductionWorkflowUsesBoundedPureExclusionAndStopsWhenCalibrationIsUnavailable`; `joint_workflow_cli_smoke_test`; `rhbm_tests_integration_command` |

The current checked-in 6Z6U artifacts remain historical inputs, not fresh
measurements for this baseline: the geometry artifact records source HEAD
`5de15aa10e6f08a3d65aa55d110d53df59520577`, while the production workflow
report records earlier implementation commits. The 6Z6U model/map inputs are
not present in this checkout, so no result has been attributed to a rerun at
the baseline HEAD. The 37,406-target, 5,938,087-observation full-fit 4 GB RSS
gate was not run and remains open.

## Baseline verification

At the baseline HEAD, `cmake --build build --parallel 4` succeeded. With
Eigen, OpenMP, OpenBLAS, and Accelerate thread counts set to one, these nine
CTest entries passed with CTest parallelism 4:

```text
joint_production_bounded_synthetic
rhbm_tests_joint_component_contract
rhbm_tests_joint_component_workflow
rhbm_tests_joint_component_numerical
joint_persistence_contract_test
rhbm_tests_integration_command
joint_workflow_cli_smoke_test
joint_component_regression
joint_component_physical_smoke
```

The CTest run reported `100% tests passed (9/9)` in 43.09 seconds. These
results establish the pre-change local test baseline; exact floating-point
parity remains governed by each test's existing tolerances.
