# Fixed-target Halo FullABC Joint Estimation — P2 implementation

Date: 2026-10-20  
Branch: `develop`  
Status: research implementation; Production remains `geometry-local-target-median-v1`.

## Contract

`EstimateHaloFullABCJoint` accepts only a qualified Pure Halo-Support Exclusion Target endpoint. It reads the immutable `JointProblem::OriginalInput()` row set and the captured Target component states. Target A/B/C are never inserted into the second solver parameter vector.

Each Halo observation count is the number of unique OriginalInput support row IDs. Halos with fewer than three rows receive the componentwise median of qualified Joint Target A/B/C and stay fixed. Every other Halo enters the FullABC solver; the classification path has no rank, condition, label, residue, or group pre-filter.

The second response starts from OriginalInput observations and subtracts the frozen Target endpoint and every fixed Halo profile. The optimizer receives every OriginalInput row, including rows outside all free Halo supports, and only sparse free-Halo support incidence. The existing OriginalInput response scale is reused, so zero-sensitivity rows retain their residual-constant contribution to objective and replay.

For free Halo widths, the implementation calls the existing calibration helper on the captured first-stage Target initialization. It requires its qualified-donor coverage and computes `Bmax = 4 * Q95`; it does not construct Halo FirstStage estimates. Free Halo A/C are profiled by the existing constrained solve and widths use the FixedNeighbor Active-set LM path, including the existing endpoint and independent replay checks.

The returned in-memory result separates each Halo's fixed/free role and observation count from its reconstructed A/B/C. A failed free-Halo solve leaves the reconstruction unavailable. An all-fixed case returns its Global Median profiles without invoking the solver. No Halo estimates are written to the Production Joint snapshot or SQLite contract.

## Peeling integration

`JointPostprocessingContext` has an explicit optional research-result argument. When supplied, it validates run provenance, complete model Halo coverage, the `n_h < 3` fixed rule, frozen Global Median values, width bound, and free-solver replay evidence. It supplies those profiles to both OriginalInput contributors and Reconstruction-only Peeling nodes. Calls that omit the research result retain the existing Peeling path and persisted snapshot semantics.

## Targeted evidence

Built with:

```text
cmake --build build/debug --target rhbm_tests --parallel 4
```

Passed:

```text
build/debug/bin/RHBM-GEM-TEST --gtest_filter='HaloFullABCJointEstimationTest.*:PureHaloSupportExclusionTest.ExplicitFullABCHaloProfilesReachReconstructionOnlyPeelingNodes'
  6 tests passed.
ctest --test-dir build/debug -R '^rhbm_tests_joint_component_contract$' --output-on-failure -j1
  passed; 26.93 seconds.
```

The synthetic solver regression used 20 frozen Targets, 32 OriginalInput rows, one free Halo with 30 distinct rows, one two-row fixed Halo, and one zero-row reconstruction-only Halo. The free and fixed Halo supports shared rows. The fitted free A/B/C matched the synthetic endpoint within `2e-3`; two zero-sensitivity rows retained their residual contribution to the reported objective. Additional cases verified all-fixed no-solver behavior, exactly-three-row classification as free, missing qualified donors as unavailable, no-Halo no-solver behavior, chemical-label invariance, and passing profiles at reconstruction-only Peeling nodes.

This is focused algorithmic evidence only. It does not measure 6Z6U estimate error, downstream Local/Group quality, large-component runtime, or process RSS; those remain P3/P4 work and no Production promotion is implied.
