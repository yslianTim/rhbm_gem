# Fixed-target Halo FullABC P0 numerical audit (2026-10-21)

## Scope and reproduction

This audit re-ran `three_observation_free_stress` and `four_observation_free_stress` from the P3 validation executable after adding failure-only component conditioning and block-stop diagnostics. Both runs used one OpenMP/BLAS thread and wrote their source data under `docs/research/results/halo_fullabc_p0/`. The result contract is `fixed-target-observation-driven-halo-fullabc-p3-v2-nearest-target-distance`.

In both cases, the first-round Target Joint endpoint qualified. Halo 103 had exactly three or four distinct `OriginalInput` rows, was classified `FreeFullABC`, and formed a one-Halo observation component with no shared rows. The other free component, rooted at Halo 101, qualified independently in 48 sweeps. The final FullABC reconstruction remained unqualified because Halo 103 did not converge.

## Failure evidence

| Measurement | 3 observations | 4 observations |
|---|---:|---:|
| Halo 103 distinct rows / local effective rows | 3 / 3 | 4 / 4 |
| Free A/C columns / numerical rank | 2 / 2 | 2 / 2 |
| A/C condition number | 14.6819 | 13.6466 |
| A/C rank threshold | 9.40e-16 | 1.25e-15 |
| Projected width rank | 1 / 1 | 1 / 1 |
| Projected width column norm (response-scaled) | 1.58297e-5 | 2.29683e-5 |
| Final width | 1.00000000191 | 1.00000000321 |
| Local A/C projected KKT | 1.73e-18 | 4.34e-18 |
| Final global width-gradient infinity norm | 8.93705e-9 | 6.56096e-9 |
| Required width-gradient threshold | 1.0e-9 | 1.0e-9 |
| First unchanged block / profile stop | sweep 171 / `bounded-tiny-step` | sweep 166 / `bounded-tiny-step` |
| Search result | 200-sweep budget exhausted | 200-sweep budget exhausted |
| Endpoint trust | `primary-bounded-width-kkt-failed` | `primary-bounded-width-kkt-failed` |

The measured A/C design is full rank in both cases. After projecting out A/C, the single width direction is also rank one, so these runs do not demonstrate exact local rank deficiency. The A/C condition numbers are moderate. The width direction is weak in the response-scaled units, and the optimizer stops proposing a usable width step while the gradient remains about 6.6–8.9 times the qualification limit. With three rows and three local parameters there are no residual degrees of freedom; with four rows only one remains. This is weak width information combined with a bounded tiny-step limitation in the current search, rather than a measured rank failure.

Local A/C feasibility is not the first failure: its projected KKT is below 5e-18 and no amplitude is active. The aggregate `ac-feasibility`, endpoint, and numerical-qualification checks fail or become unavailable only after the component has no trusted endpoint. Accepted block updates match the global objective replay to within 4.3e-18, the final unchanged blocks have zero replay delta, and the captured objective/cache replay errors are zero. Independent endpoint replay is unavailable because width stationarity was not met. There is no first failed block: the first no-update records appear at sweep 171 and 166, then repeat through sweep 200.

Both cases share this mechanism. The nonzero projected width derivative and rank-one projection leave room for a fixed-width A/C recovery experiment; this audit does not claim that such a recovery is qualified. The successful Halo 101 component demonstrates that the sparse failure does not numerically invalidate unrelated components, although the current reconstruction wrapper still reports the whole result unavailable.

## Contract audit changes

- Halo preclassification now gets its contributor eligibility from the same frozen `BuildJointPeelingContributorUniverse` used by peeling.
- Peeling uses the frozen eligible contributor identity and target role rather than reclassifying the atom from its current chemical element label. The upstream hydrogen exclusion policy remains in place.
- Regression coverage includes chemical-label invariance, frozen eligibility mismatch rejection, reconstruction-only Halo coverage, unchanged Target endpoints through the reconstruction path, and `n_h=3` with equal radii remaining free.
- P3 Halo CSV output now records the actual nearest Target distance. The output header and result contract identify that meaning explicitly.

## Files and checks

Raw stress evidence is in:

- `docs/research/results/halo_fullabc_p0/three_observation_free_stress/`
- `docs/research/results/halo_fullabc_p0/four_observation_free_stress/`

The executable's captured process peak RSS was 49,643,520 bytes for the 3-row run and 49,954,816 bytes for the 4-row run. These small stress cases are numerical diagnostics, not large-workload memory evidence.

The build and 48 focused regression tests passed:

```text
cmake --build build/debug --target rhbm_tests joint_halo_fullabc_validation --parallel 4
build/debug/bin/RHBM-GEM-TEST --gtest_filter='HaloMedianReconstructionTest.*:PureHaloSupportExclusionTest.*'
```

The earlier intermediate test run exposed a test fixture that altered an off-domain contributor not checked against the OriginalInput identity table; the final test targets an OriginalInput contributor and passes. No production Halo policy, first-round Target estimate, SQLite contract, or default reconstruction method changed in P0.
