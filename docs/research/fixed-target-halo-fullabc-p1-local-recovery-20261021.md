# Fixed-target Halo FullABC P1 local recovery (2026-10-21)

## Scope

This phase adds an explicit research-only `EstimateHaloFullABCJointWithRecovery` entry point. It always attempts FullABC first. Recovery runs only after that attempt is unqualified. The existing preclassification remains `n_h < 3` fixed at the Global Qualified Target Median and `n_h >= 3` free for FullABC.

The four policies compared are:

| Policy | Recovery behavior |
|---|---|
| A — per-Halo A/C | Hold B at the bounded FullABC initial width and solve constrained A/C for each failed Halo, holding other Halos in its component at the Global Target Median. |
| B — component median | Fix only the failed FullABC component at the Global Target Median. Preserve separately qualified components. |
| C — component A/C | Solve constrained A/C jointly over the failed component at the bounded FullABC initial width; if this solve fails, use an explicitly labeled component median fallback. |
| D — isolate | Fix a localized failed Halo at the Global Target Median and re-fit remaining Halos in its component. If the failure cannot be localized, the component median fallback is recorded explicitly. |

The fixed width is `min(Global Target Median B, Bmax4)`, the bounded initial width already used by the FullABC search. The recovery A/C solve has its own full-rank check, conditioning diagnostic, constrained KKT check, and independent replay. Each result carries a parameter source; recovered values are never labeled FullABC unless they came from a separately qualified FullABC endpoint or a qualified FullABC refit.

The compact component solve retains the parent OriginalInput normalization scale. The final mixed-source result is independently replayed across every OriginalInput row, including rows with no free-Halo sensitivity. Its full-domain normalized objective and objective replay discrepancy are recorded.

## Results

Both deterministic stress cases failed the original whole-set FullABC qualification while retaining two individually qualified FullABC Halos in a separate component. Halo 103 was the only free Halo in its failed component. The 101/102 component shared no rows with Halo 103.

| Case | Original FullABC | A and C | B and D |
|---|---|---|---|
| 3 observations | Reconstruction unavailable; 2 of 3 free Halos qualified component-wise | 2 FullABC + 1 A/C recovered; 0 unavailable | 2 FullABC + 1 median recovered; 0 unavailable |
| 4 observations | Reconstruction unavailable; 2 of 3 free Halos qualified component-wise | 2 FullABC + 1 A/C recovered; 0 unavailable | 2 FullABC + 1 median recovered; 0 unavailable |

For these cases, A and C produce the same result because the failed component contains one Halo. D has no other Halo in that component to re-fit; the other qualified component is independent by OriginalInput row coupling and remains unchanged. B and D therefore produce the same parameter sources and downstream metrics for these cases.

| Metric | 3 obs A/C | 3 obs B/D | 4 obs A/C | 4 obs B/D |
|---|---:|---:|---:|---:|
| Halo-only signal RMSE | 6.97e-5 | 5.59e-5 | 7.99e-5 | 6.67e-5 |
| Shared-support RMSE | 3.16e-4 | 3.16e-4 | 6.25e-4 | 6.25e-4 |
| Original-domain prediction RMSE | 2.366e-3 | 2.367e-3 | 2.366e-3 | 2.366e-3 |
| Peeling RMSE | 2.12448e-3 | 2.12449e-3 | 2.12210e-3 | 2.12211e-3 |
| Tail Ratio MAE | 6.6413e-3 | 6.6427e-3 | 6.5607e-3 | 6.5619e-3 |
| Local A RMSE | 0.2520 | 0.2520 | 0.2521 | 0.2521 |
| Local B RMSE | 0.0850 | 0.0850 | 0.0849 | 0.0849 |

The A/C recovery of Halo 103 had rank 2/2 in both cases. Its design condition indicator was 18.450 for three observations and 17.147 for four. Primary and independent A/C projected KKT values were `4.34e-19` and `3.47e-18`; local objective replay discrepancies were 0 and `1.69e-21`. The recovered widths were 1.10868 and 1.10787, below the 6.28850 upper bound. Full mixed-source OriginalInput prediction and objective replay discrepancies were zero in all four runs. The inherited normalization scales were 6.94012 and 7.07140.

Every method preserved the first-round Target endpoint exactly (`target_endpoint_drift_max_abs = 0`). All four recovery methods completed downstream peeling, Local, and Group steps in both cases. Group Alpha was 0 in these small stress workloads, so they do not validate Group Alpha training. These are two deterministic cases, not a multi-seed statistical comparison; the small signal and peeling differences do not establish a general accuracy winner.

## Recovery contract and limitations

- A and C retain the independently qualified FullABC component and assign `ObservationACFixedWidthRecovery` only to the recovered Halo.
- B assigns `PostFailureGlobalMedianRecovery` only to the failed component and preserves the separate qualified component.
- D records that Halo 103 was isolated at the Global Target Median. No refit was needed because it was the only Halo in the failed component; the separate 101/102 component has no shared observation rows.
- If D sees a multi-Halo failed component, it uses a localized invalid/core atom diagnostic when available, then re-fits the remaining component with the existing bounded FullABC solver. If no individual failed Halo can be identified, it records a component-level median fallback rather than inventing a failed-atom identity.
- The postprocessing gate now validates parameter source and recovery evidence before allowing recovered parameters into peeling.
- Recovery remains opt-in. The production default remains `geometry-local-target-median-v1`.

The synthetic process peak RSS was 50,266,112 bytes for the 3-observation case and 50,429,952 bytes for the 4-observation case. These small runs do not establish large-Halo memory scalability.

## Artifacts and verification

Detailed per-Halo profiles, source counts, recovery rank/KKT/replay diagnostics, downstream metrics, runtime, and peak RSS are in:

- `docs/research/results/halo_fullabc_p1/three_observation_free_stress/`
- `docs/research/results/halo_fullabc_p1/four_observation_free_stress/`

The following passed:

```text
cmake --build build/debug --target rhbm_tests joint_halo_fullabc_validation --parallel 4
build/debug/bin/RHBM-GEM-TEST --gtest_filter='HaloMedianReconstructionTest.*:HaloFullABCJointEstimationTest.*:PureHaloSupportExclusionTest.*'
```

All 53 targeted tests passed. Both stress runs executed with `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `VECLIB_MAXIMUM_THREADS=1`, and their experiment assertions checked recovery sources, rank/replay qualification, preservation of the independent FullABC component, and unchanged Target endpoints.
