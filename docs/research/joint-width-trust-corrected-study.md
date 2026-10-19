# Joint Component Width Trust：修正後的配對實驗

## 範圍與重現方式

本研究在 `develop` 上依序完成 P0 scope 修正、P1 truth normalization 修正，再使用同一組 observations、initialization、12-target forward FixedNeighbor partition、200-sweep budget、objective、qualification gates 與 independent replay，比較 A（Production baseline）、B1（conditioning-aware）及 B2（conditioning + cancellation-aware）。B1/B2 只改變研究用的 accepted-width-trial callback；沒有修改 Production 數值策略或 qualification thresholds。

實驗程式為 `tests/experiments/joint_width_trust_experiment.cpp`；配對執行、scope audit 與摘要程式為 `tests/experiments/run_joint_width_trust_paired.py`。完整 32 cases / 96 runs 的逐 case、逐策略指標、commands、輸入 hashes 與 scope audit 在 [paired results JSON](joint-width-trust-corrected-results.json)。Raw per-run trace 與 log 保留於 `/private/tmp/rhbm-joint-width-trust-corrected-20261019`；它們不包含在 Git commit。

可重跑完整批次的命令：

```sh
python3 tests/experiments/run_joint_width_trust_paired.py \
  build/bin/joint_width_trust_experiment \
  /private/tmp/rhbm-joint-width-trust-rerun \
  /Users/yslian/data/6z6u.cif \
  /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  /Users/yslian/data/emd_11103_unsharpened.map \
  --seeds 1001 1002 1003 1004 1005
```

P2 synthetic seeds 1001–1003 原樣保留；1004–1005 是額外 evaluation seeds。策略設定沒有依 evaluation 結果調整，沒有另外建立 calibration set。Synthetic 六種 profile 各跑五 seeds；兩個 6Z6U 固定輸入各跑 A/B1/B2。

外部輸入 SHA-256：

| 輸入 | SHA-256 |
| --- | --- |
| `6z6u.cif` | `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa` |
| 6Z6U simulation map | `6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a` |
| simulation `.simulation.json` | `33a3c221ec3f2fff18dcd2a0ea6484f68334971fc688f7cd100e220ad5848931` |
| 6Z6U unsharpened map | `fa9a6b15ad86d471130fc83495754cbf0a2e4084adc5327bf09f21d7e5ec60d1` |

## P0：Policy Scope 與 evidence availability

舊 runner 只追蹤 serial 1032，callback 卻會檢查其他所有 block。那些 block 的 `tracked_atoms` 為空，`Evaluate()` 產生非空 `evidence-unavailable` 字串；FixedNeighbor 將它當成 rejection reason。舊 unscoped runs 因而在 simulation 出現 27,559 次、unsharpened 出現 27,994 次 rejection；主要來源是 out-of-scope 的空 tracked evidence，並非 B1/B2 的 conditioning/cancellation 判定。

修正後，serial 1032 先經目前 Joint selection、Structural Component 與現行 12-target mutable-core partition 定位，再生成 policy scope；沒有硬編碼 block number。相同 CIF 的 simulation 與 unsharpened map 因 support partition 不同，實際定位如下：

| 輸入 | component | 1-based mutable core block | core atom ids |
| --- | --- | ---: | --- |
| normalized simulation research | `1` | 6 | `1032`, `1034` |
| unsharpened map | `1` | 5 | `1029`, `1030`, `1031`, `1032`, `1033`, `1042`, `1043`, `1044`, `1045`, `1048`, `1060`, `1061` |

32 個 paired cases 中，64 個 B1/B2 scope audits 全部通過：scope 與 A 相同、首次 rejection 位於 scope 內、`out_of_scope_interventions = 0`，且首次 intervention 前 block trace 與 A 完全一致。真實 map 上，simulation 的首個 intervention 前有 5,939 個 block events（5,896 個 out-of-scope events），unsharpened B1 有 1,264 個（1,255 個 out-of-scope events），B2 有 844 個（838 個 out-of-scope events）；這些前綴都逐項等於 baseline。首次 intervention 後的軌跡差異是該 scope 內拒絕 trial 後的 state propagation。

Policy Scope、Diagnostic Scope 與 Evidence Availability 現在分開處理。out-of-scope callback 是 no-op；in-scope evidence unavailable 也記錄原因後 no-op，不 veto trial。64 個 B1/B2 runs 的 evidence-unavailable trial 總數為 0。實驗 capture 上限為 100,000 events、32 atoms、4,096 local rows；serial 1032 block 最大 3,132 rows，沒有 event/atom/row limit exhaustion，也沒有 spectrum-unavailable event。Regression test 另以零 capture budget 確認搜尋結果與 baseline 一致。

## P1：Simulation truth 與 normalization 尺度

Simulation sidecar 由 map generator 建立，記錄原始 map 生成模型的 `element`、`charge_used`、`effective_gaussian_width` 與 `effective_charge_width`；生成 manifest 的 `normalization_applied` 為 `false`。`SINGLE_GAUS` 的 Gaussian coefficient 是元素 atomic number，charge-kernel coefficient 是 `charge_used`；兩個 kernel 和 Joint estimator 的公式、有效 width 與 2.5 Å cutoff 一致。sidecar 的 A/C 因此處於 original-domain scale，B 是物理 Å。

`MapValueArrayNormalization()` 將每個 voxel 除以 map standard deviation，沒有減 mean。若 divisor 為 `d`，normalized fit 的真值應是 `A/d`、`C/d`，而 B 不變。Runner 明確記錄 truth、estimate、comparison scale，RMSE 只在相同尺度比較。兩種 map 的 normalization 為：

| 輸入 / mode | normalization | divisor | truth |
| --- | --- | ---: | --- |
| simulation map / `normalized-simulation-research` | 已執行 `MapValueArrayNormalization()` | `0.44031545732594374` | original-domain A/C 轉 normalized-domain A/C；B 保持 Å |
| unsharpened map / `real-map-production-normalized` | 已執行 `MapValueArrayNormalization()` | `0.005236894964862669` | N/A — no truth |

這次 normalized simulation research mode 與 Production CLI `--simulation true` **不是同一種 preprocessing**：Production CLI simulation flag 會停用 map normalization；本研究 runner 仍執行 normalization，並按上述 divisor 轉換 sidecar truth。兩種 mode 沒有混在同一個 paired comparison。

修正前 simulation A 的 A/C RMSE `8.48009 / 0.626977` 將 normalized estimate 直接與 original-domain sidecar truth 比較，不能作為準確度數字。修正後 Production baseline 的 A/C RMSE 為 `9.5087230e-6 / 2.9210961e-6`；B RMSE 為 `1.4648388e-7 Å`。sidecar manifest 不正規化、AC 除 divisor、B 保持物理單位均有 regression coverage。

## P2：Synthetic paired 結果

以下是五 seeds 的中位數。A/B/C RMSE 欄的順序為 A/B/C；RMSE 對 unqualified run 是明確標記的 candidate metric，不能當成正式 Joint estimate。每個 seed 的 qualification 與 status 仍保留在 JSON。Prediction RMSE 是 fitting-domain observations 上的 candidate residual RMSE。Wall time 和 peak RSS 來自 test-instrumented research runner，不是 production benchmark。

| Profile（truth width；noise σ） | Strategy | qualified / 5 | median A/B/C RMSE | median prediction RMSE | median max accepted B (Å) | median κ / |ρAC| / cancellation | median rejections | median wall / peak RSS |
| --- | --- | ---: | --- | ---: | ---: | --- | ---: | --- |
| Well（0.72–0.83 Å；5e−4） | A | 5/5 | 0.004150 / 0.000548 / 0.000125 | 4.94384e−4 | 0.830202 | 9.810 / 0.774997 / 1.07861 | 0 | 0.113 s / 25.72 MiB |
|  | B1 | 5/5 | 0.004150 / 0.000548 / 0.000125 | 4.94384e−4 | 0.830202 | 9.861 / 0.774997 / 1.07861 | 1 | 0.062 s / 25.75 MiB |
|  | B2 | 5/5 | 0.004150 / 0.000548 / 0.000125 | 4.94384e−4 | 0.830202 | 9.861 / 0.774997 / 1.07861 | 1 | 0.062 s / 25.72 MiB |
| Weak（8.0–9.1 Å；1e−5） | A | 2/5 | 31.7882 / 1.64189 / 0.041798 | 9.92251e−6 | 9.13218 | 1106.9 / 0.999980 / 1.82076 | 0 | 2.649 s / 28.59 MiB |
|  | B1 | 5/5 | 39.8809 / 2.58453 / 0.073714 | 1.00404e−5 | 8.55727 | 745.6 / 0.999975 / 1.82076 | 1 | 1.055 s / 27.44 MiB |
|  | B2 | 5/5 | 41.4116 / 2.61291 / 0.074512 | 1.00404e−5 | 8.28098 | 714.8 / 0.999973 / 1.82076 | 1 | 1.066 s / 28.44 MiB |
| Moderate（1.8–2.24 Å；5e−4） | A | 5/5 | 1.87341 / 0.195422 / 0.007313 | 4.94521e−4 | 2.60571 | 78.69 / 0.996677 / 1.00214 | 0 | 1.236 s / 27.39 MiB |
|  | B1 | 2/5 | 1.88276 / 0.195262 / 0.007355 | 4.94521e−4 | 2.60571 | 78.69 / 0.996677 / 1.00214 | 104 | 3.330 s / 29.30 MiB |
|  | B2 | 2/5 | 1.88276 / 0.195262 / 0.007355 | 4.94521e−4 | 2.60571 | 78.69 / 0.996677 / 1.00214 | 110 | 3.300 s / 28.91 MiB |
| High noise（0.72–0.83 Å；5e−3） | A | 5/5 | 0.041690 / 0.005488 / 0.001240 | 0.00494377 | 0.833994 | 9.810 / 0.776141 / 1.07803 | 0 | 0.058 s / 25.69 MiB |
|  | B1 | 5/5 | 0.041690 / 0.005488 / 0.001240 | 0.00494377 | 0.833994 | 9.861 / 0.776141 / 1.07803 | 1 | 0.072 s / 25.75 MiB |
|  | B2 | 5/5 | 0.041690 / 0.005488 / 0.001240 | 0.00494377 | 0.833994 | 9.861 / 0.776141 / 1.07803 | 1 | 0.072 s / 25.67 MiB |
| Mixed width（0.55–2.31 Å；5e−4） | A | 5/5 | 1.44600 / 0.142638 / 0.005210 | 4.94770e−4 | 2.90744 | 78.92 / 0.997935 / 1.12303 | 0 | 0.449 s / 26.83 MiB |
|  | B1 | 5/5 | 1.44605 / 0.142638 / 0.005210 | 4.94770e−4 | 2.90744 | 78.92 / 0.997935 / 1.12303 | 1 | 0.363 s / 26.62 MiB |
|  | B2 | 5/5 | 1.44605 / 0.142638 / 0.005210 | 4.94770e−4 | 2.90744 | 78.92 / 0.997935 / 1.12303 | 1 | 0.365 s / 26.88 MiB |
| A/C correlation stress（3.2–3.53 Å；5e−4） | A | 0/5 | 8.37187e15 / 1981.57 / 3.03620e7 | 4.98444e−4 | 4872.37 | 2.878e8 / 1.000000 / 9.933e6 | 0 | 0.700 s / 26.77 MiB |
|  | B1 | 0/5 | 4.88952e13 / 556.512 / 1.34570e6 | 4.98444e−4 | 1804.87 | 3.155e7 / 1.000000 / 1.136e6 | 120 | 2.749 s / 28.42 MiB |
|  | B2 | 0/5 | 3416.65 / 4.31295 / 1.05646 | 4.99952e−4 | 14.9779 | 2257 / 1.000000 / 145.06 | 127 | 2.599 s / 29.06 MiB |

Interpretation: B1/B2 preserve the well and high-noise profiles to displayed precision. Weak-identifiability qualification rises from 2/5 to 5/5, but A/B/C RMSE worsens while prediction RMSE changes only slightly; convergence is not parameter identifiability. Moderate qualification falls from 5/5 to 2/5, with no max-width or prediction benefit. Mixed-width accuracy is slightly worse without width reduction. The A/C-correlation stress is unqualified for every strategy; B2 strongly reduces runaway candidate widths and coefficient errors but still yields no qualified estimates. These rows are candidate diagnostics only.

## P2：6Z6U paired 結果

Simulation input 有 1,418 selected targets、180,478 fitting observations；unsharpened input 有 1,418 targets、189,377 observations。B1/B2 的 width、condition 與 cancellation 指標追蹤 serial 1032 所在的 mutable core；objective、qualification 與 prediction 指標是整個 Joint component 結果。unsharpened map 沒有 truth，因此 A/B/C RMSE 全部為 N/A — no truth。

| Dataset / strategy | Qualification / replay / stop | Objective | Global AC KKT / width gradient | A/B/C RMSE | Prediction RMSE | max accepted B (Å) | max κ / cancellation | Rejections / out-of-scope | Sweeps / evaluations | Wall / peak RSS |
| --- | --- | ---: | --- | ---: | ---: | --- | --- | --- | --- |
| Simulation / A | Passed / Passed / block-stationary | 1.2333160e−13 | 9.43707e−8 / 8.77489e−10 | 9.50872e−6 / 1.46484e−7 Å / 2.92110e−6 | 9.20736e−7 | 0.500000 | 3.3273 / 1.13168 | 0 / 0 | 83 / 22,908 | 310.36 s / 532.28 MiB |
| Simulation / B1 | Passed / Passed / block-stationary | 1.2335678e−13 | 9.43707e−8 / 8.77489e−10 | 9.51391e−6 / 1.46561e−7 Å / 2.92110e−6 | 9.20830e−7 | 0.49999983 | 3.3273 / 1.13279 | 37 / 0 | 83 / 22,908 | 310.74 s / 549.43 MiB |
| Simulation / B2 | Passed / Passed / block-stationary | 1.2334538e−13 | 9.43707e−8 / 8.77489e−10 | 9.51161e−6 / 1.46526e−7 Å / 2.92109e−6 | 9.20787e−7 | 0.49999988 | 3.3273 / 1.13279 | 37 / 0 | 83 / 22,908 | 312.89 s / 472.91 MiB |
| Unsharpened / A | Unavailable / Failed / block-recovery-exhausted | 0.018417515 | 5.27587e−6 / 4.42769e−5 | N/A — no truth | 0.92603847 | 6099.81 | 2.3143e8 / 2.83166e7 | 0 / 0 | 70 / 19,705 | 277.78 s / 1325.25 MiB |
| Unsharpened / B1 | Unavailable / Failed / block-recovery-exhausted | 0.018417536 | 6.57208e−7 / 4.03181e−6 | N/A — no truth | 0.92603900 | 6099.81 | 2.3143e8 / 2.81452e7 | 79 / 0 | 97 / 27,200 | 380.69 s / 1372.72 MiB |
| Unsharpened / B2 | Unavailable / Failed / block-sweep-budget | 0.018418317 | 7.35062e−10 / 1.75400e−5 | N/A — no truth | 0.92605863 | 7.64778 | 356.92 / 42.0411 | 194 / 0 | 200 / 56,213 | 560.93 s / 1361.60 MiB |

Simulation B1/B2 each reject 37 in-scope trials and retain qualification/replay, with near-identical A/B/C and prediction RMSE to A. The cancellation multiplier adds no measurable benefit where the tracked cancellation ratio is about 1.13. The old 27,559 unscoped vetoes become 37 genuine in-scope numerical rejections.

On the unsharpened map, B1 does not reduce accepted width or improve prediction and retains the baseline `block-recovery-exhausted` / replay failure. B2 has a real cancellation-aware effect: accepted width drops from 6,099.81 Å to 7.65 Å, and observed condition/cancellation drop substantially. However, its proposed width path still reaches 1,365.30 Å; the search hits the 200-sweep budget, independent replay fails, objective regresses from 0.018417515 to 0.018418317, prediction RMSE increases from 0.92603847 to 0.92605863, and runtime doubles. This map has no parameter truth, so the width change alone cannot establish scientific accuracy. The old 27,994 unscoped vetoes become 79 B1 and 194 B2 in-scope rejections.

All real map capture runs emitted at most 400 events, with no event/atom/row-limit exhaustion and no unavailable spectra. Capture limits therefore did not silently alter any run. No global dense matrix, full-length per-block row mapping, or global covariance was added. Full-target 6Z6U policy was not attempted: the scoped real B2 case reduces width but still fails qualification and replay, so the evidence does not justify widening the intervention. Synthetic full-target tests are limited to the 12-atom controlled component.

## P2 結論與 P3 gate

The corrected evidence separates numerical convergence, parameter identifiability and scientific accuracy. B1/B2 can increase qualification in weak synthetic cases and reduce width in extreme A/C-cancellation cases, but they worsen parameter accuracy in weak cases, reduce qualification in the moderate profile, and do not resolve the real unsharpened block failure. B2's strongest real width reduction still consumes twice the runtime and fails the unchanged budget/replay gates. B2 provides an independent cancellation effect in high-cancellation cases, but not a qualified endpoint.

P3 promotion decision is recorded in the follow-up commit after reviewing this P2 result. Until then, the Production estimator remains unchanged.
