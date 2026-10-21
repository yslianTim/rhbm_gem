# Fixed-target Halo FullABC P2：真實原子幾何與 6Z6U 配對驗證

日期：2026-10-21
分支：`develop`  期別：P2

## 實驗來源與配對契約

結構使用 `/Users/yslian/data/6z6u.cif`（SHA-256 `349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa`）。局部真實幾何實驗在這些 6Z6U 原子座標及 sidecar 的 simulation truth 上合成 map 與 observation；這不是實驗量測 map。來源 map 為 `/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map`（SHA-256 `6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a`），其 truth sidecar 提供 Target/Halo A/B/C，可計算 simulation-truth RMSE。

每個案例只估計一次 qualified Target endpoint，所有 Halo methods 共用相同 map、noise realization、原子/contributor universe、OriginalInput row identity、kernel support、normalization 與該 Target endpoint。Target A/B/C drift 在每例均為 0。Noise sigma 為 0.003；短鍵組使用 deterministic seeds 4101、4102、4103。Alpha-free methods 不讀取 chemical labels；label 只用來把短接觸辨認為 CIF bond。

## 幾何稽核

| 案例 | Halo 數 | 最短 Halo-Halo 距離 (Å) | <1.5 Å bonded pairs | C–O/C–N | <1.0 Å nonbonded | 重合點 | support-overlap pairs / shared voxels |
|---|---:|---:|---:|---:|---:|---:|---:|
| seq145–147，seed 4101–4103（每 seed 同一結構） | 104 | 1.1790 | 45 | 39 | 0 | 0 | 71 / 522 |
| seq80–82，seed 4101 | 122 | 1.1782 | 57 | 40 | 0 | 0 | 111 / 1,228 |
| boundary seq5–8，seed 4102 | 85 | 1.1910 | 38 | 34 | 0 | 0 | 84 / 1,152 |

舊的 `distance > 1.5 Å` 實驗門檻已移除。短鍵 assertion 已以最新 binary 通過；合法約 1.18 Å 的 bonded contacts 保留，不合理 exact collision 或 nonbonded <1.0 Å 才作為 geometry error 拒絕。P0 的全 Chain A census 另記錄 659 supported Halo 間 1,783 組 support overlap pairs；完整距離分箱與實際 shared voxel 分布見 P0 report 與 census JSON。

## 局部多 seed 結果

五個場景共 519 Halo：343 個 `n_h<3` 事前固定、176 個 `n_h>=3` 先嘗試 FullABC。原始 aggregate FullABC 為 0/5 案例 qualified。每種 recovery 均為 5/5 reconstruction-available；per-Halo 與 component-local A/C 各把 168 個 free Halo 分配為 A/C recovery、保留 8 個 per-component FullABC來源，沒有 unavailable Halo；median recovery 將同一批 168 個失敗自由 Halo 固定為 global median。`available` 只表示來源齊全且 assembly/replay 通過，不代表 A/B/C 參數有可信精度。

下表的前四列為三個短鍵 seed 的 median；parameter/signal 數字以 sidecar truth 計算。Peeling、Local A/B 是同一 Target endpoint 的下游 paired error。Observation A/C 僅提供 Halo parameter/signal evaluation，沒有在此欄執行 Peeling。

| 方法 | Halo signal RMSE | shared-support RMSE | OriginalInput prediction RMSE | Peeling RMSE | Tail Ratio MAE | Local A / B RMSE |
|---|---:|---:|---:|---:|---:|---:|
| Production geometry kNN-4 | 0.073627 | 0.196613 | 0.073707 | 0.628933 | 24.1575 | 0.279100 / 0.021782 |
| Observation A/C，固定 kNN-4 width | 0.006804 | 0.027237 | 0.007404 | 未執行 | 未執行 | 未執行 |
| FullABC，無 Recovery | 不可用 | 不可用 | 不可用 | 不可用 | 不可用 | 不可用 |
| FullABC → per-Halo A/C recovery | 0.049004 | 0.189818 | 0.049074 | 0.628606 | 24.1461 | 0.279398 / 0.021698 |
| FullABC → failed-component Global Median | 0.067718 | 0.199878 | 0.067778 | 0.628826 | 24.1448 | 0.272838 / 0.021536 |
| FullABC → component-local A/C recovery | **0.006523** | **0.027232** | **0.007114** | 0.628634 | 24.1254 | 0.275943 / 0.021474 |
| Oracle Halo truth（診斷用） | 0 | 0 | 0 | 不適用 | 不適用 | 不適用 |

Component-local recovery 的 median signal error 最低，但 recovery-qualified Halo count 不是 accuracy 證明；FullABC raw solver 仍 0/5 aggregate success。三個短鍵 case 的 recovered Halo source count為各 2 個 FullABC、30 個 A/C recovered、72 個 pre-fixed，0 unavailable。Observation A/C 與 component-local recovery 的 A RMSE 很大（約 180–212），但 combined signal RMSE 低；這是個別參數弱可識別、合計訊號較穩的直接例子。

## 可用率之外的 worst case

最差 boundary seq5–8 case 中三種 recovery 仍標記 `available=true`，但 atom 70791 的一個 isolated component 被保留為 `FullABC` 且 `qualified=true`，A=`8.35e122`、B=`0.0934 Å`；其 6 個有效 observation rows 的徑向平方距離只落在約 5.14–6.00 Å。該 endpoint 的 signal 與 downstream Peeling 數值不可信，讓 recovery 方法的 A RMSE 約 `9.06e121`、Peeling RMSE 約 `1.91e57`。Recovery assembly replay prediction error 仍僅 `1.1e-16` 至 `4.4e-16`，objective replay error 為 0；這證明 replay 一致不足以排除窄 width/巨大振幅造成的不可辨識 endpoint。此結果是 P4 worst-case gate 的明確失敗，不可由 median 改善抵銷。

全部五個 local cases 都沒有 Group members；`group_alpha_g`、Group posterior、Group uncertainty 為 null。這些小型場景不作 Group Alpha training 通過的證據。

## 6Z6U Chain A 完整 selection

指定的 Chain A paired run 使用實際 6Z6U coordinate universe、production map normalization 及 sidecar truth，執行入口預期驗證 1,418 Targets、659 supported Halos、229,783 OriginalInput rows。相同 Target endpoint 對比 geometry kNN-4、raw FullABC、component-local Recovery C，並記錄 observable prediction、Peeling、Local/Group、source coverage、wall time 與 process peak RSS。

此完整 paired run 在 2026-10-21 約 2 小時後仍未輸出中間 checkpoint、summary或results JSON；為釋放環境而以 Ctrl-C 中止，exit code 130。沒有取得該 run 的 target census assertion結果、solver outcome、completed wall time或process peak RSS；不可把它寫成 solver failure或資源超限。P0獨立 census已記錄該輸入組合的 1,418 Targets、659 supported Halos、229,783 OriginalInput rows，但不替代這次完整 paired run。這表示 P2 的局部真實幾何與多seed驗證完成，實際 Chain A end-to-end paired selection未完成。

## 執行與產物

- Build：`cmake --build build/debug --target joint_halo_fullabc_validation --parallel 4`；成功。
- Runtime threads：`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`。
- 最新 binary 的短鍵 assertion rerun 完成；process peak RSS 256,360,448 bytes。
- 多 seed suite process peak RSS 267,124,736 bytes，5/5 cases 完成。
- Chain A完整 paired runner啟動後運行約 2 小時未寫出任何案例結果，人工中止、exit code 130；completed runtime與RSS未測得。
- 完整每案例結果、逐 Halo source 與 metrics：`docs/research/results/halo_fullabc_p2/real_geometry_multiseed/results.json` 及同目錄 CSV。
- 短鍵 assertion rerun：`docs/research/results/halo_fullabc_p2/real_geometry_short_bond_assertion/results.json`。
- Chain A中止前未產生非空案例輸出或 results JSON，故沒有保留空白CSV產物。
