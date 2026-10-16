# Joint Component Estimator — Halo-Support Exclusion Validation

| 日期 | 分支 | 基準 HEAD | 範圍 |
| --- | --- | --- | --- |
| 2026-10-16 | develop | `5de15aa10e6f08a3d65aa55d110d53df59520577` | P0 exact support-incidence contract |

本報告記錄 Halo-Support Row Exclusion 的受控驗證。Production 維持 Fixed-Neighbor Block-Coordinate Target-Only Joint Estimator；這個 mask 目前只供研究與測試使用。

## P0 — Observation row 與 support 契約

`BuildJointProblem` 先取所有 selected Target 的 2.5 Å support union 作為原始 observation domain。`row_ids` 保存 map grid index 的字串 identity；每個 contributor 以稀疏 `(row, squared_distance)` support incidence 連到該 domain。Supported non-target contributors 是 Halo。Joint support 條件為 `squared_distance <= 6.25 Å²`，邊界屬於 support。

現行 Fixed-Halo 問題以 qualified FirstStage A/B、C=0 建立固定 Halo prediction，`WithFixedHalos` 從原 observation 中扣除該 prediction。Fixed-Neighbor block、global objective、KKT/width criteria 與 independent replay 都以同一份 `JointProblemInput` 為 row domain。研究 mask compact rows 後，保存新的 row identity 與 `fit row -> original row` index；target contributor indices 也獨立 remap。

`JointPostprocessingContext` 目前也從同一份 `JointProblemInput` 建立 map-grid-to-row lookup，並讀取完整 contributor support/state 以重建鄰居 prediction。Tricubic stencil 超出這個 domain 時會回報 `outside-joint-domain`；缺 contributor state 時會回報 `missing-contributor-state`。因此 masked fit rows 不能直接取代原始 postprocessing geometry。

Dynamic Halo Mean 是獨立 research-only 實作：它以 `DynamicHaloMean` 建立 Halo ownership、mean state、augmented local design、width derivative 和 replay；目前沒有共用 Halo exclusion helper，也沒有改變其實驗程式或結果。

### Exact mask

對原始 domain `V_original`，精確 mask 由 contributor support incidence 產生：若任一 non-target contributor 在 row `v` 有 support，則排除 `v`。因此 `V_fit = V_original \ union_h support(h)`。它不依距離再次估算 support 半徑，也不使用估計參數、residual 或 simulation truth 選 row。沒有 supported Halo 時保留所有 rows；空 fitting domain 或任一 Target 沒有 retained observations 時，masked input builder 明確拒絕建立 fit input。

Deterministic tests 覆蓋無 Halo、單一與多個 Halos、Target/Halo 共用 row、多個 Target 共用 row、`2.5 Å` support 邊界、無 observation support 的 Halo、全空 mask、Target 零 coverage、component fragmentation，以及 original row identity 與 compact remapping。對每一個保留 row，測試斷言沒有 Halo support incidence。

### 6Z6U geometry-mask 對照

資料為 6Z6U chain A、non-hydrogen targets，1,418 Targets、659 supported Halos、229,783 原始 rows。Exact mask 與先前 geometry mask `nearest_halo_distance >= 2.5 Å` 的結果如下：

| 指標 | 結果 |
| --- | ---: |
| Original / retained / excluded rows | 229,783 / 180,478 / 49,305 |
| Geometry R=2.5 與 exact row-ID digest | `71e231a2313d94aa` / `71e231a2313d94aa` |
| Geometry-only / exact-only rows | 0 / 0 |
| Halo support entries 正好位於 2.5 Å | 0 |
| Retained rows 上的 Halo support violations | 0 |
| Target components / Targets without retained rows | 1 / 0 |
| 保留低於 50% 原始 rows 的 Targets | 34 |
| FirstStage-qualified Target seeds | 1,415 / 1,418 |

兩種定義在邊界的規則不同：geometry mask 只排除 `d_H < 2.5 Å`，精確 support mask 會排除 `d_H = 2.5 Å`，因為 support 使用 `r² <= 6.25 Å²`。這份 map 的 Halo support 沒有精確落在該邊界的 entries，因此 row identity sets 完全相同。逐 Target retained row counts、Target/block local rank 與 conditioning 保存於 [P0 6Z6U audit JSON](joint-halo-support-exclusion-p0-6z6u.json)。139 個可分析 blocks 均 full rank；無 rank-deficient 或 ill-conditioned block。最大 normalized condition indicator 為 AC 10.6933、ABC 12.1596、projected width 1.32441；另有 3 個 blocks 因原始 FirstStage target seeds 不合格而 unavailable，未將它們算作 mask 導致的 rank failure。

P0 audit 在 Eigen/OpenMP/OpenBLAS/Accelerate 各 1 thread 下耗時 13.196 s，process peak RSS 為 590,331,904 bytes。沒有新增 global dense design 或逐 block full-length row mask。

## 階段結論

P0 通過：support-incidence exclusion 與可追溯 row remapping 已實作並有 deterministic contract tests；6Z6U 的 exact mask 與歷史 R=2.5 geometry mask row set 完全一致。這只驗證 mask 與資料結構，不足以證明含噪聲 accuracy、model support mismatch、weighted objective 或 postprocessing/group workflow 正確。Production default 未改變。
