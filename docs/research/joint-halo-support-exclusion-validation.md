# Joint Component Estimator — Halo-Support Exclusion Validation

> Historical strategy-comparison tools have been retired after the study
> closed. This report and its result artifacts remain the research record.
> Current production regressions cover Pure Exclusion; FixedHalo execution is
> retired and only historical result decoding remains.

| 日期 | 分支 | 基準 HEAD | 範圍 |
| --- | --- | --- | --- |
| 2026-10-16 | develop | `5de15aa10e6f08a3d65aa55d110d53df59520577` | P0 exact support-incidence contract |

本報告記錄 Halo-Support Row Exclusion 的受控驗證。報告當時的
Production baseline 是 Fixed-Neighbor Block-Coordinate Target-Only Joint
Estimator；此 mask 在報告日期仍屬研究用途。後續 Pure Halo-Support
Exclusion 升格 Production 的結論見正式升格報告，本文件的策略比較
runner 已退役。

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

## P1 — Noise 與 simulation support mismatch

### 6Z6U source support 與 observation noise

原始 sidecar 已查證：`single_gaus-v1`、partial charge、outer cutoff 2.5 Å、charge-term cutoff 2.5 Å、比較條件 `squared_distance <= cutoff*cutoff`、0.5 Å voxel spacing、0.5 Å blurring。模擬端本身未 normalization；受控 fitting runner 先以完整 map SD 正規化，再將 iid Gaussian noise 加到整張 map voxel array，然後由同一張 noisy map 建立 FirstStage Halo initialization 與 Joint observations。σ 以 normalized map SD 為單位，所以 σ=0.02 是 map SD 的 2%；同一 realization 供 Fixed Halo、exact exclusion、R=2.0 geometry exclusion 共用。runner 記錄 seed、整張 noisy map 的 SHA-256、noise 實測 mean/SD 及 thread settings。σ=0.02、seed 17 的 realization hash 為 `b73f4596a9ef92da2d4599c18c8242232439096ba055830dd0286bfaa85d4335`，實測 SD 為 0.0200003911。

這次取得的來源 cutoff 與 Joint support 完全相同，故原始 6Z6U 近 Oracle 結果本身不能檢驗 model-support mismatch。逐 Target noisy-map metrics 同時記錄 point quality、parameter validity、local conditioning、parameter stability；NotAssessed 與 Unavailable 保留為各自狀態，不視為 Stable。

### 6Z6U 十組 paired noisy-map 結果

固定 `sigma=0.02`、seeds `17..26`，每個 seed 的三種方法共用同一張 noisy map。A 是 Production Fixed Halo，B 是 exact Pure Exclusion，R2 是 geometry `R=2.0 Å` sensitivity control。表中 candidate 統計只納入 finite、Global/Target convergence、Independent Replay、parameter validity 與 local conditioning 通過的候選；A 有 8/10 個 seed 可用，B/R2 是 10/10。A/B/C 的單位分別是 normalized raw units、Å、normalized raw units。

| 方法 | 可用 seeds | A bias / variance / MSE | B bias / variance / MSE (Å) | C bias / variance / MSE |
| --- | ---: | ---: | ---: | ---: |
| Fixed Halo A | 8/10 | 0.14404 / 2.813e-4 / 0.260388 | 0.003205 / 3.000e-7 / 1.367e-4 | −0.01015 / 2.552e-6 / 0.004351 |
| Pure exclusion B | 10/10 | 2.788e-4 / 2.771e-4 / 2.767e-4 | 9.943e-6 / 3.057e-7 / 3.047e-7 | −9.991e-6 / 2.515e-6 / 2.530e-6 |
| Geometry R2 | 10/10 | 0.09212 / 2.629e-4 / 0.107530 | 0.002043 / 2.906e-7 / 5.832e-5 | −0.006683 / 2.065e-6 / 0.001604 |

這些 MSE 是逐 Target 先計算十組 seed 的誤差，再對 Targets 取平均；artifact 另外保存 population/sample variance、median/P95/max absolute error。B 的 candidate RMSE 為 A `0.01663`、B `5.520e-4 Å`、C `0.001590`；相同定義下 Fixed Halo 為 `0.51028 / 0.01169 Å / 0.06596`，R2 為 `0.32792 / 0.007637 Å / 0.04005`。低 coverage（exact mask retained/original <50%，34 Targets）下，B 的 candidate MSE 是 A `3.373e-4`、B `2.547e-7 Å²`、C `5.182e-6`；Fixed Halo 分別為 `2.355`、`7.809e-4 Å²`、`0.07620`，R2 分別為 `1.146`、`4.507e-4 Å²`、`0.03216`。這支持 B 在這個 matching-support、2% map-SD noise 的 controlled realization 下改善 bias，但不是所有 noise/support regime 的保證。

Fixed Halo 在 seeds 19、23 因一個 FirstStage Halo 不合格而回報 `Unavailable`，沒有候選參數被填入統計。B 的十個 seed 都是 Global/Target Passed；每次只有 12/1,418 Targets 的 point quality/stability formal-qualified（總計 120/14,180 target-seed，low-coverage 為 3/34 targets），其餘 14,060 個 stability 狀態是 `NotAssessed`，不能稱為 Stable。A 的 8 個可用 seed 全部沒有 formal-qualified target；R2 的 10 個 seed 也沒有 formal-qualified target。每個 B fit 的 local conditioning 1,418/1,418 `Adequate`、parameter validity 1,418/1,418 `Passed`，但 stability coverage 仍受 bounded perturbation budget 限制。

完整 per-target/seed 結果在 [6Z6U noise summary JSON](joint-halo-support-exclusion-noise-sigma002-summary.json) 與 [per-target CSV](joint-halo-support-exclusion-noise-sigma002-per-target.csv)。B 相對 A 的 paired target MSE 平均差為 A `−0.2601`、B `−1.364e-4 Å²`、C `−0.004348`；low-coverage 子群分別為 `−2.354`、`−7.807e-4 Å²`、`−0.07620`。負值表示 B 的候選 MSE 較小，但 worst-target regressions 與可用 seed 數仍保留於 artifact，不能用平均值宣稱沒有個別退化。六 Z 六 U 的 B fit peak RSS 最大為 `722,452,480` bytes（約 689 MiB），每次只使用 Eigen/OpenMP/OpenBLAS/Accelerate 各 1 thread。
完整 prediction diagnostics 也保存在每個 seed 的 metrics JSON：B 的 retained fitting domain 是 180,478 rows，`actual_model_retained` RMSE 平均為 `1.74872` normalized map units；R2 的 205,895 retained rows 為 `1.68181`，而 Fixed Halo 的 229,783 original rows 為 `1.62168`。這些是 full-map model residual diagnostics，不是 target parameter MSE；它們不能繞過 P3 對 post-fit reconstruction domain 的 identity contract。

### 固定 seed synthetic paired support control

`joint_halo_support_exclusion_synthetic` 使用 1 Target、1 Halo、515 個 Target-support observation rows；Joint support 固定 2.5 Å，Halo source support 分別為 2.5、3.5、5.0 Å。每一組使用 10 paired seeds（1001–1010）與 σ={0, 0.005, 0.02, 0.05} normalized synthetic observation units。資料層是 Joint rows，與 6Z6U 的 full-map voxel noise 分開標示。四組方法是 Fixed Halo A/B、exact Pure Exclusion B、geometry R=2.0、Full Oracle ABC D。Oracle 僅作 diagnostic control。

主要結果如下；表內是 10 repetitions 的 per-target sampling RMSE，變異與 MSE 完整保留於 [synthetic noise artifact](joint-halo-support-exclusion-synthetic-noise.json)。

| R_sim | σ | 方法 | A RMSE | B RMSE (Å) | C RMSE | candidate numeric pass / 10 |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 2.5 Å | 0.02 | Fixed Halo A | 0.1989 | 0.02546 | 0.08426 | 10 |
| 2.5 Å | 0.02 | Pure exclusion B | 0.1355 | 0.02714 | 0.00307 | 10 |
| 2.5 Å | 0.02 | geometry R=2.0 | 0.1168 | 0.05197 | 0.03757 | 10 |
| 2.5 Å | 0.02 | Oracle D | 0.08721 | 0.01689 | 0.00208 | 10 |
| 2.5 Å | 0.05 | Fixed Halo A | 0.2845 | 0.04774 | 0.08405 | 10 |
| 2.5 Å | 0.05 | Pure exclusion B | 0.4007 | 0.07929 | 0.00918 | 10 |
| 3.5 Å | 0.00 | Fixed Halo A | 0.1525 | 0.02795 | 0.1264 | 10 |
| 3.5 Å | 0.00 | Pure exclusion B | 1.0253 | 0.1278 | 0.05497 | 10 |
| 3.5 Å | 0.00 | geometry R=2.0 | 0.2183 | 0.04068 | 0.08821 | 10 |
| 3.5 Å | 0.00 | Oracle D | 0.0000 | 0.0000 | 0.0000 | 10 |
| 5.0 Å | 0.00 | Fixed Halo A | 1.1398 | 0.1407 | 0.1593 | 10 |
| 5.0 Å | 0.00 | Pure exclusion B | 0.9212 | 0.1305 | 0.1215 | 10 |
| 5.0 Å | 0.00 | geometry R=2.0 | 0.9614 | 0.1263 | 0.1325 | 10 |
| 5.0 Å | 0.00 | Oracle D | 0.0000 | 0.0000 | 0.0000 | 10 |

Matching support 下，σ=0.02 的 Pure Exclusion 對 A、C 優於 Fixed Halo，但 B width RMSE 略大；σ=0.05 時 Pure Exclusion 的 A/B 誤差明顯高於 Fixed Halo，10 個 fit 都 qualified。Support 擴至 3.5 Å 而 Joint 仍是 2.5 Å 時，Pure Exclusion 雖移除了 Joint support 內 Halo rows，保留 rows 上的真實 Halo RMSE 仍為 0.04923；σ=0 時 A error bias 為 +1.0253，Oracle 可將三個 parameter errors 降至浮點誤差。R_sim=5 Å 時 Pure Exclusion 的 A bias 為 −0.9212，retained-row Halo RMSE 為 0.05658。這是明確的 support-mismatch bias 證據，不是 solver failure。σ=0.05、R_sim=3.5 Å 時 Pure Exclusion 僅 8/10 repetitions qualified，另外 2 個因數值 qualification 缺失而未併入參數統計。

此表列的是 solver、global/target convergence、independent replay、parameter validity 與 local conditioning 通過的候選 diagnostic，並非全數正式 qualified estimates。以 point quality `Available` 且 parameter stability `Stable` 作為 formal-qualified 子集時，matching support、σ=0 的 Pure Exclusion 與 Oracle 各有 10/10；σ>0 時四組皆為 0/10。support mismatch 的 σ=0 只有 Oracle 為 10/10。這些 status 與 candidate/formal 統計分開保存在 artifact。Synthetic control 的 Pure Exclusion 在含噪資料下並非普遍優於 Fixed Halo，且不匹配 support 時出現大幅 bias。10 paired seeds 可估這個單 Target controlled design 的 candidate sampling variance，但不能外推為 6Z6U target population variance，也不能取代 formal stability qualification。

## P2 — Piecewise Halo-aware weighting reference

研究 reference 固定使用 `w=lambda` 於 exact Halo-support rows、其他 rows 為 1；λ={0,0.01,0.05,0.1,0.25,1}。Halo correction 是相同 synthetic FirstStage-like A/B、C=0 prior；prior qualification 在這個縮小 reference 中是假設值，沒有呼叫專案 FirstStage qualification protocol。λ=0 與 hard mask 每個 paired seed 的 A/B/C/objective 完全相同。constrained A/C active-set、log-width grid 加 golden-section profile search 使用固定 row weights；objective 及 A/C/logB gradient 由同一組 sqrt(weight) 等價加權逐列計算。analytic gradient 對 central finite-difference 最大誤差為 `1.21e-9`；將所有零權重列改動 123 個單位，objective 與 gradient delta 都是 0；獨立重建的 weighted objective replay delta 為 0。完整結果與每個 seed 結果在 [weighting artifact](joint-halo-support-weighting-reference.json)。

這個 reference 不是 Production profile LM，亦未更改 Joint optimizer 或 Independent Replay 程式；它是因 P1 已顯示明確支援域風險而執行的縮小控制，不可用來宣稱 production weighting 已驗證。support matching、σ=0.02 時，λ=0.25 把 A RMSE 從 0.1355 降至 0.1142、B 從 0.02714 降至 0.02480，但 C RMSE 從 0.00307 升至 0.03277；λ=1 的 A/C RMSE 更升至 0.2552/0.09302。R_sim=3.5 Å、σ=0.02 時 λ=1 的 A/B RMSE（0.2545/0.04239）低於 hard mask（1.1489/0.1422），但 C RMSE 從 0.05368 升至 0.1347，target-only prediction RMSE 從 0.04361 升至 0.07130；full prediction RMSE 下降至 0.04688。加權 full residual 下降並不表示 Target estimates 更準。

因此目前沒有可跨 parameter、noise 與 support regime 的 Soft Weighting 優勢；Hard Mask 在匹配 support 的 low/moderate-noise 下保持較小 C error，但在高 noise 與 support mismatch 時仍可能退化。沒有理由把 Soft Weighting 納入 Production。它也不能消除有限 support 與資料模型不匹配的問題。

## P3 — Postprocessing integration gate

`JointPostprocessingContext` 目前以同一個 `JointProblemInput` 建立 row-ID→map-grid mapping、contributor support、prediction workspace，並要求 `snapshot.row_ids == problem.row_ids`。Workflow 接著用該 context 執行 PostFitPeelingSamples、Post-Joint local MDPDE 與 Group Fitting。遮罩 fit row IDs 無法同時代表原始 reconstruction domain。

新增的 `ExcludedFitRowsAreUnavailableToPostFitPeeling` integration contract test 使用 exact mask 建立 fit problem；Joint RuntimeConvergence 與 TargetRuntimeConvergence 均 Passed。對 target 中心 raw sampling point 執行 post-fit reconstruction 時，因 tricubic stencil 涉及被排除 row，樣本得到 `outside-joint-domain`。將 masked-fit snapshot 與 original problem 一起交給目前 context，則以 row/contributor identity mismatch 明確拒絕。現有 tests 也分別覆蓋 `outside-joint-domain` 與 `missing-contributor-state` failure reason。

因此 PostFitPeelingSamples → Local MDPDE → Group Fitting 的完整 masked workflow 尚未驗證；target 中心 sample 已不足以完成後續鏈路。方案 A 可能保留 FirstStage Halo state，但不會補回被移除的 reconstruction row geometry；方案 B 需要一份 provenance 清楚且取得 qualified Halo state 的 reconstruction domain。沒有實作、驗證雙 domain 及 bounded reconstruction metadata，故 P3 未通過，不能把 masked samples 送入 local/group fitting 當作有效資料。

## P4 — Promotion decision

**Outcome C（報告當時）：不升格 Pure Exclusion，保留當時的 Fixed-Halo Production default。** P1 matching-support noise、明確的 support-mismatch bias 以及 P3 的 post-fit domain failure 均不足以支持全面升格。未放寬 Global Practical Criteria 或 Independent Replay，未改變當時的 Production estimator、public CLI 或 Dynamic Halo Mean 的程式/結果。後續 Pure Exclusion 升格 Production，FixedHalo 估計路徑也已退役。

本報告的 6Z6U 十組 paired summary、低 coverage 統計、stability coverage、RSS/thread settings 已固定於上述 artifacts。十組是預先指定的初始 empirical sample，足以揭示本次 matching-support controlled regime 的 variance 與兩次 Fixed-Halo initialization failure，但不代表所有 geometry/noise 的 sampling variance 已收斂；support-mismatch 與 P3 workflow failure 仍使 Promotion Gate 不通過。
