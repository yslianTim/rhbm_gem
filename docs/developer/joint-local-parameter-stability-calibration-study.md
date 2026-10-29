# Joint Component Local Parameter Stability 門檻校準研究

## 決策摘要

**Partial Qualification。** 保留 Production 的 **LocalParameterChangeWarning = 1e-8** 與 **LocalPredictionChangeWarning = 1e-8**。本研究找到可重現的門檻分類效果，也隔離出 T4 solver-only 的 baseline eligibility 耦合；但高精度 synthetic chain 不能證明放寬至 1e-6 不會隱藏弱識別或不同起點造成的參數差異。不同訊號尺度的受控樣本也不足以支持新的 normalization 或拆分參數門檻。

研究只更動分類邏輯的程式結構與回歸測試；Production KKT、practical-v2、求解器、接受條件和 endpoint 行為維持原值。所有實驗性開關與量測程式都已移除。

## 範圍與方法

- Initial branch/HEAD：develop / d866e062a45206e82a33069469f9d44751d200c9；工作樹乾淨。當時本機 origin/develop tracking ref 同為此 SHA。
- Release 建置使用 AppleClang 21.0.0、arm64、CMake/Ninja、-O3 -DNDEBUG。建置與 CTest 並行上限為 4；數值研究固定 OMP_NUM_THREADS=1、OpenBLAS/MKL/VECLIB 為 1。最終 CTest 使用既有可通過的 OMP_NUM_THREADS=2，數值函式庫執行緒仍為 1。
- 基線與門檻候選使用同一份 deterministic chain synthetic 資料、partition、初值、width 上界與 solver budgets。局部門檻矩陣是對已收集的同一組合格變化量重新分類，不為九個門檻重跑 solver。
- 六個 workload-policy 條件納入門檻統計：3-target T0、30-target T0、30-target T4 solver-only、30-target T4 full-policy、300-target T0、1,000-target T0。統計單位是「被選中的 block × run/policy」，不是獨立隨機樣本；其中 16 筆完成原始量測，另有 2 筆 Unavailable、104 筆未選中的 NotAssessed。不得把此數量解讀為獨立樣本量。
- T0 指 Production search criteria {A/C 1e-7, Width 1e-9}；T4 solver-only 指 search 使用 {1e-7, 1e-5}，local assessment 仍使用 practical criteria；T4 full-policy 則在臨時實驗編譯中讓 local assessment 也接受 T4 endpoint。所有臨時 policy 變更在結束後復原。
- Unavailable 和 NotAssessed 不參與 warning 門檻重分類。Stable/Sensitive/Unstable 僅對資格足夠的量測套用分類函式。

## A. Git 狀態

| 項目 | 結果 |
|---|---|
| Initial HEAD | d866e062a45206e82a33069469f9d44751d200c9 |
| Branch | develop |
| 最終 HEAD | 本報告與 README 索引同一提交後的 develop HEAD（本報告提交即最終提交；完整 SHA 見該分支 git log） |
| Working tree | 提交前工作樹乾淨；報告提交後以 git status 確認乾淨 |
| origin/develop | 本機 tracking ref 起始於 Initial HEAD。最終 develop 預計領先此本機 ref 3 commits；即時遠端查詢因執行環境無法解析 github.com 而失敗，故不能聲稱已確認遠端最新狀態 |
| Commits | 88b2e865 refactor: isolate local stability classification；a01cd7cf test: cover local stability threshold matrix；本報告提交列於最終 git log |
| Push | 否 |
| Branch/worktree | 未建立 branch 或 worktree。初始檢查看到另一個先前已存在、與本研究無關的 detached worktree；未碰觸 |
| Notion | 未讀取、修改或新增任何 Notion 文件 |

## B. S0–S5 Roadmap 驗收

| Stage | 目標 | 結果 | Commit |
|---|---|---|---|
| S0 | 原始數值與分類審計 | **部分完成。** 重現 30-target 三種分類與 3/30/300/1,000 T0 工作負載；保留所有被選 block 的彙總變化量、主要 KKT、conditioning、求解/評估次數及資源估計。造成 T4 full-policy 分類的方向與端點數值完整；其餘案例未保留每個正負 seed 的完整端點診斷，詳見限制。 | 無獨立 S0 程式變更 |
| S1 | Qualification／Classification 隔離 | **程式與回歸測試通過。** 新增純 classifier，將已完成量測與閾值分類分開；同一組 T4 原始量測可重分九組門檻，且不重跑 solver；無有效量測或方向未齊不會被轉成 Stable。replay 和 KKT 仍由上游 qualification 負責。 | 88b2e865, a01cd7cf |
| S2 | 二維門檻矩陣 | **完成於已取得的 16 筆量測。** 九組門檻共用相同資料；同時補上等於門檻、略高於門檻、非有限值、缺少量測等純分類測試。 | a01cd7cf |
| S3 | Truth／Noise／Scale／Difficult Cases | **部分完成；不足以升格。** 有 chain 真值誤差、exact-fit 尺度探索與一個 noise 探索，也有既有 weak/active 單元 fixture；多起點、完整噪聲梯度、near-boundary、rank-deficient/near-singular 的 stability matrix、不同真實 width/A/C 幅度及真實 cryo-EM map 均未完成或無可用 stability 證據。 | 無 |
| S4 | Diagnostic Promotion／Cleanup | **Partial Qualification；保留原門檻並完成清理。** G3 的程式隔離與完整回歸通過；G1 及 G2 因 coverage、qualifier coupling 與 noisy-scale replay 缺口未完整通過。 | 88b2e865, a01cd7cf 及本報告提交 |
| S5 | Future Width KKT Recommendation | **只提出獨立後續 gate；不變更 Production。** 本次 Stability 校準不足以支援 Width KKT 1e-8 promotion。 | 本報告提交 |

## C. 30-target 異常分析

### T0 與 T4 分解

| Run | Block | Baseline A/C KKT | Baseline Width KKT | Eligibility / replay | D_B | D_AC | D_pred | Classification at 1e-8 / 1e-8 |
|---|---|---:|---:|---|---:|---:|---:|---|
| T0, 5 sweeps | block-1 | 1.6023889e-9 | 1.8806623e-11 | Eligible / pass | 8.8314940e-9 | 3.0084309e-9 | 3.2812498e-9 | Stable |
| T0, 5 sweeps | block-2 | 1.4894172e-9 | 1.7431612e-11 | Eligible / pass | 8.2140479e-9 | 2.7974584e-9 | 3.0503010e-9 | Stable |
| T0, 5 sweeps | block-3 | 1.5266e-16 | 1.9984e-15 | Eligible / pass | 3.7825298e-13 | 7.8737449e-14 | 1.1533822e-13 | Stable |
| T4 solver-only, 4 sweeps | block-1 | 9.4629833e-8 | 1.1122324e-9 | **Ineligible / replay failed** | — | — | — | Unavailable |
| T4 solver-only, 4 sweeps | block-2 | 8.7952211e-8 | 1.0307811e-9 | **Ineligible / replay failed** | — | — | — | Unavailable |
| T4 solver-only, 4 sweeps | block-3 | — | — | Eligible / pass | 1.0040e-11 | 2.1113e-12 | 2.3525e-12 | Stable |
| T4 full-policy, 4 sweeps | block-1 | 9.4629833e-8 | 1.1122324e-9 | Eligible / pass | 5.2135492e-7 | 1.7762220e-7 | 1.9376222e-7 | **Unstable** |
| T4 full-policy, 4 sweeps | block-2 | 8.7952211e-8 | 1.0307811e-9 | Eligible / pass | 4.8488116e-7 | 1.6515706e-7 | 1.8011287e-7 | **Unstable** |
| T4 full-policy, 4 sweeps | block-3 | — | — | Eligible / pass | 1.0040e-11 | 2.1113e-12 | 2.3525e-12 | Stable |

T4 full-policy 的兩筆 Unstable 來自不同 block/direction：

- block-1：coordinate 方向的 minus seed；requested/effective infinity displacement 1e-3。D_B 為 5.2135492e-7（約為 parameter 門檻的 **52.14×**），D_AC 為 1.7762220e-7，D_pred 為 1.9376222e-7（約為 prediction 門檻的 **19.38×**）。造成最大值的 perturbed endpoint A/C KKT 1.3878e-16、Width KKT 0、Independent Replay pass、objective delta 7.4871e-15。
- block-2：common-mode 方向的 plus seed；requested/effective infinity displacement 2.8867513e-4。D_B 為 4.8488116e-7（**48.49×**），D_AC 為 1.6515706e-7，D_pred 為 1.8011287e-7（**18.01×**）。造成最大值的 perturbed endpoint A/C KKT 1.9429e-16、Width KKT 5.5511e-16、Independent Replay pass、objective delta 6.4693e-15。
- 所有 3 個 block 在 T4 full-policy 都是 WellConditioned；block-1/2 free A/C rank 為 24/24、projected Width rank 為 12/12、width bound active 數為 0。T0/T4 選中的是同一組 block，方向標籤集合也相同（weak-singular、coordinate、common-mode）。目前資料沒有顯示這次改變由不同的 width active-bound 或 rank 分類造成；但沒有保存所有 seed 的完整 active-set trace，不能排除更細的 A/C 邊界行為。
- T4 solver-only 的 block-1/2 不是 perturbation solver 或 perturbed replay 失敗：兩個 block 在開始擾動前，baseline Width KKT 分別高於 local practical 1e-9 約 11.2% 與 3.1%。A/C KKT 已低於 1e-7；記錄的 replay failure reason 是 primary-bounded-width-kkt-failed，且兩個 block 都是 0 次 perturb solve。這證明 local baseline eligibility 仍與 practical Width KKT 門檻耦合；單純改 warning threshold 無法解除 Unavailable。
- 沒有其他已選 block 在 30-target T4 full-policy 下為 Sensitive。T4 full-policy 的 probe endpoint 確實通過該實驗 policy 的 eligibility/replay，且量測超過 warning 門檻；這不是只因更改文字分類造成的假變化。另一方面，T4 停在 4 sweeps，T0 停在 5 sweeps，端點並不相同。故數值顯示該較早端點沿已測方向有實際 local prediction response，但目前無法判斷這是科學上重要的估計敏感度，還是提前停止留下的端點差異。低 truth error 也不能單獨證明 warning 過嚴。
- T4 full-policy 的 conditional reprofile 另有一筆 StrictlyQualified 結果，但它自己的 Independent Replay 標記為 Failed；失敗原因未保存。這不是上表中兩個已完成且 replay 通過的 stability probe endpoint，故不改寫其測量分類；它仍是 G1 的未解證據缺口。

## D. Threshold Sensitivity Matrix

下表涵蓋六個 run-policy 條件中所有 16 筆已選取且可重分類量測；Unavailable 與 NotAssessed 保持固定。計數不是獨立抽樣統計。

| Parameter Threshold | Prediction Threshold | Stable | Sensitive | Unstable | Unavailable | NotAssessed |
|---:|---:|---:|---:|---:|---:|---:|
| 1e-8 | 1e-8 | 7 | 0 | 7 | 2 | 104 |
| 1e-8 | 1e-7 | 7 | 0 | 7 | 2 | 104 |
| 1e-8 | 1e-6 | 7 | 7 | 0 | 2 | 104 |
| 1e-7 | 1e-8 | 7 | 0 | 7 | 2 | 104 |
| 1e-7 | 1e-7 | 7 | 0 | 7 | 2 | 104 |
| 1e-7 | 1e-6 | 7 | 7 | 0 | 2 | 104 |
| 1e-6 | 1e-8 | 7 | 0 | 7 | 2 | 104 |
| 1e-6 | 1e-7 | 7 | 0 | 7 | 2 | 104 |
| 1e-6 | 1e-6 | 14 | 0 | 0 | 2 | 104 |

### 按 workload 的結果與轉移

- 3-target T0 與 30-target T0：共 4 筆，九組門檻全為 Stable。
- 30-target T4 solver-only：1 筆 Stable、2 筆 Unavailable；九組門檻不改變結果。
- 30-target T4 full-policy：1 Stable + 2 Unstable at prediction 1e-8 或 1e-7；prediction 1e-6 時，parameter ≤1e-7 為 1 Stable + 2 Sensitive，parameter 1e-6 時為 3 Stable。
- 300-target T0：所選 block-1、block-7 為 Unstable，block-21 為 Stable；prediction 1e-7 仍有兩筆 Unstable，prediction 1e-6 改為兩筆 Sensitive（parameter ≤1e-7）或 Stable（parameter 1e-6）。另外 23 個 block 因代表性 probe 上限未選，為 NotAssessed。
- 1,000-target T0：所選 block-1、block-3、block-40 均為 Unstable；prediction 1e-7 仍為 Unstable，prediction 1e-6 時三筆變為 Sensitive（parameter ≤1e-7）或 Stable（parameter 1e-6）。另外 81 個 block 為 NotAssessed。
- 相對 baseline 1e-8 / 1e-8：單獨將 parameter 門檻提高到 1e-7 或 1e-6，16 筆分類都不變；將 prediction 提到 1e-7 也都不變；只有 prediction 1e-6 會讓 7 筆 Unstable → Sensitive。再把 parameter 門檻提高到 1e-6，這 7 筆再由 Unstable（相對 baseline）成為 Stable；相對前一列則是 Sensitive → Stable。未觀察到 Stable → Unstable 或 Unavailable/NotAssessed → Stable。
- 主要分類槓桿是 Prediction threshold；本資料沒有任何只因 parameter threshold 放寬而改變分類的 block。這是該批測量的結果，不能外推成 parameter threshold 一般不重要。
- 以 T4 full block-1/2 為例，D_pred 約 1.8–1.94e-7，只把 prediction threshold 放寬到 1e-7 仍會保留 Unstable；到 1e-6 才降為 Sensitive。300/1,000-target 的 D_pred 約 1.78–1.94e-7，也相同。放寬到 1e-6 / 1e-6 可把目前已選的所有有效量測判為 Stable；這種轉換缺少困難案例支持，不能作為 promotion 理由。
- 已保留的純 classifier regression 覆蓋固定 T4 block-1 metrics 的九種組合、等於門檻、nextafter 略高於門檻、缺少 qualified perturbation、不完整方向與非有限輸入。原始 per-seed 參數向量未全部保存，故不能從報告重建每一個方向的分類；九格重分類則可由已記錄的 max D_B/D_AC/D_pred 獨立重建。

### 分類量測與成本

| Run | Selected blocks / total blocks | Stability cost | Perturbation work | Workspace peak estimate |
|---|---|---:|---:|---:|
| 3-target T0 | 1 / 3 | block CPU 約 10.8 ms | 6 solves / 30 evaluations | 未留獨立 peak 記錄 |
| 30-target T0 | 3 / 15 | 約 245 ms wall；block CPU 約 247 ms | 共 18 solves；每 block 28–29 evaluations | 約 23.7 MB |
| 30-target T4 solver-only | 3 / 15 | 未完成的兩個 block 無擾動成本 | block-1/2 為 0 solve；block-3 有 6 solves | 約 23.7 MB |
| 30-target T4 full-policy | 3 / 15 | 約 250 ms wall | 成功 block 每個 6 solves / 約 28–29 evaluations | 約 23.7 MB |
| 300-target T0 | 3 / 26 | 0.2153 s wall / 0.2141 s CPU | 18 solves / 78 evaluations | 53,378,250 B |
| 1,000-target T0 | 3 / 84 | 0.3003 s wall / 0.2990 s CPU | 18 solves / 72 evaluations | 130,208,180 B |

Memory 欄是 Parameter Stability workspace estimate，不是整個 process 的 peak RSS。300/1,000-target 的代表性選擇分別只覆蓋 3/26、3/84 個 block；其餘是未選取，而不是未通過或已穩定。

## E. Scientific Accuracy

### Synthetic chain 精度

| Targets | Policy | Max A truth error | Max C truth error | Max B truth error | Endpoint |
|---:|---|---:|---:|---:|---|
| 3 | T0 | 6.5353e-5 | 1.4822e-5 | 3.3070e-6 | 3 sweeps; certification/replay pass |
| 30 | T0 | 6.7098e-5 | 1.4793e-5 | 3.6832e-6 | 5 sweeps; certification/replay pass |
| 30 | T4 full | 6.7098e-5 | 1.4787e-5 | 3.6832e-6 | 4 sweeps; certification/replay pass |
| 300 | T0 | 6.5510e-5 | 1.5113e-5 | 3.2833e-6 | 4 sweeps; certification/replay pass |
| 1,000 | T0 | 6.5510e-5 | 1.5113e-5 | 3.2833e-6 | 4 sweeps; certification/replay pass |

30-target T0/T4 真值誤差相對既有應用標準 0.005 很小。先前研究記錄兩個端點間 13,710 個 prediction rows 的 max/RMS/P95/P99 絕對 prediction 差為 3.2855e-7 / 3.6282e-8 / 1.2596e-7 / 1.5385e-7；objective 從 1.51707603e-10 到 1.51721796e-10，相對差約 9.4e-5。此組比較支持 T4 端點整體 fit 精度仍高，但不否定被局部擾動測得的敏感方向。

300-target T0 residual prediction RMSE 為 3.7780e-6；1,000-target 為 3.7789e-6。這是 fit residual RMSE，不是 local D_pred。兩者的 selected blocks 都是 WellConditioned 且滿秩；其 D_pred 約 1.78–1.94e-7，真值誤差仍低。這表示低真值誤差與 Stable/Unstable 不是同一種品質量測，不能把低 error 當成必須 Stable 的證明。

### 困難情況與覆蓋率

| Case | Stability evidence | 狀態與限制 |
|---|---|---|
| Exact synthetic，9 observations | 近 machine-precision 變化；見 F | 已跑尺度探索；幾乎零殘差，不代表有噪聲場合 |
| 3/30/300/1,000-target chain | 上表與 D 節 | 已跑；大型 workload 只選 3 個代表 block |
| High-overlap cube（8 targets） | 先前 Width study 有合格 fit，但沒有可用的同一 endpoint stability matrix | Stability NotAssessed |
| Weak-width / weak singular | 正式兩原子單元 fixture 有小型單側/不可行擾動，結果 Stable；另有 weak-width workload 的 Stability NotAssessed | fixture 規模太小，不能代表弱識別族群 |
| Active-width upper bound | 正式單原子 fixture 的 inward one-sided probe 為 Stable；width bound active 行為有測試 | 有限單元案例，不是多種上界位置與幅度矩陣 |
| Noise-only | 過去合格 fit 的 stability 為 NotAssessed；本次 noisy-scale 較小探索只一組有 baseline/replay 合格 | 部分 NotAssessed；不可用於完整 noise sensitivity 結論 |
| Near-boundary | 無適用 stability 原始量測 | Not Run |
| Rank-deficient / near-singular | 先前 Width study 的 fit 未通過 qualification，故不應分類為 Stable | Stability threshold matrix Unavailable |
| Different initializations / multi-start | 未執行 | Not Run |
| Noise levels（noiseless/low/moderate/high 完整梯度） | 未執行完整矩陣 | Not Run |
| Different true widths、A/C amplitudes | 未執行 | Not Run |
| Real cryo-EM map | repo 沒有適用 fixture，未執行 | Unavailable |

既有 weak-singular 兩原子 fixture：D_B 9.87e-14、D_AC 1.25e-13、D_pred 2.73e-16，3 個方向/seed 評估中含單側或 bound-infeasible 方向。既有 active upper-bound 單原子 fixture：D_B 0、D_AC 1.27e-16、D_pred 4.81e-17，只有 inward minus 方向有效，另一側為零位移。這兩個小型測試只驗證單側語義，不提供不同寬度、overlap 或噪聲狀況下的錯誤分類率。

受控合成工作的 endpoint 認證與 Independent Replay 在列明的 chain workload 均通過，真值誤差也低；但 multi-start、弱識別矩陣、noise-level/scale 梯度與真實 map 都不足。因此放寬到 1e-6 確會把此批約 2e-7 的局部 prediction response 隱藏為 Stable，現有證據無法證明這種反應普遍沒有科學意義。

## F. Normalization Analysis

目前：
- D_B = ||Δ log B||∞，小變化下近似 B 的相對變化。
- D_AC = max_j |Δβ_j| ||X_j||₂ / ObservationScale，將 A/C 係數變化按特徵及觀測尺度縮放。
- D_pred = ||Δŷ||₂ / max(1, ||ŷ_base||₂)。

D_pred 在 ||ŷ_base||₂ < 1 時是未做相對尺度化的絕對 L2 差；此區域會隨訊號幅度和重複觀測數增加而變化。在 baseline norm 大於 1 時則近似 relative prediction L2；若每筆觀測完全複製，分子及分母同乘 sqrt(N)，比值大致不變。分母中的 1 使 map/觀測 normalization 穿越 norm=1 時改變 metric 的意義。以 residual noise 或每觀測 RMSE 正規化可能改善可解釋性，但要先明訂 noise estimator 與零噪聲處理，再獨立驗證；本次未選定替代定義。

精確合成 9-row scale probe：

| Signal multiplier / case | ||ŷ||₂ | context.scale | D_pred | D_B / D_AC | Status |
|---|---:|---:|---:|---:|---|
| 0.001× | 7.8594e-4 | 1 | 1.7780e-18 | 2.60e-14 / 3.72e-17 | Stable |
| 1× | 0.785943 | 1 | 1.7701e-15 | 2.57e-14 / 3.66e-14 | Stable |
| 1,000× | 785.943 | 785.943 | 2.2624e-15 | 2.61e-14 / 4.79e-14 | Stable |
| 每筆觀測複製兩次 | 1.11149 | 1.11149 | 2.2391e-15 | 2.58e-14 / 4.70e-14 | Stable |

這些近乎 exact-fit 的 D_pred 位於浮點殘差量級，不足以判定 noisy case 的分類是否 scale invariant。一個加噪探索（0.03 sin + 0.015 cos）中，0.001× 資料 baseline/replay 合格、分類 Sensitive，D_B 2.0471e-8、D_AC 1.5096e-10、D_pred 1.2992e-12；尺度 1、1,000 與 duplicate rows 的 replay 失敗，故這些 case 為 Unavailable，不能與前者組成尺度結論。相應的 conditional-reprofile replay failure reason 未被保留。

D_B 無量綱但衡量相對 width；scaled A/C 用 observation/feature 尺度，與 D_B 的數值含義不同。相同 1e-8 不等於兩者有相同科學意義。應先累積不同 width、A/C amplitudes、feature scales 及 noise 的證據，再考慮把 log-width 與 scaled A/C 警告分開；本研究不改正式 API。

## G. Promotion Decision

決策：**Partial Qualification**。

| Gate | 結果 | 理由 |
|---|---|---|
| G1 — Diagnostic Numerical Validity | **部分通過，整體未完成** | T0/T4 full-policy 實際合格 probes 的 endpoints 有 KKT/replay 證據；classifier 不重跑 solver，invalid/incomplete 及缺少方向結果不能變 Stable。solver-only local eligibility 確實耦合 practical Width KKT；noisy-scale 多數 case replay failure，失敗 reason 未保存；所有 probe 的每 seed qualification/replay 欄位未完整留存。 |
| G2 — Scientific Classification Quality | **未通過 promotion 標準** | T0/大 workload truth error 低，但不同 endpoint、multi-start、noise/scale、弱/active/近奇異覆蓋不足。1e-6 / 1e-6 將本批 7 個 Unstable 重分類為 Stable，尚無證據證明這不會掩蓋有意義的敏感度。 |
| G3 — Engineering and Resource Safety | **通過程式面檢查** | classifier 是純函式，變更不呼叫 solver/replay；門檻只改分類。完整 Release build/CTest 通過；Global KKT、practical-v2、求解器與 persistence schema 未改。 |

Production 繼續使用：

- LocalParameterChangeWarning = 1e-8
- LocalPredictionChangeWarning = 1e-8
- Global A/C KKT 1e-7
- Global projected Width KKT 1e-9
- practical-v2

本研究沒有證明原本的 1e-8 是最佳門檻；結論只是目前沒有足夠證據採用候選值。下一輪需先改善完整 endpoint/seed 證據留存，並解耦 local assessment eligibility 的研究規格，而不是直接放寬它。

### S5：未來 Width KKT 1e-8 的獨立 gate

本次 stability evidence 尚不足以啟動 promotion。未來獨立研究至少應：

1. 從同一個已接受 baseline endpoint 同時保存並延續所有 η=log B 與 β=(A,C) 狀態；確認 T0/T4 各自的 certification、replay、objective、KKT 與 sweep，不以只複製 width 的 shadow run 代替 full-state continuation。
2. 重複 30-target 的 sweep-saving 對照，並在 3、300、1,000 targets 測量 paired end-to-end wall/CPU、stability overhead、solve/evaluation 數及 process RSS。A/C/Width qualification 與 warning classification 分別報告。
3. 覆蓋 multi-start、不同噪聲、訊號/觀測尺度、true widths、A/C amplitudes、high-overlap、weak-width、active/near-boundary、rank-deficient/near-singular、noise-only；未通過 qualification 的案例只標 Unavailable。
4. 使用同一批已 replay 的 perturbation raw metrics 比較 3×3 matrix，保證 Unavailable/NotAssessed 不變；將 local stability 評估放在相同認可 endpoint 或解釋 endpoint 差異。
5. Promotion 前要求：兩政策皆通過 KKT/feasibility/certification/replay；完整狀態 continuation 無隱藏參數或 prediction 差異；多 workload 效益可重複且超過量測變異；分類不把已識別的弱方向/不同起點差異轉為 Stable；point-quality 與 JSON contract regression 通過。Width KKT 1e-8 不得因本次 warning study 自動採用。

## H. Cleanup Audit

- 移除臨時 T4 criterion 編譯、研究 harness、測試期 console instrumentation、per-seed diagnostic prints，以及臨時加進 ParameterStability_test.cpp 的大型 workload case；正式測試檔已復原至提交版本。
- 恢復 ConvergenceCriteria.hpp 中 practical criteria {1e-7, 1e-9}；production warning constants 保持 {1e-8, 1e-8}。
- 保留正式用途的 pure classifier、boundary/non-finite/incomplete 輸入測試，以及以固定 T4 原始 metrics 重分九格的 regression test。
- 沒有 research-only runtime switch、threshold policy branch、JSON exporter、第三方 dependency 或 persistence schema 變更。
- 未發現本次新增而未使用的 helper/include。局部 classifier 只讀已計算量測，不影響 A/C/B accepted state。
- Global KKT、practical-v2、solver、endpoint certification/replay 和研究外的 production 數值契約皆未改。純分類重構不改 fit inputs、search 或 candidate acceptance；以同一既有 fixture 的 regression 和完整 CTest 驗證無數值流程回歸。

## I. Final Verification

- Compiler / configuration：AppleClang 21.0.0，arm64，CMake/Ninja Release，-O3 -DNDEBUG。
- Build command：cmake --build build/release --parallel 4；最終清理後重新建置通過。
- Full CTest：env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 ctest --test-dir build/release --parallel 4 --output-on-failure；**19/19 通過，16.73 秒**。
- Focused stability tests：JointParameterStabilityTest.* 為 8/8 通過。
- Numerical qualification：表列 T0 chain 3/30/300/1,000 與 30-target T4 full-policy endpoints 通過 certification / Independent Replay；T4 solver-only 有 2 個 local baseline width-gate failure，故 Unavailable。近奇異與 rank-deficient 不視為可分類的成功 endpoints。
- JSON round-trip：完整 CTest 中 persistence/JSON test 通過；schema 未變。
- Performance / memory：T0 stability assessment 約 11 ms（3 target）、245 ms（30）、215 ms（300）、300 ms（1,000）；300/1,000 workspace estimates 約 53.4 MB / 130.2 MB。因 Production 數值行為未改，沒有新增求解成本；上述數據是局部 stability diagnostic 成本，不是完整 fitting wall time。
- 工作樹與最終提交：本報告和索引加入後的 develop HEAD；最終提交與 clean status 由完成訊息所列 Git 檢查確認。未 push。遠端即時狀態因 DNS 無法驗證；本機 tracking ref 仍指向初始 HEAD。
- Known limitations：沒有 multi-start、完整 noise × scale × amplitude × true-width 矩陣、real cryo-EM map；大量 block 未被代表性 probe 選中；每個 probe seed 的完整 endpoint trace 與 noisy-scale replay failure reason 不齊。以上均為 Not Run 或 Unavailable，不列為 Pass。
