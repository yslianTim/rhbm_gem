# Joint Component Local KKT 精度與局部擾動敏感度研究

## 決策摘要

**Partial Qualification。** Local baseline eligibility 已使用與 Global practical-v2 相同的 A/C 1e-7、Width 1e-9。Local perturbed endpoint 則另外要求 Strict A/C、Strict Independent Replay，並檢查 practical Width。成功量測的端點遠低於 Strict 門檻。30-target T4 的兩個 unavailable baseline 即使使用 Global-equivalent Q3，仍因 Width KKT 略高於 1e-9 而不合格。

目前證據不足以把 Strict Local endpoint qualification 改成 practical-v2，也不足以改動 1e-8 / 1e-8 Warning Threshold。T4 已量到局部 prediction response，但 T0、T4 使用不同 sweep endpoint；缺少同端點 null perturbation、完整逐 seed KKT/replay 原始值，以及安全隔離的 Local Solver Precision 控制。因此維持所有 Production 門檻與 solver 行為，不升格。

研究新增本報告，沒有修改 Production 程式碼、測試、persistence schema、Global KKT 或 Local Warning Threshold。

## 範圍與基線

- 初始 branch：develop。
- 初始 HEAD 與本機 origin/develop tracking ref：463e06206acf5b614b8c1078e937d439e5c81fab。
- 初始工作樹乾淨；未建立 branch/worktree。檢查到一個先前存在的 detached worktree，未使用、未修改。
- 3/30/300/1,000-target 與 T0/T4 stability 原始研究摘要來自同一個初始 HEAD 的既有報告 joint-local-parameter-stability-calibration-study.md；Projected Width baseline 參照 joint-width-kkt-promotion-study.md。本輪重新審核程式與測試，但不把舊報告未保存的逐 seed 數值補成新量測。
- 本輪 Release build：cmake --build build/release --parallel 4，成功，Ninja 無待建工作。
- 本輪完整 CTest：19/19 通過，CTest 並行 4，最終執行 15.87 秒。Local stability focused tests：8/8 通過，17 ms。
- 本輪重跑 Production closed-loop benchmark：3、30、300、1,000 targets；數值函式庫固定單執行緒，每個 workload 先 warmup 一次，再量測五次。

同一 Initial HEAD 的既有 Local Stability evidence 摘要：

| Workload / policy | selected blocks | Stability summary | perturbation solves / evaluations | stability CPU / wall | workspace peak estimate |
|---|---:|---|---:|---:|---:|
| 3-target T0 | 1/3 | 1 Stable | 6 / 30 | 約 10.8 ms CPU | 未記錄 |
| 30-target T0 | 3/15 | 3 Stable | 18 / 每 block 28–29 | 約 247 ms CPU / 245 ms wall | 約 23.7 MB |
| 30-target T4 solver-only | 3/15 | 1 Stable、2 Unavailable | block-1/2 為 0 solve；block-3 有 6 solves | block-1/2 無擾動成本 | 約 23.7 MB |
| 30-target T4 full-policy | 3/15 | 1 Stable、2 Unstable | 成功 block 各 6 / 每 block 約 28–29 | 約 250 ms wall | 約 23.7 MB |
| 300-target T0 | 3/26 | 1 Stable、2 Unstable | 18 / 78 | 0.2141 s CPU / 0.2153 s wall | 53,378,250 B |
| 1,000-target T0 | 3/84 | 3 Unstable | 18 / 72 | 0.2990 s CPU / 0.3003 s wall | 130,208,180 B |

整體是 16 個 selected block summaries：7 Stable、7 Unstable、2 Unavailable；另有 104 NotAssessed。這些數字來自同 HEAD 的既有研究文件；workspace peak estimate 不是 process peak RSS。30-target global endpoints 中，T0 是 5 sweeps、A/C KKT 1.60e-9、Width KKT 1.88e-11、objective 1.51707603e-10；T4 是 4 sweeps、A/C KKT 9.46e-8、Width KKT 1.1122e-9、objective 1.51721796e-10。既有研究記錄 max truth error A/C/B：T0 為 6.7098e-5 / 1.4793e-5 / 3.6832e-6，T4 為 6.7098e-5 / 1.4787e-5 / 3.6832e-6；完整 accepted parameter vectors 沒有附在 developer report。T0/T4 global endpoint certification 與 replay 通過；T4 solver-only 的兩個 local stability baseline 則未通過 local Width KKT eligibility。

本輪 Production benchmark：

| Targets | Search median (range) | Wall median (range) | Sweeps / accepted blocks | Global A/C KKT | Projected Width KKT | Objective | Process peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2.830 ms (2.819–2.956) | 3.711 ms (3.679–3.837) | 3 / 3 | 2.78e-16 | 6.7065e-10 | 2.01991e-10 | 24.1 MB |
| 30 | 100.937 ms (94.772–209.267) | 123.940 ms (114.509–234.093) | 5 / 15 | 1.6024e-9 | 1.8807e-11 | 1.51707603e-10 | 38.8 MB |
| 300 | 776.210 ms (757.128–865.205) | 977.934 ms (959.496–1,078.930) | 4 / 104 | 2.9920e-8 | 1.1108e-10 | 1.47922971e-10 | 91.9 MB |
| 1,000 | 2.599 s (2.546–2.619) | 3.302 s (3.238–3.445) | 4 / 336 | 1.6387e-8 | 3.3321e-11 | 1.47923161e-10 | 240.5 MB |

四個 endpoint 都是一輪 certification、Independent Replay 通過、0 failed blocks。Process peak RSS 是 benchmark 程序 high-water mark，與上面的 stability workspace estimate 分開。30-target 某次量測偏慢；該 workload 的 median 仍約 124 ms。

## A. Git 狀態

| 項目 | 結果 |
|---|---|
| Initial HEAD | 463e06206acf5b614b8c1078e937d439e5c81fab |
| Branch | develop |
| Working tree at start | Clean |
| origin/develop | 本機 tracking ref 與 Initial HEAD 相同；未 fetch，故不宣稱已核對即時遠端 |
| Final HEAD | 本報告文件 commit 後的 develop HEAD；完整 SHA 見最終訊息 |
| Commit | 本報告文件一個 commit；沒有 Production 程式或測試變更 |
| Push | 否 |
| Branch/worktree | 未建立；既存 detached worktree 未觸碰 |
| Notion | 未讀取或修改 |

## B. L0–L4 執行結果

| Stage | 工作 | 結果 | Commit |
|---|---|---|---|
| L0 | Local KKT contract 與 baseline 審計 | 完成契約審計；數值分布沿用同 HEAD 的既有報告。已追蹤 Global endpoint、Fixed-neighbor response、Local baseline、profile search、perturb endpoint、Certificate、Width KKT、Independent Replay、classification。舊報告沒有保存完整逐 seed KKT/replay，故不補造分布。 | 本報告 |
| L1 | Q0–Q3 qualification matrix | 部分完成。以可重建的 30-target T0/T4 baseline 數值獨立判斷 Q0–Q3；已保存的成功 perturb endpoints 由現行 Strict gate 保證通過 Q0。完整 workload/seed 矩陣、失敗 seed 的 raw endpoint 欄位不在舊 evidence 中，不能完整計數。 | 本報告 |
| L2 | Full-state reference、null control、solver precision | 未完成 promotion 所需實驗。程式已有同 Width A/C conditional reprofile 指標，但不是 Local Search 的 full-state null continuation；舊報告的 reprofile replay 曾失敗且未留完整原因。沒有 Local-only precision injection；不改共享 solver 來補做。 | 本報告 |
| L3 | Local/Global KKT 尺度一致性 | 完成程式路徑與正規化審計。Local 與 Global 使用同一 context.scale，Local objective 條件於固定鄰居；相同 KKT 數值適合表示同一尺度下的 per-block 一階 stationarity，不代表兩個不同變數集合的 max KKT 數值必然相同。 | 本報告 |
| L4 | KKT × Warning Threshold | 沿用並重驗同 HEAD 已保存的 16 個 block summary。九格 Warning 結果可重建；strictly qualified response 在所存資料中沒有因放寬 KKT 而新增，Warning 分類變化主要由 Prediction Threshold 驅動。 | 本報告 |

## C. Local KKT 契約與責任

### 數值路徑

| 階段 | Objective／scale | KKT／停止條件 | 職責 |
|---|---|---|---|
| Global joint search / endpoint | Global observation rows；context.scale = max(1, ||y||₂) | practical-v2：A/C 1e-7、bound-projected Width 1e-9；Fixed-neighbor stationarity 使用呼叫者傳入的 criteria | Global sweep 停止與 endpoint qualification |
| Local fixed-neighbor workspace | core atoms 影響的 affected rows；y_local = old_core_prediction - global_residual；沿用同一 Global context.scale | 不另以 local row count 重設 scale | 把鄰居視為固定，建立條件式 local objective |
| Local baseline raw certificate | 條件式 local objective；A/C column norm 與 residual 均以 context.scale 正規化 | CertifyLinear() 一律計算 1e-10 的 certificate.kkt_passed；baseline eligibility 本身明確比較 raw KKT 與 practical {1e-7,1e-9} | Strict bit 是 certificate 欄位；baseline eligibility 採 Global practical-v2 |
| Local baseline Independent Replay | 獨立 scalar 重算同一 local design、effective response、residual 與 objective | 呼叫時傳入 practical {1e-7,1e-9}；數值一致性界線不隨 criteria 改變 | 核對 baseline KKT qualification 與重算一致性 |
| Conditional reprofile | 固定原 eta，重新解 A/C | 用 certificate.kkt_passed 標記 Strict A/C；Replay criteria 為 practical | 診斷 baseline 參數經 A/C re-fit 的變化，不是 full-state continuation |
| Perturbed profile search | 同一條件式 objective，從擾動 Width 初始化；每次 profile evaluation 重新求 A/C | LM ftol=1e-14、xtol=1e-12、gtol=1e-12；Local stability budget 最多 11 次 profile evaluation、6 accepted updates | Local profile optimizer；實際可因 budget 提前停止 |
| Perturbed endpoint primary gate | 同一條件式 objective | 要求 candidate.certificate.kkt_passed，即 A/C 1e-10，並且 feasible | Local endpoint Strict A/C qualification |
| Perturbed endpoint Replay | 獨立重算 | CheckTrust() 預設 Strict {A/C 1e-10, Width 1e-12}；另外再檢查 practical Width 1e-9 | Strict KKT qualification 與獨立數值一致性 |
| Warning classification | 已完成的 D_B、D_AC、D_pred | parameter/prediction warning 各 1e-8 | 只分類合格量測，不應替代 solver 或 replay qualification |

其中 A/C 只限制 Gaussian coefficient 為非負，Charge coefficient 無下界；A/C KKT 是經 observation scale 與各 column norm 正規化後的 projected infinity norm。Width gradient 除以 scale²，並用 eta - min(eta-gradient, eta_upper_bound) 計算 bound-aware projected KKT。local objective 是條件式 residual 的平方和，不會依 affected-row 數重設 observation scale。

CheckIndependentReplay() 將 optimality criteria 作為參數；prediction/residual、objective、gradient 和 KKT-agreement 的一致性界線獨立於 KKT threshold：prediction/residual 2e-12 + 2e-13 × |primary|、objective 1e-12 + 2e-12 × |objective|、A/C KKT agreement 1e-13、Width gradient agreement 1e-13 + 2e-9 × |gradient|。放寬 optimality threshold 不會自動放寬這些一致性檢查。

CertifyLinear() 的 strict kkt_passed 是共用欄位，LinearSolve refinement、profile accepted-trial replay 等路徑也依賴 strict certificate。EvaluateProfile() 遇到 strict A/C certificate 不合格時會回傳無效 profile，因此也無法只靠改 Local endpoint 的事後 qualification 來模擬較寬鬆的 solver stopping。不能把共用 certificate 改成 practical 來當作局部研究開關。Local perturbed endpoint 還會經過共用 SearchProfile() 的 strict CheckReplay()；其 LM gtol/ftol/xtol 是 profile solver tolerance，不等同於 KKT criteria，也沒有 Local-only precision parameter。

### KKT 尺度可比性

Global context 將 observation scale 設為 max(1, ||y||₂)；準備 local workspace 時沿用該值。local effective response 以固定鄰居後的 residual 還原 core atoms 的條件式 response。Core 的 affected rows 包含其支援區，因此 A/C 梯度仍使用相同 response scale 與 column-norm normalization；Width 梯度則用相同 scale²，並在上界使用 bound-aware projected KKT。

因此 Global 與 Local KKT 可在相同 observation scale 下比較同一個 block 變數的一階 stationarity。它們不是相同的最佳化問題：Global KKT 是整個 component 所有變數的最大值；Local KKT 是固定其餘鄰居後的條件式 block 指標。相同門檻有合理一致的單變數尺度意義，但不保證兩個 max statistic 數值相同。

Stability metrics 是 D_B = ||Δ log B||∞、D_AC = max_j |Δβ_j| ||X_j||₂ / context.scale、D_pred = ||Δŷ||₂ / max(1, ||ŷ_base||₂)。D_pred 在 baseline prediction norm 小於 1 時是絕對 L2 差；在大於 1 時才接近相對 L2 差。Local Warning 是分類門檻，不是 KKT 或 replay tolerance。

## D. L1 Qualification Matrix

Q0–Q3 定義同任務規格。下表只列舊研究中能從原始 baseline KKT 數值重建的 30-target 三個 T0 block 與兩個 T4 block；其餘 workload 的 per-block baseline 原始 KKT 未保存。

| 配置 | 30T T0 baseline eligible | 30T T4 block-1 / block-2 eligible | 30T T4 solver-only unavailable | 可重建原因 |
|---|---:|---:|---:|---|
| Q0: A/C 1e-10, Width 1e-12 | 1/3（僅 block-3） | 0/2 | 2 | T0 block-1/2 A/C 分別 1.6024e-9、1.4894e-9；T4 A/C 9.4630e-8、8.7952e-8，皆超過 Q0 |
| Q1: A/C 1e-10, Width 1e-9 | 1/3 | 0/2 | 2 | 已知 T0 block-1/2 仍因 A/C 超標；T4 block-1/2 也因 A/C 超標 |
| Q2: A/C 1e-8, Width 1e-9 | 3/3 | 0/2 | 2 | T0 三個 baseline 都合格；T4 block-1/2 A/C 仍超標 |
| Q3: A/C 1e-7, Width 1e-9 | 3/3 | 0/2 | 2 | T4 block-1/2 A/C 合格，但 Width 1.1122e-9、1.0308e-9 仍超標 |

依任務要求另列 Qualified／Unavailable 分類欄位。Eligible 與 Unavailable 是同一組已知的五個 baseline block 子集；其他 workload 和逐 seed endpoint 的分母/原始記錄未保存，故以「未保存」表示，不以 0 代替。

| Local A/C | Local Width | Eligible baseline（30T identified subset） | Perturbed Qualified（逐 seed） | Newly accepted | Unavailable baseline（identified subset） | Replay numerical failure |
|---:|---:|---|---|---|---|---|
| 1e-10 | 1e-12 | T0 1/3；T4 block-1/2 0/2 | 未保存；兩個已記錄 T4 最大反應端點通過 | 未保存；無已記錄的 KKT-only 新接受端點 | T0 block-1/2、T4 block-1/2 | 未保存；T4 記錄的原因是 KKT qualification |
| 1e-10 | 1e-9 | T0 1/3；T4 block-1/2 0/2 | 未保存；兩個已記錄 T4 最大反應端點通過 | 未保存；無已記錄的 KKT-only 新接受端點 | T0 block-1/2、T4 block-1/2 | 未保存；T4 記錄的原因是 KKT qualification |
| 1e-8 | 1e-9 | T0 3/3；T4 block-1/2 0/2 | 未保存；T0 30-target 原 summary 為 3 Stable blocks | 未保存；T4 block-1/2 A/C 仍未合格 | T4 block-1/2 | 未保存；已知 T4 failure reason 不指向 replay numeric disagreement |
| 1e-7 | 1e-9 | T0 3/3；T4 block-1/2 0/2 | 未保存；T0 30-target 原 summary 為 3 Stable blocks | 0 個已知新增 Qualified T4 endpoint | T4 block-1/2 | 未保存；已知 T4 failure reason 不指向 replay numeric disagreement |

此表不是完整 workload matrix：它列出能以已保存 raw baseline KKT 獨立重建的部分。對成功 endpoint，Strict gate 意味著其 KKT 已落在 Q0；對失敗 endpoint，原始 primary/independent KKT 沒有逐 seed 留存，所以不能判斷放寬能否新增接受。Numerical replay disagreement 與 replay 的 KKT rejection 必須分開計數。

T0 baseline 已記錄值：

| Block | A/C KKT | Projected Width KKT |
|---|---:|---:|
| 1 | 1.6023889e-9 | 1.8806623e-11 |
| 2 | 1.4894172e-9 | 1.7431612e-11 |
| 3 | 1.5266e-16 | 1.9984e-15 |

T4 solver-only baseline 已記錄值：

| Block | A/C KKT | Projected Width KKT | 原因 |
|---|---:|---:|---|
| 1 | 9.4629833e-8 | 1.1122324e-9 | Q0–Q2 A/C 不合格；Q3 Width 不合格 |
| 2 | 8.7952211e-8 | 1.0307811e-9 | Q0–Q2 A/C 不合格；Q3 Width 不合格 |
| 3 | 未保存 | 未保存 | Production practical 下完成 perturbation 並為 Stable；不能推斷其 Q0/Q1 狀態 |

### Perturbed endpoint 重新資格化

- 現行實作的已完成 perturbation endpoint 必須通過 Strict primary A/C certificate、Strict Independent Replay KKT 與 replay numerical consistency。因此已保存的成功 endpoint 在同一原始 state 下不會因 Q0→Q3 而新增 acceptance。
- 舊報告記錄 T4 full-policy block-1/2 造成最大反應的 endpoint A/C KKT 為 1.3878e-16、1.9429e-16，projected Width KKT 為 0、5.5511e-16；兩者 Independent Replay pass。它們通過 Q0–Q3。
- 此結論只涵蓋已完成且保存數值的 endpoint。完整逐 seed primary/independent KKT、failed endpoint raw state、各 Q 的 Qualified endpoint 總數未保存；不得將其寫成完整的 endpoint matrix。
- T4 solver-only 的兩個 unavailable block 在 baseline eligibility 階段即停止，perturbation solve 為 0。其 replay failure reason 是 primary-bounded-width-kkt-failed；按程式檢查順序，失敗點是 Width optimality qualification，不是 prediction/residual/objective/gradient replay disagreement。Q3 仍失敗，故放寬 A/C 不能解開這兩個 Unavailable。
- 先前研究記錄的 16 個選中 block summary 為 7 Stable、7 Unstable、2 Unavailable、104 NotAssessed；這是 block summary，不是獨立樣本數。只憑該表不能重算 Q0–Q3 對所有 workload 的 Eligible、Qualified、Replay-failure 總數。

## E. L2 Solver Precision、Null 與 Perturbation

| 實驗 | 本輪狀態 | 可用結果 |
|---|---|---|
| P0 Strict Local solver | 有既有 T0/T4 synthetic chain evidence；本輪 8 個 stability regression tests 通過 | 已完成 endpoints 使用 strict solver/replay；T4 full-policy 兩個方向量到 D_pred 1.80e-7、1.94e-7，endpoint replay 通過 |
| P1 Qualification-only Q0–Q3 | 部分可重建 | T0 30-target baseline block-1/2 僅在 Q2/Q3 eligible；T4 solver-only block-1/2 即使 Q3 仍因 Width KKT 不合格；已記錄 strict qualified perturb endpoint 的 acceptance 不變 |
| Null perturbation | 未完成 | 現行 conditional_reprofile 是固定同一 eta 後重新解 A/C，不是使用與正負 Width perturb 完全相同 solver/budget 的 full-state null continuation。舊研究有一個 conditional reprofile Independent Replay failure，未保存其完整原始欄位 |
| Strict full-state Local reference | 未完成 | SearchProfile() 只接受初始 Width，profile evaluation 會重新計算 A/C；不能把 Width-only restart 說成 (eta,beta) full-state continuation |
| P2 Relaxed Local Solver Precision | 未執行 | 沒有 Local-only stopping tolerance；改 CertifyLinear() 或 CheckReplay() 會觸及共用 profile / candidate acceptance。要安全測 P2 需先建立有界的 Local Assessment solver/qualification seam；本次沒有依據冒險改寫 |

30-target T0 與 T4 的舊報告比較使用不同 sweep endpoint：T0 5 sweeps、T4 4 sweeps。T4 相較 T0 的全域 objective 是 1.51721796e-10 對 1.51707603e-10，相對差約 9.4e-5；13,710 個 prediction rows 的 max/RMS absolute difference 是 3.2855e-7 / 3.6282e-8。這是 endpoint 間差異，不是 null drift estimate。T4 擾動 endpoints 的 objective 變化約 6.5e-15–7.5e-15，而 D_pred 約 1.8e-7–1.94e-7；低 objective difference 不代表 prediction response 小。

因此本研究無法將該 D_pred 拆成真正 perturbation response 與 baseline 重新最佳化 drift，也無法估計 relaxed solver 的 prediction error 是否低於 1e-8。現有資料只能說：在已認證的 T4 端點，已測方向的 response 超過既有 warning；其跨 endpoint 科學意義仍未定。

## F. L3 Local／Global 一致性

1. Baseline eligibility 的 {1e-7,1e-9} 已與 Global practical-v2 相同；Local Width KKT 也使用同一 bound-aware projection。
2. Perturbed endpoint 的 A/C primary certificate 和 Independent Replay 仍要求 1e-10；Independent Replay Width KKT 預設 1e-12，並有額外的 practical 1e-9 Width check。
3. CheckIndependentReplay(criteria) 可分離 optimality threshold 與 prediction/residual/objective/gradient consistency；但 CheckTrust() 無 criteria overload，預設 strict，且 profile accepted trial 的 CheckReplay() 固定 strict。
4. Local 與 Global 的 observation/gradient normalization 可比，Local objective 則是固定鄰居下的條件式 objective。
5. Q3 對 30-target T4 的 A/C baseline 足夠，但 Width KKT 不足。只改 A/C KKT 或把 replay consistency 與 optimality 混在一起，不會解除這個 unavailable。
6. 成功的 T4 perturbed endpoints 已遠低於 Q0；本批沒有證據顯示 strict endpoint optimality 門檻造成 false Unavailable。同時缺少 failed candidate raw KKT，故不能排除其他 workload 有此情況。

## G. L4 Warning Threshold 交互作用

下表重列已存的同一批 16 個 selected block summaries。Unavailable／NotAssessed 固定不參與重分類。

| Parameter threshold | Prediction 1e-8 | Prediction 1e-7 | Prediction 1e-6 |
|---:|---|---|---|
| 1e-8 | 7 Stable, 7 Unstable, 2 Unavailable, 104 NotAssessed | 7 Stable, 7 Unstable, 2 Unavailable, 104 NotAssessed | 7 Stable, 7 Sensitive, 2 Unavailable, 104 NotAssessed |
| 1e-7 | 7 Stable, 7 Unstable, 2 Unavailable, 104 NotAssessed | 7 Stable, 7 Unstable, 2 Unavailable, 104 NotAssessed | 7 Stable, 7 Sensitive, 2 Unavailable, 104 NotAssessed |
| 1e-6 | 7 Stable, 7 Unstable, 2 Unavailable, 104 NotAssessed | 7 Stable, 7 Unstable, 2 Unavailable, 104 NotAssessed | 14 Stable, 2 Unavailable, 104 NotAssessed |

- 對已保存資料，單獨把 parameter warning 提高到 1e-7 或 1e-6 不改變分類；把 prediction warning 提到 1e-6 才使 7 個 Unstable 轉成 Sensitive。parameter warning 再到 1e-6，它們才轉成 Stable。
- Q0–Q3 只影響 KKT qualification；warning thresholds 只改分類。已完成且通過 Q0 的 perturb endpoint 不會因 qualification 放寬增加；Q3 也沒有解除 30-target T4 的兩個 baseline unavailable。
- 本資料顯示分類轉移主要由 Prediction Warning 門檻驅動。把 1e-6 / 1e-6 套在這批高精度 chain evidence 會把全部 14 個有效 summaries 歸為 Stable，但沒有 multi-start、weak-width、noise-level、scale、幅度、near-singular 與真實 map coverage 支持這個升格。

Scientific coverage audit:

| Case | Evidence | 結果 |
|---|---|---|
| Exact synthetic / 9 observations | 舊報告的 0.001×、1×、1,000×、duplicate-row cases | 數值變化接近 machine precision；不足以校準 noisy-data warning |
| Noise / signal scale | 一組 noisy exploration：0.001× baseline/replay pass；scale 1、1,000 與 duplicate-row cases replay failed | 不能形成 scale-invariance 結論 |
| 3/30/300/1,000-target chain | 同 HEAD 的 Local Stability summary 與本輪 3/30/300/1,000-target Production benchmark | 選中 block 的 KKT/replay 與分類見前表；大型 workload 只代表 3 個 blocks |
| High-overlap | 8-target fixture 有舊 Width fit，但沒有同 endpoint Stability matrix | NotAssessed |
| Weak-width / active-width | 小型單元 fixtures 有 weak direction 與 active upper-bound one-sided probe | 只驗證語義，不代表多尺度科學分類矩陣 |
| Noise-only | 既有 fit stability 為 NotAssessed | 不納入 warning 分類 |
| Near-boundary | 無適用 Stability raw measurement | Not Run |
| Rank-deficient / near-singular | 既有 fit 未通過 qualification | Unavailable；不可分類成 Stable |
| Multi-start、不同真實 Width、A/C amplitude 梯度 | 未有完整 evidence | Not Run |
| Real cryo-EM map | repository 無適用 fixture | Unavailable |

## H. Production Promotion Decision

決策：**Partial Qualification**。

| Gate | 狀態 | 理由 |
|---|---|---|
| G1 Numerical Validity | 部分通過，整體未完成 | 已保存成功端點維持 feasibility、strict endpoint gates 與 replay；基線 T4 失格可歸因於 Width KKT。完整 per-seed KKT/replay raw data 未保存，不能完成全部 Q matrix。 |
| G2 Sensitivity Reliability | 未通過 promotion 標準 | 無完整狀態 continuation、null control、P2 precision matrix；不能量化 baseline reoptimization drift 或 relaxed solver 的 1e-8 warning error budget。不同 sweep endpoint 已知會有 measurable prediction difference。 |
| G3 Engineering Safety | 通過本輪修改範圍 | Production solver、criteria、acceptance、replay、Global sweep、JSON 與 point-quality 未變；Release build 和完整 CTest 通過。 |

維持：

- Global practical-v2 A/C KKT 1e-7、Width KKT 1e-9。
- Local baseline eligibility A/C 1e-7、Width 1e-9。
- Local perturbed endpoint strict A/C certificate、Strict Independent Replay 及 practical Width check。
- LocalParameterChangeWarning = 1e-8。
- LocalPredictionChangeWarning = 1e-8。

本研究不支持把所有 Local KKT 門檻簡化成一組 Global 數字。未來可獨立實作 Local Assessment optimality criteria 與 Replay numerical-consistency tolerances 分離的內部介面；promotion 前須先用 same-baseline、full-state null、paired P0/P2、per-seed raw evidence 驗證，不得藉此更動 Global solver。

## I. Cleanup 與最終驗證

- 新增：本研究報告。
- Production／測試程式：未修改；沒有 temporary policy switch、experimental adapter、debug print、exporter、test hook 或 orphan helper。
- 保留：既有純 classifier、KKT boundary tests、Independent Replay failure、Width bound、Local Conditioning、Parameter Stability、Target point-quality 與 JSON tests。
- Benchmark／memory：已用 Release Production closed-loop benchmark 測 3/30/300/1,000 targets，每組 warmup 1 次、量測 5 次、數值函式庫單執行緒；wall median 為 3.711 ms、123.940 ms、977.934 ms、3.302 s，process peak RSS 為 24.1、38.8、91.9、240.5 MB。既有 stability workspace estimate 仍是獨立數據。
- Release build：成功；--parallel 4。
- 完整 CTest：19/19 通過；最終執行 15.87 秒，CTest --parallel 4。
- Local Parameter Stability tests：8/8 通過。
- Numerical equivalence：四種 Production benchmark endpoint 都通過 certification / Independent Replay；因本輪未執行不同 policy 的 paired equivalence，沒有宣稱 T0/T4 endpoint 逐參數等價。既有 truth error 數值見同 HEAD 的 Local Stability 校準報告。
- JSON round-trip：完整 CTest 中既有 persistence tests 通過；schema 未變。
- 未解問題：缺逐 seed qualification raw evidence、同端點 null/full-state continuation、local-only solver precision interface、完整噪聲/尺度/幅度/多起點矩陣、near-boundary 與真實 cryo-EM map。

## 完成條件與解讀

本研究能回答「是否有證據把 Local endpoint KKT 放寬至 Global practical-v2，並同時保證 1e-8 stability warning 的數值可信度？」——目前沒有足夠證據。已確認 Global-equivalent baseline criteria 不會解決 30-target T4 Width KKT 不合格；已確認所有保存的 T4 成功 endpoint 都有很大的 Strict KKT 餘裕。缺少的 null 與 precision 實驗是判斷 solver error 是否接近 Warning 門檻所必需，不能用分類比例替代。因此 Production 維持原契約，結論為 Partial Qualification。
