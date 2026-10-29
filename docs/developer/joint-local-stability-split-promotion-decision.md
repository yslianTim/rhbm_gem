# Joint Local Stability：科學品質、分類政策與升格決策

日期：2026-10-29。承接 [P0 baseline](joint-local-stability-split-baseline.md) 與 [P1/P2 drift/seed 驗證](joint-local-stability-drift-seed-validation.md)。本報告評估 P3–P6。凡未經本輪直接量測者明確標為未完成。

## P3：困難案例與科學品質

**狀態：Partial。** 現有 production regression fixtures 覆蓋若干重要防線，但沒有把新的 drift/seed metrics 套用到完整 synthetic matrix，故不能據此認證新的 point-quality 政策。

| 案例 | 現有證據 | 本輪 split 指標驗證 |
|---|---|---|
| Exact／小誤差 | `ParameterStability_test.cpp` 的 exact-local-fit deterministic regression；raw run 有 exact/low-noise workloads | 未對整套 ground-truth matrix 重新跑 split metrics |
| Noise-only／weak-width | stability 與 conditioning tests；raw N2 耗盡 update budget，strict reference 與一般 N2 的 `η` 差 `2.398e-5` | 明確保持 seed robustness unavailable；small seed spread 不作識別性證據 |
| Active-width／near-boundary | 既有 one-sided inward probe、zero-displacement rejection 與 active-width endpoint tests；raw noise-only 有 active width bounds | 沒有完成全矩陣的 paired split 量測 |
| Rank deficient／near singular | `Numerics_test.cpp` rank-23 fixture、conditioning exact singular、block-coordinate rank/near-singular recovery failure fixtures | 新 drift/seed classifier 未逐 fixture 驗證 |
| High overlap | 現有 cube fixture condition 約 `1.07`，不能代表病態高重疊 | 未建立/量測本輪要求的高 condition high-overlap case |
| Multi-start | raw evidence 中 50 個 low/high width endpoints 全數 strict replay；prediction 差最大 `3.92e-12`，scaled A/C `5.70e-12` | 沒有發現不同 prediction basin；不代表唯一解或全域可識別 |
| 尺度、A/C amplitudes、真實 map | raw evidence 有多個 synthetic scale/幅度情境；本輪沒有真實 cryo-EM map fixture | 要求的完整尺度矩陣與真實 map 驗證未完成；無可用真實 map 證據 |

現有測試能分辨 rank failure、weak/singular conditioning、active face、Replay failure 與 point-quality Warning；這些防線不得被單一小 `D_seed` 覆寫。已知 multi-start 和 conditioning 數據不足以證明 ground-truth A/C/B 精度、物理唯一性或所有 local minima 的排除，因此 P3 promotion gate 未通過。

## P4：分類政策與 threshold matrix

候選矩陣只對 22 個具有合格 N2 reference 的 blocks 計數；3 個 split diagnostic incomplete blocks 不在可比較分母。舊 C0 status 仍以全部 25 blocks 計數。

| 政策 | 可比較 blocks | Legacy status | Drift warning annotation | Seed warning annotation | Split unavailable/incomplete | 本輪決定 |
|---|---:|---|---:|---:|---:|---|
| C0 Production baseline | 25 | 15 Stable / 8 Unstable / 2 Unavailable；0 NotAssessed | — | — | — | 保留 |
| C1 Split-only | 22 | C0 狀態完全不變 | 7 | 0 | 3 | 僅改善解釋，不改 status |
| C2 Split + calibrated thresholds | 22 | 視候選門檻把舊 Unstable 分為 drift/parameter warning 或 Stable candidate | 0 或 7 | 0 | 3 | 不升格 |
| C3 Quality-aware mapping | 未產生完整 matrix | 未計算 | 未計算 | 未計算 | KKT/Replay/Rank 不合格仍須不可用 | 不升格 |

C2 門檻的實測分類（22 個 qualified split blocks）：

| `T_drift,pred` | `T_parameter` | Stable candidates | Drift warning | Parameter sensitive | Seed warning |
|---:|---:|---:|---:|---:|---:|
| `1e-8` 或 `1e-7` | `1e-8` 或 `1e-7` | 15 | 7 | 0 | 0 |
| `2e-7`、`3e-7` 或 `1e-6` | `1e-8` 或 `1e-7` | 15 | 0 | 7 | 0 |
| `2e-7`、`3e-7` 或 `1e-6` | `1e-6` | 22 | 0 | 0 | 0 |

`T_seed,pred` 在 `1e-8`、`1e-7`、`1e-6` 三種候選皆為零個 warning；本輪 22 個 valid blocks 的最大 seed prediction difference 為 `4.32e-12`。把 drift 門檻設為 `2e-7` 會讓觀測的 `~1.94e-7` drift warning 消失；這是候選門檻的結果，沒有構成採用該門檻的科學理由。較高 `T_parameter=1e-6` 也會把 7 個舊 Unstable 轉成 Stable candidate，證據不足以支持。

因此 C1 最適合作為解釋報告：分開說明 reoptimization drift、seed endpoint robustness、qualification/conditioning，但舊分類與 point-quality 原樣保守保留。Split 不完整仍須標記為 Unavailable/incomplete，不可把已完成 seeds 的小差異當 Stable。

## P5：工程成本與數值等價

**狀態：Partial。** 本輪有新鮮 B0 stability-on benchmark；沒有實作或量測 B1（新增 split diagnostics）與 B2（安全重用）。數字與 fingerprint 在 [P0 benchmark log](data/joint-local-split-production-benchmark-2026-10-29.txt) 及 [P0 report](joint-local-stability-split-baseline.md)。

| Targets | B0 search wall | B0 end-to-end wall | Stability CPU | solves / evaluations | Stability workspace estimate | Peak process RSS |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2.803 ms | 14.179 ms | 10.441 ms | 6 / 30 | 5.44 MB | 24.48 MB |
| 30 | 85.811 ms | 340.115 ms | 235.071 ms | 18 / 85 | 23.72 MB | 35.49 MB |
| 300 | 737.431 ms | 1.142 s | 205.922 ms | 18 / 78 | 53.38 MB | 88.39 MB |
| 1,000 | 2.609 s | 3.622 s | 304.927 ms | 18 / 72 | 130.21 MB | 236.45 MB |

各 workload 為 Release、單數值執行緒、1 warmup + 5 runs；同一 workload 的 accepted η/β/prediction/residual/objective fingerprint 五次相同。B0 之外的早期 stability-off runner wall 約為 3.663 ms、104.409 ms、930.285 ms、3.141 s；和 B0 相較約為 +287%、+226%、+22.8%、+15.3%。此為現有 Stability 的估計增量，非 B1/B2 的增量；benchmark runner 的 claim policy 有記錄，跨 runner 比較仍有少量量測差異。

Workspace estimate 與 process RSS 是不同量：前者是 stability workspace preflight/peak estimate，後者是整個 benchmark process 的 peak RSS。沒有 B1/B2 數據，故新分離診斷是否符合 10–15% wall、10% RSS 候選門檻未知。不能以省略 KKT、Rank 或 Replay 降成本。本輪未修改 Production solver；既有 15/15 Stability on/off synthetic evidence 和本輪 5/5 accepted-state fingerprints 均逐值一致，但沒有新增完整的 B1/B2 equivalence matrix。

## P6：promotion、cleanup 與 gates

**決策：Partial Qualification。** 本輪不 promote Production diagnostics/classification；保留現行 Local KKT、`1e-8` parameter/prediction warning constants、C0 status 語意和 Target point-quality contract。

| Gate / 決策項目 | 結果 | 理由 |
|---|---|---|
| Numerical validity G1 | Partial | 414/414 perturbation endpoints 通過 endpoint KKT/Replay/feasibility；但 split reference 只有 22/25 blocks 合格，2 個 baseline eligibility 不符，1 個 noise-only N2 budget exhausted |
| Diagnostic reliability G2 | Partial / not passed | vector decomposition、同 response/scale/bounds 下 paired comparisons 支持本資料的 drift/seed 分離；P3 困難矩陣、識別性及全方向/真實 map 覆蓋不足 |
| Engineering G3 | Partial | Release build/full CTest 最終結果列於最終驗收；Production solver 未改，B1/B2 成本/equivalence 未跑，不能宣稱新候選已達工程 Gate |
| Drift/seed 欄位升格 | No | 僅研究報告保留，未新增 Production persistence/schema |
| Warning threshold 升格 | No | threshold grid 不足以支持改值；不挑門檻去消除約 `1.94e-7` 的現有 warning |
| Target point-quality 變更 | No | 保留所有現有 Unstable/Sensitive/Unavailable Warning 與 NotAssessed 語意 |
| Local KKT / Global solver | 保留 | 未改任何 Production tolerance、practical-v2、求解/接受/Candidate、certification 或 Replay policy |

### Cleanup audit

- 沒有留下臨時 JSON/CSV exporter、one-off runner、runtime switch、temporary null observer、precision control 或逐 seed console instrumentation。
- 保留 benchmark-only accepted-state fingerprint 與 stability work counters，供後續 deterministic closed-loop benchmark 比較；Production fitting path 未使用此 instrumentation。
- 未變更 Persistence schema，也未刪除原始 gzip 證據；原始端點資料支援報告中的 vectors、qualification、seed endpoint 重建。
- 保留本輪 P0/P1–P2/P3–P6 報告及必要 benchmark log。既有測試包含 null identity、one-sided active-width、zero displacement、Replay/AC failure、conditioning/rank、point-quality Warning、JSON/target evidence regression；本輪未新增 split-metric C++ regression，因未將 split policy promote。
- 沒有新增 Git branches 或 worktrees；沒有未提交研究 runtime hook。原有 detached worktree untouched。

## 實際限制

本輪無法連到 GitHub 核實 live `origin/develop`：fetch 受 `.git/FETCH_HEAD` 唯讀限制，`ls-remote` 因 DNS 無法解析 GitHub 失敗。本地 `origin/develop` ref 保持初始 SHA。未 push。P3–P5 完整研究矩陣（困難 synthetic truth cases、high-overlap conditioning、B1/B2 benchmarking）未執行；因此最終推薦保守維持 Production，而非推論新政策已可靠。
