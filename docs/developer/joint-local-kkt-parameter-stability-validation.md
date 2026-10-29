# Joint Component Local KKT 與 Parameter Stability 數值驗證

## 研究範圍與基準

本研究從 `develop` 的 `8de693f9f6b2e7ddddcf2de1f22a619b3df6b2f8` 開始；`origin/develop` 當時指向相同 commit。初始 working tree 乾淨。沒有建立 branch 或 worktree，也沒有 push 或修改 Notion。唯一既存的 detached worktree 未觸碰。

程式核對確認以下正式數值契約：Global A/C Projected KKT `1e-7`、Global Projected Width KKT `1e-9`、`practical-v2`；Local baseline eligibility `1e-7 / 1e-9`；Strict A/C certificate `1e-10`；Strict Independent Replay `1e-10 / 1e-12`；Parameter／Prediction warning 均為 `1e-8`。本研究不調整這些正式值。

| 階段 | 結果 | 說明 |
|---|---|---|
| P0 原始證據 | 完成 | 25 個 selected block 快照、138 個完成的擾動 seed，含逐端點 KKT 與 Replay 欄位 |
| P1 Null／Reference | 部分完成 | N0、N1、baseline-η-seeded N2 與 30/300/1,000 target reference；SearchProfile 會重新 profile β，未直接注入完整 `(η₀, β₀)` 狀態 |
| P2 Solver Precision | 部分完成 | 實際改變 Local profile LM `gtol`；A/C constrained linear solve 沒有可安全注入的 KKT stopping tolerance |
| P3 Qualification | 完成 | 對 25 個 baseline、25 個 conditional reprofile 與 138 個成功 perturbation endpoint 以同一端點重算 Q0–Q3 |
| P4 困難案例 | 部分完成 | 合成 exact/noise/scale/amplitude/width、cube 與 3–1,000 target；尚無 real map、rank-deficient/near-singular workload |
| P5 Promotion／Cleanup | 完成（決策 Partial Qualification） | 不升格候選門檻；研究 hooks／runner 已清除，保留兩個 null／replay regression tests；Release build 與 19/19 CTest 通過 |

## P0 — Per-seed 原始數值證據

研究工具只連結測試專用的 `rhbm_gem_test`；追蹤結構與觀察回呼在 `RHBM_GEM_TEST_INSTRUMENTATION` 下編譯。Production library 沒有此定義。輸出保留 baseline snapshot、fixed-neighbor effective response、支援列與距離、η/β、方向、requested/effective displacement、端點、目標函數軌跡、solver work、primary certificate、Independent Replay 和 Q0–Q3 重評。無效／不可用值保留為 null，沒有填成零。

原始 JSONL gzip 檔含 15 個 deterministic synthetic 執行案例：3、30、300、1,000 targets；30-target T0、T4 solver-only、T4 full-policy；8-target cube；exact、noise-only、兩種 noise、兩種 signal scale、兩種幅度及不同真實 width。每個 global summary 記錄 baseline HEAD、fixture hash、policy、Release build、單執行緒設定與資源資料。重新執行 T4 full-policy 後，排除時間和 RSS 欄位的 canonical JSONL 逐行完全一致。

原始資料：[joint-local-stability-2026-10-29.jsonl.gz](data/joint-local-stability-2026-10-29.jsonl.gz)。此檔案內每一行都是 JSON record；先 gzip 解壓即可用 JSONL 工具讀取。

四個重點 Global endpoint 均完成 endpoint certification 與 Independent Replay：

| Workload | Sweeps | Global A/C KKT | Global Width KKT | Endpoint / Replay |
|---|---:|---:|---:|---|
| 30-target T0 | 5 | `1.6024e-9` | `1.8807e-11` | 通過 / 通過 |
| 30-target T4 | 4 | `9.4630e-8` | `1.1122e-9` | 通過 / 通過 |
| 300-target T0 | 4 | `2.9920e-8` | `1.1108e-10` | 通過 / 通過 |
| 1,000-target T0 | 4 | `1.6387e-8` | `3.3321e-11` | 通過 / 通過 |

T4 solver-only 設定是 Global `{1e-7,1e-5}`、Local eligibility `{1e-7,1e-9}`；T4 full-policy 僅把測試用 Local eligibility 改成 `{1e-7,1e-5}`。兩種 T4 run 的 accepted η、β 和 objective 逐值相同；eligibility 只改診斷判定。

### P0.5 Conditional reprofile Replay failure

受控重現於 30-target T4 solver-only block 1/2。Baseline A/C KKT 為 `9.4630e-8`／`8.7952e-8`，Width KKT 為 `1.1122e-9`／`1.0308e-9`，因此 Local practical eligibility 在 Width 條件拒絕兩個 block。固定 η 後的 conditional reprofile A/C KKT 降至約 `8.33e-17`／`1.67e-16`，但 Width KKT 為 `5.0148e-9`／`4.6642e-9`。

`conditional_reprofile_status` 回報 `StrictlyQualified`，因該欄位只表示 A/C certificate 通過 `1e-10`。Independent Replay 的完整結果為 Failed，原因 `primary-bounded-width-kkt-failed`；Primary／Independent Width KKT 相同至數值精度。Prediction、residual、objective、gradient agreement 與數值一致性 enclosure 均通過，Width feasible，沒有 active width bound。這是 Width optimality gate 失敗，不是 replay 數值不一致；兩個結果欄位應分開解讀。

## P1 — Null Control 與局部參考

本研究的 N0 直接以完整 baseline `(η₀, β₀)` 呼叫 `EvaluateState`，重算 prediction、residual、objective 與 KKT；N1 固定 `η₀` 並重新 profile A/C；N2 從 `η₀` 開始 `SearchProfile`；Strict Reference 對同一 local domain / response / scale / bound 使用 profile budget 200、update budget 100。`SearchProfile` 只接受 width seed，並在每次 profile evaluation 重新求解 β；raw record 明記 `initial_beta_passed_to_search=false`。因此 N2 是 baseline-η-seeded variable-projection continuation，不是直接把 β₀ 傳給 optimizer 的 full-state optimizer。N1 對 N2 的差值明確量出首輪 A/C reprofile 後的 width-search 增量。

所有下列 prediction change 均相對原 baseline、以 `||Δŷ||₂ / max(1, ||ŷ₀||₂)` 正規化。N2 與同 η 起點的 Strict Reference 在這些 block 收斂到逐值相同端點。

| Workload / Block | N1 固定寬度重 profile `D_pred` | N2 / Strict Reference `D_pred` | N1→N2 `D_pred` | N2 objective change |
|---|---:|---:|---:|---:|
| T0 30 / block 1 | `2.9764e-9` | `3.2813e-9` | `1.3813e-9` | `2.1471e-18` |
| T0 30 / block 2 | `2.7665e-9` | `3.0503e-9` | `1.2847e-9` | `1.8555e-18` |
| T4 30 / block 1 | `1.7577e-7` | `1.9376e-7` | `8.1542e-8` | `7.4871e-15` |
| T4 30 / block 2 | `1.6337e-7` | `1.8011e-7` | `7.5840e-8` | `6.4693e-15` |
| T0 300 / block 7 | `1.7586e-7` | `1.9387e-7` | `8.1592e-8` | `7.4822e-16` |
| T0 1,000 / block 40 | `1.7580e-7` | `1.9380e-7` | `8.1576e-8` | `2.2429e-16` |

T4 full-policy block 1 的 N0 identity prediction/objective difference 是 0；N1 把 A/C KKT 降到 `8.33e-17`，同時造成 `1.7577e-7` prediction change；N2 再把 `η` 移動 `5.2135e-7`，終點 A/C KKT `1.39e-16`、Width KKT 0。N2 與 Strict Reference 完全相同，objective change `7.4871e-15`。block 2 的 N1/N2 prediction change 分別為 `1.6337e-7`、`1.8011e-7`。T0 30-target 相應漂移約 `3.0e-9`，低於目前 `1e-8` warning。

在 T4 block 1，正負 perturbation 的 `D_common` 對 δ=`1e-4, 1e-3, 1e-2` 約為 `1.9376222e-7`，而 coordinate-direction 的 `D_odd` 分別約 `6.94e-16, 6.01e-14, 9.80e-17`。block 2 的 `D_common` 約 `1.8011287e-7`，`D_odd` 最高約 `4.19e-15`。T0 30-target 的 `D_common` 約 `3.05e-9–3.28e-9`。兩側 perturbation seeds 的有效位移相同。這表示本組對稱擾動的反應主要是共同的 baseline reoptimization 位移；odd response 很小且不隨 δ 呈一階成長。此分解是該方向、該 solver 與已測 δ 的觀測，不是唯一的模型敏感度／solver drift 證明。

300-target block 7 與 1,000-target block 40 也量到約 `1.94e-7` 的 N2 baseline-to-reference drift；N1 已解釋約 `1.76e-7`。因此 T0 global policy 本身不足以保證固定鄰居 local baseline 已經 profile 收斂。

N0 identity 的 prediction/objective change 在所有 selected baseline 均為精確 0；它不會掩蓋 baseline 的資格問題。例：T0 30 block 1 的 Primary A/C／Width KKT 為 `1.6024e-9 / 1.8807e-11`，Strict Replay 因 A/C KKT 未達 `1e-10` 而 fail；T4 30 block 1 為 `9.4630e-8 / 1.1122e-9`，Strict Replay 也 fail，但 prediction、residual、objective 與 gradient 一致性通過。N2 Strict Reference 對這兩個 block 的 A/C／Width KKT 分別降至約 `1.25e-16 / 5.55e-16`、`1.39e-16 / 0`，Replay 通過；N2 各使用 1–2 accepted updates、2–3 profile evaluations，未耗盡 budget。

T0/T4 30-target block 1 的 Strict Reference 終點也很接近：兩端 η 的最大差 `8.83e-9`、β 最大絕對差 `6.44e-8`；同一 local design 上兩組預測的 normalized L2 差為 `3.28e-9`，RMS row 差 `7.15e-10`。兩組 effective response 本身仍相差 `9.21e-9`（normalized L2），所以這是相近資料下的端點比較，不是完全相同 response 的同一最佳化問題。

## P2 — Local Solver Precision 對照

Q 與 S 分開執行。Q0–Q3 只對已保存端點重算資格，不變動 solver。S 矩陣則在同一 T4 block 1 baseline-η-seeded N2 null、同一 local domain/response/scale/bounds、profile budget 11、update budget 6 上，只注入 `SearchProfile` 的 LM `gtol`；`ftol=1e-14`、`xtol=1e-12` 固定。Strict Reference 使用 production 預設值，profile/update budget 200/100。下表數字比較各 S 端點與 Strict Reference；objective 是 normalized profile objective。

| 設定 | 實際 LM gtol | A/C KKT | Width KKT | `D_B` vs ref | `D_AC` vs ref | `D_pred` vs ref | `ΔF` vs ref | Updates / evaluations | CPU / reference wall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Strict Reference | `1e-12` | `1.39e-16` | `0` | `0` | `0` | `0` | `0` | 2 / 3 | CPU 未獨立記錄（wall `14.22 ms`） |
| S0 | `1e-12` | `1.39e-16` | `0` | `0` | `0` | `0` | `0` | 2 / 3 | `16.19 ms` |
| S1 | `1e-9` | `1.39e-16` | `0` | `0` | `0` | `0` | `0` | 2 / 3 | `16.45 ms` |
| S2 | `1e-8` | `1.39e-16` | `0` | `0` | `0` | `0` | `0` | 2 / 3 | `16.29 ms` |
| S3 | `1e-7` | `8.33e-17` | `3.00e-14` | `3.12e-12` | `6.54e-13` | `4.88e-13` | 約 `5.8e-24` | 1 / 2 | `11.69 ms` |

S0–S2 的 endpoint 逐值相同；S3 提早一個 accepted update 停止，但其 prediction solver error 約為 warning `1e-8` 的 `4.9e-5` 倍。S3 相對 Strict Reference 的 max row prediction difference `4.07e-12`、RMS `1.06e-13`、P95 `1.82e-13`、P99 `5.79e-13`；Independent Replay 仍通過。四種設定 active A/C constraints 與 active width-bound targets 均為 0。已測案例中，LM stopping precision 沒有解釋 T4 的 `~1.9e-7` common-mode 位移，也沒有接近 `1e-8` warning。

這是實際的 Width profile solver stopping precision 對照，但不是 A/C linear solver tolerance 對照。程式審計發現 `CertifyLinear()` 的 A/C KKT certificate 固定以 `1e-10` 設定 `kkt_passed`；`SolveLinear()` 的 rank policy、active-set 與 refinement controls 不提供獨立的 local KKT target。把這個 certificate 門檻改小／改大會改變共同使用的 evaluation 與 candidate qualification，不會單純改變局部線性求解精度。故本輪不把 A/C KKT threshold 冒充 A/C solver precision，也沒有放寬 Certificate、`CheckReplay()`、`CheckTrust()` 或 Candidate Acceptance。

Precision 個別 run 沒有配置獨立 RSS recorder；T4 30-target process peak RSS 為 `41.2 MB`，local-stability workspace estimate 為 `23.7 MB`。這兩個是整個診斷執行的資源資料，不是每一 gtol 設定的增量記憶體。P2 因此只對 Width LM precision 完成受控矩陣，A/C solver precision 與 per-setting memory 仍未定量。

## P3 — 同端點 Q0–Q3 Qualification Matrix

下表以同一份 raw endpoint 的 Primary／Independent A/C KKT、Primary／Independent Width KKT、bound feasibility 與 replay consistency 計算。Q 放寬不重跑 solver。Qualification 必須同時通過 optimality 與 replay consistency；prediction、residual、objective、gradient 一致性不會隨 Q 門檻放寬。

| Endpoint 類型 | 端點數 | Q0 qualified | Q1 qualified | Q2 qualified | Q3 qualified | Q3 仍失敗原因 |
|---|---:|---:|---:|---:|---:|---|
| Baseline | 25 | 6 | 14 | 16 | 21 | 4 個 T4 block 1/2 baseline 的 Width KKT 略高於 `1e-9` |
| Conditional reprofile | 25 | 6 | 21 | 21 | 21 | 同 4 個 T4 block 1/2 reprofile Width KKT 在 `4.66e-9–5.01e-9` |
| 完成的 perturbation seed | 138 | 138 | 138 | 138 | 138 | 無 |

`Q0={1e-10,1e-12}`、`Q1={1e-10,1e-9}`、`Q2={1e-8,1e-9}`、`Q3={1e-7,1e-9}`，數值依序為 A/C、Width。25 個 baseline/reprofile 的 replay numerical consistency 欄位全部通過；baseline run-criteria replay 中僅 T4 solver-only block 1/2 因 Width KKT fail。完成的 138 個 seeds 全部 linear solver valid、feasible、width feasible、Independent Replay 通過，並在最嚴格 Q0 通過。132 個 seed 沒有 active A/C 或 width bound；另 6 個 noise-only seed 有 2 個 active A/C constraints 與 2 個 active width bounds，也通過各自 replay consistency 和 Q0 optimality。Raw endpoint 中沒有 rank invalid、nonfinite、budget-exhausted 或 replay-disagreement 的 perturbation seed。

Q3 相對 Q0 多讓 15/25 個 baseline 和 15/25 個 conditional-reprofile 快照通過；它沒有讓任何新的 perturbation endpoint 合格，因為已完成 endpoints 原先都已通過 Q0。T4 full-policy 雖用 Local eligibility width `1e-5` 完成 block 1/2 的診斷，但 Q3 仍拒絕它們：兩個 baseline Width KKT 是 `1.1122e-9`／`1.0308e-9`，兩個 N1 Width KKT 是 `5.0148e-9`／`4.6642e-9`。因此放寬 A/C qualification 或 T4 Local eligibility 都沒有理由改動 Global `1e-9` Width criterion。

## P4 — 科學品質、尺度與資源覆蓋

本輪 15 個 synthetic run variant 包含 3/30/300/1,000 targets、T0/T4 solver-only/full-policy、8-target cube topology、無噪聲資料、兩個 noise amplitude、noise-only、兩種 signal scale、不同 A/C amplitude 與 true width `0.8`。每個 variant 的 accepted state 在「有／無 stability assessment」對照中，η、β、prediction、residual、objective 與 sweep count 都逐值相同（15/15）。

| Fixture | Truth／擾動 | 最大 width error | 最大 A error | 最大 C error | Prediction RMSE | 局部結果／觀察 |
|---|---|---:|---:|---:|---:|---|
| exact 3 targets | 無噪聲，B=.5, A=2, C=.2 | `5.45e-9` | `2.42e-8` | `7.00e-10` | `5.65e-10` | Well-conditioned，Stable |
| low noise | noise amplitude `1e-6` | `3.77e-7` | `1.93e-6` | `9.09e-8` | `3.22e-8` | Well-conditioned，Stable |
| higher noise | noise amplitude `1e-5` | `3.72e-6` | `1.91e-5` | `9.03e-7` | `3.20e-7` | Well-conditioned，Stable |
| noise-only | true A=C=0，noise `1e-5` | `0.720`（不可解讀為真實 width 誤差） | `5.58e-6` | `5.23e-7` | `1.67e-7` | ActiveFace；6 seeds 均有兩個 width bound 與兩個 A/C active；Unstable、parameter change `.763` |
| small scale | signal scale `.001` | `5.45e-9` | `2.42e-11` | `7.00e-13` | `5.65e-13` | Width recovery 與大尺度相同，絕對 A/C 誤差按尺度縮放 |
| large scale | signal scale `1000` | `5.45e-9` | `2.42e-5` | `7.00e-7` | `5.65e-7` | 絕對誤差按尺度放大，condition estimate 幾乎不變 |
| wider true width | B=`.8` | `2.64e-12` | `1.82e-11` | `3.75e-13` | `8.16e-14` | 11 sweeps，Well-conditioned、Stable |
| lower amplitudes | A=`.2`, C=`.02` | `5.45e-9` | `2.42e-9` | `7.00e-11` | `5.65e-11` | Well-conditioned、Stable |
| cube 8 targets | 中等空間重疊 | `8.92e-6` | `6.96e-5` | `1.36e-5` | `6.15e-6` | full rank；projected width condition 約 `1.073`，不足以代表病態 overlap |

Noise-only 的 Width 不具真值識別性，upper-bound active 且出現 Unstable，說明 Parameter Stability warning 不能單獨解讀成可識別參數的物理敏感度；應搭配 target point-quality、active face 與 rank/conditioning。low/high noise 與尺度測試皆有合理的誤差縮放，但不構成多 map 統計驗證。

### 擾動大小與不同初始值

T4 full-policy block 1 在 coordinate direction 的 `D_common` 對 δ=`1e-4 / 1e-3 / 1e-2` 分別是 `1.9376222241e-7 / 1.9376222238e-7 / 1.9376222240e-7`；`D_odd` 分別是 `6.94e-16 / 6.01e-14 / 9.80e-17`。block 2 同樣呈現 `D_common≈1.80112867e-7`，`D_odd<4.2e-15`。在所測範圍內 common response 幾乎不依賴 δ，符合 baseline reoptimization offset；沒有觀察到一階成長的方向性 odd response。upper-bound projection 與非線性仍限制此結論。

另外對 T4 3 個 selected blocks 分別從 baseline η 的 ±`.25` log-width offsets 啟動 multi-start；六個 endpoint 都有 Strict Replay。相對 strict reference，`D_B` 約 `2.2e-16–3.3e-16`、`D_AC` 約 `7.5e-17–1.8e-16`、`D_pred` 約 `2.5e-16–4.7e-16`。在這六個起點未見不同局部極小值，但不能據此排除其他 component 或更遠初值的多模態。

### 覆蓋率與 production benchmark

跨 15 個 synthetic variants 的 25 個 selected-block assessment 為 15 Stable、8 Unstable、2 Unavailable；沒有 Sensitive/NotAssessed。兩個 Unavailable 是 T4 solver-only block 1/2。這些是跨 scenario 的 selected-block 結果，不是獨立科學樣本。300 targets 選 3/26 個 local blocks（11.5%），1,000 targets 選 3/84（3.6%），不能當成整個 component 的 stability 證明。

Local perturbation diagnostic 的 measured CPU 約為 T0 30 targets `0.266 s`、300 targets `0.219 s`、1,000 targets `0.311 s`；對應 workspace estimates 為 `23.7 MB`、`53.4 MB`、`130.2 MB`。Workspace 是估算值；research runner 的 process RSS 高水位包含 stability-on/off 兩次 fit，不能解讀成診斷增量。另以 Release、單執行緒、warmup 1 + measured 5 runs 執行 production closed-loop benchmark，完整結果在 [joint-local-production-benchmark-2026-10-29.txt](data/joint-local-production-benchmark-2026-10-29.txt)：

| Targets | Median search | Median wall | Peak RSS | Memory preflight estimate | Sweeps | Endpoint / Replay |
|---:|---:|---:|---:|---:|---:|---|
| 3 | `2.83 ms` | `3.69 ms` | `24.1 MB` | `0.45 MB` | 3 | 通過 / 通過 |
| 30 | `89.5 ms` | `108.6 ms` | `35.5 MB` | `3.93 MB` | 5 | 通過 / 通過 |
| 300 | `0.753 s` | `0.953 s` | `89.7 MB` | `36.7 MB` | 4 | 通過 / 通過 |
| 1,000 | `2.569 s` | `3.255 s` | `242.7 MB` | `121.6 MB` | 4 | 通過 / 通過 |

Production Release benchmark 的 accepted objective/KKT、sweep count、candidate counts 均符合 T0 synthetic baseline，cache replay error 為 0。Production `.git` artifact 中沒有 cryo-EM `.mrc/.map/.ccp4` fixture；real-map validation 未執行。獨立 Rank-deficient、near-singular 與 high-overlap hard workload 也未完成；cube fixture condition 約 1.07，不能充當這些案例。

## P5 — Promotion、Cleanup 與最終驗收

### C. 30-target 核心異常

T4 solver-only block 1/2 baseline Width KKT 為 `1.1122e-9`／`1.0308e-9`，略高於 Local baseline eligibility `1e-9`，因此是 `Unavailable`；這不是 perturbation endpoint 的 replay disagreement。T4 full-policy 讓這兩個 block 可被量測，但其 N2/reference 相對 accepted baseline 的 prediction change 是 `1.9376e-7`／`1.8011e-7`。N0 identity 的 prediction/objective change 精確為 0；N1 fixed-width reprofile 是 `1.7577e-7`／`1.6337e-7`；N1 到 N2 的額外 change 是 `8.1542e-8`／`7.5840e-8`。T0 相應 N2 change 約 `3.05e-9–3.28e-9`。

T4 block 1 的 N2/reference 相對 accepted baseline `η` 最大 change 約 `5.2135e-7`、objective change `7.4871e-15`，A/C KKT 約 `1.39e-16`，Width KKT 為 0；block 2 同樣收斂到該條件 local profile reference。N2 與 Strict Reference 的已測 endpoint 相同，且未耗盡 solver budget。正負 pair 在三個 δ 上的 `D_common` 約 `1.9376e-7`，`D_odd` 最大約 `6.01e-14`；這表示所測回應主要是共同 baseline reoptimization drift，但不能視為唯一的模型敏感度分解。六個 `±0.25` log-width multi-start endpoint 通過 Strict Replay，且彼此相符；沒有在這些起點觀察到另一個局部極小值。這些測試未顯示 T4 perturbation 的 active-set 跳變。

### D. Local KKT Qualification

同端點 Q0–Q3 重評維持 replay consistency 獨立檢查：25 個 baseline 中 6/14/16/21 個在 Q0/Q1/Q2/Q3 通過；25 個 conditional reprofile 中為 6/21/21/21；138 個成功 perturbation seed 在全部 Q0–Q3 都通過。四個 T4 block 1/2 baseline/reprofile 仍被 Q3 的 Width `1e-9` 拒絕，其 N1 Width KKT 約 `4.66e-9–5.01e-9`。已完成 perturbation endpoints 的 Primary／Independent KKT 與 Replay consistency 都通過；沒有因真正 replay disagreement、rank invalid 或 budget exhaustion 而失敗的 seed。Q3 只新增 15 個 baseline 與 15 個 reprofile classification，沒有新增 perturbation endpoint qualification。

### E. Solver Precision

已完成的隔離 precision 對照只控制 Width profile LM `gtol`；沒有把 AC certificate threshold 誤稱為 A/C linear-solver precision。S0/S1/S2 與 Strict Reference 完全相同；S3 (`gtol=1e-7`) 相對 reference 的 `D_B=3.12e-12`、`D_AC=6.54e-13`、`D_pred=4.88e-13`，CPU `11.69 ms`，仍通過 Independent Replay。這比原始 `1e-8` prediction warning 小約五個數量級。S0 CPU `16.19 ms`，S1 `16.40 ms`，S2 `16.27 ms`。Precision 各設定沒有獨立 memory recorder，故不報每設定記憶體；T4 整體 process peak RSS `41.2 MB`、workspace estimate `23.7 MB`。

### F–G. Null Perturbation 與 Stability 分布

N0、N1、N2、正負 perturbation 與 Strict Reference 的定義、逐端點 KKT、objective、solver work 及 replay 狀態詳見 P0–P3；N2 是 baseline-η-seeded variable-projection continuation，每次 profile evaluation 都重求 β。因 API 不接受 `(η₀,β₀)` 完整 optimizer seed，本報告不把它稱為直接 full-state 初始化。N1 隔離量出先重新 profile A/C 的 prediction contribution；對稱正負 seeds 的 common/odd 分解再量出所測 δ 下的共同 solver reoptimization 位移與 odd response。

25 個 selected-block assessment 為 15 Stable、8 Unstable、2 Unavailable；沒有 Sensitive 或 NotAssessed。兩個 Unavailable 是 T4 solver-only block 1/2。Noise-only fixture 有弱 Width 識別、active A/C/Width bounds 與 Unstable 分類，表明 warning 要搭配 point-quality、active face 與 conditioning 解讀。S3 precision 沒有接近 `1e-8` warning；Q3 沒有新增任何 perturbation endpoint；因此本輪不調整 warning 門檻，也不把 solver drift 宣稱為真實模型敏感度。

### H. Promotion 決策

| Gate | 結果 | 依據 |
|---|---|---|
| G1 Numerical Validity | 在已測合成 endpoint 通過 | KKT、feasibility、replay consistency、Independent Replay 與 failure reason 都有保存；錯誤的 Width optimality endpoint 被明確標為 failed，未當成 qualified perturbation |
| G2 Sensitivity Reliability | 未通過 | 直接完整 `(η₀,β₀)` optimizer initialization 未建立；A/C solver precision 無獨立 stopping-tolerance seam；尚缺 real map、high-overlap、rank-deficient 與 near-singular 驗證 |
| G3 Engineering Safety | 通過 | Global accepted state 在 15/15 對照逐值相同；Global sweep count、正式 KKT、acceptance、certification 與 replay contract 未變；清理後 Release build 與全套 CTest 通過 |

整體結論為 **Partial Qualification**。維持 Global A/C KKT `1e-7`、Global Width KKT `1e-9`、`practical-v2`、Local strict endpoint qualification 與原 Warning `1e-8 / 1e-8`。本研究不能回答所有資料情境下的 drift/sensitivity 比例，也不支持升格新的 Local KKT 或 Warning policy。

### I. Cleanup 與 Verification

Temporary research hooks、precision injection seam、observer、one-off runner 與 CMake target 均已移除；Production 沒有研究 runtime switch。保留兩個長期 regression tests：N0/N1/profiled-reference exact-fit identity，以及「Conditional reprofile AC strict status 不蓋過 Width replay failure」。Release 全量建置成功，CTest 19/19 通過；本次只在既有 `develop` 建立本地 commits，沒有 push、建立 branch/worktree 或操作 Notion。Production benchmark 四個規模均完成 endpoint certification 與 Independent Replay；1,000-target 最終量測 peak RSS `242.7 MB`、memory preflight estimate `121.6 MB`。尚未解決的項目是完整 β-seeded local optimizer 接口、A/C precision 控制，以及 real-map／病態與高重疊資料驗證。
