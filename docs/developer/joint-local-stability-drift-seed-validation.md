# Joint Local Stability：Baseline Drift 與 Seed Robustness 驗證

日期：2026-10-29。此報告承接 [P0 baseline](joint-local-stability-split-baseline.md)，重建既有 raw endpoints 的 prediction vectors，完成 P1/P2 數值分離。資料來源為 [`joint-local-stability-2026-10-29.jsonl.gz`](data/joint-local-stability-2026-10-29.jsonl.gz)；Raw schema 與 fresh B0 measurements 見 P0。

## 定義與端點配對

`s₀=(η₀,β₀)` 是 accepted global endpoint 在 fixed-neighbor local domain 上的 accepted state；`η=log B`、`β=(A,C)`。N1 固定 `η₀` 並 profile A/C。N2 以 `η₀` 做 variable-projection `SearchProfile`，每次 evaluation 都重新 profile β；raw evidence 明記 `initial_beta_passed_to_search=false`，所以 N2 不是從完整 `(η₀,β₀)` 注入的 full-state search。perturbation seed 也是以寬度 seed 啟動並重新 profile A/C。

每一比較沿用相同 local support/rows、effective response、observation scale、A/C constraints、width upper bound 及 solver/replay policy。T4 full-policy 是研究用 local baseline eligibility 寬限到 `1e-5` 的診斷 run；其 accepted global state 與 T4 solver-only 逐值一致。正式 Local eligibility 仍是 A/C `1e-7`、Width `1e-9`，因此 T4 solver-only blocks 1/2 仍是 Unavailable。

從保存的 squared-distance supports 與 Production `EvaluateKernel()` 重建每一端點 prediction vector。以 `M=max(1,||ŷ₀||₂)` 計算 prediction changes；`D_reopt^AC` 與 `D_seed^AC` 使用同一個 `X(η_N)` reference column norm，`D_reopt^B`/`D_seed^B` 使用 log-width infinity difference。這與 P0 的舊 baseline-column A/C normalization 分開。

## Qualification 與涵蓋率

舊分類由 low-level endpoints 重建後仍是 15 Stable、8 Unstable、2 Unavailable（25/25 一致）。138 個 `δ=1e-3` perturbation endpoints 均通過其 endpoint eligibility Replay、Strict Independent Replay 與 Q0 KKT。另有 `δ=1e-4`、`1e-2` 的 276 個 endpoints；三個 displacement 的 414/414 endpoints 都通過上述 endpoint gates、A/C/width feasibility 與 Independent Replay。

Seed Robustness 還要求 N2 Null Reference 本身 resource-qualified。22 個 selected blocks 的 N2 通過 KKT、Replay 且未耗盡 budget，共 132 個 main seeds（每 block 3 directions × 正負）。兩個 T4 solver-only blocks baseline 不合格，無 seed runs；noise-only block 的 default N2 hit accepted-update budget，雖然端點 Replay/KKT 通過，仍不作為正式 `D_seed` reference。該 case 的高 budget Strict Reference 不耗盡 budget，但 η 與 default N2 相差 `2.398e-5`；故這 6 個原本合格的 perturbation seeds 不得用來宣稱 seed robustness。split diagnostics 的完整分母是 22/25，不是 25/25。

## P1：Drift 與 Robustness 數值

`D_profile=||ŷ₁−ŷ₀||/M`，`D_reopt=||ŷ_N−ŷ₀||/M`，`D_width-reopt=||ŷ_N−ŷ₁||/M`；`D_original` 是合格 perturbation endpoints 對 `s₀` 的最大 prediction change；`D_seed` 是同一組端點對合格 N2 reference 的最大 prediction difference。下表使用 `δ=1e-3`；`D_seed^B`/`D_seed^AC` 對 132 個有效 main seeds 取 max。

| Run / block | D_profile | D_reopt | D_width-reopt | D_original | D_seed | D_reopt^B / D_reopt^AC | D_seed^B / D_seed^AC |
|---|---:|---:|---:|---:|---:|---:|---:|
| T0 30 / 1 | 2.9764e-9 | 3.2813e-9 | 1.3813e-9 | 3.2813e-9 | 5.47e-14 | 8.83e-9 / 3.01e-9 | 2.29e-13 / 4.75e-14 |
| T0 30 / 2 | 2.7665e-9 | 3.0503e-9 | 1.2847e-9 | 3.0503e-9 | 2.31e-14 | 8.21e-9 / 2.80e-9 | 8.92e-14 / 1.85e-14 |
| T4 full 30 / 1 | 1.7577e-7 | 1.9376e-7 | 8.1542e-8 | 1.9376e-7 | 1.20e-13 | 5.21e-7 / 1.78e-7 | 7.73e-13 / 1.62e-13 |
| T4 full 30 / 2 | 1.6337e-7 | 1.8011e-7 | 7.5840e-8 | 1.8011e-7 | 2.03e-14 | 4.85e-7 / 1.65e-7 | 6.16e-14 / 1.28e-14 |
| T0 300 / 7 | 1.7586e-7 | 1.9387e-7 | 8.1592e-8 | 1.9387e-7 | 3.84e-13 | 5.22e-7 / 5.62e-8 | 2.38e-12 / 1.57e-13 |
| T0 1000 / 40 | 1.7580e-7 | 1.9380e-7 | 8.1576e-8 | 1.9380e-7 | 5.40e-14 | 5.22e-7 / 3.07e-8 | 2.29e-13 / 8.23e-15 |

全部表列 N2 endpoints 通過 Strict Independent Replay；A/C KKT 約 `3e-17–1.4e-16`、Width KKT `0–3e-15`。N2 normalized objective difference：T0 30 block1/2 `2.1e-18 / 1.9e-18`；T4 full 30 block1/2 `7.5e-15 / 6.5e-15`；300 block7 `7.5e-16`；1,000 block40 `2.2e-16`。N2 與高 budget Strict Reference 在除 noise-only 外的 24/25 block η、β 完全相同；noise-only 兩 endpoint 的 normalized prediction 差只有 `1.15e-10`，objective 近同，但 width 仍不具識別性。

結果回答 P2 核心比較：T4 30 block1/2 的 N2 drift 仍是約 `1.94e-7 / 1.80e-7`，而 Seed Robustness 只有約 `1.2e-13 / 2.0e-14`；T0 30 block1/2 drift 約 `3.28e-9 / 3.05e-9`，seed 差約 `5.5e-14 / 2.3e-14`。300-target block7 與 1,000-target block40 有相同型態：reopt 約 `1.94e-7`，seed 約 `3.8e-13 / 5.4e-14`。因此該批 Unstable 並非由不同 perturbation seeds 落到不同 prediction endpoints 所主導，而是 accepted baseline 到 local reoptimization reference 的共同位移；這仍不證明全域可識別性。

### Vector 分解與正負配對

對每一 perturbation vector 逐列驗證：

`(ŷ_k−ŷ₀) = (ŷ_N−ŷ₀) + (ŷ_k−ŷ_N)`。

在重建 vectors 上，最大 normalized decomposition residual `2.6e-23`。逐個向量滿足 triangle inequality；其 scalar norms 不可相加等同，因為 drift 與 seed residual 可同向或抵消，而且最大值可能來自不同 endpoint。代表性座標方向正負配對如下：

| Run / block | D_common | max D_odd | max D_common|null | δ values |
|---|---:|---:|---:|---|
| T0 30 / 1 | 3.28125e-9 | 6.92e-16 | 8.25e-15 | 1e-4, 1e-3, 1e-2 |
| T4 full 30 / 1 | 1.93762222e-7 | 6.01e-14 | 6.01e-14 | 1e-4, 1e-3, 1e-2 |
| T4 full 30 / 2 | 1.80112867e-7 | 2.83e-16 | 7.88e-16 | 1e-4, 1e-3, 1e-2 |
| T0 300 / 7 | 1.93867348e-7 | 5.70e-14 | 3.72e-13 | 1e-4, 1e-3, 1e-2 |
| T0 1000 / 40 | 1.93804020e-7 | 2.93e-16 | 2.34e-14 | 1e-4, 1e-3, 1e-2 |

T4 block1 的 common shift 幾乎不隨 δ 變動；相對 N2 的 common residual 仍只有 `~6e-14`。這是此方向與這三種位移的觀察，不能推廣成所有模型敏感度都為零。其餘兩個 direction families 也沒有顯示比代表性 seed robustness 更大的方向差異。

### δ、方向、active face 與精度

T4 block1 的 D_seed max 對 δ=`1e-4 / 1e-3 / 1e-2` 為 `1.96e-15 / 1.20e-13 / 1.96e-12`；block2 為 `1.04e-15 / 2.03e-14 / 3.82e-16`。T0 30 block1 為 `8.38e-15 / 5.47e-14 / 1.96e-12`。以上遠低於候選 `1e-8` Seed threshold，且沒有隨 δ 單調放大或改變定性結論。22 個 valid N2 blocks 的最大 `D_seed` 在所有 direction/sign/δ 均不超過 `4.32e-12`。

414/414 perturbation endpoints 的 primary/independent AC KKT 最大 `4.45e-16`，Width KKT 最大 `4.30e-13`，Objective replay difference 最大 `3.69e-21`；所有 Replay、finite、feasibility 與 width-bound checks 通過。每個 δ 有 138 endpoints：132 個沒有 active AC/width constraint；6 個 noise-only endpoints 有 2 個 active AC constraints 及 2 個 active width bounds。T4 表列 endpoints 沒有 active-set changes。main δ 有 3 個、`1e-4` 有 5 個 endpoint 在 accepted-update budget stop，但最終仍通過 KKT/Replay；這些 stop reason 保留在 evidence，沒有把「達 budget」改寫成自然收斂。

高 budget Strict Reference (profile/update 200/100) 與 LM `gtol=1e-12, 1e-9, 1e-8, 1e-7` 候選對照已在 raw evidence 中。T4 precision run 的 S0–S2 endpoint 逐值一致；S3 (`gtol=1e-7`) 相對 reference 最大 `D_pred=4.88e-13`、`D_AC=6.54e-13`、`D_B=3.12e-12`，Replay 通過，遠小於目前 warning。完整 A/C solver tolerance 並未被獨立調整，也未放寬共用 `CertifyLinear()` 或 `CheckTrust()`。

50 個 low/high width multi-start endpoints 均通過 Strict Replay。相對 raw reference，最大 prediction difference `3.92e-12`、scaled A/C difference `5.70e-12`；全組最大 log-width difference `4.47e-8`，出現在 weak/noise 情境，未見明顯不同 prediction basin，但這不能排除其他起始值或真實 map 上的多模態。

## P1/P2 Gate 結論

- **P1：Partial。** 30-target T0/T4（full-policy block1/2）、300-target block7、1,000-target block40/1 完成 N0/N1/N2、vector drift/seed 分解、KKT/Replay、objective、work。兩個 T4 solver-only block因 baseline width KKT 不合格為 Unavailable；noise-only N2 預設 search 耗盡 update budget，因此 6 個 seed 不納入 seed robustness。
- **P2：Complete（qualified paired matrix）。** 3 個 δ、正負方向、N2 reference、strict/high-budget 與 LM precision comparisons 均有實際端點資料；無效/null 不合格情形保留為 Unavailable。此完成狀態只表示已測 qualified paired matrix，不等於所有 synthetic 或 physical cases 已驗證。
