# Local Parameter Stability Drift/Seed 研究 — P0 基準

日期：2026-10-29。此文件記錄本輪 P0；P1/P2 endpoint 分析及 P3–P6 gate 評估分別見後續研究報告。

## Git 與執行環境

- 初始 HEAD：`97198fd96ccae3d4d8015468ad1f5d6dbe734401`；分支 `develop`；開始時工作樹乾淨，`git diff --check` 通過。
- 本機 `origin/develop` 與 HEAD 相同。`git fetch origin develop` 無法更新 `.git/FETCH_HEAD`（唯讀）；唯讀 `git ls-remote` 因執行環境無法解析 `github.com` 失敗，所以本輪無法核實 live remote 的最新 SHA。
- `git worktree list` 顯示一個原有 detached worktree `/Users/yslian/.codex/worktrees/d178/em_map_analysis`；未讀寫或變更它。未建立 branch/worktree，未 push，未操作 Notion。
- `origin` URL 指向 `https://github.com/yslianTim/rhbm_gem.git`，與指定專案相符。

## Production 資料流審計

Production `AssessLocalParameterStability()` 以 accepted `(η₀, β₀)` 建立 local baseline、做 Independent Replay，並以固定 `η₀` 重新 profile A/C（N1）。它再由寬度方向產生 perturbation seeds，重新 profile A/C，最後以相對 accepted local baseline 的最大 parameter/prediction change 分類。該 Production 函式目前沒有 baseline-η-seeded Null Search（N2），也不保留逐端點 prediction vector。Target Evidence 將 `Unstable`、`Sensitive`、`Unavailable` 維持為 point-quality Warning；`NotAssessed` 仍是 NotAssessed。

現有分類 normalization 已核對：prediction denominator 為 `max(1, ||ŷ₀||₂)`；A/C 參數差用 baseline design 欄位 norm `||X_j(η₀)||₂/context.scale`；width 差為 log-width infinity norm。新分析採 N2 共同 reference design 的 A/C normalization，與舊參數指標分開報告，沒有混用。

## Raw evidence 重建

來源：[`joint-local-stability-2026-10-29.jsonl.gz`](data/joint-local-stability-2026-10-29.jsonl.gz)。以標準函式庫 gzip/JSONL 解讀，共 1,701 records、15 個 global summaries、25 個 selected-block snapshots/classifications，以及 138 個 main (`δ=1e-3`) perturbation endpoints。摘要分類由 25 個 block-classification rows 重算為 15 Stable、8 Unstable、2 Unavailable。

Raw JSON 沒有直接儲存 prediction vectors，但 baseline snapshot 保存完整 effective response、local rows、每 atom 的 row/squared-distance support、scale、η、β、width bound；Null 與 perturbation records 保存完整 endpoint η/β、directions、requested/effective displacement、solver stop/work、KKT、feasibility 與 Replay evidence。因此可按 Production `EvaluateKernel()` 從相同 local rows 重建每個 prediction vector，再重建 residual。endpoint records 未逐筆重複 `scenario`；重建時依 JSONL 中 snapshot/global-summary 的 run 邊界及 block ID 對應。資料有足夠的 null/endpoint state 與 normalization 欄位，無須以先前報告摘要替代。

用重建 vectors/columns 重新計算原分類輸入，25/25 status 與 raw classification 一致；對可對照的舊 prediction/parameter scalar 指標，最大 absolute difference `2.66e-17`。N0/N1/N2、Strict local reference 各有 25 筆；每個 baseline snapshot 都有有效 row geometry/response，η/β 無缺。主要擾動 endpoints 含 Replay 結果；沒有直接的 prediction vector 欄位，因此 vector-level恆等式的獨立重建依賴上述 kernel/geometry，而不是 byte-for-byte 對照已存 vectors。

## Fresh Release B0 baseline

為符合 Production claim budget，benchmark runner 現在在 benchmark-only policy 設定每 fit 最多 8 個 stability block claim（保留 1 個 active-bound 優先名額）。另外輸出 FNV-1a 64 accepted-state fingerprint，涵蓋完整 η、β、prediction、residual、objective；每個 workload 的五次正式量測 fingerprints 完全一致。Production solver、acceptance、KKT、certification、Replay 都沒有修改。

執行命令：

```sh
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  cmake --build build/release --parallel 4
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  build/release/bin/rhbm_fixed_neighbor_closed_loop_benchmark <targets>
```

每個 workload 均 1 warmup + 5 runs、Release、單執行緒。完整輸出及 hashes 保存於 [`joint-local-split-production-benchmark-2026-10-29.txt`](data/joint-local-split-production-benchmark-2026-10-29.txt)。摘要如下：

| Targets | Search wall median | End-to-end wall median | Stability CPU median | Stability solves/evals | Stability workspace estimate | Process peak RSS | Sweeps / accepted blocks | Endpoint / Replay |
|---:|---:|---:|---:|---:|---:|---:|---|
| 3 | 2.803 ms | 14.179 ms | 10.441 ms | 6 / 30 | 5.44 MB | 24.48 MB | 3 / 3；Pass / Pass |
| 30 | 85.811 ms | 340.115 ms | 235.071 ms | 18 / 85 | 23.72 MB | 5 / 15；Pass / Pass |
| 300 | 737.431 ms | 1.142 s | 205.922 ms | 18 / 78 | 53.38 MB | 88.39 MB | 4 / 104；Pass / Pass |
| 1,000 | 2.609 s | 3.622 s | 304.927 ms | 18 / 72 | 130.21 MB | 4 / 336；Pass / Pass |

這是含現行 Stability 的 B0。當前 benchmark log 沒有輸出全域 η/β/ŷ/residual 數值陣列；它保存覆蓋完整 accepted state 的 deterministic fingerprint，以及 objective/KKT/sweeps/certification/replay。前次 JSONL 的 15/15 on/off synthetic comparisons 另有 accepted η/β，prediction、residual、objective、sweeps 逐值相同。

## P0 結果

**Complete。** 程式資料流與 normalization 已審計；原始分類已從 raw records 獨立重建；Release B0 benchmark 已重新執行四種規模，並保存 accepted-state fingerprints、CPU/wall、stability work、workspace/RSS、KKT、objective、sweeps、certification 與 replay 結果。Live `origin/develop` 仍無法查驗是外部網路限制，已標示為未確認。
