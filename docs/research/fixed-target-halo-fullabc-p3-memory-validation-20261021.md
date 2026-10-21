# Fixed-target Halo FullABC P3：記憶體界線與擴展性稽核

日期：2026-10-21
分支：`develop`  期別：P3

## 程式記憶體路徑稽核

| 階段 | 目前配置/估算 | 稽核結果 |
|---|---|---|
| OriginalInput snapshot | observations、support vectors、Atom/row IDs、selection provenance；`EstimateSnapshotBytes` 按容量逐項 checked arithmetic | 保留完整輸入以維持 identity、normalization 與 replay；Halo path 另建 effective-response buffer，這是目前可見的完整 rows 拷貝 |
| Problem/structural construction | `ProblemConstruction`、`StructuralConstruction` 估算 row incidence、identity/local maps、union-find、core traversal state | `PlanMemory` 在 construction 前檢查 3.2 GB planning budget，overflow 與超 budget 分別返回明確原因 |
| Halo preclassification/global median | 逐 contributor 計算 OriginalInput distinct rows，唯一事前規則仍是 `n_h<3` | n=0 Halos仍保留於重建 contributor universe；n>=3 均先進 FreeFullABC，不新增 rank/chemistry pruning |
| Target/fixed-Halo subtraction | 依去重 support rows 扣除 frozen Target 與 pre-fixed median prediction | 暫存有效 response 與原始 row identity 同維；符合 global scaling/replay，遠端 n=0 halos不加入 free solver support |
| Free support layout/FixedNeighbor | `FixedNeighborPreparation` 依 affected rows/support entries估算；每 component用 compact parent/local row mapping | 不建立 observation×Halo dense matrix或每個 block 的全域 row map；Component domain使用 `DomainAtoms::Select` 的 sparse local support index |
| Local constrained A/C / Active-set LM | `LocalConditioning` 估算最大 block rows/entries、neighbor atoms、free columns、width columns與 8192 row tiles | dense factor 受最大局部 block限制；無全域 covariance。Search每個 local block仍會配置小型 Eigen vectors/matrices |
| Endpoint replay | `EndpointVerification` 估算 full residual/compensated prediction buffers與 O(P) vectors | replay掃描 OriginalInput 全域以保留獨立語意；不建立 N×P design |
| Recovery | 每 Halo、component-local、或固定 median/refit | recovery workspace仍依既有各階段 planner；本輪未增加跨 component的同時 resident work |
| Snapshot / Peeling / Local / Group | postprocessing本身有獨立 payload估算與既有 profiling | 本輪沒有新增儲存契約；production snapshot仍是 `geometry-local-target-median-v1` |

`MemoryShape` 使用 checked add/multiply；budget target為十進位 3.2 GB，hard limit為十進位 4.0 GB。錯誤理由含 `size-overflow` 與 `estimated-peak-exceeds-budget`。`ProcessPeakRssBytes()` 使用 process high-water RSS；目前是 process級最後觀測，不提供每個階段的 actual RSS time series。

## Active rows 與 objective 等價性

Halo FullABC仍保留完整 `OriginalInput` row IDs/observations作為 replay domain。Support-connected components以 sparse row membership建立 local component views；`FixedNeighborPreparation`、`LocalConditioning`使用實際 affected rows與block最大工作集，不掃描 N×P dense design。Halo diagnostics另記錄 `free_support_entries` 和 `zero_sensitivity_rows`。端點 replay在完整 domain計算 prediction、objective及相同原始 normalization scale，Recovery assembly也獨立比較完整 response。

基於以上已存在的 sparse component path，本輪沒有再複製一份 active-row observation/row-ID global snapshot。這會多出一份長期 resident buffer及 parent→active row mapping，且需重做 fixed residual constant、OriginalInput identity、no-Halo/all-fixed及global replay數學證明；目前沒有測量證據支持這項重構能降低 peak RSS。此決策是保留既有已測的component-local row reduction，不是宣稱全域 active-row API已完成。

## 大型 workload 證據

| Scale | 工作量 | 結果 | 限制 |
|---|---:|---|---|
| A：小型 synthetic / P2 真實局部 geometry | 85–122 Halo、32–46 free Halo | P2 suite process peak RSS 267,124,736 bytes；5 cases均完成 | 不可外推到大案例 |
| B：6Z6U Chain A | P0 census：1,418 Targets、659 supported Halos、229,783 OriginalInput rows，35,329 candidate Halos為n=0 | 完整 paired run運行約2小時後中止，exit 130；沒有summary或results JSON | elapsed等待時間>2小時；completed workflow wall time與process peak RSS未量測 |
| C：更大的含 Halo Partial Selection | 尚未執行 | Not measured | Chain A完整run未完成前不啟動更大selection；不以 Scale D替代 |
| D：5,938,087 observations / 37,406 contributors | 無 Halo；37,406皆Target | 先前FixedNeighbor preparation/zero-sweep endpoint preflight process peak RSS 1,471,479,808 bytes | 只驗證 no-Halo prepare/preflight，不是 Halo FullABC fitting或Recovery |

既有 Scale D 記錄見 [zero-Halo giant-component report](joint-bounded-giant-component-resource-validation-20261020.md)。它不能作為本 Halo estimator 的大型驗證。Small synthetic P3的約 50 MB RSS也不能支持 4 GB gate。

## P3 結論

- 確認目前沒有 global dense observation×Halo design；solver以 connected component sparse rows和局部最大workspace進行。
- Existing planner包含 checked size arithmetic、per-stage estimated peak及 hard-limit planning budget；本輪沒有取得完整含 Halo chain A process的actual RSS。
- Scale B完整 workflow未完成，Scale C未執行；因此 `<=4 GB`、大型 Halo search/runtime及Recovery peak workspace未通過驗收。
- 在完整Scale B未完成前不進行高風險 active-row representation重構；現有 component-local reduction保留，完整OriginalInput normalization與replay契約不變。

## 執行追蹤

- Build：`cmake --build build/debug --target joint_halo_fullabc_validation --parallel 4`；成功。
- Local P2 suite：單執行緒 process peak RSS `267,124,736` bytes。
- 6Z6U Chain A：單執行緒 runner運行約2小時未寫出中間/最終結果；以Ctrl-C中止，exit 130。這是long-run operational evidence，不是solver failure；completed peak RSS未測得。空輸出目錄為 `docs/research/results/halo_fullabc_p2/6z6u_chain_a_paired/`。
- 本報告中的No-Halo Scale D數據沿用既有research report，未冒充本輪實測。
