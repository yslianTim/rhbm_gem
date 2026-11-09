# Joint Component Estimator 測試診斷精簡驗收

## Git 狀態

- Branch 為既有的 `develop`。初始 HEAD、`origin/develop` 與唯讀查核的遠端 ref 均為 `58c5e020f9e37a3247be346e8d0ffa5c1672306e`。
- 初始工作區乾淨。修改與提交皆在原本的 checkout，原有的其他 worktree 未動用。作業期間未建立 branch 或 worktree。
- 最後一筆程式碼 HEAD 為 `478735488fda9ab59a7e8abc4d87d695a2ec5bde`。交付 HEAD 為加入本報告的 commit，可用 `git log -1 --format=%H -- docs/developer/joint-test-instrumentation-simplification.md` 查詢。
- 四筆程式碼提交已用一般 `git push origin develop` 快轉同步，並核對 `origin/develop` 與程式碼 HEAD 一致。提交報告前工作區乾淨。

| 階段 | Commit | 內容 |
| --- | --- | --- |
| P1 | `1c8016aa0d80d4cfc003fdca04c31087e88a24e1` | 移除 bounded LM trace |
| P2 | `f29a11cb9956ec388f8ea6b04b2ff5ac95e858e0` | 精簡 Profile LM trial |
| P3 | `921ebcd4ed475e8fb260275dc65cb58645d6630e` | 將候選觀察縮限到增量 replay 審核 |
| P4 | `478735488fda9ab59a7e8abc4d87d695a2ec5bde` | 移除 component observer 並精簡診斷收集 |
| P5 | 加入本報告的 commit | 完整驗收紀錄 |

## P0–P5 結果

| 階段 | 狀態 | 修改範圍與驗收 |
| --- | --- | --- |
| P0 | Completed | 修改前盤點巨集、欄位、寫入者與讀取者。五組相關 CTest 通過，保存 Profile 軌跡、BCD 認證、回退、復原及 CLI 持久化基準。以 checkout 的程式碼為準，唯讀參考 Notion 與 repository 內三份 Joint 開發文件。 |
| P1 | Completed | `ActiveSetBoundedLM.{hpp,cpp}`、`ActiveSetBoundedLM_test.cpp`。刪除 thread-local trace、逐次序列及 trace 參數。用正式步長、gradient、multiplier、KKT 與工作計數檢查放開 working bound、rank recovery、邊界與跨呼叫獨立性。15 項替代斷言所在的測試先在原實作通過，刪除後 LM 與 Profile 共 34 項通過。 |
| P2 | Completed | `Numerics.hpp`、`ProfileLmSearch.cpp`、`GaussSeidelBcd.{hpp,cpp}`、`Search_test.cpp`、`JointRuntimeJson.hpp`。移除子問題複本、trial 詳細欄位及 test-only 狀態序列化。保留接受軌跡、trial replay、Width 限制、stop reason、工作計數、正式 failure log 與 driver 狀態機。19 項 Profile 測試先在原診斷實作通過，刪除後五組相關 CTest 通過。 |
| P3 | Partial | `GaussSeidelBcd.{hpp,cpp}`、`GaussSeidelBcd_test.cpp`。移除路徑 enum 和未使用事件欄位，回退、fast path、故障和復原測試改讀現有結果。五組相關 CTest 通過。仍保留增量候選的最小逐筆觀察。 |
| P4 | Completed | `BcdPolicy.hpp`、`ComponentOrchestration.cpp`、`GaussSeidelBcd.{hpp,cpp}`、`ComponentSweepBudget_test.cpp`、`GaussSeidelBcd_test.cpp`、`PureHaloSupportExclusion_test.cpp`、`Runtime_test.cpp`、`JointBcdDiagnostic.cpp`。移除整份搜尋結果 observer、solver 內 sweep 歷史及 telemetry toggle。保留策略覆寫、block 記錄、conditioning、正式進度及 `final_sweep`。五組相關 CTest 通過。 |
| P5 | Completed | 本報告。掃描殘留符號，驗證 Release 正式與測試目標、完整 CTest、確定性數值和正式 JSON／SQLite。 |

### P0 依賴矩陣

| 原項目 | 定義與寫入 | 讀取與測試 | 行為與處置 |
| --- | --- | --- | --- |
| `BoundedLmDampingTrace`、`BoundedLmTestWork`、`BoundedLmWorkForTesting`、`capture_trace` | `ActiveSetBoundedLM.hpp/.cpp`，在 working-set loop、QR 和 damping search 擷取 | `ActiveSetBoundedLM_test.cpp` | 測試專用純觀察。改由 `ActiveSetStep`、解析解、KKT 與既有計數驗證。 |
| `last_bounded_subproblem`、`LmTrial` 詳細欄位 | `ProfileLmSearch.cpp` 記錄求解子問題並複製到 `Numerics.hpp` 的 trial | `Search_test.cpp`。子問題狀態另流向 BCD failure formatter，正式建構一向輸出 `unavailable` | 刪除測試複本與詳細欄位。保留 logger 使用的步長及 reduction 資料，原正式 log 欄位與值照常輸出。 |
| `SearchResult::lm_status` | `Numerics.hpp`、`SearchProfile` | `Search_test.cpp`、`JointRuntimeJson.hpp` | 測試專用副本。改驗證 `stop_reason`，driver 自身 status 仍有獨立測試。 |
| `BcdCandidateEvent`、路徑 enum 與 `candidate_observer` | `GaussSeidelBcd.hpp/.cpp` 候選提交前組裝 | `GaussSeidelBcd_test.cpp` | 測試專用純觀察。多數測試改用結果與獨立 replay，逐筆增量審核仍依賴原 callback 的資料。 |
| `observe_search_for_testing` | `BcdPolicy.hpp` 宣告，`SolveBcdComponent` 呼叫 | Sweep Budget、Runtime、Pure Halo 測試 | 測試專用純觀察。改用直接回傳結果、component evidence、正式進度與持久化 sweep 計數。 |
| `collect_diagnostics`、`blocks`、`sweeps`、`final_sweep` | `SearchBcdComponent` 寫入記錄。開關另讓有 certification 的失敗終點執行 Local Conditioning | BCD 診斷輔助程式與 block failure 測試讀歷史。正式進度與 component evidence 使用最後一個 sweep | 歷史只在測試收集。成功 sweep 歷史移到測試輔助程式，失敗 sweep 從 `final_sweep` 驗證。保留 `blocks` 和 conditioning gate。 |
| `collect_telemetry` | `BcdSearchPolicy` 的測試開關控制 `BcdWork` 既有計數更新 | BCD 測試使用 replay、fast path 和 scratch 計數 | 正式建構一向啟用。取消測試開關，沿用原計數語意。 |

P3 的逐筆審核比較增量 prediction／residual、未受影響的列、objective、local／global delta 與獨立接受判定。搜尋終態只反映後續更新後的 cache，成功 sweep 又會用完整 replay 重設 objective。這些終態無法取得每個提交前候選的原始列與 objective。現有 replay 實作是檔案內部函式，另外公開測試入口仍無法觀察實際拒絕的候選。保留範圍限於受影響列重算，沿用原 callback，無新增測試 solver。

## 變更量

| 範圍 | 新增行 | 刪除行 |
| --- | ---: | ---: |
| `src/` 正式來源及其測試條件分支 | 30 | 207 |
| `tests/` 測試及輔助程式 | 316 | 401 |

- 移除 2 個 struct、1 個 enum、1 個 component-result callback。候選 callback 保留，已移除完整 replay 與 unchanged 路徑的事件組裝。
- `RHBM_GEM_TEST_INSTRUMENTATION` 的 `#ifdef` 由 50 處減至 30 處。移除欄位共 47 個，含 LM trace 型別內 21 個、Profile LM 與其診斷橋接 15 個、候選事件 8 個、BCD 診斷與策略 3 個。
- GTest 測試案例新增 0、移除 0。既有故障測試擴展為七種故障各跑受影響列與完整 replay 路徑。

## 測試與數值驗收

- Release 正式 CLI、正式函式庫及測試版建構通過。檢查編譯設定確認 114 個正式函式庫與 1 個 CLI 編譯命令均未定義巨集，207 個含巨集的命令全部屬於測試目標。
- 受影響的五組 CTest 通過。完整 CTest 21／21 通過，含 GTest 註冊稽核、Pure Halo Support Exclusion、下游 Local／Group fitting、持久化及 CLI smoke。
- 完整 CTest 首次使用單執行緒 OpenMP 限制，Map Simulation 的 2／4 執行緒測試及 HRL 的 Eigen 執行緒復原測試因此失敗。重跑使用 4 執行緒上限、CTest 序列執行。既有空間查詢並行測試把 worker 數寫死為 16，完整驗收以暫存編譯副本調為 4，測試邏輯與斷言照常執行，Git 來源與正式建構產物保持原樣。
- P3 兩項替代斷言在原實作暴露了比較前提錯誤。獨立 replay 的 objective 採不同加總順序，因此沿用現有 enclosure。非有限 objective 仍可能進行 block solve，因此改用候選 replay 計數與未提交終態驗證。重跑 BCD 38 項全數通過，未調整容差。
- 使用同一組固定 `chain(3)` 資料、Release 設定與單執行緒，前後以十六進位浮點格式記錄數值。Profile 的每次接受更新、eta、A/C、objective 與 replay 結果相同。BCD 的最終 A/C/τ、objective、兩種 projected KKT、stop reason、integrity、endpoint verification、Independent Replay、復原及 work 計數逐位元相同。前後快照 SHA-256 皆為 `e859d215f94e606b4b3dea29d1e57b5ae659f485648a7c4252b62053669eaccd`。

| BCD 固定情境 | Stop reason | Sweeps | 復原嘗試／成功 | 增量／完整候選 replay | 終點認證／Independent Replay |
| --- | --- | ---: | --- | --- | --- |
| 單 block | `block-stationary` | 3 | 0／0 | 3／0 | 通過／通過 |
| 雙 block 增量 | `block-stationary` | 5 | 0／0 | 10／0 | 通過／通過 |
| 雙 block 記憶體回退 | `block-stationary` | 5 | 0／0 | 0／10 | 通過／通過 |
| 注入一次 Width LM 失敗 | `block-stationary` | 5 | 1／1 | 9／0 | 通過／通過 |

- BCD 回退測試亦比對兩條路徑的 eta、beta、prediction、residual 及最終 objective，皆為位元相同。每個 block 的接受判斷與理由一致，objective 使用原有 enclosure 比對。
- CLI 同一固定輸入的正式 JSON 匯出、SQLite schema v21、四筆 Joint payload 及四筆 stage payload 通過比較。只正規化耗時、暫存路徑、變動的建置雜湊和每次執行生成的 run UUID。初次比較定位到 stage payload 的隨機 UUID。依 UUID 格式與同一 run 在各階段的引用對應正規化後，四筆 stage payload 均一致。
- 正式建構原本已啟用計數，這次直接保留原更新位置。保留的候選審核仍只編入測試目標。沒有新增正式配置或數值運算。未進行可比較的 wall time／peak RSS benchmark。

## 保留的測試能力

| 巨集位置 | 區塊數 | 保留理由 |
| --- | ---: | --- |
| `PotentialFittingWorkflow.hpp/.cpp` | 3 | 宣告、存放並執行 workflow JointFitResult fault hook。 |
| `JointComponentEstimator.hpp/.cpp` | 4 | Problem Construction 記憶體預算的宣告、儲存、讀取與設定。 |
| `Numerics.hpp`、`ProfileLmSearch.cpp` | 2 | Width LM failure injection 與失敗工作計數。 |
| `BcdPolicy.hpp` | 3 | `configure_for_testing` 的相依 include、forward declaration 和策略注入。 |
| `GaussSeidelBcd.hpp` | 3 | 候選故障與增量審核的型別、測試策略和故障選項、逐 block Width／untrusted 注入。 |
| `GaussSeidelBcd.cpp` | 15 | 搜尋／conditioning／認證設定 2 處、Global State 預算 1 處、Width／untrusted 控制 2 處、候選故障及避免跳過注入 3 處、增量審核 5 處、sweep replay corruption 1 處、component 策略覆寫 1 處。 |
| `tests/CMakeLists.txt` | 建構設定與隔離檢查 | 只為測試目標定義巨集，並檢查正式目標的 compile definitions。 |

## 結論

純觀察與詳細診斷來源已縮減。現有測試維持 Active-set 與 Profile 的數值契約、BCD 接受與 rollback、故障復原、認證與下游資料驗證。正式演算法、容差、Width 校準、資源預算與資料契約依原實作執行。候選逐筆審核及必要的 block 記錄留在現有介面，只有在能等效驗證提交前結果與失敗原因時才適合另案移除。
