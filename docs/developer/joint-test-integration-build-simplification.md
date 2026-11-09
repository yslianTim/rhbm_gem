# Joint Component Estimator 測試介面與建構隔離驗收

## A. Git 狀態

- 分支為既有的 `develop`。初始 HEAD、`origin/develop` 與 GitHub 遠端 ref 均為 `b8b0e88249ca0b19f3ebcd1c9cb8fab082ca55f0`。附件的參考提交為 `7df11391`。
- 初始工作區乾淨。全部變更在原 checkout 完成，現有的另一個 worktree 僅列出位置。未建立分支或 worktree。
- 最後的程式碼 HEAD 為 `e1c287c865ada9d146e0f28b321e0ad321cf45d8`。交付 HEAD 是加入本報告的提交，可用 `git log -1 --format=%H -- docs/developer/joint-test-integration-build-simplification.md` 查詢。寫入本報告前，追蹤工作區乾淨。
- P1–P4 已用一般 `git push origin develop` 快轉同步，並以 `git ls-remote origin refs/heads/develop` 核對到程式碼 HEAD。報告提交的最終同步結果另以相同 ref 查核。

| 階段 | Commit |
| --- | --- |
| P0 | 基準與依賴盤點合併記錄於本報告 |
| P1 | `9eac5fe7` — `refactor(test): simplify joint BCD diagnostic facade` |
| P2 | `1ffcbfcb` — `refactor(joint): narrow component test policy integration` |
| P3 | `b383d5da` — `refactor(test): simplify joint workflow integration hooks` |
| P4 | `e1c287c8` — `refactor(build): clarify joint test helper ownership and isolation` |
| P5 | 本報告 — `test(joint): validate integration interface cleanup` |

## B. 階段驗收與依賴盤點

| 階段 | 狀態 | 修改範圍與驗收 |
| --- | --- | --- |
| P0 | Completed | 查核 `b8b0e882` 的依賴、編譯命令、GTest／CTest 註冊、六種搜尋情境、Profile LM、正式組裝與 CLI 持久化。完整 CTest 21／21 通過。 |
| P1 | Completed | `JointBcdDiagnostic.{hpp,cpp}`、`GaussSeidelBcd_test.cpp`、`ComponentSweepBudget_test.cpp`。刪除自製 fit 與搜尋結果副本。既有正式 fit 呼叫直接負責組裝驗證。三組 Joint CTest 與兩組 GTest 註冊稽核通過。 |
| P2 | Partial | `BcdPolicy.hpp`、`GaussSeidelBcd.cpp`、`ComponentSweepBudget_test.cpp`。六個使用點均需執行真實 component 流程，保留借用策略並限制七個欄位。三組 Joint CTest 與持久化測試通過。 |
| P3 | Partial | `CommandApiPipeline_test.cpp`。API outcome、程序內 CLI exit mapping、兩份 SQLite reload 及發布 gate 各有明確斷言。兩種故障仍用原 snapshot hook。命令、Joint workflow、持久化及正式 CLI smoke 共四組通過。 |
| P4 | Completed | `tests/CMakeLists.txt`、`tests/cmake/check_build_isolation.py`。純測試 helper 由 `rhbm_tests` 編譯。新增獨立 compile-command 檢查，保留配置期保護。 |
| P5 | Completed | 本報告。Release 正式與測試建構、完整 CTest、數值逐位元比對、正式 CLI JSON／SQLite 對照與殘留介面掃描完成。 |

正式 Component Orchestration 建立各元件 evidence 和 qualified state。Component Assembly 負責映射、available rows、constant terms、parent replay 與發布。Workflow 使用捕捉後的 snapshot 決定後處理。盤點以實際程式碼、兩輪 Git 歷史驗收文件及以下唯讀資料為依據。[前輪 instrumentation 驗收](https://github.com/yslianTim/rhbm_gem/blob/daecf691/docs/developer/joint-test-instrumentation-simplification.md)、[前輪 fault injection 驗收](https://github.com/yslianTim/rhbm_gem/blob/7df11391/docs/developer/joint-test-fault-injection-simplification.md)、[Joint Component Estimator 流程](https://app.notion.com/p/3e27cca52580813f9b21fe12a2c599ef)、[Runtime](joint-component-runtime.md)、[Numerics](joint-component-numerics.md)、[Stage Integration](joint-stage-integration.md)。

| 分類與介面 | 定義／使用 | 邊界與等效性 |
| --- | --- | --- |
| A `AssembleResult` | `tests/support/JointBcdDiagnostic.cpp` 內部。BCD Search、Sweep Budget 透過 `BcdDiagnosticResult::fit` 間接讀取。 | 測試檔仿製 component state、qualification evidence、identity mapping、mask、prediction、objective 和 consistency。搜尋改讀原始證據。發布由既有 `FitBcdComponents`（同 fixture 已呼叫）及 Runtime、Sweep Budget 的正式入口驗證。 |
| A `BcdDiagnosticResult` | `tests/support/JointBcdDiagnostic.hpp/.cpp`。上述兩組測試。 | 多數欄位、assessment、trust、attempted 與 verified 是原始搜尋副本。以 `BcdSearchResult::certification` 取代。無正式資料契約變更。 |
| B `RunBcdDiagnostic` | 同上。 | 保留測資取用、block records、按需接收 sweep progress。需 instrumentation。移除人造 validation 與組裝，正式搜尋本身驗證 domain、layout、bound 和 accepted state。 |
| B／C `search_policy_for_testing` | `BcdPolicy.hpp` → `GaussSeidelBcd.cpp`。Component Sweep Budget、Runtime、Pure Halo。 | 六處使用均經正式元件入口驗證組裝、拒絕發布、progress／capture 或下游。搜尋預算與故障需要可控，保留同步借用。可改成只讀這些控制，不轉接 observer 或診斷。 |
| C `JointPostprocessingSnapshotHookForTesting` | `PotentialFittingWorkflow.hpp/.cpp`。`CommandApiPipeline_test.cpp` 的兩種故障，各經 API 與 in-process CLI 一次。 | direct postprocessing 或 repository 無法覆蓋實際 command outcome、exit 3／4、命令寫入 DB。獨立 CLI smoke 覆蓋真實 CLI 成功路徑。保留必要的整合 seam。 |
| C `RHBM_GEM_TEST_INSTRUMENTATION` | `tests/CMakeLists.txt` 啟用。型別定義 `EvaluationContext`、`BcdExecutionPolicy`、`BcdSearchPolicy` 和正規程式流程均受影響。 | 正式與 instrumented library 須分開編譯並在 binary 間分離。測試 helper 原放入 instrumented library，可移到測試 executable。 |

## C. Facade 與整合入口

- `SearchBcdFixture` 讀取既有 `JointProblem` 的 domain、layout、observations 與 context，開啟 block 記錄，按需用既有 progress callback 收集成功發出的 sweep records，直接回傳 `BcdSearchResult`。根 Fitting Domain 的 `BuildParameterLayout` 使用連續 identity indices，因此既有 fixture 可直接借用 observations。正式搜尋負責尺寸、有限性、Width bound 和 layout 驗證。
- `BcdEndpointCertification` 直接提供 assessment、trust、verification 與 qualification。舊 `AssembleResult` 的 component state/evidence、mask、identity mapping、prediction、objective 及假造的 consistency 均已移除。`SmallSingleBlockFixturesUseExistingEndpointCertification` 仍執行原有的搜尋與正式 fit 各一次，直接比較正式 fit 的 A/C、log τ、prediction、objective 和 stationarity／certification sweep。其他純搜尋測試只執行原本的搜尋。
- 元件策略指標借用測試呼叫端的區域變數，僅於同步 `FitBcdComponents`／`SolveBcdComponent` 讀取。白名單為 core size、nominal／maximum sweeps、Global State memory override、Width failure、candidate fault、sweep replay fault。正式 criteria、Width bound、component identity 和 progress 仍從 production policy 或正式 observer 建立。搜尋的 candidate observer 與 block 記錄開關只由直接搜尋測試設定。
- `JointPostprocessingSnapshotHookForTesting` 仍在 capture 後、postprocessing 前以 thread-local 指標提供兩種固定故障。可直接驗證 postprocessing 的測試已使用 `RunJointPostprocessing`。單獨呼叫它或 `DataRepository` 可驗證 gate／codec，但無法確認命令對已持久化失敗結果的 outcome、CLI mapping 與實際資料庫寫入。現有正常測資會收斂，小選區會先走 calibration failure。保留 hook 以固定觸發「初始化有效但 Joint 不合格」和「Joint 合格但沒有 Halo donor」兩個命令情境。

## D. 建構隔離

| Release 目標 | P0 編譯來源 | P5 編譯來源 | P5 巨集 | P5 純測試 helper |
| --- | ---: | ---: | --- | --- |
| `rhbm_gem` | 114 | 114 | 無 | 0 |
| `rhbm_gem_cli` | 1 | 1 | 無 | 0 |
| `rhbm_gem_test` | 115 | 114 | 有 | 0 |
| `rhbm_tests` | 92 | 93 | 有 | 10 |

- `rhbm_gem_test` 從 `rhbm_gem` 取得正式 source list，沿用 C++20、compile properties、include paths 與 PIC。兩版連結相同的 Eigen／Boost／OpenMP 依賴，Release 的 UMAP 功能與連結仍啟用。UMAP 巨集只由正式屬性複製一次。
- 正式 CLI 的連結輸入為正式函式庫，測試 executable 的連結輸入為 instrumented 靜態函式庫。兩套編譯保留 `EvaluationContext`、`BcdSearchPolicy`、`BcdExecutionPolicy` 的各自定義，避免同一個 executable 混入兩種編譯契約。正式 install 匯出只有正式目標。
- `test_build_isolation` 使用 `compile_commands.json` 檢查巨集、目標歸屬及兩版 library 的來源集合。Ninja／Makefiles 開啟 compile commands 時註冊此檢查。其他 generator 仍有 CMake 配置期 guard。此檢查在 P0 資料中抓到原本的 helper 歸屬，在 P5 通過。額外使用修改過的命令副本驗證可拒絕 CLI 巨集洩漏與 test library 巨集缺失。114 組相對應的 Release 編譯旗標集合僅有 shared-library export define 與測試巨集差異。另完成 UMAP 關閉的 Debug 及純正式 Release 配置查核。完整建構及執行驗收在 macOS arm64 Release 執行。
- GTest 保持 1,097 案。CTest 21 → 22 組，新增的檢查不啟動 GTest。兩個註冊 coverage 測試通過，原有測試分組與 labels 保留。

## E. 覆蓋等效矩陣

| 測試或分支 | 變更後驗證 | 判定 |
| --- | --- | --- |
| `SmallSingleBlockFixturesUseExistingEndpointCertification` | 四個 chain／cube 測資仍比對搜尋與正式 fit。真正的 component 和 assembly 結果提供 publication、convergence、映射後參數、prediction、objective 與 sweeps。 | Equivalent |
| `ReplayWorkspacePreservesExactSearchResultsAcrossRuns` | 直接比較兩次搜尋的 accepted state、block records、每次 sweep 的判斷與 objective 位元。 | Equivalent |
| `CandidateObserverDoesNotChangeSearch` | 直接比較相同搜尋契約、認證有無及 verification，以及 accepted state 的位元。 | Equivalent |
| `IncrementalObjectiveMatchesIndependentFullReplay` | 保留逐筆獨立 replay、接受決策與 sweep reconciliation。由正式 certification 取得結果。 | Equivalent |
| `IncrementalObjectiveQualifiesAtActiveWidthUpperBound` | 原 projected KKT、Width face、replay 與 reconciliation 照常檢查。 | Equivalent |
| `EverySweepReconcilesObjectiveAfterReplayAudit` | 用正式 certification 與 accepted-state integrity 代替模擬 RuntimeConvergence。保留逐 sweep 與獨立 objective。 | Equivalent |
| `SweepReplayMismatchFailsBeforeObjectiveReconciliation` | 原故障與 failed final sweep、無成功 sweep callback、objective 未刷新均保留。 | Retained Integration Hook |
| `IncrementalObjectiveRemainsQualifiedAcrossInitialWidths` | 三個初值照常比較。以正式 certification、integrity 與 replay 驗證搜尋。 | Equivalent |
| `SignedCancellationQualifiesWhenOnlyKktGateIsSatisfied` | 仍驗證 signed cancellation、實際候選判定、兩種 KKT、integrity 與合格終點。 | Equivalent |
| `HardCaseIncrementalObjectiveAndReplayQualification` | 保留所有測資及獨立 replay。失敗分支直接檢查無 certification。low-signal 分支檢查正式 qualification 和 integrity。 | Equivalent |
| `UnchangedCandidateSkipsOnlyCandidateReplay` | 仍以實際發出的 sweeps 和 block records 檢查 replay fast path。 | Equivalent |
| `ProductionBoundedSearchReleasesTheStartingUpperFace` | 保留 bounded LM release、全域 replay 和投影梯度對照。 | Equivalent |
| `MultiBlockBudgetFailureRemainsDistinctFromRuntimeConvergence` → `MultiBlockZeroSweepBudgetCannotCertify` | 兩種 topology 保留零預算停止與零認證。移除對空白測試組裝的 RuntimeConvergence 斷言。正式不可發布的契約仍由 Sweep Budget 與 Runtime 驗證。 | Equivalent（搜尋契約） |
| `CheapStationarityNeedsCompleteSweepConfirmation` | 仍比較一次與完成確認的 sweeps，保留實際 callback 計數。 | Equivalent |
| `BcdAcceptsAtMostOneTrustedLocalUpdatePerBlockVisit` | 使用同一正式 block records 驗證每次 visit 的接受上限。 | Equivalent |
| `EasyQualifiedSearchStopsWithoutExtraWorkOrParameterChanges` | 比較 fixed／adaptive 搜尋的參數、prediction、objective、sweeps，並直接驗證認證與 integrity。 | Equivalent |
| `ExtensionPreservesTheOriginalPathAndQualifiesProductionAssembly` | 保留延長後的正式 assembly 及證據。新增正式 Width bound、criteria、identity 和 callback 作用域斷言。 | Retained Integration Hook |
| `UnqualifiedHardLimitAndMemoryRejectionCannotPublishComponentState` | 保留硬限額及記憶體拒絕到正式 assembly。 | Retained Integration Hook |
| `MultipleComponentsExtendAndStopIndependently` | 保留各 component 的不同 sweep 停止與共同組裝。 | Retained Integration Hook |
| `TransientWidthFailureRecoversAndAssemblesQualifiedComponent` | 原 Width 注入、recovery counters、進度及合格 assembly。 | Retained Integration Hook |
| `WidthRecoveryExhaustionAndUnqualifiedTerminalStatesRejectAssembly` | 原 persistent failure、early budget、candidate fault 與實際發布拒絕。 | Retained Integration Hook |
| `CanonicalReasonReachesComponentProgressAndCapturedResult` | 保留七種 stop reason 經正式 Solve、Fit、Progress 與 Capture 的驗證。 | Retained Integration Hook |
| `ExtendedSearchQualifiesAssemblyAndFeedsEveryDownstreamStage` | 保留延長 sweep 後的正式組裝與直接執行 postprocessing，含 Halo、Peeling、Local 和 Group。 | Retained Integration Hook |
| `JointFailureAndUnavailableOutcomesPersistThroughApiCliAndReload` | 保留兩種注入及 outcome、exit 3／4。兩條命令路徑現在都 reload 並檢查完整 gate 與 failure／donor evidence。 | Retained Integration Hook |
| `joint_workflow_cli_smoke_test`、`ProductionBoundedJointDefaultPersistsReloadsAndFeedsPostprocessing` | 前者仍啟動正式 CLI 程序。後者仍驗證命令中的合格結果與正式獨立求解的一致性。 | Equivalent |

## F. 數值與持久化

基準和最終探測都使用 Release、Apple Clang 21、arm64 和單執行緒。固定 `chain(3)` 使用初始 τ = 0.55、上界 τ = 2.2、practical-v2。單 block 使用預設 core，雙 block 使用 core 2 與 5 sweep 上限。記錄搜尋 state、objective、兩種 KKT、integrity、stop reason、sweeps、recovery、certification、replay 與工作計數。另記錄同資料的 Profile LM 接受軌跡，以及正式組裝的 A/C/τ、row mask、prediction、evidence 與 objective。十六進位浮點快照以 `cmp` 通過，兩份 SHA-256 均為 `01e6c5c02a94470dc18fa323435365dc1e747981fcea7a1a83e0783bde3366f9`。

| 項目 | 單 block | 雙 block 增量／Full Replay |
| --- | --- | --- |
| A₁ | `0x1.0000cb9e2047fp+1` | `0x1.0000cba58766dp+1` |
| C₁ | `0x1.99988cf009dbfp-3` | `0x1.99988ce62d1d7p-3` |
| τ₁ | `0x1.00003822cec87p-1` | `0x1.00003824d5411p-1` |
| A₂ | `0x1.fffbb76f1c1cbp+0` | `0x1.fffbb6fae35b3p+0` |
| C₂ | `0x1.99a15f046ca59p-3` | `0x1.99a15f51a128ep-3` |
| τ₂ | `0x1.ffff2203093a3p-2` | `0x1.ffff21e3270cep-2` |
| A₃ | `0x1.ffff5e2bf6584p+0` | `0x1.ffff5e3ad3569p+0` |
| C₃ | `0x1.999c5ffaa65bbp-3` | `0x1.999c5ff0c61d7p-3` |
| τ₃ | `0x1.fffff9a43eb45p-2` | `0x1.fffff9a852442p-2` |
| Objective | `0x1.bc2eb370d637fp-33` | `0x1.bc2eb3fd335ep-33` |
| Global A/C KKT | `0x1.cp-52` | `0x1.22c1bfp-29` |
| Projected Width KKT | `0x1.cdp-42` | `0x1.5c87cp-34` |

| 情境 | Stop reason | Sweeps | Recovery 嘗試／成功 | Integrity | 增量／Full 候選 replay | Certification／Independent Replay |
| --- | --- | ---: | --- | --- | --- | --- |
| 單 block | `block-stationary` | 3 | 0／0 | 通過 | 3／0 | 通過／通過 |
| 雙 block | `block-stationary` | 5 | 0／0 | 通過 | 10／0 | 通過／通過 |
| 雙 block 記憶體回退 | `block-stationary` | 5 | 0／0 | 通過 | 0／10 | 通過／通過 |
| 暫時 Width 失敗 | `block-stationary` | 5 | 1／1 | 通過 | 9／0 | 通過／通過 |
| 持續 Width 失敗 | `block-recovery-exhausted` | 1 | 1／0 | 失敗 | 1／0 | 未執行 |
| Global State 超預算 | `resource-limit:global-state-evaluation:estimated-peak-exceeds-budget` | 0 | 0／0 | 未建立 | 0／0 | 未執行 |

- 正式組裝的 mapping 與 consistency 通過，objective 為 `0x1.bc2eb370d637ep-33`，available-row mask 與逐列 prediction 位元相同。帶正式 Target selection 的資格與下游另由 Adaptive Sweep、Pure Halo 及 CLI 測試驗證。
- 真實 CLI 用固定 25 個原子的模型和同樣的 map simulation、joint arguments、SQLite 寫入與 JSON 匯出流程。比對 Joint JSON v24、SQLite v21 schema、四筆 Joint payload 及四筆 Stage v5 payload。只正規化耗時、暫存路徑、建置／來源雜湊及 run UUID 與其引用。Stage payload 內的 Peeling、Local 與 Group 證據均納入比對。
- 完整測試包含 Active-set QR、Bounded LM、Profile LM、BCD、candidate acceptance、incremental／full replay、sweep extension、recovery、endpoint certification、independent replay、assembly、Pure Halo、Post-Joint Local、Group、Command API、CLI mapping、JSON／SQLite 與正式 CLI smoke。建構使用 `--parallel 4`。CTest 序列執行，OpenMP 限制 4、巢狀層級 1，沿用各組的單／雙執行緒設定。基準首次強制 OpenMP 1 使 HRL thread restore 測試失敗，改回 4 上限並保留其測試契約後 21／21 通過。

## G. 簡化量與成本

- `JointBcdDiagnostic.{hpp,cpp}` 由 202 行減為 30 行。移除六個私有 helper（`Select`、`IndicesOf`、`Values`、`SelectBeta`、`State`、`AssembleResult`）和一個平行結果型別。刪除的 `AssembleResult` 有 89 行。搜尋入口由一個 facade 換成一個測資 helper，重複 Assembly／Qualification 實作由一份減至零。
- 不含本報告，程式與建構設定新增 281 行、刪除 363 行，淨減 82 行。保留元件策略與 workflow hook，兩套 library target 繼續分離。
- 最終完整 CTest 22／22 通過，耗時 41.48 秒，基準 21／21、41.42 秒。Joint numerical 分組 14.57 秒，基準 14.64 秒。計時為個別執行的 wall time。P0 使用已建構的 Release。P4 重新配置後重建耗時 4 分 36.53 秒。這兩次建構的工作量不同，無法比較 clean-build 速度。
- 測試 executable 為 15,439,008 → 15,424,688 bytes。測試靜態函式庫為 12,161,120 → 12,094,952 bytes。正式 CLI 和正式函式庫大小分別維持 33,496 與 6,321,376 bytes。

## H. 保留的測試介面

| 介面 | 範圍與用途 |
| --- | --- |
| `BcdSearchPolicy::candidate_observer` | 在候選提交前觀察受影響列、預測、殘差、objective 和接受判定。終態無法還原每筆被拒絕的候選。 |
| `force_width_failure_for_testing`／`force_width_subproblem_failure_for_testing` | 逐 block visit 注入 Width 子問題失敗，覆蓋暫時復原、retry exhaustion 和 rollback。 |
| `candidate_fault_for_testing` | 仍覆蓋有限 objective mismatch、NaN objective／parameter、Inf prediction／residual 的拒絕與狀態完整性。 |
| `global_state_memory_budget_for_testing` | 低記憶體時完整 replay 回退與 resource preflight 拒絕。 |
| `perturb_sweep_replay_objective_for_testing` | Sweep replay mismatch、objective reconciliation 與不發布驗證。 |
| `search_policy_for_testing` | 同步借用搜尋預算和四個故障控制，支援正式 component assembly、capture、progress 與下游整合。 |
| `JointPostprocessingSnapshotHookForTesting` | 僅用於命令的未收斂及 donor 缺失，保留 exit 3／4 和 SQLite 整合路徑。 |
| `collect_diagnostics`、`SearchBcdFixture` | 直接搜尋的 block records、選用的 sweep 證據與既有數值測試。 |
| `RHBM_GEM_TEST_INSTRUMENTATION` | 只由測試靜態函式庫公開傳給 `rhbm_tests`。正式 build 另行編譯。 |
