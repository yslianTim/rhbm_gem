# Joint Component Estimator 故障注入與策略覆寫驗收

## A. Git 狀態

- 作業分支為既有的 `develop`。初始 HEAD、本機 `origin/develop` 與遠端 `refs/heads/develop` 均為 `1904a44418d4442dd5d6a72e8de6254d72d829f7`。初始工作區乾淨。
- 作業在原 checkout 進行。未建立分支或 worktree，既有其他 worktree 保持原狀。未使用 detached HEAD 進行作業。
- 最後一筆 JCE 程式碼提交為 `2d07408f803d71275a494bc782d4ca1b7bcf2e9e`。Final HEAD 為納入本報告及 CLI 巨集隔離檢查的 P6 提交，可用 `git log -1 --format=%H -- docs/developer/joint-test-fault-injection-simplification.md` 查詢。P6 提交前 `git diff --check` 通過。
- P0–P5 已以一般 `git push origin develop` 快轉同步，並查核遠端 ref 等於 `2d07408f803d71275a494bc782d4ca1b7bcf2e9e`。P6 報告提交的同步狀態以交付時的本機 HEAD 與遠端 `refs/heads/develop` 核對。

| 階段 | Commit | 內容 |
| --- | --- | --- |
| P0 | `34bc1c9d` | test(joint): establish fault injection cleanup baseline |
| P1 | `7dc9b545` | refactor(joint): remove BCD search policy bypasses |
| P2 | `0af981e6` | refactor(joint): remove forced rejection of trusted profiles |
| P3 | `36e53eb4` | refactor(joint): narrow candidate fault injection to replay checks |
| P4 | `7dbc6a20` | refactor(joint): simplify memory budget test overrides |
| P5 | `2d07408f` | refactor(joint): narrow test policy and workflow hooks |
| P6 | 本報告所在提交 | test(joint): validate fault injection cleanup |

## B. P0–P6

| 階段 | 結果 | 檔案、遷移與驗收 |
| --- | --- | --- |
| P0 | Completed | 唯讀檢視 Notion「Joint Component Estimator 流程」與 repository 內 runtime、numerics、stage integration。舊報告在初始 HEAD 已刪除，從 `daecf691` 唯讀取得 fixture 定義後重新建構並實測。將 `DataObjectAssemblySpatialQuery_test.cpp` 的並行查詢從 16 個工作執行緒調為 3 個，連同主執行緒共 4 個，保留啟動同步及所有斷言。完整 CTest 21／21，保存固定數值與正式持久化基準。 |
| P1 | Completed | `GaussSeidelBcd.{hpp,cpp}`、`JointBcdDiagnostic.cpp`、`GaussSeidelBcd_test.cpp`、`LocalConditioning_test.cpp`、`Runtime_test.cpp`。移除 stationarity、certification、conditioning 三個開關及診斷端補認證。測試版與正式版共用停機與認證；僅已收斂且合格的 endpoint 執行 Local Conditioning。零觀測 fixture 在正式停機前覆蓋 unchanged fast path。三組 Joint CTest 通過，快照逐位元一致。 |
| P2 | Partial | `GaussSeidelBcd.{hpp,cpp}`、`GaussSeidelBcd_test.cpp`、`Runtime_test.cpp`。移除強制丟棄可信 profile 的 callback。以現有 Profile LM 預算耗盡覆蓋沒有可信 trial 的停機與認證門檻，以原有 scoped Width 注入覆蓋有鄰居更新的復原與 rollback。保留直接 Profile LM 及逐 block 的 Width 注入。三組 Joint CTest 通過，快照逐位元一致。 |
| P3 | Partial | `GaussSeidelBcd.{hpp,cpp}`、`GaussSeidelBcd_test.cpp`。以 `chain(1)` 自然出現的局部 objective 上升取代人工覆寫，對照單 sweep 前綴的完整 accepted state。移除 Inf objective／parameter 的重複分支及注入對 unchanged fast path 的依賴。五種各自覆蓋不同驗證位置的候選 fault 與 sweep-end mismatch 注入保留。三組 Joint CTest 通過，快照逐位元一致。 |
| P4 | Partial | `JointComponentEstimator.{hpp,cpp}`、`PartialSelection_test.cpp`、`Numerics_test.cpp`。移除 Problem Construction thread-local 預算、setter 及 getter。固定 `MemoryShape` 驗證預估位元組、臨界預算、stage、rejection、overflow。另以正式預算實跑同一個 builder preflight。保留 Global State 的單一預算覆寫，使 Full Replay 回退及兩種 workspace 均超預算仍可在小 fixture 內驗證。三組 Joint CTest 通過，快照逐位元一致。 |
| P5 | Partial | `BcdPolicy.hpp`、`GaussSeidelBcd.cpp`、`PotentialFittingWorkflow.{hpp,cpp}`、`ComponentSweepBudget_test.cpp`、`Runtime_test.cpp`、`PureHaloSupportExclusion_test.cpp`、`CommandApiPipeline_test.cpp`。移除任意執行的 component configure callback，整合測試借用唯讀 `BcdSearchPolicy`，正式 criteria、Width bound、component identity 和 progress 由 solver 設定。延長 sweep 的合格 snapshot 直接交給正式 postprocessing。命令失敗測試僅在 capture 後注入 snapshot，成功命令以同資料獨立正式求解對照。Joint、persistence 與 command 的五組 CTest 通過，快照逐位元一致。 |
| P6 | Completed | `tests/CMakeLists.txt` 與本報告。增列正式 CLI 巨集隔離檢查，重新盤點所有剩餘區塊。完整 Release CTest 21／21 通過，含 Active-set Bounded LM、Profile LM、BCD、incremental／full replay、memory、recovery、Adaptive Sweep Budget、認證、Independent Replay、assembly、Pure Halo、Post-Joint Local／Group、JSON／SQLite、CLI smoke 與註冊稽核。 |

## C. P0 依賴與 Coverage Equivalence Matrix

分類 A 表示能維持相同契約。B 表示遷移後有可描述的覆蓋差異。C 表示保留驗證特定失敗或整合路徑所需的最小入口。

| 控制與定義 | 設定、消費位置及正式路徑 | 分類與處置 |
| --- | --- | --- |
| `stop_after_stationarity`、`certify_endpoint`、`assess_local_conditioning`，`BcdSearchPolicy` 的三個 bool | BCD、Local Conditioning 與 Runtime 測試設定。`SearchBcdComponent` 消費；診斷端也讀 stationarity／certification 開關。影響停機、重試認證、最終認證和 conditioning。 | A。固定正式路徑，診斷端只讀搜尋認證。刻意跳過認證的測試狀態已移除。 |
| `force_width_subproblem_failure_for_testing`，`EvaluationContext` bool | `Search_test`、`GaussSeidelBcd_test` 可直接設定；`ProfileLmSearch.cpp::SolveBoundedStep` 回傳明確故障並記入 work。 | C。現有缺少 Width 上界、非有限及不合格輸入涵蓋各自 guard；缺少穩定 fixture 同時通過初始信任並在 Width 子問題內失敗。 |
| `force_width_failure_for_testing`，`BcdSearchPolicy` 的 zero-based block/sweep callback | BCD 與 Runtime 測試設定；每次 block visit 設定該 block 的測試 context，再由 Profile LM 使用。 | C。復原、鄰居更新後重試、耗盡 rollback 與合格組裝需要指定某次 Width 求解失敗。移除每次 block visit 的測試 context 複製。 |
| `force_no_trusted_profile_for_testing`，`BcdSearchPolicy` callback | BCD 與 Runtime 測試設定；原 `SearchBcdComponent` 在選擇可信 trial 前攔截，進入復原 bookkeeping。 | B。移除。零 profile 預算產生沒有可信 trial 的正常結果；scoped Width 失敗重用共用復原分支，維持鄰居更新與發布驗證。 |
| `BcdCandidateFaultForTesting`、`candidate_fault_for_testing`，`GaussSeidelBcd.hpp` | BCD 與 Runtime 測試以 1-based sweep/block 設定；原本在局部 objective 比較前選定 fault，後續覆寫 replay 候選。 | A／B／C。局部增加已有自然 fixture，兩種重複 Inf 分支移除，五個各有獨立判定作用的 fault 保留。只對已建立的 replay 候選注入。 |
| `perturb_sweep_replay_objective_for_testing`，`BcdSearchPolicy` bool | BCD 與 Runtime 測試設定；`SearchBcdComponent` 在 sweep audit 前改變重算 objective。 | C。有限 mismatch 的 cache 拒絕、未刷新 objective 與禁止發布仍需注入。邊界判定測試直接驗證共用 enclosure，非有限初值另有測試。 |
| `problem_memory_budget_for_testing`、`SetJointProblemMemoryBudgetForTesting`，`JointComponentEstimator.{hpp,cpp}` | Partial Selection 測試設定，原 builder 原始 row buffer preflight、後續規劃及拒絕 log 消費。 | B。持續測試改為直接 planner 檢查。正式門檻的真實 builder 拒絕另做一次性驗證。 |
| `global_state_memory_budget_for_testing`，`BcdSearchPolicy` 的 optional bytes | BCD、Sweep Budget、Runtime 設定；`SearchBcdComponent` 用同一預算依序規劃 affected-row 和 full replay。 | C。小型輸入可驗證真正切換路徑及拒絕後未建立 accepted state。接近正式預算的 Global State 工作負載會先建構大型拓樸與支援資料。 |
| `BcdExecutionPolicy::configure_for_testing`，`BcdPolicy.hpp` 的 callback | Sweep Budget、Runtime、Pure Halo 設定；`SolveBcdComponent` 原本在填好正式策略後允許 callback 改寫整份策略。 | C。移除 callback，保留借用唯讀 search policy，正式門檻、寬度及元件身分在採用該配置後設定。 |
| `JointWorkflowFitHook`、`JointWorkflowFitHookForTesting()`，`PotentialFittingWorkflow.{hpp,cpp}` 的 thread-local function pointer | Pure Halo 替換整個 fit；Command API 測試修改 search qualification、donors 或擷取 fit。原本在估計之後、capture 之前執行。 | B／C。移除 fit hook。Pure Halo 直接傳入延長 sweep snapshot，正常 command 改為同輸入獨立求解對照。保留 post-capture snapshot hook 專供 CLI exit 3／4、資料庫與發布門檻整合測試。 |

| 原測試或明確分支 | 對應結果 | 判定 |
| --- | --- | --- |
| `ReplayWorkspacePreservesExactSearchResultsAcrossRuns`、`CandidateObserverDoesNotChangeSearch` | 有限 sweep budget 下比較同一正式流程；observer 測試另比對認證與搜尋狀態。 | Equivalent |
| `UnchangedCandidateSkipsOnlyCandidateReplay` | 零觀測使候選在正式停機前 bitwise unchanged，仍驗證候選 replay 略過、完整 sweep replay 和工作計數。 | Equivalent |
| `ProductionBoundedSearchReleasesTheStartingUpperFace`、`BcdAcceptsAtMostOneTrustedLocalUpdatePerBlockVisit` | 原 sweep 限制、獨立投影 KKT 及每次 visit 的斷言，使用正式認證。 | Equivalent |
| `EverySweepReconcilesObjectiveAfterReplayAudit`、`SignedCancellationQualifiesWhenOnlyKktGateIsSatisfied`、`HardCaseIncrementalObjectiveAndReplayQualification` | 直接依正式認證及 conditioning 規則取得證據。 | Equivalent |
| `PreparedConditioningUsesValidatedContextWidthBound`、`MultiBlockDiagnosticsLeaveTheAcceptedFitUnchanged` | 比較繼承與明示 Width bound 的合格狀態；開關只控制 block 記錄，兩種執行皆有正式 conditioning。 | Equivalent |
| `CanonicalReasonReachesComponentProgressAndCapturedResult` 的 bypass variant | 移除跳過認證卻標記搜尋成功的人工狀態。保留合格結果及真實未認證／未合格的 sweep budget 與失敗終態。 | Not Equivalent |
| `InitialAcAcceptanceDoesNotHideAWidthSubproblemFailure`、`TrustedInitialProfileCannotCompleteFailedWidthRecovery`、`WidthFailureIsScopedToOneBlockVisit` | 原注入與斷言保留，驗證初始信任、Width failure、visit 作用域與復原。 | Retained Injection |
| `TransientWidthFailureRecoversAndAssemblesQualifiedComponent`、`WidthRecoveryExhaustionAndUnqualifiedTerminalStatesRejectAssembly` | 同一實際 component、搜尋及組裝入口，保留原注入與發布斷言。 | Retained Injection |
| `FailedBlockRetriesAfterLaterBlocksUpdateTheNeighborState`、`ExhaustedBlockRecoveryLeavesTheLastAcceptedComponentStateUntouched` | 改由 Width failure 進入共用復原與 rollback。鄰居更新、失敗 log、accepted state 位元及 replay 照原契約檢查。原觸發原因為人工拒絕可信 trial。 | Partially Equivalent |
| `PracticalCriteriaDoNotQualifySweepsWithFailedBlocks` | 現有 Profile LM 預算為零，沒有有效或可信的初始 trial，且沒有 Width 子問題故障。搜尋照正式規則耗盡且未認證。 | Partially Equivalent |
| `CanonicalReasonReachesComponentProgressAndCapturedResult` 的 recovery budget／exhausted 變體 | 用 scoped Width 失敗取代人工丟棄 profile，原 terminal reason、progress 及 capture 斷言保留。 | Partially Equivalent |
| `CandidateLocalObjectiveIncreasePreservesAcceptedBeforeState` | 自然第二輪拒絕局部 objective 上升，對照第一輪的 eta、beta、prediction、residual、objective 全數逐位元。 | Equivalent |
| `CandidateConsistencyAndNonfiniteFaultsPreserveAcceptedBeforeState` 的 finite delta mismatch、NaN objective／parameter、Inf prediction／residual | 在 affected-row／Full Replay 兩條路上保留五個注入，直接驗證拒絕前後狀態、計數、認證與 integrity。Runtime 的 NaN parameter 發布防護也保留。 | Retained Injection |
| 同一候選測試的 Inf objective／parameter | 這兩個值原先與對應 NaN 通過同一 `isfinite`／`allFinite` 拒絕；保留 NaN 候選的完整路徑，以及 NaN／Inf 的 enclosure 邊界。Inf 不再於這兩個欄位中途注入。 | Partially Equivalent |
| `SweepReplayMismatchFailsBeforeObjectiveReconciliation`、Runtime 的 cache mismatch 變體 | 保留同一 sweep-end mismatch、objective 未刷新、integrity 及發布 gate。 | Retained Injection |
| `ProblemBuilderMemoryPreflightRejectsBeforeSupportAllocation` | 自動測試改為 `ProblemConstructionPreflightUsesTheEstimatedPeak`。另外將 763 個選取原子置於同一原點，搭配原點偏移 (−1, −1, −1)、格距 0.03125 Å 的 64³ 格點，產生 200,015,872 個重疊 support incidences 實跑 builder，確認原始 row buffer 在 reserve 前回傳 `resource-limit:build-joint-problem:raw-target-row-buffer-exceeds-budget`。該次探測的最大 RSS 約 23 MB，峰值記憶體 footprint 約 7 MB。 | Partially Equivalent |
| `AffectedReplayFallsBackBeforeScratchExceedsMemoryBudget`、`E1MemoryPreflightRejectsBeforeCreatingAcceptedState`、Sweep Budget 及 Runtime 的記憶體拒絕變體 | 同一 global budget override 實際執行 replay fallback、兩條路數值對照及拒絕發布。 | Retained Injection |
| `ExtensionPreservesTheOriginalPathAndQualifiesProductionAssembly`、`UnqualifiedHardLimitAndMemoryRejectionCannotPublishComponentState`、`MultipleComponentsExtendAndStopIndependently` | 以唯讀 search policy 提供相同 core、預算及 fault 設定，保留正式組裝、獨立元件進度與發布斷言。 | Retained Injection |
| `ComponentValidationProvidesReasonWithoutEnteringSearch`、`ComponentViewPreflightPreservesReasonAndMemoryDiagnostics` | 以零 block solve／零 certification attempt 及既有 progress、state、reason 斷言驗證預先拒絕，移除專為觀察 configure callback 的旗標。 | Equivalent |
| `ExtendedSearchQualifiesAssemblyAndFeedsEveryDownstreamStage` | 原延長 sweep、合格組裝與相同 snapshot 的 Halo／Peeling／Local MDPDE／Group 測試保留，直接呼叫 `RunJointPostprocessing`。全正式 workflow 另由 `RunsFullPureWorkflowWithGeometryLocalHaloSources` 驗證。失去單一 workflow 呼叫中途替換搜尋結果的路徑。 | Partially Equivalent |
| `JointFailureAndUnavailableOutcomesPersistThroughApiCliAndReload` | 在 capture 後注入 search 不合格或沒有 donor，保留原 API outcome、CLI exit 3／4、SQLite reload 與 stage gate。 | Retained Injection |
| `ProductionBoundedJointDefaultPersistsReloadsAndFeedsPostprocessing` | 以同模型、selection 和正規化 map 的獨立正式求解比對完整數值記錄。保留正式 CLI、SQLite 與各下游檢查；原先從同一次命令內擷取 pre-postprocessing fit 的能力移除。 | Partially Equivalent |

## D. 程式碼簡化量

- 移除 5 個測試控制欄位或狀態，包含 3 個 BCD 開關、1 個 forced-profile callback 欄位及 1 個建構預算 thread-local 狀態。另以借用指標取代 1 個 component configure callback 欄位。
- 原 callback 介面移除 3 個，分別為 forced-profile、component configure 及 workflow fit。保留 1 個作用域較窄的 postprocessing snapshot callback。候選 enum 保留，移除 3 個枚舉值。
- `RHBM_GEM_TEST_INSTRUMENTATION` 的來源 `#ifdef` 由 30 處降至 21 處。
- 不含本報告，`src/` 新增 55 行、刪除 124 行，`tests/` 新增 185 行、刪除 194 行，合計新增 240 行、刪除 318 行。
- GTest 註冊案例 1097 → 1097。移除 1 個 builder 人工預算案例，加入 1 個 Problem Construction planner 案例。移除和改寫的參數分支另列於矩陣。

## E. 數值回歸

以原 Release preset 的既有建構目錄與 Apple Clang 21 arm64 執行。Build 最多 4 個工作，CTest 序列，OpenMP 上限 4、巢狀層級 1。數值快照和 CLI 持久化對照用單執行緒。初始完整測試 21／21，未見基準失敗。P1 初版測試 fixture 在正式停機前沒有 unchanged candidate，改用零觀測輸入。P2 初版預期錯誤的 stop reason，改以無 Profile 評估、無可信初值、無 Width failure 及原復原狀態驗證。P5 的獨立求解對照只出現元件 conditioning 的實測耗時差異，已在數值比較中正規化該欄位。上述驗證均在各自提交前重跑，容差未調整。

固定 `chain(3)` 使用初始 τ = 0.55、上界 τ = 2.2 及 practical-v2。單 block 使用 12-atom core 上限，雙 block 使用 2-atom core 上限與 5 sweep 上限。P0 與 P6 的完整十六進位數值快照使用 `cmp` 對照，SHA-256 均為 `039c51a516df4870e4a768214394c31e582392e3552500c81b46bd5a0611ae8b`。表列共同數值。

| 參數 | 單 block | 雙 block affected-row／Full Replay |
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

| 情境 | Stop reason | Sweeps | 復原嘗試／成功 | Accepted-state integrity | Affected／Full 候選 replay | 認證／Independent Replay |
| --- | --- | ---: | --- | --- | --- | --- |
| 單 block | `block-stationary` | 3 | 0／0 | 通過 | 3／0 | 通過／通過 |
| 雙 block | `block-stationary` | 5 | 0／0 | 通過 | 10／0 | 通過／通過 |
| 雙 block 記憶體回退 | `block-stationary` | 5 | 0／0 | 通過 | 0／10 | 通過／通過 |
| 暫時 Width 失敗 | `block-stationary` | 5 | 1／1 | 通過 | 9／0 | 通過／通過 |
| 持續 Width 失敗 | `block-recovery-exhausted` | 1 | 1／0 | 失敗 | 1／0 | 未執行 |
| Global State 超預算 | `resource-limit:global-state-evaluation:estimated-peak-exceeds-budget` | 0 | 0／0 | 未建立 | 0／0 | 未執行 |

- Profile LM 接受軌跡、sweep objective 對帳、affected-row scratch、Width subproblem solve／QR／failure、profile evaluation 和接受更新計數均與 P0 相同。正式單 block 組裝、mapping 及 assembled-consistency 通過，組裝 objective 兩版皆為 `0x1.bc2eb370d637ep-33`。
- 正式 CLI 的 Joint JSON、SQLite schema v21、4 筆 Joint payload 及 4 筆 stage payload 均精確相同。僅正規化牆鐘耗時、暫存輸入路徑、來源／建置雜湊及同一 run 的 UUID 與引用。Joint JSON schema v24，stage payload v5。
- 114 個正式 library 與 1 個 CLI 編譯命令均無測試巨集。含巨集的 207 個命令只屬於測試 library／executable。兩個正式 binary 的 symbol 與字串查核也排除測試 hook 與注入 Width LM 狀態。

## F. 剩餘 Instrumentation

| 位置 | 區塊數 | 用途 |
| --- | ---: | --- |
| `PotentialFittingWorkflow.hpp/.cpp` | 3 | 宣告、thread-local 儲存及呼叫 `JointPostprocessingSnapshotHookForTesting`。只用於命令整合的非合格結果和 Halo donor 缺失。 |
| `Numerics.hpp`、`ProfileLmSearch.cpp` | 2 | Width 子問題故障的 context 布林值與既有 LM work 計數。 |
| `BcdPolicy.hpp` | 2 | 測試 search policy 的 forward declaration 與 borrowed pointer。 |
| `GaussSeidelBcd.hpp` | 3 | 候選事件／fault enum、區塊診斷／記憶體預算／候選 observer／replay faults 欄位、逐 visit Width fault。 |
| `GaussSeidelBcd.cpp` | 11 | Block 記錄 1、Global State 預算 1、Width visit 1、候選 fault 1、增量 observer 5、sweep objective mismatch 1、測試 search policy 1。 |
| `tests/CMakeLists.txt` | 建構設定 | 只在測試 target 定義巨集，檢查正式 library 與 CLI 的巨集隔離。 |

`GaussSeidelBcd.cpp` 的五個 observer 區塊分別定義事件組裝、記錄非有限候選、記錄 delta mismatch、記錄全域 objective 增加拒絕，以及在狀態提交前記錄最後接受判定。

`candidate_observer` 保留受影響列、候選 objective 與接受判定的逐筆稽核。成功 sweep 會重設 objective cache，終態無法還原提交前的候選。測試在每筆觀察發生時與獨立 Full Replay 比對。
