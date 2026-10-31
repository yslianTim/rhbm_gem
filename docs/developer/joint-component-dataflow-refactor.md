# Joint Component Estimator 資料流重構驗收報告

日期：2026-10-31。審計與數值比較基準：`9ba33414130535e1b93e9fddf1967ed516a16eb9`。

## A. Git 與執行範圍

- 初始 branch：`develop`；初始 working tree clean；初始本地與遠端實際 HEAD 都等於上述基準。
- 所有修改、階段 commits 都直接在既有 develop 完成；沒有建立 branch、worktree，也沒有切換分支。
- P1：`c28e02ba` — `refactor(joint): remove redundant numerical arguments and context copies`。
- P2：`e9e1be74` — `refactor(joint): reduce persistent prepared block state`。
- 數值驗證的 code HEAD：`4e0e4827e894a06263dac9940871a0f3e63b3429`。
- P3：`4e0e4827` — `refactor(joint): simplify width policy propagation`。
- 本報告隨最後的驗收文件 commit 保存；最終完整 HEAD 與 push 後三方同步狀態列於交付訊息。
- 沒有新增第三方依賴、修改 Notion、JSON/SQLite schema 或無關 production 模組。
- 所有 CMake build 使用 `--parallel 4`、CTest 使用 `--parallel 4`，OpenMP 上限 4。正式 CLI `-j 4`；保留 Joint 群組原有 `OMP_NUM_THREADS=1` 與 solver Eigen serial scope。

## B. P0 參數與欄位審計

審計涵蓋 Joint 目錄 41 個 C++ source/header，向整個 `src/`、`include/`、`tests/` 追蹤呼叫／欄位讀取；另檢查 workflow、command、serialization、progress 及 downstream。方法包括 `rg`、header/definition 比對、相關 Git history，以及 Release compile commands 的 `-fsyntax-only -Wunused-parameter -Wunused-private-field -Wunused-function`。

修改前、修改後各對 32 個相關 translation units 的 instrumentation 開／關版本編譯，共各 64 次；均無錯誤或 unused diagnostics。匿名未讀取參數另經人工讀取定義確認，未僅依 compiler silence 判定有用途。測試 support / dense reference / failure injection 的呼叫也包含在內。

分類：A 明確無用途；B 安全可推導或只需暫存；C 有條件重複；D 獨立必要責任。

| 類別 | 盤點項目 | Removed | Simplified | Kept | Deferred |
|---|---:|---:|---:|---:|---:|
| Function parameters | 5 | 3 | 1 | 1 | 0 |
| Context fields / copies | 14 | 3 | 0 | 11 | 0 |
| Policy fields / copies | 7 | 0 | 1 | 6 | 0 |
| Prepared Block fields | 9 | 3 | 0 | 6 | 0 |
| Intermediate results | 13 | 0 | 0 | 13 | 0 |
| Progress parameters | 2 | 0 | 1 | 1 | 0 |
| Helpers / overloads | 7 | 1 | 1 | 4 | 1 |

數量以以下可追蹤審計項目計算，包含明列的 Context/Policy copies；同一 helper 的簽章與個別參數分屬不同審計面向，不代表互斥的刪除行數。位置以初始 HEAD `9ba33414` 為準。未刪除 production serialization schema 欄位。

| Symbol | 定義位置 | 呼叫／讀取位置 | 分類 | 決策 | 風險／依據 |
|---|---|---|---|---|---|
| CheckTrust(replay)::Domain | Numerics.hpp:159 / Numerics.cpp:446 | 定義完全未讀取；CheckTrust 四參數入口、CertifyFixedNeighborEndpoint、Runtime_test 兩處 | A | Removed | 低；只刪簽章與轉傳，TrustEvidence 指派／條件不變 |
| CheckTrust(replay)::y | Numerics.hpp:159 / Numerics.cpp:446 | 定義完全未讀取；CheckTrust 四參數入口、CertifyFixedNeighborEndpoint、Runtime_test 兩處 | A | Removed | 低；只刪簽章與轉傳，TrustEvidence 指派／條件不變 |
| CheckTrust(replay)::EvaluationContext | Numerics.hpp:159 / Numerics.cpp:446 | 定義完全未讀取；CheckTrust 四參數入口、CertifyFixedNeighborEndpoint、Runtime_test 兩處 | A | Removed | 低；只刪簽章與轉傳，TrustEvidence 指派／條件不變 |
| PrepareLocalEndpointWorkspace::eta_upper_bound | FixedNeighborBlockCoordinate.cpp:260 | 唯一 caller 為 SearchFixedNeighborComponent 的 conditioning；準備用 Context 已通過入口 validation | B | Simplified | 低；保留 missing/nonfinite/mismatch 拒絕 |
| PrepareDerivative / PrepareProjectedWidthDerivative::scale | TiledDerivative.hpp:31 | Search、LocalConditioning、dense reference；可配合 null Context 或顯式 normalization | D | Kept | 獨立數值工具契約，不能假定等於任意 Context scale |
| EvaluationContext::snapshot_hash | Numerics.hpp:30 | CreateContext；tests/support/JointTestContext.cpp:56、JointComponentChecks.cpp:61 | D | Kept | 測試 snapshot/provenance evidence；Production 為空也不代表無讀取者 |
| EvaluationContext::observations | Numerics.hpp:30 | CreateContext / Observe；SolveFixedNeighborComponent 擷取 local_y | D | Kept | 父 observation view 與不可變 snapshot ownership，不能由局部 effective response 取代 |
| EvaluationContext::atom_ids | Numerics.hpp:30 | Partition、Search 的尺寸驗證、FailureDiagnostics、Independent Replay 最大 gradient atom ID | D | Kept | identity、mapping 與失敗診斷 |
| EvaluationContext::row_ids | Numerics.hpp:30 | Search 尺寸驗證、MakeComponentProblemView、prepared row selection、測試 evidence | D | Kept | 局部列身份與 mapping validation |
| EvaluationContext::scale | Numerics.hpp:30 | EvaluateProfile、streaming evaluation、Replay、Search、ComponentState、Assembly | D | Kept | Original parent objective normalization；不能改用 local norm |
| EvaluationContext::rank | Numerics.hpp:30 | TiledDerivative、Search weak-width assessment、LocalConditioning、tests ContextEvidence | D | Kept | 父列數與 family dimensions 的 threshold 責任 |
| EvaluationContext::linear | Numerics.hpp:30 | EvaluateProfile / SolveLinear；CreateContext、reconstruction domain 設定 release_response_norm | D | Kept | strict local A/C correctness、parent response release threshold |
| EvaluationContext::profile_budget | Numerics.hpp:30 | SearchProfile、Profile::Get/retry、Progress callback | D | Kept | evaluation budget 與 progress |
| EvaluationContext::update_budget | Numerics.hpp:30 | SearchProfile；prepared 設定為 1；Progress | D | Kept | 每次 block visit 最多一次 accepted update |
| EvaluationContext::eta_upper_bound | Numerics.hpp:30 | Search 入口與投影、EvaluateState、Replay、LocalConditioning | D | Kept | 驗證後 effective bound，與入口政策責任不同 |
| EvaluationContext::force_width_subproblem_failure_for_testing | Numerics.hpp:30 | Search::Profile::solve_bounded_step；FixedNeighbor test-only visit Context | D | Kept | instrumentation 開啟時的故障注入；關閉時不編譯 |
| SearchProfile::search_context copy | Search.cpp:202 | Profile 僅同步讀取；InstrumentedLM 與所有 references 均於函式內結束 | A | Removed | 低；借用 caller const Context，無 callback/async 逃逸 |
| Profile::scale | Search.cpp:42 | 唯一初始化等於 search_context.scale；4 個數值讀點改讀 context.scale | B | Removed | 低；值、公式及浮點運算順序不變 |
| LocalEndpointWorkspace::context | FixedNeighborBlockCoordinate.cpp:258 | 只複製 prepared.context 並覆寫相同 bound，再立即傳給 AssessFixedNeighborLocalConditioning | B | Removed | 低；直接借用仍存活的 prepared.context |
| FixedNeighborSearchPolicy::eta_upper_bound | FixedNeighborPolicy.hpp:24 | FitJointComponentsInitialized → Orchestrate → Solve → Search；測試可直接建立低層 policy | C | Kept | 入口設定與 effective Context 間的 mismatch validation 必須保留 |
| FixedNeighborSearchPolicy::component_index | FixedNeighborPolicy.hpp:24 | SolveFixedNeighborComponentView 覆寫為 view identity；FailureDiagnostics、block ID 與 test diagnostic facade | C | Kept | 移除需新增 identity 入參／adapter；保留失敗格式及獨立低層測試 |
| FixedNeighborSearchPolicy::component_id | FixedNeighborPolicy.hpp:24 | SolveFixedNeighborComponentView 覆寫為 view identity；FailureDiagnostics、block ID 與 test diagnostic facade | C | Kept | 移除需新增 identity 入參／adapter；保留失敗格式及獨立低層測試 |
| FixedNeighborPolicy::eta_upper_bound | FixedNeighborBlockCoordinate.hpp:121 | FitJointComponentsInitialized → Orchestrate → Solve → Search；測試可直接建立低層 policy | C | Kept | 入口設定與 effective Context 間的 mismatch validation 必須保留 |
| FixedNeighborPolicy::component_index | FixedNeighborBlockCoordinate.hpp:121 | SolveFixedNeighborComponentView 覆寫為 view identity；FailureDiagnostics、block ID 與 test diagnostic facade | C | Kept | 移除需新增 identity 入參／adapter；保留失敗格式及獨立低層測試 |
| FixedNeighborPolicy::component_id | FixedNeighborBlockCoordinate.hpp:121 | SolveFixedNeighborComponentView 覆寫為 view identity；FailureDiagnostics、block ID 與 test diagnostic facade | C | Kept | 移除需新增 identity 入參／adapter；保留失敗格式及獨立低層測試 |
| FitFixedNeighborComponentsImpl::bounded_policy copy | JointComponentEstimator.cpp:755 | 複製後從未修改，只讀 bound 並傳入 Orchestrate | A | Simplified | 低；直接 const-reference production_policy |
| PreparedFixedNeighborBlock::block | FixedNeighborBlockCoordinate.hpp:112 | conditioning block ID | D | Kept | 必要 block identity |
| PreparedFixedNeighborBlock::core_atoms | FixedNeighborBlockCoordinate.hpp:112 | 每次 visit 擷取／更新 accepted core、candidate observer、conditioning evidence | D | Kept | 保持 core ordering |
| PreparedFixedNeighborBlock::affected_rows | FixedNeighborBlockCoordinate.hpp:112 | effective response、incremental objective、candidate prediction/residual、replay mapping | D | Kept | 保持 affected-row ordering |
| PreparedFixedNeighborBlock::profile_atoms | FixedNeighborBlockCoordinate.hpp:112 | 只在 preparation Select、atom identities、local dimension 讀取 | B | Removed | 改為 per-block 區域暫存；順序不變 |
| PreparedFixedNeighborBlock::profile_rows | FixedNeighborBlockCoordinate.hpp:112 | 只在 preparation mapping、Domain row count、row identities 讀取 | B | Removed | 改為 per-block 區域暫存；順序不變 |
| PreparedFixedNeighborBlock::row_mapping | FixedNeighborBlockCoordinate.hpp:112 | 只在 preparation Select 與 mapping 統計讀取；DomainAtoms::Select 已共享持有 | B | Removed | 只移除額外 owner；mapping data、Find 與 shared ownership 保留 |
| PreparedFixedNeighborBlock::domain | FixedNeighborBlockCoordinate.hpp:112 | SearchProfile、local Design、LocalConditioning；DomainAtoms owns snapshot/parent/mapping/supports | D | Kept | 不可用失效 reference 或大型 Problem copy 取代 |
| PreparedFixedNeighborBlock::context | FixedNeighborBlockCoordinate.hpp:112 | local SearchProfile、FailureDiagnostics、LocalConditioning；test-only visit copy 可修改故障旗標 | D | Kept | 局部 identity/rank/budget 與父 scale/linear policy 有各自責任 |
| PreparedFixedNeighborBlock::local_atoms | FixedNeighborBlockCoordinate.hpp:112 | local Eigen vector dimensions、rank 初始化、maximum_block_columns | B | Kept | 保留一次建立的 Eigen::Index dimension，避免反覆轉型 |
| ComponentProblemView::parent_row_count | Problem.hpp:32 / Components.cpp:66 | MakeComponentProblemView 保存完整 view.rows.size；Solve 設定 rank.rows，不能改成 compact informative-row count | D | Kept | 保留 parent/local/fitting/constant row 與 identity 責任邊界 |
| ComponentProblemView::local_to_parent_atoms | Problem.hpp:32 / Components.cpp:66 | initial widths selection、root atom mapping、identity validation | D | Kept | 保留 parent/local/fitting/constant row 與 identity 責任邊界 |
| ComponentProblemView::local_to_parent_rows | Problem.hpp:32 / Components.cpp:66 | 從 parent observations 擷取 response；local row identity／support mapping | D | Kept | 保留 parent/local/fitting/constant row 與 identity 責任邊界 |
| ComponentProblemView::local_layout | Problem.hpp:32 / Components.cpp:66 | compact local parameter indices 與 informative rows；不是 parent layout | D | Kept | 保留 parent/local/fitting/constant row 與 identity 責任邊界 |
| ComponentProblemView::domain | Problem.hpp:32 / Components.cpp:66 | nested compact support view，供 Search 使用 | D | Kept | 保留 parent/local/fitting/constant row 與 identity 責任邊界 |
| ComponentProblemView::atom_ids | Problem.hpp:32 / Components.cpp:66 | 原 snapshot atom identities 的 indexed view，供 local Context／diagnostics | D | Kept | 保留 parent/local/fitting/constant row 與 identity 責任邊界 |
| ComponentProblemView::row_ids | Problem.hpp:32 / Components.cpp:66 | original parent row identities 的 indexed view，供 local Context／validation | D | Kept | 保留 parent/local/fitting/constant row 與 identity 責任邊界 |
| FixedNeighborComponentResult::search | FixedNeighborBlockCoordinate.hpp:170 | 候選狀態、KKT、certification、accepted integrity、work evidence；Components 組裝正式輸出，tests diagnostic facade 亦讀取 | C | Kept | 移除單一欄位需擴大修改多個早期失敗路徑；維持既有語意 |
| FixedNeighborComponentResult::succeeded | FixedNeighborBlockCoordinate.hpp:170 | search_converged 且 accepted-state integrity/certification qualified 的一次性結果；Components 組裝正式輸出，tests diagnostic facade 亦讀取 | C | Kept | 移除單一欄位需擴大修改多個早期失敗路徑；維持既有語意 |
| FixedNeighborComponentResult::stop_reason | FixedNeighborBlockCoordinate.hpp:170 | Solve 的 unobserved/invalid-initial-widths 可早於 Search；preflight、Search、Certification/recovery 路徑另有 reason；Components 組裝正式輸出，tests diagnostic facade 亦讀取 | C | Kept | 移除單一欄位需擴大修改多個早期失敗路徑；維持既有語意 |
| LocalEndpointWorkspace::eta | FixedNeighborBlockCoordinate.cpp:255 | PrepareLocalEndpointWorkspace 建立、AssessFixedNeighborLocalConditioning 同步讀取 | D | Kept | 不同 local numerical operands，非純轉傳 wrapper |
| LocalEndpointWorkspace::beta | FixedNeighborBlockCoordinate.cpp:255 | PrepareLocalEndpointWorkspace 建立、AssessFixedNeighborLocalConditioning 同步讀取 | D | Kept | 不同 local numerical operands，非純轉傳 wrapper |
| LocalEndpointWorkspace::effective_response | FixedNeighborBlockCoordinate.cpp:255 | PrepareLocalEndpointWorkspace 建立、AssessFixedNeighborLocalConditioning 同步讀取 | D | Kept | 不同 local numerical operands，非純轉傳 wrapper |
| NotifyJointProgress positional helper | JointProgress.hpp:89 | ComponentStarted、ComponentCompleted、CertificationStarted；沒有 Event 拆散再重建模式 | D | Kept | 保留 phase/counts/objective/KKT/reason/trusted-state 以及延遲建構 Event |
| production_progress Context capture | FixedNeighborBlockCoordinate.cpp:1267 | lambda 唯一讀取 context.profile_budget | B | Simplified | 只 capture profile_budget 值；無 ownership 或 reference 逃逸問題 |
| CheckTrust replay-evidence overload | Numerics.hpp:159 | 四參數 replay 入口、endpoint certification、Runtime_test | A | Simplified | 5→2 參數；獨立 replay qualification 不變 |
| CheckTrust executing-replay overload | Numerics.hpp:158 | tests/support/JointTestNumerics.cpp:69；本身呼叫 CheckIndependentReplay | D | Kept | 真正執行獨立驗證的入口 |
| NotifyJointProgress(observer, event) | JointProgress.hpp:83 | 全 repository 僅 declaration/definition，無呼叫者；既有 Event 直接 observer(event) | A | Removed | 不新增 Event builder／conversion |
| CheckIndependentReplay assembled arguments | Numerics.hpp:145 / Numerics.cpp:544 | 現有 callers 使用 default/nullptr；實作包含 assembled mapping、prediction、objective validation | C | Deferred | 新增一般 endpoint 具名 overload 會增加 API；移除 assembled 分支會碰驗證算法，不在本次範圍 |
| TiledDifferential::Rows convenience overload | TiledDerivative.hpp:14 | JointDenseReference.cpp:118、LocalConditioning_test.cpp:61 | D | Kept | 測試／獨立 dense reference 仍有 caller |
| FitFixedNeighborComponents wrapper | Problem.hpp:53 / TargetEvidence.cpp:265 | 正式 Initialized fit 與測試；執行 Impl 後 AddTargetEvidence | D | Kept | 具有 target qualification 與 Eigen thread scope 責任 |
| SolveFixedNeighborComponentView wrapper | Problem.hpp:41 / Components.cpp:311 | Orchestrate；memory preflight、parent/local mapping validation、identity adaptation | D | Kept | 不是單純轉傳，不能移除前置檢查 |


## C. P1 結果

Replay-evidence `CheckTrust` 從 5 參數降至 2，只接收 endpoint 與 replay；移除的 Domain、response、Context 在舊定義中均未被讀取。保留另一個真正執行 Independent Replay 的四參數入口；TrustEvidence 所有欄位指派、判定條件與 qualification 未變。

SearchProfile 直接借用 caller 的 const Context。Profile、InstrumentedLM 及其 lambda 均在函式返回前銷毀，無 async 或 escaping callback，測試注入也只讀布林旗標。Profile 的 scale 欄位唯一初始化為該 Context 的 scale，因此四個數值使用點直接改讀 `context.scale`，公式與運算順序不變，未重新計算 local response norm。

P1 Release build、7 個 Joint/persistence/coverage CTest entries 通過。六組 deterministic synthetic traces 的三次執行自身一致，且與初始基準逐 bit 一致；沒有另外留下已確認而未處理的低風險候選。

## D. P2 結果與記憶體持有權

Prepared Block 由 9 欄位減至 6 欄位，保留 `block`、`core_atoms`、`affected_rows`、`domain`、`context`、`local_atoms`。profile atom/row selectors 及 mapping handle 僅存活於每個 preparation iteration；Domain selection、identities selection、索引排序、mapping 統計與 support 計算順序未改變。

`DomainAtoms::Select()` 仍以 shared ownership 保存 row mapping、parent domain、snapshot 及 local supports；只刪除了 Prepared 上的額外 mapping owner。沒有以裸指標或借用 reference 替換 immutable snapshots。`DomainAtoms::selected` 與 `Identities::Select()` 的 index maps 各自保存必要索引，不會借用已釋放的 profile vectors。新增生命週期測試確認暫存 handles 與 parent 物件銷毀後，巢狀 mapping 各由 Domain 唯一持有且 support row/squared distance 正確。另有 standalone AddressSanitizer 檢查，100 次巢狀 mapping / Observe 生命週期全部通過；這不是對整個 CLI 的 ASan 宣稱。

Production ABI `sizeof(PreparedFixedNeighborBlock)` 為 **416 → 352 bytes**；instrumentation ABI 為 **424 → 360 bytes**。每個 mapping 在 preparation scope 結束後的直接 owners 為 2 → 1，mapping 配置本身保留。

Memory Preflight 同步調整：

- Global resident estimate 不再計入已釋放的 profile selectors，且 `sizeof(PreparedFixedNeighborBlock)` 自動反映較小物件。
- Preparation estimate 保留一組最大 block 的 selector scratch；尚無 census 的入口使用 component dimensions 作保守上限。
- `PartitionStorageBytes()`、`LocalSupportIndexBytes()`、`CompactRowMapping::StorageBytes()`、mapping/support storage 的計數與 checked arithmetic 保留。
- 3.2 GB planning budget、4.0 GB hard limit、overflow 拒絕與 full replay memory fallback 不變。
- 新增 admission 與 overflow 檢查，避免將所有 selectors 繼續視為同時常駐；原有 fallback/preflight failure tests 通過。

6CVM Chain A，同一固定輸入、同一 configuration、相同 thread limits；只執行 preparation 的 production benchmark 前後各三次：

| 項目 | 修改前 | 修改後 |
|---|---:|---:|
| Prepared blocks | 742 | 742 |
| Core atoms | 8,193 | 8,193 |
| Fitting rows | 508,449 | 508,449 |
| Affected-row entries | 972,329 | 972,329 |
| Prepared mapping entries | 4,149,010 | 4,149,010 |
| Prepared mapping bytes | 57,625,624 | 57,625,624 |
| Global resident estimate bytes | 205,184,753 | 197,293,089 |
| Preparation peak estimate bytes | 212,858,777 | 204,998,553 |
| Preparation median seconds | 0.296651250 | 0.297791333 |
| Process Peak RSS range bytes | 2,543,943,680–2,566,471,680 | 2,192,719,872–2,554,544,128 |

Resident estimate 下降 **7,891,664 bytes**，恰等於 `742*64 + (8193+972329)*sizeof(Eigen::Index)`。這是移除的物件／selector payload accounting，不等同作業系統 RSS 降幅。Mapping bytes 沒有下降，因為仍然必要。RSS 包含載入模型與 map normalization；波動明顯，不能據此宣稱穩定加速或固定 RSS 改善。

P2 Release build、7 個相關 CTest entries、ownership/memory checks 及六組 baseline bitwise traces 通過。

## E. P3 條件式精簡

Width 流程維持 `calibration → production policy → component Context / numerical policy → Search validation → prepared Context`。移除 LocalEndpointWorkspace 的 bound 入參與 Context 複本後，Local Conditioning 直接使用 prepared Context。入口原有 `invalid-global-width-bound`、`invalid-width-upper-bound`、initial projection 與 active-bound KKT 判定不變。新增測試驗證 inherited/explicit bound 的 accepted state、sweeps、reason、conditioning width evidence 相同，mismatch/missing/NaN/infinity 在 preparation 前拒絕。

Progress callback 只捕獲實際讀取的 `profile_budget` 值；移除從未修改的 bounded policy copy，以及完全無呼叫者的 Event helper overload。保留少量 callers 使用的 positional helper，因為沒有先拆 Event 再重建的資料流。

Component identity、Result stop reason、ComponentProblemView 與 EvaluationContext 的最後決策見審計表。保留的 Context copies 會設定本地 identities/rank/budget/bound、保持 assembly policy，或供 instrumentation 逐 visit 注入；不是未使用的只讀副本。Root observations / identities 都保留 snapshot ownership，未切分出新的通用 Context。

Replay API 未加新 overload；assembled mapping/prediction/objective validation 分支也未刪除。此候選 Deferred，原因是一般 endpoint wrapper 會增加 API，而刪除 assembled 分支會擴大到驗證算法。

P3 Release build、7 個相關 CTest entries、新增 Width validation test、六組 baseline bitwise traces 全部通過。沒有修改 production numerical branches 或新增 test-only production branch。

## F. P4 回歸證據

命令：

```sh
cmake --preset release -DBUILD_TESTING=ON
cmake --build build/release --parallel 4
OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  ctest --test-dir build/release --parallel 4 --output-on-failure
```

Release configure/build 成功；完整 CTest **21/21** 通過。Coverage Guard：**1,056 registered/enabled/covered unique、missing=0、intentional duplicates=18、disabled=0**。基準為 1,054；增加的兩個測試為 `NestedDomainsOwnMappingsAfterPreparationHandlesExpire` 與 `PreparedConditioningUsesValidatedContextWidthBound`。原有測試未刪除或改成寬鬆 tolerance。

第一次 baseline 執行把全域 `OMP_NUM_THREADS` 設為 1，造成既有 HRL thread restoration 測試要求的 2 threads 被 runtime 限為 1；原樣重跑也失敗。改為全域 4 並保留 Joint CTest 的 serial 設定後，修改前與修改後都 21/21 通過，沒有為此改動 HRL 程式或測試。Build 重編譯時有三個既有 persistence fixture missing-field-initializer warnings，本次未修改該檔案；相關 unused compiler audit 無診斷。

Synthetic 比較涵蓋 chain/cube、active upper bound、transient recovery、recovery exhaustion 及 nonfinite bound。比較內容包含 A/B/C、eta、accepted trials、LM status/reason、candidate prediction/residual、incremental/local objective、accepted/rejected decisions、core/affected ordering、逐 sweep objective/KKT、mapping statistics、TrustEvidence、replay 與 certification；浮點值以 64-bit 表示比較。P1/P2/P3 每階段與同一初始 trace 完全相同，未採用 tolerance 接受差異。四份 trace 的 SHA-256 均為 `88c0590a34f035e5b86b411f43c32a31b0e1cdee417b47f1a621fca493129fa1`。

6CVM Chain A 使用相同 model/map SHA-256、normalization、selection、solver provenance、configuration hash 與 threads。基準與修改後各重新跑正式 CLI，SQLite quick_check 均為 ok、SQLite schema 20、Joint schema 22。比較時僅排除執行時間、對應下降的 memory estimates、隨變更而改變的 software source/build hashes，以及每次執行的 run UUID；scientific floats 轉成 IEEE-754 64-bit 表示逐項比較。

| 正式數值／契約 | 修改前與修改後 |
|---|---|
| Search | 250 sweeps；block-stationary |
| Objective | 0.019770945620759765 |
| A/C projected KKT | 4.269076906392166e-8 |
| Projected Width KKT | 9.699454295297727e-10 |
| Initialization、component/assembled A/B/C、eta、gradient | 逐 bit 一致 |
| 全部 250 次已輸出的 progress sweep fields | 一致，僅 elapsed time 不比較 |
| Certification / Independent Replay / Assembly | Passed；所有非 timing/memory evidence 一致 |
| Local conditioning | 742 blocks；數值 evidence 逐 bit 一致 |
| Targets / Halos / fitting rows | 8,193 / 1,645 / 508,449 |
| Halo donor provenance / reconstruction | 正規化 run UUID 後逐 bit 一致 |
| Original-domain peeling samples | 1,638,600；結果逐 bit 一致 |
| Local MDPDE | 8,013 available；180 unavailable；輸入與結果逐 bit 一致 |
| Local unavailable 原因 | invalid-local-parameters 176、local-mdpde-unqualified 3、insufficient-positive-responses 1 |
| Group inputs / qualifications | 167 groups；166 qualified；member/eligible/excluded/reason 一致，posterior source/有限值驗證通過 |

Group alpha 使用原有隨機訓練，沒有要求不同執行之間的 Group posterior 逐 bit 一致；沒有為此固定 random_device、修改算法或修改正式預期值。

前後各執行兩次完整正式 CLI，四次都通過資格、250 sweeps，initialization 與 component/assembled states 逐 bit 相同。兩次補充量測依序執行；最初基準與修改版執行期間曾有建構／審計活動，不能將 wall time 差異全歸因於重構。

| 執行 | Full CLI wall seconds | Search seconds | Peak RSS bytes |
|---|---:|---:|---:|
| baseline | 570.640423209 | 495.005441663 | 3,271,655,424 |
| baseline-repeat | 543.663221583 | 472.454804332 | 3,408,920,576 |
| final | 550.041574834 | 475.440803873 | 3,353,198,592 |
| final-repeat | 541.980304792 | 470.425689883 | 3,083,698,176 |

第一次前後 RSS 由 3,271,655,424 增至 3,353,198,592 bytes（+2.49%），因此補做一對完整量測；第二對由 3,408,920,576 降至 3,083,698,176 bytes。Baseline 觀測範圍 3.272–3.409 GB，修改版 3.084–3.353 GB，範圍重疊。沒有觀察到明顯持續的記憶體或時間回歸；兩次樣本與 RSS 波動仍不足以宣稱穩定改善，也不保證每次峰值均下降。結構性的 resident estimate 降幅則有明確 allocation/ownership 依據。

完整 CTest 另外涵蓋 initialization、width calibration/projection、partition、Constrained A/C、bounded active-set LM、Adaptive Sweep Budget、Global KKT、accepted-state integrity、recovery/rollback、Full Replay、certification、assembly、Halo reconstruction、Original-domain peeling、Local MDPDE、Group fitting、progress 及 persistence。正式 core size 12、nominal/max sweeps 200/300、practical A/C KKT 1e-7、projected width KKT 1e-9、strict local correctness 1e-10、one accepted update per block visit 均未修改。

`git diff --check` 通過。原始 logs、compiler audit、baseline binaries、synthetic traces、CLI databases、位元比較與量測 scripts 保存在本機 `/private/tmp/rhbm-joint-dataflow/`；沒有把大型測試資料或研究路徑加入 repository。

## G. 精簡效益

| 項目 | 修改前 | 修改後 | 效益 |
|---|---|---|---|
| Replay CheckTrust 入參 | 5 | 2 | 移除 3 個完全未讀取值 |
| PrepareLocalEndpointWorkspace 入參 | 3 | 2 | 不再重複傳遞 bound |
| SearchProfile Context copy | 每次 search 1 份 | 0 | 保留 const borrowing 生命週期 |
| Local endpoint workspace Context copy | 每 block conditioning 1 份 | 0 | 直接使用 prepared immutable Context |
| Progress Context capture | 每 component observer 1 份 | 1 個 int | 不再攜帶 snapshot/identity/bound 等無關狀態 |
| Unchanged bounded policy copy | 每 fit 1 份 | 0 | 減少 policy/std::function copy |
| Profile scale | 1 個重複 double | 0 | normalization 唯一來源 |
| Prepared persistent fields | 9 | 6 | 不保存 preparation-only selectors/owner |
| Production Prepared object | 416 bytes | 352 bytes | 每 block 減少 64 bytes，另釋放 selectors |
| Mapping direct owners/block | 2 | 1 | Domain 仍擁有 mapping，資料配置不變 |
| 6CVM global resident estimate | 205,184,753 bytes | 197,293,089 bytes | 模型對應實際釋放的 payload |
| Uncalled Event overload | 1 | 0 | 無新增 adapter/framework |
| Preparation time | median 0.296651250 s | median 0.297791333 s | 無明顯回歸；不宣稱加速 |

## H. 結論

1. 確認完全無用途的是 replay-evidence CheckTrust 的 Domain、response、Context 三個入參。
2. 已移除 Profile scale、三種不必要的 Context 保存／複製、未修改的 policy copy，以及三個 Prepared persistent fields。
3. Parent scale/local response、Original/Fitting domains、Primary/Independent Replay、strict/practical KKT、identity/provenance、parent-row rank、snapshot ownership 都有必要責任，保留。
4. 無呼叫者的 Event overload 已移除；剩餘 Fit/Solve wrappers 含 target qualification、preflight 或 mapping validation，並非單純轉傳。
5. Prepared 常駐資料確實減少，mapping/support storage 與排序不變。
6. Local Conditioning 使用 Prepared Context 的 bound；Progress 只保存 budget；Context/Policy 邊界更清楚，未引入新的大型 Context。
7. 正式數值結果一致：synthetic traces 與 6CVM scientific state/evidence、Halo/Peeling/Local 及 Group inputs 均通過逐 bit 比較；Group 隨機 posterior 依正式資格檢查。
8. Memory Preflight、overflow、hard limit、fallback、recovery、replay/certification/assembly 失敗安全性保留；所有相關測試通過。
9. 可另案評估 Context 型別切分、diagnostic identity API、單一 result reason 與 assembled Replay API；本次不為少量欄位引入更多 adapters 或改動失敗語意。Deferred 不代表未完成本次條件式 P3。
10. 可以正式結案：P0–P4 的結構、數值、memory preflight、downstream 與測試驗收均通過；最後文件 commit、working tree clean 及一般 fast-forward push 的核對結果列於交付訊息。
