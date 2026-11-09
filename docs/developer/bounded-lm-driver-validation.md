# Bounded LM Driver 重構驗收

## Git

在既有 `develop` 上執行。起始 HEAD 為 `c28c797b35d892ef1c0c7a298f9ac83364e5b820`，與本機 `origin/develop` 及 GitHub 唯讀查詢的分支 SHA 相同。起始工作目錄乾淨，沒有未提交或未追蹤的檔案。

| 階段 | Commit |
|---|---|
| P0 | `c42bae8d` `test(joint): establish bounded LM behavioral baseline` |
| P1–P2 | `59b1c6ee` `refactor(joint): implement project-owned bounded LM driver` |
| P3–P4 | `2557ab9f553609703cd245314530e421d04d31d7` `refactor(joint): migrate profile search and retire Eigen-derived wrapper` |
| 驗收紀錄 | 本檔的新增 commit 為交付 HEAD。完整 SHA 可由 `git log -1 --format=%H -- docs/developer/bounded-lm-driver-validation.md` 取得。 |

交付後工作目錄乾淨。GitHub 最後一次唯讀確認仍為起始 SHA，`develop` 相對於 `origin/develop` 超前 4 個 commit，落後 0 個。所有 commit 皆保留在本機。

## 責任與相依性

| 舊介面或狀態 | 現在的責任 |
|---|---|
| `minimizeInit` 的維度、容許值、初始殘差、阻尼與計數器 | Driver。以 `BoundedLmOptions`、`BoundedLmState`、`BoundedLmStatus` 和 `Initialize` 表達。Accepted iterations 從零開始，函數評估計數包含初始嘗試及成功評估的候選。 |
| `minimizeOneStep` 的 pivoted QR、縮放梯度與尺度更新 | Driver。Profile 提供由 tiled QR 約化的模型及 column norms。Driver 以同一個 pivoted QR 計算梯度及每次 retry 的模型預測。 |
| `solve_bounded_step`、可行步長與阻尼 | Profile 驗證 Width 上界與目前位置，`SolveBoundedTrustRegionStep` 負責 Active-set、damped QR、半徑及子問題 KKT。Driver 檢查回傳值並用同一對步長與阻尼評估候選。 |
| Actual / predicted reduction、ratio、信賴域更新 | Driver。Actual 使用完整殘差 norm，predicted 使用約化模型。接受門檻維持 `ratio >= 1e-4`。 |
| `Trial` 與 `retry` | Profile 的 `AssessTrial` 和 `CanRetry` 負責候選的 Replay、weak-direction safeguard、trace、Profile 評估預算與失敗原因。Driver 分開處理不可信候選與低 ratio 候選。 |
| Stationarity、tiny step、浮點表示與停止 | Driver 檢查小步長、非有限值、不可表示更新、容許值與精度極限。Profile 提供 accepted Width 的投影 KKT 資格及搜尋結果。 |
| Eigen 專用狀態、可選 `Scalar`、數值差分回傳碼、1-based iteration 與常駐暫存 | 已移除。正式 Profile 使用 double 與解析微分。Driver 從 QR 直接取用 permutation 與三角因子，並在接受時交換兩個既有殘差緩衝區。 |

舊程式的輸入檢查、梯度與預測的 QR 運算、ratio 半徑更新及容許值終止沿自 Eigen LM。原專案已另外加入約化模型、`factor * ||D||` 初始半徑、約束步長、候選可信度、實體步長與 Width KKT。重構保留這些數值規則及浮點運算順序。Header-only 模板讓 Profile 和固定測試問題直接使用同一控制流程。Profile 專屬型別留在 `ProfileLmSearch.cpp`，QR 只存活於單次 `Advance`。

`SearchResult::lm_status` 的讀取端是測試診斷 JSON（`tests/support/JointRuntimeJson.hpp`）。`DiagnosticLmStatus` 明確對應舊代碼 `-2` 到 `9`，包含 tolerance、precision、budget 與 problem stop。正式 BCD 讀取 `initial_accepted`、可信候選與 `width_failed`。子問題失敗、非定常 tiny step，以及 accepted Width 未達 KKT 的不可表示步長仍會阻止該 Width block 通過。BCD 的重試與獨立 Endpoint Certification 保留原入口及資格條件。公開結果和 SQLite 格式未變更。

原本 retry 之後的 stationarity 查詢可能讀到最新候選的 certificate。Profile 現在只在 `Accept` 時更新 accepted Width 的 stationarity，讓 tiny-step 與不可表示步長都按 accepted eta 判定。另一個重現案例從 `tau=.98` 朝 `eta_max=-1e-18` 前進，受限步長相加捨入為零，候選會在建立 trace row 前因超界而被拒絕。修正前這會改寫 initial accepted row 的步長、半徑及接受提案欄位。修正後保留原 row 與相同的拒絕重試路徑。

## 檔案與授權

- 新增 `src/core/detail/joint_component/BoundedLmDriver.hpp`。介面使用 `Initialize`、`Advance`、`AcceptedIterations` 與 `FunctionEvaluations`。
- 修改 `ProfileLmSearch.cpp` 以接上單一 Driver。修改 `Search_test.cpp` 以保留並擴充固定輸入、縮小半徑及 accepted-state 迴歸測試。`Numerics_test.cpp` 移除未使用的 unsupported LM include。盤點亦發現該檔案的既有 `Direct` 測試 helper 無呼叫端。
- 刪除正式 `InstrumentedLM.hpp` 與 `tests/support/InstrumentedLM.hpp`。Joint Component 改用 Eigen Core/QR。HRL 的獨立 MDPDE 求根使用 `Eigen::HybridNonLinearSolver`，因此該處仍需要 Eigen 的 `NonLinearOptimization` 模組。
- `BoundedLmDriver.hpp` 內仍有改寫自 Eigen 的縮放、reduction、信賴域更新與終止流程，保留 Eigen 5.0.1 來源、Thomas Capricelli 著作權及 MPL 2.0 通知，並在 `THIRD_PARTY_NOTICES.md` 標記該檔案。Eigen 的模組也標示 MINPACK 來源，因此一併收錄 Eigen 隨附的 MINPACK 授權通知及 University of Chicago / Argonne National Laboratory 致謝。參見 [MPL 2.0 §1.10、§3.1 與 §3.4](https://www.mozilla.org/MPL/2.0/) 與 [Eigen 5.0.1 的 MINPACK 通知](https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.MINPACK)。

## 數值證據

固定輸入使用 `SyntheticJointWorkload("chain",3)`、`tau=(.46,.49,.52)`、`eta_max=log(.53)` 及既有 Profile 預算。初始 eta 是 `(-0.776528789499, -0.713349887877, -0.653926467407)`。表內每列在重構前後一致。完整 17 位有效數字的 trace 以 `diff -u` 比對，涵蓋 accepted eta、候選序列、每輪 objective 與 KKT、Replay、計數器及停止狀態。

| Profile 評估 | 步長（log-width） | 信賴半徑 | 阻尼 | Actual | Predicted | Ratio | Accepted update |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `(0, 0, 0)` | 0 | 0 | 0 | 0 | 0 | 0（初始） |
| 2 | `(.0826690046727, .0203831597713, -.0389109997399)` | .0290531450193 | 0 | .999926678309 | .982515000917 | 1.01772153848 | 1 |
| 3 | `(.000716010571, -.000187053645, -.000309894205)` | .0581062900386 | 0 | .993623671690 | .993618073232 | 1.00000563442 | 2 |
| 4 | `(-6.03496516e-8, -1.45892533e-8, -8.70510040e-9)` | .000503430554233 | 0 | 9.64869631e-7 | 9.64864918e-7 | 1.00000488430 | 3 |
| 5 | `(-3.50838635e-13, 1.21172471e-13, -1.55675733e-14)` | 3.92921562e-8 | 0 | -1.10578213e-13 | 3.30595660e-17 | -3344.81744581 | 拒絕 |

最終 eta 為 `(-0.69314383460502904, -0.69315379634036556, -0.69314737005669791)`，objective 為 `2.019907813286478e-10`，primary 與 independent projected Width KKT 均為 `3.4194869158454821e-14`。獨立 Replay 與 endpoint qualification 通過。共 5 次 Profile 評估、4 次微分、3 次 accepted update、4 次 bounded solve 及 4 次 damping QR。子問題最大 KKT 為 `5.4905352408500466e-19`，無失敗、bracket solve、working-set release 或 active bound。停止碼為 `2`，原因為 `native-lm-stop`，`width_failed=false`。永久測試 `ProfileSearchPreservesAcceptedTrajectory` 將 update budget 固定為 3，驗證前三次更新、計數器及 endpoint 證據，讓跨平台測試避開最後一輪接近機器精度的 actual reduction。

其他對照包含以下結果。

- `BoundedLmUsesTheConstrainedStepForTrustRatio` 從 `eta=(10,10)` 前進，第一維步長為 `0.2`。依實際 constrained step 與回傳阻尼獨立計算 actual、predicted 及 ratio。上界定常測試在 `(10.2,0)` 得到 `bounded-stationary`。
- `TrustedLowRatioAndUntrustedTrialUseDifferentRetries` 的兩次求解讀取同一個 accepted eta 與 model residual。第一輪低 ratio 的可信候選使半徑縮為原本的 `0.1`，Replay 類不可信候選使半徑縮為 `0.25`。兩者都在第二輪才接受。
- 完整殘差與約化模型的測試以同一個 7×2 問題比較 dense 與 tiled 模型的實際步長、預測與 ratio。近單位 Width、極小尺度 Jacobian、非定常 tiny step、rank deficiency、上界可行性及 Replay disagreement 由重點測試覆蓋。
- 注入一次 Width 子問題失敗後，BCD 從 `retry-pending` 進入第二次 sweep 的成功重試。測試確認 1 次 recovery attempt、1 次成功、0 個 unresolved failure，global A/C KKT `<=1e-7`、projected Width KKT `<=1e-9`、上界違規為 0，且 endpoint、獨立 Replay 及組裝皆通過。持續失敗及不合格 endpoint 的組裝拒絕也通過測試。

## 建構與測試

Apple M1（16 GiB）、macOS 26.6.2 arm64、Apple Clang 21、Eigen 5.0.1、Ninja。使用既有 Release 與 Debug CMake 目錄，SYSTEM dependencies，`BUILD_TESTING=ON`、UMAP 與 experimental feature 開啟、Python bindings 關閉。OpenMP 5.1 可用。正式 `rhbm_gem` 與 CLI 使用未加測試 instrumentation 的編譯設定，`rhbm_gem_test` 與重點測試使用 `RHBM_GEM_TEST_INSTRUMENTATION`。CTest 保留個別測試的執行緒限制。

```sh
export OMP_NUM_THREADS=4
cmake --build build/release --parallel 4
ctest --test-dir build/release --output-on-failure --parallel 4
build/release/bin/RHBM-GEM-TEST --gtest_filter='JointComponentSearchTest.*:ActiveSetBoundedLMTest.*:JointComponentRuntimeTest.*:JointBcdSearchTest.*:JointLocalConditioningTest.*:JointComponentNumericsTest.*:JointTestNumericsTest.*'
cmake --build build/debug --parallel 4
build/debug/bin/RHBM-GEM-TEST --gtest_filter='JointComponentSearchTest.*:ActiveSetBoundedLMTest.*:JointComponentRuntimeTest.*:JointLocalConditioningTest.*:JointBcdSearchTest.EndpointCertificationRetryScheduleIsPreserved:JointBcdSearchTest.ProductionBoundedSearchReleasesTheStartingUpperFace:JointBcdSearchTest.WidthFailureIsScopedToOneBlockVisit'
printf '#include "BoundedLmDriver.hpp"\n' | /usr/bin/c++ -std=c++20 -I src/core/detail/joint_component -isystem /opt/homebrew/include/eigen3 -x c++ -fsyntax-only -
```

| 檢查 | Passed | Failed | Skipped |
|---|---:|---:|---:|
| P0 Release CTest | 21 | 0 | 0 |
| P0 新增案例後的 Search gtest | 14 | 0 | 0 |
| P1–P2 Driver 遷移後的 Search gtest | 15 | 0 | 0 |
| P3 捨入超界案例，修正前 | 0 | 1 | 0 |
| 最終 Release CTest | 21 | 0 | 0 |
| 最終 Release 重點 gtest | 121 | 0 | 0 |
| 最終 Debug 重點 gtest | 58 | 0 | 0 |

Release 與 Debug 均完成建構。完整 CTest 使用 Release 執行，Debug 另執行表中的重點測試。Debug 建構在既有資料層測試出現 aggregate initializer 警告。

## 效能與記憶體

以同一個 Release instrumentation 設定連結舊、新靜態測試函式庫，順序交替，各規模各執行五次，每次啟動新程序。重用 `SyntheticJointWorkload("chain", targets)`、全體 `tau=.55`、`eta_max=log(2.2)`、`PracticalJointConvergenceCriteria` 及預設 12 targets/block。以 `SearchBcdComponent` 執行兩次 sweeps。時間只涵蓋 BCD 呼叫，峰值 RSS 取程序的 `getrusage(RUSAGE_SELF)`。評估及更新取 BCD 計數，QR 欄為 bounded damping QR 次數。各規模的 sweep、工作計數、最後 objective 與 projected KKT 在兩個版本一致。固定 sweep 預算用於量測，獨立 endpoint 資格另由上方測試驗證。

| Targets | BCD wall 中位數（範圍）前 → 後，ms | 峰值 RSS 中位數 前 → 後，MiB | Profile 評估 | Accepted 更新 | Bounded solves | Damping QR |
|---:|---|---:|---:|---:|---:|---:|
| 3 | 4.880（4.690–6.460）→ 4.643（4.593–4.869） | 21.88 → 21.83 | 4 | 2 | 2 | 2 |
| 30 | 42.053（41.748–42.611）→ 42.123（41.786–43.017） | 31.66 → 31.62 | 12 | 6 | 6 | 6 |
| 300 | 417.399（414.732–423.833）→ 416.545（413.669–427.512） | 79.17 → 79.06 | 104 | 52 | 52 | 110 |
| 1,000 | 1,406.598（1,393.871–1,409.882）→ 1,397.235（1,382.061–1,411.774） | 210.80 → 210.62 | 336 | 168 | 168 | 226 |

四組的前後耗時範圍重疊，峰值 RSS 中位數差距小於 0.2 MiB。這些量測未見可確認的效能或記憶體退化。

可在既有建構目錄以 `python3 build/release/build-bounded-lm-evidence.py` 重連結量測程式。以 `OMP_NUM_THREADS=4 build/release/bounded-lm-evidence 3 profile` 輸出 Profile 軌跡。移除 `profile` 參數後，第一參數為 BCD synthetic target 數量。量測程式及原始執行紀錄位於本機忽略的 `build/release/bounded-lm-*`。上述測試輸入、計數與摘要納入 Git。這次驗收採用合成資料，6Z6U / EMD-11103 未納入效能數據。

## 結論

單一正式 `BoundedLmDriver` 取代 `InstrumentedLM`，保留測試診斷狀態碼與 Eigen 衍生授權。測得的候選與 BCD 工作量等價，Width KKT、Replay、recovery、Certification 及組裝測試維持原資格。唯一有意調整的行為是讓拒絕後的定常性與初始接受證據繼續對應 accepted state。
