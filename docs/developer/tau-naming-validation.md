# JCE τ 命名與數值契約驗收

## 範圍與建構設定

起始 HEAD 為 `4e8be0ce95b446654bb157345d12dc928fa6419b`，分支為 `develop`，
起始工作目錄乾淨。依 P0 → P1 → P2 → P3 → P4 → P5 執行。
建構使用 `build/release`、Ninja、AppleClang 21.0.0、arm64、
`Release`（`-O3 -DNDEBUG`），SYSTEM dependencies，
UMAP、ROOT、OpenMP、experimental features、testing 與 Python bindings 均啟用。

所有建構均為 `--parallel 4`；CTest 採 `--parallel 1`，設定
`OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`。
Joint 測試沿用自身的一執行緒設定，Eigen math 測試沿用二執行緒設定。
各建構、測試與探針程序循序執行，未同時執行多個建構／測試。

P0 重新建構成功；完整 CTest 在上述環境為 22/23 組通過。唯一既有失敗是
runtime 文件連向不存在的 `jce-terminology-validation.md`。
該處目前指向本驗收文件。初次以 `OMP_NUM_THREADS=1` 執行時另有既有 Eigen
thread-restoration assertion 失敗；同一個未修改測試在四執行緒環境下通過。
未變更求解器或放寬容差。

## P0 符號與契約盤點

盤點涵蓋 `include/`、`src/`、`tests/`、`docs/`、`resources/examples/`、
工具、JSON fixtures、SQLite、CLI 與 Python bindings；專案沒有獨立 `examples/` 目錄。

| Original | Final | File / consumers | Reason / contract impact |
| --- | --- | --- | --- |
| `JointInitialization::b` | `tau` | `JointEstimationTypes.hpp`、FirstStage、全部初始化呼叫端 | MDPDE 索引 1（第二座標）的 width，單位 Å；同步 JSON |
| `original_b`, `used_b` | `original_tau`, `used_tau` | 核心型別、FirstStage、JSON codec | 分別為原始與實際使用的初始化 τ |
| `JointState::b` | `tau` | 核心型別、assembly、Halo、peeling、stage adapter、exports | `exp(eta)`；欄位位置與型別不變 |
| `JointState::log_b` | `log_tau` | 核心型別與 state consumers | 儲存 η，並非物理 τ |
| `JointState::width_gradient` | `eta_gradient` | 核心型別、Numerics、JSON | 實際為 `dF/deta`，未加入或移除 chain-rule factor |
| `candidate_width_gradient_inf_norm` | `candidate_eta_gradient_inf_norm` | candidate、component、JSON | raw η gradient norm |
| `global_width_gradient_threshold` | `global_eta_gradient_threshold` | candidate、JSON | bound-projected η stationarity threshold |
| `global_width_gradient` | `global_eta_gradient` | `ConvergenceCriteria.hpp` 與求解器／測試 | 門檻數值與 practical-v2 / strict-v2 定義不變 |
| `*_width_gradient*` | `*_eta_gradient*` | BCD sweep、progress、independent replay、diagnostics | 明確區分 raw、projected、stationarity η gradient |
| `local_widths` | `local_tau` | `JointBlockConditioning`、`LocalConditioning.cpp` | 內部 `exp(eta)` 向量；原欄位並未序列化 |
| `bmax4`, `bmax` | `tau_max4`, `tau_max` | `WidthCalibration.hpp/.cpp`、呼叫端 | `4 * Q0.95`，未更改公式 |
| `BoundedInitialization::widths` | `tau_values` | calibration、estimator、tests | 投影後的 τ 向量 |
| `InitialWidthProjection` | `InitialTauProjection` | `WidthCalibration.hpp` | 僅用於 Joint τ 初始化 |
| `original_first_stage_width` | `original_first_stage_tau` | projection record | 原始 FirstStage τ |
| `original_initial_width` | `original_initial_tau` | projection record | 投影前 τ |
| `bounded_initial_width` | `bounded_initial_tau` | projection record | 投影後 τ |
| `gaussian_log_width` | `d_gaussian_d_log_tau` | `Numerics.hpp/.cpp`、test basis | `dG/deta` |
| `charge_log_width` | `d_charge_d_log_tau` | `Numerics.hpp/.cpp`、test basis | `dK/deta` |
| `initial_b`, `parent_initial_b`, physical parameter `b` | `initial_tau`, `parent_initial_tau`, `tau` | Search、Components、BCD、estimator、adapter、test helpers | 追蹤為正的物理 width 後更名 |
| physical locals `width`, `widths`, `width_squared` | `tau`, `tau_values`, `tau_squared` | kernels、replay、FirstStage、validation | 保留全部浮點運算順序 |
| `gaussian_width_dot`, `charge_width_dot`, `width_compensation` | `gaussian_eta_dot`, `charge_eta_dot`, `eta_compensation` | `Numerics.cpp` | η derivative reduction |
| `maximum_target_width`, `max_B` | `maximum_target_tau`, `max_tau` | BCD、progress、Logger | 最大 `exp(eta)` |
| `width_upper_bound` | `tau_upper_bound` | BCD failure diagnostics | `exp(eta_upper_bound)` |
| `width_step_inf` | `eta_step_inf` | BCD failure diagnostics | LM step 的 η 座標 norm |
| `width-grad`, `gradient-inf` | `eta-grad` | CLI progress | η gradient norm |
| `width-kkt` | `eta-projected-kkt` | CLI progress | projected η KKT |
| `b_gradient`, `b_gradient_inf` | `eta_gradient`, `eta_gradient_inf` | test diagnostic snapshots | 原資料已是 η gradient |
| JSON `b`, `log_b`, `original_b`, `used_b`, `width_gradient` | `tau`, `log_tau`, `original_tau`, `used_tau`, `eta_gradient` | `JointResultJson.cpp` | schema 24 僅接受新 key |
| units `B`, `joint-kernel-map-units-v1` | `tau`, `joint-kernel-map-units-v2` | JSON units | 單位仍為 Å；嚴格拒絕舊 metadata |
| CSV `A,B,C` | `A,tau,C` | `JointAnalysisFileIO.cpp`、CLI consumers | 只變更表頭，欄位順序不變 |
| Joint summary `A/B/C`, `B mean` | `A/tau/C`, `tau mean` | `StageSummary.cpp` | 僅 Joint 模式標籤更新 |
| `joint_abc`, `joint-abc-profile` | `joint_component`, `joint-tau-profile` | test namespace、diagnostic label | 清理測試中的舊三參數記號 |
| test design `abc`, `log_width` | `full_design`, `eta_derivative` | `JointHaloObservationExclusion.*` | 前者是完整設計矩陣 spectrum，後者是 η derivative matrix |

刻意保留：`eta`、`eta_max4`、`eta_upper_bound`、`candidate_eta`、`accepted_eta`；
`WidthCalibration`、`CalculateWidthBoundKkt`、`width_bound_active`、
`active_width_bound_count`、`projected_width_kkt`、projected-width spectra；
`JointHaloMedianReconstruction::width` 及其 JSON `width`；
`GaussianModel3D::GetWidth()` / `SetWidth()`、通用 Width uncertainty 與 local/group API。
這些名稱表達對數座標、物理性質、限制或通用 Gaussian 模型。

`LinearSolve::b` 是 A/C coefficient workspace，`TiledQR::b` 是右側向量；
比較運算的 `b`、顏色的 `b`、A/C regression fixture 的 `b`、CIF B-factor、chain B、
BCD、Bounded、Block、Budget、Build、Byte 及 Benchmark 不屬於原子 τ，均保留。
文件連結中的既有 `joint_abc_*` 檔名作為穩定文件路徑保留，內容採現行 τ 記號。

## P1–P3 數值驗證

`JointState` 保留原始 aggregate 順序 `(ac, tau, log_tau, eta_gradient, objective)`。
FirstStage OLS／MDPDE 陣列與 Gaussian `(A, width, offset)` 的位置不變。
原始資料、median fallback、donor qualification、Hyndman–Fan Type 7、20 donor 最小值、
95% coverage、overflow 判斷、投影與 calibration-unavailable reasons 均維持原規則。

固定 calibration fixture 有 25 個 qualified donors、coverage 1，
`p95=0.628`、`tau_max=2.512`、`eta_max=0.9210792486059522`，恰有一個投影。
每階段比較相同輸入與相同建構設定。

數值探針包含 15 組 kernel/derivative samples，chain／cube 各 4 atoms，
完整 design、kernel derivative、profile residual／objective／A/C／η gradient、
projected derivative、Profile Jacobian、每一筆 LM trial 的 step／damping／radius／ratio／
active indices／multipliers，以及兩個 block 的 BCD update sequence、sweep KKT、
accepted values、prediction／residual 與 endpoint / independent replay。
兩組 BCD 都在 5 sweeps 通過 endpoint qualification。

P1、P2、P3、P4 每階段各比較 **100,680 個浮點值與 313 個其他 scalar**，
浮點以 IEEE-754 binary64 bit patterns 比較，全部差異為 0。
44 個數值與下游 source/header 在移除註解、套用已核准名稱對照後，
運算 token 完全相同；沒有改變任何 arithmetic、常數、solver、迴圈或接受條件。
Stage summary 的額外條件只選擇 Joint 顯示標籤。

## P4 持久化與輸出契約

| Contract | Version |
| --- | --- |
| Joint Result JSON | **24** |
| Units metadata | **joint-kernel-map-units-v2** |
| Stage Result JSON | **4，保留** |
| SQLite schema | **20，保留** |
| Solver provenance | **joint-solver-provenance-v6，保留** |

Reader、writer、validator、required-key lists、fixtures、snapshots 與 consumers 同步採用
新 key。舊 schema、舊 key、混合 key 與舊 units contract 皆明確拒絕；沒有 alias、wrapper、
雙名稱輸出、migration、fallback 或舊 reader。
SQLite 的 `result_json` 仍由共同 strict Joint codec 解碼，不變更資料表或資料庫版本。
拒絕 schema 23 的 SQLite load 測試同時驗證資料庫位元組不被改寫。
Python command API 沒有公開需更名的 JCE `b` binding；整合測試同步驗證 schema 24。

## P5 驗證結果

完整 Release CTest **23/23 組通過、0 失敗，40.02 秒**。
註冊檢查確認 **1,063 個啟用的 GoogleTest 全部覆蓋、0 遺漏**。
P1 的五組相關測試、P2/P3 各三組 Joint 測試、P4 的七組契約／整合測試亦均通過。

| 驗證項目 | 結果 |
| --- | --- |
| Numerics、WidthCalibration、ActiveSetBoundedLM、Search | 通過 |
| GaussSeidelBcd、LocalConditioning、Runtime、Progress | 通過 |
| HaloMedianReconstruction、TargetEvidence、PureHaloSupportExclusion | 通過 |
| 新增物理 τ finite difference / η chain-rule regression | 通過 |
| JSON schema 24 round-trip、舊 key／混合 key／舊 units 拒絕 | 通過 |
| schema 23 拒絕、SQLite 不改寫、Stage JSON v4 regression | 通過 |
| 真實 CLI Joint 流程與 CSV／SQLite exports | 通過 |
| Python bindings runtime 與 examples smoke | 通過 |
| 安裝至本機 build 目錄後，外部 C++ consumer 編譯／執行／JSON／SQLite round-trip | 通過 |
| Markdown link validation、GTest registration coverage | 通過 |
| `git diff --check` | 通過 |

正式 `src/`、`include/` 已無退休的 B 參數符號。
精確名稱掃描僅在 `DataObjectPersistence_test.cpp` 找到四行舊名稱拒絕資料；
candidate 舊 gradient key 與舊 units 同樣僅作負向測試資料。
既有文件檔名與 BCD 驗收文件中的歷史 schema 版本不代表現行可接受的契約。
通用 Width、BCD／Bounded／Block 等字樣及無關的 B 均保留。

本機證據保存在 `build/release/jce-tau/`：
`inventory.md`、完整搜尋清單、`baseline-ctest-threads4.log`、`final-ctest.log`、
各階段建構／測試／比對結果、`probe.cpp`、`compile_probe.py`、`compare.py`、
`baseline-numerics.json`、`final-numerics.json`、baseline/final CLI 的 SQLite 與 CSV、
source audit、安裝與 consumer logs。這些為本機驗證產物，不編入 production。

重現最終數值比較可執行：

```sh
python3 build/release/jce-tau/compare.py final
```

重建數值探針與執行時，沿用上列環境與原建構設定：

```sh
python3 build/release/jce-tau/compile_probe.py
OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 VECLIB_MAXIMUM_THREADS=1 \
  build/release/jce-tau/probe > build/release/jce-tau/final-numerics.json
```

## 數值一致性與邊界

最終比較涵蓋探針及四筆 CLI SQLite record 的 Joint / Stage snapshots，
共 **233,080 個 binary64 浮點值與 321,537 個其他 scalar，非預期差異 0**。
所有 CSV 在僅正規化 `A,B,C` → `A,tau,C` 表頭後也完全相同。
比較只正規化核准的 key／schema／units 變更；排除耗時、輸入檔案位置、build/source
fingerprint 與兩次執行各自產生的 run UUID。未對任何科學數值套用容差。
輸入 model/map 的內容雜湊、科學值、狀態、qualification 與離散更新結果皆參與比較。

| 代表值 | P0 與最終值（相同） |
| --- | --- |
| chain 首原子 A | `2.000023927212772` |
| chain 首原子 C | `0.19999802701496863` |
| chain 首原子 τ | `0.5000016492454843` |
| chain 首原子 η | `-0.6931438820744168` |
| chain 最終 sweep objective | `1.788801924022418e-10` |
| chain 最終 A/C KKT | `1.9239990156627584e-09` |
| chain 最終 raw η gradient norm | `6.225088632265865e-11` |
| chain 最終 projected η KKT | `6.22508711245473e-11` |
| cube 最終 sweep objective | `1.3161114143953504e-10` |
| CLI 25 Target objective | `8.880905620094166e-17` |

44 個相關數值 source/header 的運算 token 稽核通過；
另外對全部 **27 個有變更的 production C++ source/header** 逐一與初始 HEAD 比較，
扣除核准更名、註解、JSON/version 與 Joint 顯示標籤後沒有差異。
所有既有 `EXPECT_NEAR` / `ASSERT_NEAR` 的變動行保留原容差；新增測試沒有放寬舊測試。
Halo fixture 的幾何 donor ordering、kNN-4 與 medians 之既有精確斷言，
以及 peeling／charge correction／Local MDPDE／Group 結果與完整 Stage snapshots 均通過。

**在已驗證輸入與固定建構設定下，沒有科學數值結果差異；無待修復的數值差異。**

科學運算的核准變更僅為名稱。未引入第三方依賴，未變更 kernel、normalization、
sampling、Pure Halo-Support Exclusion、Constrained A/C QR、Variable Projection、
Profile Jacobian、Bounded Active-set LM、trust region、damping、KKT、replay、certification、
Halo donor/kNN-4/median、peeling、charge correction、Local MDPDE 或 Group Fitting。
通用 Width API 與 Halo `width` persistence 仍保留。

本次直接在既有 `develop` 修改與建立本機 commits，未建立 branch／worktree，
未使用 detached HEAD、force push 或破壞性 reset。未 push，未新增依賴，
未讀取或更新 Notion。起始沒有使用者未提交修改。
