# RHBM-GEM — Post-Joint Local Domain Correctness Validation

| 項目 | 結果 |
| --- | --- |
| 日期 | 2026-10-19 |
| 工作分支 | `develop` |
| 程式驗證基準 | `7cc9a35e8ae2e4ad497cebf092c4cdaedd0f1852` |
| support radius | 2.5 Å |
| Joint schema | 17 |
| 最終判定 | **B — Production Functional, Scientific Limitations Remaining** |

本報告記錄 Pure Halo-Support Exclusion 的 Post-Joint Local charge-correction correctness 修正、P0–P4 驗證，以及修正後 6Z6U A/B/C/D 研究結果。歷史研究 artifact 保留不覆寫；本次重跑結果另存於
[`joint-halo-support-exclusion-validation-6z6u-post-joint-local-fix.json`](joint-halo-support-exclusion-validation-6z6u-post-joint-local-fix.json)。

## 1. 結論摘要

已確認的 root cause 是 `BuildPostJointLocalSamples()` 將 Original contributor index 用在 compact Fitting `support` 上，並連帶把 fitting row key 與 Original row lookup 混用。這不是假設性風險：P0 fixture 在修正前產生數值錯誤與 target endpoint failure。

修正只將 Post-Joint Local charge correction 的 support source 統一為 `m_problem.OriginalInput()`。Joint objective、FixedNeighbor solver、Halo median policy、Peeling neighbor subtraction、Local MDPDE、Group Fitting、Independent Replay 與 schema 17 均未改寫。

修正後：

- P0 independent tricubic charge reference 通過；正、負、零 charge control 及 Halo-independence control 通過。
- 6Z6U B/C/D 的 Joint A/B/C、Halo field、Peeling 完全保持原值；Local/Group 明顯改善，但 Oracle D 仍有 Local `A=0.709793`、`B=0.021414` RMSE，故不能宣稱科學精度完全恢復。
- 含 Halo 的 3-target synthetic full workflow 通過：同 GroupKey Halo 使用 GroupKey median，跨組 Halo 使用 Global fallback，local 與 multi-member group summary 均通過明確數值 bounds。
- 最新程式 commit 的 Release build 與 CTest 為 43/43 Passed、0 Failed、0 Skipped。

## 2. Original/Fitting domain contract

本次確認的 index ownership 如下：

| 資料 | index domain | 用途 |
| --- | --- | --- |
| `OriginalInput().atom_ids` | Original contributor | 完整 Target/Halo ordering |
| `OriginalInput().support` | Original contributor / Original row | raw geometry、Peeling、own-charge correction |
| `Input().atom_ids` | compact Fitting contributor | Target-only Joint objective |
| `Input().support` | compact Fitting contributor / Fitting row | exclusion 後的 solver support |
| `FittingToOriginalContributors()` | Fitting → Original | contributor identity conversion |
| `FittingToOriginalRows()` | Fitting → Original | row identity conversion |
| `m_contributor_index_by_id` | Original contributor | `atom_id` lookup |
| `m_row_index_by_grid` | Original row | Original map-grid lookup |
| `m_fitting_row_index_by_grid` | Fitting row | compact fitting lookup |
| `JointAnalysisResult.components` | captured contract | component atoms/rows 另由 snapshot identity validation 驗證 |

`BuildPostFitPeelingSamples()` 與 raw sampling geometry 使用 Original domain。修正前唯一漏掉的是 `BuildPostJointLocalSamples()` 建立 own-charge profile 時仍取 `m_problem.Input()`；因此 Original contributor index 被解讀為 Fitting contributor index，且 `support.row` 也被解讀為另一個 row space。

修正後的資料流為：

```text
Pure Exclusion Joint
        ↓
Qualified Target A/B/C
        ↓
GroupKey / Global Halo Median
        ↓
Original-domain PostFitPeelingSamples
        ↓
Original-domain Target own-charge correction
        ↓
Gaussian-only Local Samples
        ↓
Post-Joint Local MDPDE
        ↓
Group Fitting
```

## 3. P0 — Root cause reproduction

P0 在 production fix 前先加入，commit `8a1b244d`。fixture 包含 Original ordering `H1,T1,T2`，Fitting ordering `T1,T2`，非連續 exclusion rows `{1,4,7}`，Fitting rows `{0,2,3,5,6,8,9}`，contributor map `{1,2}`。T1/T2 的 width 分別為 `.65/1.15 Å`，charge 分別為 `+0.35/-0.45`，另有 `C=0` control。

test 使用獨立手寫 tricubic reference，沒有用 production interpolation helper 生成 expectation。它也改變 Halo charge 並確認 Target own-charge correction 不變。

修正前 focused run 的實際證據：

```text
target T1 fixed_charge_response = 0.02888530633150601
independent expected             = 0.20404523326715163
absolute difference              = 0.17515992693564561
```

Gaussian-only response 也有同樣的差值；T2 另觸發 `Post-Joint samples require a selected Joint target.`。這證明不是只缺少 `optional` 值，而是 contributor/row domain 混用造成的數值與 endpoint 錯誤。

修正後 P0 test 的 assertions 包含：

- `FittingToOriginalContributors()` 與 `FittingToOriginalRows()` 的非 identity mapping；
- 非連續 exclusion rows；
- 正、負、零 Target charge；
- T1/T2 不同 width；
- independent `fixed_charge_response` 與 `gaussian_response`，tolerance `1e-12`；
- Halo charge mutation 不影響 Target own-charge correction；
- no-Halo、Fixed Halo、boundary、invalid identity 與 schema roundtrip regressions。

## 4. P1 — 最小 production fix

修正 commit `a9684293`，檔案 `src/core/detail/PostFitPeeling.cpp`：

```diff
- const auto & input = m_problem.Input();
+ const auto & input = m_problem.OriginalInput();
```

修正後 `atom_index`、`input.support[atom_index]`、`support.row`、`m_row_index_by_grid` 都屬於 Original domain。own-charge 仍只取該 Target 的 Joint stage `C/B`，沒有使用 Halo median、FirstStage、Group estimate 或 Oracle truth。

沒有修改：

- Pure Exclusion `R=2.5 Å` 與 Target-only Joint objective；
- FixedNeighbor block-coordinate solver 與 global practical criteria；
- Halo GroupKey median / Global Target Median fallback；
- Peeling neighbor subtraction、fitting row mask、Independent Replay；
- Local MDPDE 或 Group Fitting statistical objective；
- Fixed Halo historical API 與 schema 17。

## 5. P2 — 6Z6U scientific revalidation

### 5.1 固定條件

使用既有 `tests/experiments/joint_halo_support_exclusion_validation.cpp`，輸入為：

- `/Users/yslian/data/6z6u.cif`；
- `/Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map`；
- 對應 simulation truth sidecar；
- normalization divisor `0.44031545732594374`，truth scale `2.2710990117699854`；
- Eigen/OpenMP/OpenBLAS/Accelerate 均 1 thread，workflow `thread_size=1`；
- B/C/D 共用同一個 Pure Exclusion Joint endpoint 與同一組 Target estimates，只替換 Halo reconstruction。

資料規模與 audit：

| 指標 | 值 |
| --- | ---: |
| Original observations | 229,783 |
| Original contributors | 2,077 |
| Targets | 1,418 |
| retained Halos | 659 |
| Joint groups | 173 |
| Same GroupKey Halo | 453 |
| Global fallback Halo | 206 |
| fallback support rows | 40,286 |
| raw local/peeling samples | 283,600 |
| available samples after fix | 262,260 |
| unavailable samples | 21,340 |
| original-support boundary samples | 21,340 |
| mapping mismatch / missing contributor / other failure | 0 / 0 / 0 |

Joint endpoint 在四組中一致：global/target runtime convergence 為 Passed，Independent Replay 為 `true`，Target A/B/C RMSE 為 `9.508723e-6 / 1.464839e-7 / 2.921096e-6`。修正只發生在 Post-Joint Local 之後，因此這個不變性是必要的 regression。

### 5.2 Before/After：Joint、Halo、Peeling

表中 `before → after`，before 來自保留的 historical artifact `joint-halo-support-exclusion-validation-6z6u.json`，after 來自本次新 artifact。各欄均為 RMSE；Peeling 另列 bias。

| 組別 | Joint A/B/C RMSE | Halo field RMSE | Peeling RMSE / bias |
| --- | --- | ---: | ---: |
| A Fixed Halo | `1.158249/0.0116784/0.149694 → same` | `0.417081 → same` | `0.139140/0.004615 → same` |
| B GroupKey + fallback | `9.50872e-6/1.46484e-7/2.92110e-6 → same` | `0.368676 → same` | `0.0396845/0.00714560 → same` |
| C Global median | `9.50872e-6/1.46484e-7/2.92110e-6 → same` | `0.522220 → same` | `0.0538315/0.00463867 → same` |
| D Oracle Halo | `9.50872e-6/1.46484e-7/2.92110e-6 → same` | `0 → same` | `3.02584e-6/8.37847e-9 → same` |

B/C/D 的 Halo field 與 Peeling 不因 own-charge correction 改變，符合修正只影響 `Gaussian-only local response` 的預期。

### 5.3 Before/After：Local MDPDE 與 Group Fitting

Group `member` 是每個有效 local member 的 truth error；`group mean` 是 Group output 對 member truth mean 的 error。這裡列出 A/B RMSE 及差值 `Δ = after - before`。

| 組別 | Local A/B before → after（Δ） | Group member A/B before → after（Δ） | Group mean A/B before → after（Δ） |
| --- | --- | --- | --- |
| A Fixed Halo | `1.385687/0.024496 → 1.385687/0.024496` (`0/0`) | `1.377811/0.024373 → 1.374747/0.024355` (`-0.003064/-0.000018`) | `0.888268/0.022343 → 0.842092/0.022066` (`-0.046176/-0.000277`) |
| B GroupKey + fallback | `8.630362/0.058381 → 0.786846/0.021981` (`-7.843516/-0.036400`) | `8.669362/0.058478 → 0.765446/0.021821` (`-7.903916/-0.036657`) | `8.529955/0.056782 → 0.618822/0.020783` (`-7.911133/-0.035999`) |
| C Global median | `8.624608/0.058349 → 0.783207/0.021912` (`-7.841401/-0.036437`) | `8.658182/0.058304 → 0.759989/0.021736` (`-7.898193/-0.036568`) | `8.536495/0.056875 → 0.610129/0.020702` (`-7.926365/-0.036173`) |
| D Oracle Halo | `8.650290/0.058503 → 0.709793/0.021414` (`-7.940497/-0.037089`) | `8.690001/0.058604 → 0.685234/0.021243` (`-8.004768/-0.037361`) | `8.547503/0.056896 → 0.562252/0.020331` (`-7.985251/-0.036565`) |

修正後 B/C/D 的 Local effective sample count 都為 79,974，Local status 都是 `Available: 1,418`；Group 結果為 173 groups、147 `Available`、26 `Unavailable`，其中 21 `single-member`、5 `singular-group-covariance`。這些 unavailable group 是既有統計資格結果，沒有因本次修正增加。

### 5.4 GroupKey Median versus Global Median

在 Halo field 層級，B 的 GroupKey policy 明顯優於只用 Global median 的 C：

| 指標 | B GroupKey + fallback | C Global only |
| --- | ---: | ---: |
| aggregate Halo RMSE | 0.368676 | 0.522220 |
| same-group Halo RMSE | 0.00975849 | 不適用 |
| global-fallback Halo RMSE | 0.480629 | 0.522220 |
| Peeling RMSE | 0.0396845 | 0.0538315 |
| 修正後 Local A/B RMSE | 0.786846 / 0.021981 | 0.783207 / 0.021912 |
| 修正後 Group mean A/B RMSE | 0.618822 / 0.020783 | 0.610129 / 0.020702 |

Halo reconstruction policy 仍只影響 Halo field 與後續 neighbor subtraction；B/C 的 Target Joint endpoint 完全相同。Global-only 在此資料上 Halo/Peeling 較差，但 Local/Group 差距相對小，說明 Local residual 不是單由 median policy 決定。

### 5.5 Oracle isolation 與剩餘科學限制

D 使用 exact/oracle Halo reconstruction，Halo field RMSE 為 0，Peeling RMSE 為 `3.025840e-6`，但修正後仍有：

```text
Local MDPDE A/B RMSE  = 0.7097926369 / 0.0214136341
Group mean A/B RMSE   = 0.5622518251 / 0.0203311722
```

B/C 相對 D 的 post-local output deviation 分別為：

```text
B versus D: A/B RMSE = 0.2248122206 / 0.0020053566
C versus D: A/B RMSE = 0.2367413813 / 0.0022192200
```

D 對自己的 deviation 為 0 是比較定義，不是宣稱 D 已達 truth。這個 controlled chain 已排除：

- B/C/D normalization 不一致；
- Target Joint estimates 不一致；
- Oracle Halo field error；
- P0 所示的 Target own-charge domain/index error；
- Local/Group 結果只是 unavailable 或 missing 的假象。

現有 Local estimator 仍有明確的 statistical contract：使用 `0–1 Å` local design range、只保留正的 adjusted response、對 response 取 log，再執行 Local MDPDE；Group Fitting 只使用通過 `AssessPostJointLocalResult` 的 local MDPDE outputs。這些規則可在 `PreparedLocalGaussianFit.cpp` 與 `PotentialFittingWorkflow.cpp` 中直接核對。因而本次證據支持「domain correctness 已修正，剩餘偏差屬 local finite-sample/model/statistical limitation」；但沒有足夠證據把約 `0.71/0.021` 的 D residual 歸因到單一統計因素，所以本報告保留 B，而不是 A，也沒有為降低 RMSE 擅自修改 MDPDE 或 Group Fitting。

## 6. P3 — Halo-aware multi-member Group Fitting regression

新增 test `PureHaloSupportExclusionTest.RunsFullPureWorkflowWithMultipleTargetsAndHaloGroupSources`，位於 `tests/core/joint_component/PureHaloSupportExclusion_test.cpp`。fixture 使用正式 production kernel 與 full workflow：

- targets 1/2/3，共用 GroupKey `(1,7)`，真值 `(A,B,C)` 為 `(2.0,.65,-.2)`、`(1.6,.75,.15)`、`(1.2,.85,-.1)`；
- Halo 4 與 targets 同 GroupKey，真值 `(.5,.8,.3)`；
- Halo 5 為不同 GroupKey，真值 `(.3,1.0,-.25)`，因此沒有同組 donors；
- map 中 Halo support 與 target sampling region 重疊，且 charges mixed-sign；
- production workflow 產生同組 `GroupKeyMedian`、跨組 `GlobalTargetMedian`，不是手工替換結果。

測試不是只檢查 `optional` 或 finite：

- target indices 精確為 `{0,1,2}`；
- Halo 4/5 source 與 donor count 都精確為 3；
- 三個 local 結果 `Available`，並以 synthetic truth 檢查 A `±0.25`、B `±0.08`；
- Group summary 精確包含 member IDs `{1,2,3}`，`point_count=eligible_count=3`、`excluded_count=0`、`Available`；
- descriptive group mean 以三個 local outputs 的明確算術平均 `1e-12` 比對；
- group MDPDE A/B 均 finite 且為有效結果。

P0 已提供更嚴格的 independent numerical charge/tricubic assertions；P3 則驗證完整 `Pure Exclusion → Halo median → Peeling → charge correction → Local MDPDE → Group Fitting` production data flow。

## 7. P4 — Production、schema 與 compatibility

Production CLI 以同一 6Z6U chain-A 輸入完成完整 run：

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
build/bin/RHBM-GEM potential_analysis \
  --estimator joint-components --jobs 4 --verbose 3 \
  --folder <output-folder> --database <output-folder>/production.sqlite \
  --model /Users/yslian/data/6z6u.cif \
  --map /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  --select-chain A --save-key post-joint-local-fix
```

實際 run 輸出 `Joint fitting qualified`，保存約 76 MB SQLite，1,418 selected targets、659 retained halos，並通過 production Pure Exclusion default。`EstimateJointComponentsFixedHalo()` 仍是明確 Fixed Halo research entry；schema 17、GroupKey classification、Global fallback、Independent Replay 與 bindings/CLI smoke tests 均在 CTest 中通過。

### RSS

本次 corrected CLI 已完整執行，但 sandbox 不允許 `/usr/bin/time -l` 讀取 `kern.clockrate`，也不允許 `ps` 取樣 process，因此**本次 corrected CLI 的 external peak RSS 未量測**。不能把它寫成新測量值。

repository 既有的同規模 measured reference 是 `986,136,576 bytes = 940.45 MiB = 0.918 GiB`，低於 4 GiB architecture target；這是 reference，不冒充本次 corrected run 的 peak RSS。此次 source diff 只有一行 domain source 替換，沒有建立 `O(N_observations × N_targets)` 額外資料結構。

## 8. Release build / CTest

最新程式 commit `7cc9a35e` 執行：

```text
cmake --build build --parallel 4       Passed
ctest --test-dir build --output-on-failure -j 4
100% tests passed out of 43
0 Failed, 0 Skipped
Total Test time (real) = 82.69 sec
```

Build configuration：

```text
CMAKE_BUILD_TYPE=Release
BUILD_TESTING=ON
CMAKE_CXX_COMPILER=/usr/bin/c++
Apple clang 21.0.0 (clang-2100.3.34.2), arm64-apple-darwin25.6.0
```

## 9. Commits

所有修改直接提交在 `develop`，沒有建立 branch/worktree，也沒有更新 Notion：

| SHA | Commit |
| --- | --- |
| `8a1b244d` | `test(joint): reproduce post-joint local domain index mismatch` |
| `a9684293` | `fix(joint): use original-domain support for post-joint charge correction` |
| `846e98f6` | `experiment(joint): revalidate pure exclusion local and group accuracy` |
| `1451095c` | `test(joint): cover multi-target halo source workflow` |
| `7cc9a35e` | `test(joint): tighten synthetic local accuracy assertions` |

本文件的 docs commit 會在上述 validation commit 之後追加；最終 `develop` HEAD 以交付時的 `git rev-parse HEAD` 為準。

## 10. Gate 判定與 remaining risks

| Gate | 判定 | 證據 |
| --- | --- | --- |
| Indexing correctness | Pass | P0 failing-before / passing-after；Original/Fitting mapping 明確分離 |
| Charge correction correctness | Pass | independent tricubic fixed-charge 與 Gaussian-only `1e-12` |
| Local scientific accuracy | Partial | D Oracle Peeling 幾乎 exact，但 Local A/B 仍 `0.709793/0.021414` |
| Group scientific accuracy | Partial/Pass for data flow | P3 Halo-aware multi-member regression pass；6Z6U truth RMSE 仍有 residual |
| Production regression | Pass | Pure default、Fixed Halo API、schema 17、Replay、CLI/bindings smoke pass |
| Final validation | Pass | Release build、43/43 CTest |

### 最終結論：B — Production Functional, Scientific Limitations Remaining

本次已解除實際存在的 Original/Fitting domain correctness defect，並以 independent numerical regression、6Z6U controlled variants、Oracle isolation 與 Halo-aware multi-member synthetic workflow 證明修正有效。Production data contract 沒有未解決的 index/charge correctness failure。

但 Oracle Halo 下 Local/Group truth error 尚未降至可稱為 fully closed 的程度；現有證據支持它已不再是本次 domain indexing bug 的直接結果，仍應視為 Local MDPDE 的有限樣本／模型限制與後續研究風險。沒有執行或宣稱 37,406-target 完整 fitting，也沒有為了通過測試而放寬任何 qualification threshold。
