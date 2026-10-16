# RHBM-GEM — Pure Halo-Support Exclusion Production Validation

| 項目 | 結果 |
| --- | --- |
| 日期 | 2026-10-17 |
| 工作分支 | `develop` |
| implementation HEAD | `521b2d74`（P5 implementation，文件 commit 另計） |
| support radius | 2.5 Å |
| estimator | Pure Halo-Support Exclusion Target-Only Joint Estimator |
| Joint schema | 17 |

本報告記錄 Pure Halo-Support Exclusion 正式整合、Production default promotion，以及 6Z6U full workflow smoke 的實際結果。未在本次執行中量測的數值明確標示為 `未量測`；歷史研究結果不當作本次 Production CLI 的重新測量。

## 1. Production Architecture

正式資料流為：

```text
Original observation domain
        │
        ├─ Halo support incidence (r² <= 6.25 Å²)
        │          └─ exact row exclusion
        │
        └─ Target-only A/B/C initialization
                   ↓
        Fixed-neighbor block-coordinate profile LM
                   ↓
        Global qualification + independent replay
                   ↓
        Qualified Target A/B/C
             ┌─────┴─────┐
             │           │
       GroupKey medians  Global Target median
             └─────┬─────┘
                   ↓
        Immutable Halo peeling reconstruction
                   ↓
        Original-domain PostFitPeelingSamples
                   ↓
        Post-Joint Local MDPDE → Group Fitting
```

Joint fitting 使用：

\[
V_{fit}=V_{original}\setminus\bigcup_{h\in H}support(h),
\qquad r_{vh}^{2}\le 6.25.
\]

Halo 在 Joint 中只有幾何角色：它決定被排除的 rows，不提供 optimizer parameter，也不需要 Halo FirstStage estimate。Joint 只估計 Targets 的 `A >= 0`、`B > 0`、`C ∈ R`。

Postprocessing 不把 `V_fit` 當成原始 geometry。`V_original` 保留原始 row identity、contributor support、map-grid lookup 與 raw sample geometry；fit rows 另外保存 compact row identity、fit-to-original row mapping 與 available-row mask。Snapshot 對兩個 domain 都做 identity validation。

## 2. GroupKey Contract

Halo 與 Target 的分類直接重用現有 Group Fitting contract：

```cpp
data_internal::GetGroupKey(atom)
```

其底層為：

```cpp
KeyPackerComponentAtomClass::Pack(
    atom->GetComponentKey(),
    atom->GetAtomKey())
```

這不是 Residue＋Element 分組，也沒有建立新的 Halo-specific classification。相同 `ModelObject` 中 GroupKey 的意義與 Group Fitting 完全一致。

## 3. Median Reconstruction

有效 donor 必須來自同一個 Joint run，且符合既有必要 qualification：Joint/global/target convergence passed、point state 存在、parameter validity passed、local conditioning `Adequate`、A/B/C finite、A 非負、B 為正。Parameter stability 的 `NotAssessed` 保留原狀，並不因此把 donor 宣稱為 `Stable`；bounded perturbation 尚未涵蓋的 Target 不會被自動排除。

對每個參數分別使用 deterministic `array_helper::ComputeMedian()`：

\[
\widetilde\theta_h=
\begin{cases}
M_{GetGroupKey(h)},&\text{同 GroupKey 有有效 Target donors},\\
M_{global},&\text{否則}.
\end{cases}
\]

空 donor 集合在呼叫 median 前先檢查；全體沒有有效 Target 時明確拒絕，錯誤為 `halo-median-no-valid-targets`，不產生零值或 FirstStage fallback。完成後的 Halo medians 是 immutable peeling-only snapshot，不會回饋 Joint objective、block update、residual 或 Group Fitting。

Snapshot provenance 使用 `group-key-target-median-v1`，每個 Halo 保存 GroupKey、source、donor count 與 Joint run ID。Pure Exclusion 使用新的 `target-only-pure-halo-support-exclusion-v1` parameterization，不偽裝成舊的 `target-only-fixed-halo-v1`。

## 4. 6Z6U Coverage Results

資料為 chain A、1,418 Targets；本次 production CLI 使用 4-job request，程式回報 Joint initialization/fitting 實際使用一個 worker。

| 指標 | 實測值 |
| --- | ---: |
| Original rows | 229,783 |
| Fitting rows | 180,478 |
| Excluded rows | 49,305 |
| Retained ratio | 78.5428% |
| Candidate supported Halos | 35,988 |
| Retained supported Halos | 659 |
| Zero-support/excluded Halos | 35,329 |
| Valid Target donors | 1,418 |
| Same GroupKey Halo medians | 453 |
| Global fallback Halo medians | 206 |
| Global fallback ratio among retained Halos | 31.2595% |
| Fallback support rows | 40,286 |

`fallback_support_rows` 是每個 fallback Halo 的 support row 集合先去重、再跨 fallback Halos 加總的 audit metric，不是 fallback Halos support union 的 row 數。它量測了 fallback Halo contribution 進入 peeling support 的影響範圍。

本次 6Z6U 全體有效 Target median 為：

```text
A = 13.626595080023293
B =  0.499999989778770
C = -0.501912921753690
```

Target assembled state 的 min / median / max 為：

| 參數 | min | median | max |
| --- | ---: | ---: | ---: |
| A | 13.6264653486 | 13.6265950800 | 36.3375884661 |
| B | 0.3999999123 | 0.4999999898 | 0.5000019311 |
| C | -2.0689712068 | -0.5019129218 | 2.0485313336 |

GroupKey coverage audit 與 fallback mapping 的 deterministic unit tests 也覆蓋單一/偶數/奇數 donor、mixed-sign C、invalid donor、traversal order independence、沒有同組 donor、沒有有效 Target 與 run identity mismatch。

## 5. Scientific Validation

### 5.1 本次 Production CLI

本次 6Z6U CLI 完成了 Joint、peeling、Post-Joint Local 與 Group Fitting 的持久化流程，但該 CLI 輸入沒有附帶可供此結果直接比對的 simulation truth table。因此下列本次 Production RMSE 尚未量測：

* Target A/B/C RMSE；
* Halo field RMSE；
* Peeling RMSE；
* Local MDPDE RMSE；
* Group Fitting parameter RMSE。

這不影響數值 qualification 或資料契約驗證；它只表示本次 production persistence smoke 沒有捏造 scientific-truth accuracy。

### 5.2 已有研究 baseline

先前 matching-support 的 6Z6U Pure Exclusion 研究結果（不是本次 production CLI 重新計算）為 Target A/B/C RMSE `4.19e-6 / 1.46e-7 Å / 1.29e-6`；Oracle FullABC Halo 對照為 `3.27e-6 / 1.16e-7 Å / 8.20e-7`。這些結果支持 solver 與 exact support mask 的研究基準，但不應被解讀成 GroupKey median peeling 的本次 RMSE。

P4 的 synthetic full-workflow fixture 已實際通過 Pure Joint → GroupKey/Global median → original-domain peeling → Local MDPDE → Group Fitting；它驗證資料流與數值契約，不取代有 truth 的 6Z6U scientific run。

## 6. Numerical Validation

6Z6U persisted Joint result：

| qualification | 實測值 |
| --- | ---: |
| Joint runtime convergence | Passed |
| Target runtime convergence | Passed |
| Candidate objective | 1.2333160419894865e-13 |
| Global A/C KKT | 9.4370721605e-8 ≤ 1e-7 |
| Width gradient infinity norm | 8.7748879855e-10 ≤ 1e-9 |
| Confirmed log-width change | 5.5563206691e-7 ≤ 1e-6 |
| Independent replay | Passed / verified |
| Endpoint verification | Passed |
| Numerical qualification | Passed |
| Assembly consistency | Passed; 8.88e-15 ≤ 2e-12 |
| Local conditioning | 138/138 blocks `WellConditioned` |
| Parameter stability | trajectory `observed`; perturbation `NotAssessed` |

未執行的 offline high-precision audits（two-step derivative、Richardson、precision 50/100、boundary audit、regular certificate）仍保存為 `not-run`，沒有被誤標成 Passed。既有 Global Practical Criteria、Local A/C KKT、profile LM、Independent Replay、parameter validity 與 conditioning rejection contract 均保留。

## 7. Workflow Integration Results

本次第三次 6Z6U production run 的 persisted stage result：

| stage | 實測結果 |
| --- | --- |
| Joint result | schema 17；1,418 Target；659 Halo median records |
| PostFitPeelingSamples | 1,418 Target outputs；283,600 raw samples；grid-consistent |
| Original-domain interpolation | fitting mask 不會改寫 original row lookup；identity checks 啟用 |
| Post-Joint Local MDPDE | 1,418/1,418 `Available`；77,380 effective samples |
| Group Fitting | 173 groups；147 `Available`；26 `Unavailable` |

26 個 unavailable group summaries 中，21 個是 `single-member`，5 個是 `singular-group-covariance`；point count 1,418、eligible count 1,418、excluded count 0。這是既有 Group Fitting 對不足統計資訊的結果，不是新增的 Halo reconstruction fallback。

Peeling 有 21,340 個 raw samples 的 tricubic stencil 超出原始 map domain，並記錄為 `outside-joint-domain`；這是 map 邊界的明確 sample-level unavailable state，不是 fitting/reconstruction row identity mismatch，也不是補零或靜默跳過。其餘 samples 依 original-domain geometry 重建；P3 regression tests 另驗證 mask boundary、Target own contribution、Halo Gaussian/charge contribution、missing state 與 identity mismatch rejection。

## 8. Resource Usage

```text
Build/test parallelism: maximum 4
6Z6U CLI job request: 4
Joint search: 83 sweeps, 22,908 profile evaluations
Joint cost fields total: 299.192900417 s
Original observations: 229,783
Retained fitting observations: 180,478
Contributors in fit: 1,418 Targets
Retained supported Halos: 659
```

Persisted numerical resource evidence reports endpoint verification memory `96,849,862` bytes under the `3.2 GB` verification budget；largest recorded local workspace estimate was `106,644,934` bytes。此 workflow 沒有建立 global dense observation-by-atom matrix。

本次 production CLI 沒有以 `/usr/bin/time` 或 equivalent external RSS sampler 包住 process，因此「本次這一個 CLI run 的 peak RSS」為 `未量測`。同一 6Z6U Pure Exclusion research/noise run 曾實測 peak RSS `722,452,480` bytes（約 689 MiB，Eigen/OpenMP/OpenBLAS/Accelerate 各 1 thread）；該數字是 measured reference，不冒充本次 CLI 的 peak RSS。

## 9. Regression and Reproduction

實作 commit 依序為：

1. `d969932b` — `experiment(joint): assess group-key halo median coverage`
2. `643b6d05` — `feat(joint): integrate halo-support exclusion fitting`
3. `6393fa78` — `feat(joint): reconstruct peeling halos from group-key medians`
4. `3de17ddf` — `refactor(joint): separate fitting and peeling domains`
5. `b91c3188` — `test(joint): validate halo-exclusion production workflow`
6. `521b2d74` — `refactor(joint): promote pure halo exclusion production`
7. 本報告 — `docs(research): validate pure halo exclusion production promotion`

主要 reproduction command：

```bash
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure -j 4

build/bin/RHBM-GEM potential_analysis \
  --estimator joint-components --jobs 4 --verbose 3 \
  --folder <output-folder> --database <output-folder>/production.sqlite \
  --model /Users/yslian/data/6z6u.cif \
  --map /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map \
  --select-chain A --save-key production
```

Production log 會輸出 `valid_targets`、`same_group_halos`、`global_fallback_halos` 與 `fallback_support_rows`。Result JSON schema 17 保存 dual-domain mapping、Pure parameterization、median source、donor count 與 Joint run identity；CSV export 將 Halo rows 標為 `HaloMedianReconstructed`。

## 10. Production Status and Remaining Risks

**已正式升格。** `RunPotentialFittingWorkflow()` 的 `JOINT_COMPONENTS` production entry 目前明確使用 Pure Halo-Support Exclusion；Halo FirstStage 不再是 Joint fitting prerequisite。PostFitPeelingSamples 使用同 GroupKey median，無同組 donor 時使用本次同一 Joint run 的 Global Target median；Group Fitting estimation algorithm 未改寫。

舊功能仍保留：

* Fixed-Halo schema 16 結果仍可讀取；
* Fixed-Halo parameterization 與研究入口保留供歷史結果、研究工具與 regression tests 使用；
* legacy detail workflow 的 explicit non-Pure entry 保留，不是 Production default。

已知且被本次決策接受的風險：

* simulation support radius 超過 2.5 Å 時的 finite-support mismatch bias 沒有被補償；
* 6Z6U 本次 Production CLI 的 target/halo/local/group truth-based RMSE 尚未量測；
* 參數穩定性 perturbation coverage 仍是 `NotAssessed`，不能宣稱所有 Target stable；
* 26 個 Group Fitting summaries 因 single-member 或 singular covariance unavailable；
* 本次 CLI peak RSS 未由外部 sampler 測量，只有 persisted memory estimate 與既有 6Z6U reference RSS。

上述項目沒有造成 Pure Exclusion 的 fitting/reconstruction identity mismatch，也沒有要求放寬 numerical criteria；因此不構成此次 Production default promotion 的 correctness blocker。
