# Joint Component Estimator — Halo-Aware Observation Exclusion Feasibility and Controlled Experiments

| 日期 | 分支 | 資料 |
| --- | --- | --- |
| 2026-10-16 | develop | 6Z6U chain A；non-hydrogen targets |

結果文件：[per-target coverage CSV](joint-halo-observation-exclusion-per-target-coverage.csv)；[sidecar accuracy JSON](joint-halo-observation-exclusion-metrics.json)

## 1. Executive Summary

**結論：Outcome E — 證據不足以升格 production。** 本研究證明 geometry mask 可以和既有 FixedNeighbor estimator 相容，A/B/C 七個完整擬合均通過原有 practical convergence 與獨立 replay。對這份幾乎無雜訊、由 sidecar 可重建的單一 simulation map，R=2.0 Å 的 Pure Exclusion B 將 A/B/C 參數 RMSE 相對 production baseline 降低 34.8–39.3%；R=2.5 Å 降至 10⁻⁶ 至 10⁻⁷ 級。

R=2.5 Å 移除了 49,305 / 229,783 個 observations（21.46%）。P1 沒有發現 block 或 target rank deficiency、ill-conditioning 或 component fragmentation；但有 34 個 target 剩不到一半原本 rows，且兩個 target 失去最外層 radial shell。P2 的擾動式 stability assessment 對多數 targets 未評估，單次 deterministic map 也不能估計 sampling variance，因此這些結果不能支持 production 推廣。

另有一項重要預測限制：R=2.5 Å 時，保留的 observation domain 內 RMSE 約 9.21×10⁻⁷ normalized map units，但排除域 RMSE 為 0.793；使用本次 fitting 實際保留的 Halo 設定，在完整原始域的 RMSE 為 0.367。該半徑下 C 組沒有任何保留 Halo，因而等同 B。精確的 Target A/B/C recovery 不代表完整物理域可由沒有 Halo 訊號的模型預測。

目前不建議改 production default，也不建議宣稱 Soft Weighting 已有優勢。R=2.0–2.5 Å 可作為後續含噪聲重複實驗的候選範圍。

## 2. Implementation Audit

### Original domain、contributors 與 blocks

Production 的 BuildJointProblem 先建立 selected targets 的 2.5 Å support union；每個 observation row 保留原 map voxel index 作為 identity。Targets 可以共用 rows：本例 178,360 / 229,783 rows 至少由兩個 targets 支持，最多一列由 11 個 targets 支持。Original target graph 是一個 component。

Halo 是 target support union 中有非空 contributor support 的 eligible non-target；本例 659 個 supported halos。Halo identification 以這個原始 domain 為準。Production 先建立完整問題與 FirstStage initialization，再建立 fixed-halo prediction、從 observations 扣除該 prediction，並以 target-only contributors 建立 FixedNeighbor blocks。Halo snapshot 使用 FirstStage A/B、C=0，整段 Joint fit 中固定。Objective、全域 practical convergence、local conditioning、parameter stability 與 independent replay 都沿用同一個 Joint problem row domain 和 row mapping。

### Mask integration 與 normalization

實驗在 original Joint problem 和 FirstStage initialization 完成後，依每個原 row 到最近 supported Halo 的幾何距離產生 mask；不依 FirstStage 參數、residual、current target estimate 或 Oracle truth 選 row。研究 helper 只 compactly remap 保留 rows，並裁切既有 sparse support/incidence；不建 global dense design 或 observation-by-halo distance matrix。若 Halo 在 mask 後沒有 target-row support，就從該 masked problem 的 contributor metadata 移除。C 組只對仍有保留 support 且 FirstStage-qualified 的 Halo 固定 A/B、C=0。

B/C 使用相同半徑、相同 row identities（每個半徑 row-ID digest 相同），相同 target initial widths、12-atom block core、sweeps、numerical gates 與 replay policy。A/B/C 全部固定沿用 original observation scale：

- 原始 observation count：229,783
- s = 787.5806867375935
- map normalization divisor：0.44031545732594374
- production practical thresholds：global AC KKT 1e-7、log-width gradient 1e-9、confirmed log-width change 1e-6

R=0 不排除任何 row；A 組直接呼叫正式 initialized fit，B/C 的 R=0 contract 亦由 deterministic controls 覆蓋。新增邏輯位於 research executable/helper 與 detail test access，不改 production estimator 預設行為，也沒有新增 CLI option。

## 3. Geometry and Coverage Results

正式半徑在讀取 simulation truth 前依 P1 geometry/identifiability 固定為 0、0.5、1.0、1.5、2.0、2.5 Å。Observation domain 是 production 的 target support union。

| 排除半徑 | 保留 rows | 保留比例 | 每 target 保留比例：min / p10 / median | Targets <50% rows | 零 row targets |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0 Å | 229,783 | 100.0000% | 100.00 / 100.00 / 100.00% | 0 | 0 |
| 0.5 Å | 229,669 | 99.9504% | 98.83 / 100.00 / 100.00% | 0 | 0 |
| 1.0 Å | 228,172 | 99.2989% | 91.81 / 99.04 / 100.00% | 0 | 0 |
| 1.5 Å | 222,078 | 96.6468% | 74.71 / 95.01 / 100.00% | 0 | 0 |
| 2.0 Å | 205,895 | 89.6041% | 46.60 / 85.47 / 97.71% | 3 | 0 |
| 2.5 Å | 180,478 | 78.5428% | 17.32 / 67.36 / 90.68% | 34 | 0 |

每一列被排除的 row 原本同時由至少一個 target 與一個 Halo 支持；R=2.5 的 49,305 個排除 rows 正好等於全部 target/Halo mixed rows。保留 rows 中由多個 targets 共享者仍有 152,184 列。各半徑 component 數維持 1，沒有 observation hole 導致 component fragmentation。

Target-Halo 最近距離分布（<1、1–1.5、1.5–2、2–2.5、≥2.5 Å）為 0、0、1、7、1,410；沒有 target 距 Halo 小於 1.5 Å。Observation 最近 Halo 距離分布（<0.5、0.5–1、1–1.5、1.5–2、2–2.5、≥2.5 Å）為 114、1,497、6,094、16,183、25,417、180,478 列。

每個 target 的 radial shell membership retained 比例如下；各列是依既有 target support 半徑統計的 memberships，不是互斥的 global rows。

| R | 0–0.5 Å | 0.5–1.0 Å | 1.0–1.5 Å | 1.5–2.0 Å | 2.0–2.5 Å |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0 Å | 100.000% | 100.000% | 100.000% | 100.000% | 100.000% |
| 0.5 Å | 100.000% | 100.000% | 100.000% | 99.998% | 99.967% |
| 1.0 Å | 100.000% | 100.000% | 99.994% | 99.903% | 99.567% |
| 1.5 Å | 100.000% | 99.981% | 99.777% | 99.069% | 97.857% |
| 2.0 Å | 99.917% | 99.488% | 98.158% | 96.082% | 92.267% |
| 2.5 Å | 98.093% | 96.215% | 93.185% | 87.728% | 83.075% |

因此 2.5 Å mask 沒有移除大多數 target 中心 rows：第一個 shell 仍保留 98.09%。2.5 Å 下 1,416 targets 保有全部五個 shells，2 個少一個 shell；所有 1,418 targets 都至少保有一列。

## 4. Local Identifiability

P1 使用既有 FixedNeighbor core=12 的 local designs；沒有建立 global dense matrix。最大 local row 數由 R=0 時的 4,667 降至 R=2.5 時的 4,132；對 12-atom block，AC/ABC/projected-width columns 分別為 24/36/12。

每個半徑有 142 個預定 blocks；其中 3 個因 target FirstStage seed 無效而 unavailable，不是 mask 造成。其餘 139 個可分析 blocks 在 AC、ABC 與 projected width sensitivity 都 full rank；所有半徑的 rank-deficient blocks、ill-conditioned blocks、rank-deficient targets、不可識別 singleton targets 均為 0。3 個 seed-unqualified target IDs 為 1264、1871、2552，應與 mask-induced loss 分開看待。

| R | 可分析 blocks | max κ(AC) | max κ(ABC) | max κ(projected width) | rank deficient / ill-conditioned |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0 Å | 139/142 | 7.033 | 9.214 | 1.255 | 0 / 0 |
| 0.5 Å | 139/142 | 7.033 | 9.215 | 1.255 | 0 / 0 |
| 1.0 Å | 139/142 | 7.053 | 9.221 | 1.255 | 0 / 0 |
| 1.5 Å | 139/142 | 7.204 | 9.246 | 1.255 | 0 / 0 |
| 2.0 Å | 139/142 | 7.718 | 9.347 | 1.257 | 0 / 0 |
| 2.5 Å | 139/142 | 10.693 | 12.160 | 1.324 | 0 / 0 |

所有 target singleton 的固定寬度 AC、完整 ABC 與 projected-width 設計均 full rank。R=2.5 雖有 34 個 targets 保留不到一半 rows，block 數值條件仍溫和；觀察到的是資料量減少，不是本例的 rank/conditioning failure。

## 5. Controlled Joint Results

A 是正式 production baseline。B 是純 geometry exclusion，不作 Halo subtraction。C 是相同 mask 加上剩餘 supported Halo 的 FirstStage A/B、C=0 correction。三個正式 mask 半徑的 B/C row hashes 完全一致；原始 observation scale固定不變。

| 組別 | R | rows | C 使用 Halo 數 | Target RMSE A / B / C | 秒 | sweeps |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A production | 0 | 229,783 | 659 | 0.509995 / 0.0116784 / 0.0659125 | 346.3 | 79 |
| B pure exclusion | 1.0 | 228,172 | 0 | 0.489748 / 0.0112578 / 0.0629519 | 340.0 | 79 |
| C exclusion + correction | 1.0 | 228,172 | 657 | 0.513016 / 0.0117284 / 0.0667129 | 336.7 | 79 |
| B pure exclusion | 2.0 | 205,895 | 0 | 0.327261 / 0.00760917 / 0.0400226 | 333.4 | 81 |
| C exclusion + correction | 2.0 | 205,895 | 606 | 0.327513 / 0.00761517 / 0.0400450 | 337.7 | 81 |
| B pure exclusion | 2.5 | 180,478 | 0 | 4.18684e-6 / 1.46484e-7 / 1.28620e-6 | 303.5 | 83 |
| C exclusion + correction | 2.5 | 180,478 | 0 | 4.18684e-6 / 1.46484e-7 / 1.28620e-6 | 308.1 | 83 |

RMSE 順序為 A（raw map units）、B（Å）、C（raw map units）。七個 fits 的 runtime convergence、target runtime convergence、數值 qualification、獨立 replay 全部 Passed；每案 1,418 個 target parameter validity 都 Passed。A 的 objective 為 0.004102082117643847，和已記錄 production practical run 0.0041020821176438495 相差約 2.6e-18。不同 mask 的 objective 使用不同 rows，不作跨域優劣判定。

R=1 時 B 的三項 RMSE 比 A 改善 3.6–4.5%；C 三項略差於 A。R=2 時 B 改善 34.8–39.3%，C 幾乎相同但略差。R=2.5 時 B/C 完全相同，因 C 沒有 Halo 在保留 rows 上還有 support。C correction 對完整域 prediction 有幫助（例如 R=2 時），但沒有改善目標參數 RMSE。

各 B/C retained row set 的 FNV-1a 64-bit identity digest：R=1.0 為 a1e1e6793e14d67b，R=2.0 為 a8b6d81e5aff5114，R=2.5 為 71e231a2313d94aa。A 原始 row set 為 17ee67ce8da3432c。

## 6. Oracle Controls

P4 的 D 是 simulation-only full-domain Oracle Halo control，不是 production fitting或mask selection的一部分。重用了既有 full-oracle-ab-c practical run；與本研究相同的 model/map checksums、normalization divisor、229,783 original rows、1,418 targets、659 halos、12-atom core、forward block order、max 200 sweeps，以及 unchanged practical thresholds 均相符。

| D Oracle Halo control | 結果 |
| --- | ---: |
| halos | truth A/B/C，659 個 |
| 收斂 | 79 sweeps，converged-and-verified |
| objective | 7.956799179355278e-14 |
| target A / B / C RMSE | 3.27291e-6 / 1.15828e-7 Å / 8.20093e-7 |
| target / full prediction RMSE | 6.54480e-7 / 6.55421e-7 normalized map units |
| independent 與 full-observation replay | Passed / Passed |
| time / process peak RSS | 343.343 s / 665,600,000 bytes |

Oracle control 顯示原有 Joint solver 在 halo contribution 正確時可精確恢復 target；production baseline 的主要差距與固定 Halo 模型誤差相容。Oracle 只在 P4 truth comparison 中使用，沒有進入 mask 選擇、initialization 或 A/B/C fitting。

## 7. Statistical Accuracy

Simulation sidecar 的 model/map SHA-256 分別為 349ca6da5546b7d8bae53e6edaca9a65985b69f8e6c24316378c4ea04404a4aa 與 6d7771a9fd584e562550da473d624307b22d20396f11a2a79a064bff2551682a；sidecar SHA-256 為 33a3c221ec3f2fff18dcd2a0ea6484f68334971fc688f7cd100e220ad5848931。由 sidecar 重建 normalized map 的 RMSE 為 4.16420e-8 map units。Accuracy 使用 1,418 targets；A/C 誤差轉回 raw map units，B 為 Å。每格依序為 RMSE、median absolute error、95th percentile absolute error、maximum absolute error。

### A parameter error

| Fit | RMSE | Median | P95 | Max |
| --- | ---: | ---: | ---: | ---: |
| A baseline | 0.509995 | 0.0313734 | 0.983242 | 6.55308 |
| B 1.0 Å | 0.489748 | 0.0419823 | 0.946379 | 6.28886 |
| C 1.0 Å | 0.513016 | 0.0279610 | 0.994129 | 6.87191 |
| B 2.0 Å | 0.327261 | 0.0161121 | 0.591376 | 5.32004 |
| C 2.0 Å | 0.327513 | 0.0161321 | 0.591384 | 5.32105 |
| B/C 2.5 Å | 4.18684e-6 | 1.96432e-7 | 4.67919e-6 | 5.79475e-5 |

### B parameter error（Å）

| Fit | RMSE | Median | P95 | Max |
| --- | ---: | ---: | ---: | ---: |
| A baseline | 0.0116784 | 0.000807360 | 0.0203483 | 0.139283 |
| B 1.0 Å | 0.0112578 | 0.000924119 | 0.0205790 | 0.134684 |
| C 1.0 Å | 0.0117284 | 0.000682404 | 0.0204035 | 0.145363 |
| B 2.0 Å | 0.00760917 | 0.000421213 | 0.0132779 | 0.114921 |
| C 2.0 Å | 0.00761517 | 0.000422411 | 0.0132788 | 0.114940 |
| B/C 2.5 Å | 1.46484e-7 | 5.74163e-9 | 1.51164e-7 | 2.13308e-6 |

### C parameter error

| Fit | RMSE | Median | P95 | Max |
| --- | ---: | ---: | ---: | ---: |
| A baseline | 0.0659125 | 0.00540476 | 0.131171 | 0.611235 |
| B 1.0 Å | 0.0629519 | 0.00764558 | 0.125482 | 0.577438 |
| C 1.0 Å | 0.0667129 | 0.00456892 | 0.132249 | 0.607730 |
| B 2.0 Å | 0.0400226 | 0.00263705 | 0.0762433 | 0.378995 |
| C 2.0 Å | 0.0400450 | 0.00263777 | 0.0762505 | 0.379071 |
| B/C 2.5 Å | 1.28620e-6 | 3.53829e-8 | 1.17015e-6 | 2.09363e-5 |

### Target-Halo distance strata

預先固定分組為 <2.5 Å、2.5–5 Å、≥5 Å，樣本數分別為 8、1,095、315。以下每格依序是 A/B/C RMSE；這些是分組中不同 targets 的 errors，不能解讀成重複抽樣 variance。

| Fit | <2.5 Å A / B / C | 2.5–5 Å A / B / C | ≥5 Å A / B / C |
| --- | --- | --- | --- |
| A baseline | 1.16074 / 0.0190940 / 0.160241 | 0.571010 / 0.0131620 / 0.0735736 | 0.0565725 / 0.00157391 / 0.00937880 |
| B 1.0 Å | 0.943119 / 0.0151829 / 0.154598 | 0.550679 / 0.0127195 / 0.0702442 | 0.0546328 / 0.00150428 / 0.00895154 |
| C 1.0 Å | 1.18020 / 0.0197390 / 0.166559 | 0.574232 / 0.0132131 / 0.0744125 | 0.0559551 / 0.00155594 / 0.00904583 |
| B 2.0 Å | 0.576914 / 0.0113693 / 0.0783521 | 0.368644 / 0.00858920 / 0.0449365 | 0.0354594 / 0.000950288 / 0.00594384 |
| C 2.0 Å | 0.577233 / 0.0113746 / 0.0783995 | 0.368930 / 0.00859603 / 0.0449617 | 0.0354657 / 0.000950457 / 0.00594507 |
| B/C 2.5 Å | 4.79894e-7 / 2.37517e-8 / 1.19398e-7 | 4.35318e-6 / 1.52075e-7 / 1.35590e-6 | 3.60980e-6 / 1.27222e-7 / 1.02756e-6 |

B 組 R=1 的改善在最靠近 Halo 的 8 個 targets 較明顯，遠組僅小幅改善；R=2 的三個距離組都改善約三至五成，R=2.5 則三組皆接近 machine precision。改善並非只出現在 Halo-near targets，尤其 R=2 與 2.5；這限制了將整體差異單獨歸因於 Halo proximity 的解讀。C 組在 R=1 的三個 strata 都略差於 baseline。

### Prediction error domains

以下均為 normalized map units。Actual model 欄使用各 fit 實際採用的 Halo policy；B 沒有 Halo correction，C 僅含保留 support 的 fixed halos。Oracle-assisted 欄只用於隔離 target prediction error。不同域樣本數不同，retained RMSE 不能單獨代表整體估計較準。

| Fit | Target-only / original | Actual / retained | Actual / excluded | Actual / original full domain | Oracle-Halo-assisted / original |
| --- | ---: | ---: | ---: | ---: | ---: |
| A baseline | 0.0984444 | 0.148817 | not applicable | 0.148817 | 0.0984444 |
| B 1.0 Å | 0.0943982 | 0.145325 | 3.90915 | 0.357923 | 0.0943982 |
| C 1.0 Å | 0.0989927 | 0.142971 | 0.522546 | 0.149036 | 0.0989927 |
| B 2.0 Å | 0.0638622 | 0.103604 | 1.07169 | 0.359188 | 0.0638622 |
| C 2.0 Å | 0.0638980 | 0.103637 | 0.368354 | 0.154044 | 0.0638980 |
| B/C 2.5 Å | 1.12666e-6 | 9.20736e-7 | 0.792632 | 0.367163 | 1.12719e-6 |

## 8. Numerical Reliability and Failure Categories

七個 A/B/C runs 的 global runtime convergence 與 target runtime convergence 均 Passed；所有 qualified FirstStage halos 與 target parameters valid。Independent replay 和 row identity mapping 皆 Passed。local conditioning 每個 target 都是 Adequate；P2 R≤2 的 142/142 blocks、R=2.5 的 138/138 blocks 均為 WellConditioned。

| Diagnostic | A / B,C 1 Å / B,C 2 Å / B,C 2.5 Å |
| --- | --- |
| 最大 AC KKT | 9.20e-8 / 9.33e-8 / 9.54e-8 / 9.44e-8（門檻 1e-7） |
| 最大 width gradient | 4.28e-10 / 4.28e-10 / 4.38e-10 / 8.77e-10（門檻 1e-9） |
| 最大最近三個 sweeps 的 log-width change | 7.44e-7 / 7.52e-7 / 7.11e-7 / 7.10e-7（門檻 1e-6） |
| 最大 scaled A/C trajectory change | 3.75e-6 / 3.79e-6 / 4.04e-6 / 4.28e-6 |

P2 的 trajectory-stability summary 有實際觀測；擾動式 stability check 則沒有完整執行。A、B/C 1 Å、B/C 2 Å 的 perturbation status 為 Unavailable，0 perturbation directions；各有 1,406 個 target stability NotAssessed、12 個 Unavailable。R=2.5 為 1 個 block Stable、137 個 blocks NotAssessed；1,406 個 targets NotAssessed、12 個 targets Stable。Point quality 在 R≤2 為 1,406 NotAssessed、12 Warning；R=2.5 為 1,406 NotAssessed、12 Available。故 run-level convergence 通過不等於所有 target 的 perturbation stability 已證明。

失敗類別：A/B/C 沒有 runtime convergence、replay、rank 或 conditioning failure；P1 的 3 個 unavailable blocks 來自原本 3 個 target FirstStage seed failures。單次 map 的主要限制是 statistical variance 與 per-target perturbation stability 未得到估計，不是已觀察到的 local rank failure。

## 9. Computational Resources and Reproduction

P1 幾何/identifiability 共用 nearest-Halo distances；單執行緒 runtime 13.394 s，process peak RSS 486,408,192 bytes。P2 七個 fits 各自執行，單案 303.527–346.283 s，process peak RSS 574,078,976 bytes；最多 229,783 observations、1 component、12-atom core。Oracle reuse 為 343.343 s、665,600,000 bytes。這些峰值均低於 4 GB；未建立大型 dense global design 或 observation-by-halo matrix。數值 experiments 設定 Eigen/OpenMP/OpenBLAS/Accelerate 各 1 thread。

重現命令（依序執行；output files 可放在 /private/tmp）：

    cmake --build build --parallel 4
    ctest --test-dir build --output-on-failure -j 4
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 ./build/bin/joint_halo_observation_exclusion /Users/yslian/data/6z6u.cif /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map /private/tmp/feasibility.json
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 ./build/bin/joint_halo_observation_exclusion_fit /Users/yslian/data/6z6u.cif /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map /private/tmp/p2-fits.json
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 ./build/bin/joint_halo_observation_exclusion_metrics /Users/yslian/data/6z6u.cif /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map /Users/yslian/Documents/simulation/sim_map_gaus_grid0.50_charge1_width_6Z6U_bw0.50.map.simulation.json /private/tmp/p2-fits.json /private/tmp/metrics.json

P1 final 與 repeat 的 coverage/identifiability records 完全相同（忽略 elapsed time 與 RSS）。P2 初次結果與 reliability diagnostic rerun 的七組 target A/B/C values 和 objective 全部逐值相同，最大參數差為 0。P4 metrics 第二次輸出與第一次 JSON 完全相同。

全 Release build 通過，使用 parallel 4。完整 CTest 42 項中 41 項通過，-j 4；joint component contract、workflow、numerical、runtime/regression、physical smoke、offline、mapping-memory scalability 與新增 deterministic controls 均通過。新增的 experiment sources 已加入 cleanup contract 白名單。唯一失敗為既有 joint_statistical_experiment_smoke：tests/integration/joint_statistical_experiment.py 第 81 行將 size=2 的 state.ac reshape 為 2×2 而失敗；單獨重跑可重現，與本次 source 無關，未修改該程式。Release build 有 5 個既有 joint_fixed_width_oracle.cpp conversion/sign-compare warnings。

## 10. Final Recommendation

本研究不修改 production defaults。R=2.0 與 2.5 Å 的 Hard Mask 在此 simulation 上是值得後續研究的 accuracy 候選；R=1.0 效果小，C 組在 target RMSE 上沒有勝過 B。P1 顯示 2.5 Å 雖移除 21.46% rows，卻沒有觀察到 rank/condition loss；因此沒有證據顯示這裡需要 Soft Weighting 才能補救 local identifiability。P3 未執行，亦未實作或建立其 commit。

採 Outcome E 的原因：

1. sidecar truth 可將 map 重建至 4.16e-8 RMSE，這不是含有可估 sampling noise 的重複觀測實驗。跨 targets 的 error spread 不能當成同一 target 的 sampling variance。
2. 多數 targets 的 perturbation stability 未評估；R≤2 的 12 個 warning 明確標為 Unavailable。
3. R=2.5 的 target parameter recovery 雖接近精確，實際 Halo-free full-domain prediction RMSE 仍是 0.367，excluded-domain RMSE 是 0.793。其 retained-domain success 不代表完整域物理預測成功。
4. 2.5 Å 下 34 targets 的 row coverage 低於 50%，且沒有含噪聲 repetitions 可測量這些 targets 的 variance/robustness。

Q6 初步可行性：本次 target-Halo geometry 中沒有 <1.5 Å 的 pairs，僅 8/1,418 targets 的最近 Halo 距離 <2.5 Å。現有 Joint input 與 sidecar 沒有提供可用的 symmetry-orbit/tied-contributor metadata，因此無法由本資料驗證特定 symmetry-tied groups；以對稱性連結 contributors 可能比刪 rows 更直接，但需另行建模，不在本次實作範圍。

建議後續先以固定 seeds 加入多組可控 noise realizations，保留 2.0 與 2.5 Å 兩個 mask；再補足 target/block perturbation stability，並為完整域 prediction 明確評估未保留 Halo 的處理方式。取得該證據前，不升 production、不宣稱已降低 estimator variance。
