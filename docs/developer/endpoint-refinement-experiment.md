# Endpoint refinement：分支、完整 operator 與 recovery 實驗

基準為 `0785a395`；本輪所有數值政策均限於 `BUILD_TESTING`。
Production 預設、公開 C++／CLI 介面、MDPDE 方程、weight floor、covariance 公式、offset IRLS、
外層接受政策與 `1e-4` 收斂門檻均不變。Observation matching 不在本輪範圍。

## 方法與判準

候選從同一 dataset 的原 OLS／固定點路徑終點啟動 Powell hybrid，座標為
`(beta0, log(beta1), log(v))`，沿用中央差分 Jacobian 與 `xtol=1e-12`。
採用前重新檢查有限正 variance、有效 Gaussian、完整 weighted rank、正 denominator，
並要求含 floor 的 fresh scaled equations 最大範數 `<=1e-8`。
Solver 原生停止狀態不等於 qualification。

每次採用前，另從**原終點**繼續原 `weights → beta → variance` 更新，
要求參考解 residual `<=1e-10`，總固定點更新不超過 10,000 次（含原求解迭代）。
三個座標逐項差異除以 `max(1, |u_candidate|, |u_reference|)` 後均不得超過 `1e-6`；
weights 最大絕對差異不得超過 `1e-6`，逐樣本 floor activation 必須一致。
這是數值分支一致性判準，不是同分支的數學證明。已知第二根由此判準拒絕。

接受時一併更新 beta、variance、fresh weights 及既有 covariance；保留 OLS 與原迭代診斷。
拒絕時保留原終點數值及其 weights／covariance，但有效 status 為未合格；原 status 另存。
參考解不作 fallback。Exact-fit、variance boundary、rank deficiency、invalid denominator、
budget exhaustion、reference failure 與 branch mismatch 分別記錄。

Testing-only 政策在第二階段進入前安裝、workers 結束後復原，所有 OpenMP workers 共用該政策。
第一階段及 alpha training 不受影響。同一候選函式供 offline replay、proposal、boundary／polish refit、
nominal evaluation、recovery 與 final certification 使用。

## A：離線候選與預算

重新編譯原始基準並捕捉 337 個 shape datasets：首次失敗 1 份、final 168 份、
recovery-current 168 份；另外精確重播 offset fixture。三個版本化 shape fixture 與新 capture 逐位元一致。
CIF、map、manifest 的雜湊與前次報告一致，基準仍為 32 attempts、31 accepted、best iteration 28、
`recovery-failed`；final nominal p99 仍是
`[0.006143670517417728, 0.008207428922010910, 0.001048905842297531]`。

337 份真實輸入加上四個不讀模擬真值的有限正 variance fixture，全部通過方程與分支核對。
先使用 2,000 次探索預算，再依既定規則選擇 `{32,64,128,256}` 中涵蓋最大觀察成本兩倍的最小值。

| 每例額外方程評估（含起點檢查、Jacobian 與終點驗證） | 數值 |
| --- | ---: |
| 最小／中位數／p99／最大 | 18／19／28／35 |
| 341 例合計 | 6,522 |
| 凍結候選預算 | **128** |
| 凍結預算後通過數 | **341/341** |

分支參考延續另需 7,306 次固定點更新，單例最大 287 次。
原本求解器的迭代、候選的方程評估／線性求解、參考延續成本分開記錄。
這是同一 fold 的相關子問題及少量 fixtures，不能視為 341 個獨立資料集，也不構成 production 效能保證。

## B：同 state 的完整 operator

在原流程取得 `final`（best-28，搭配最後 certification 的 context）及 `recovery-current`。
前者不是 iteration 28 當時的背景／partition。各自在完全相同的 state、samples、鄰居資料、
frozen background、partition、alpha、range、ridge 與設定下，重新執行完整 joint offsets → shape refits。
每個政策使用獨立 workspace；不拼接離線 shape 結果、不串接前一政策的輸出。

| 政策 | Final qualification | 精化 atom 數 | Final recovery residual mean square |
| --- | --- | ---: | ---: |
| legacy | 未合格 | 0 | unavailable |
| failed-only | 合格 | 1 | 198.58450023990505 |
| fresh-residual | 合格 | 140 | 198.58324901721994 |

`failed-only` 保留原 SUCCESS 的資格語義，並未將全部原成功結果升級到 `1e-8`；
實驗另外保存這些終點的 fresh residual，不把政策下的 qualified 誤稱為所有方程都通過新標準。

`failed-only` 的三項 p99／max 完全不變；`fresh-residual` 的 shape p99 亦完全不變，
第三個 transformed-coordinate p99 僅由 `0.001048905842297531` 變為 `0.001048903784391063`。
相對 failed-only 的 residual mean square 下降約 `6.30e-6`，不是實質消除外層 residual。
Recovery-current 的兩種候選亦皆合格，mean square 分別為
`198.5887387543184`／`198.58748794407884`，結果一致。

兩個 checkpoint 的 joint-offset 子解皆逐欄相同，physical `ΔC=0`。
第三個 nominal 座標包含 peak normalization，因此即使 C 不變，shape 精化仍可使它略變。
Final 的最大 `|ΔA|/|ΔB|`：failed-only 為 `6.23e-8 / 4.37e-8`，fresh-residual 為 `1.87e-6 / 7.09e-8`。

Legacy 重評估與原 operator 完全一致；完整輸入 snapshot 前後雜湊相同。
額外 probe 隔離 capture 編號與工作量計數。168 筆完整 DB records、品質、軌跡摘要與 final certificate
都與原始基準一致。兩種候選皆通過 B，可進入閉迴路；未要求 B 預先滿足外層收斂門檻。

## C：閉迴路結果

兩種政策皆從與基準完全相同的第二階段初始 snapshot 開始，全部新增採用均通過分支核對，
未發生精化拒絕或參考解 fallback。
它們改變了 proposal 軌跡，均在 attempt 10 進入 recovery；不能把這些狀態當作 B 的 recovery-current。

| 結果 | legacy | failed-only | fresh-residual |
| --- | ---: | ---: | ---: |
| Attempts／accepted | 32／31 | 14／13 | 14／13 |
| 保存 best iteration | 28 | 6 | 6 |
| 精化採用次數 | 0 | 20 | 2,660 |
| 候選＋fresh 診斷方程評估 | — | 3,691 | 52,188 |
| 參考路徑額外更新 | — | 4,439 | 67,599 |
| Recovery operator 評估 | 1 | 9 | 9 |
| Recovery 接受次數 | 0 | 4 | 4 |
| 最終 certificate qualified／complete | false／true | true／true | true／true |
| 停止原因 | recovery-failed | recovery-failed | recovery-failed |
| Amplitude RMSE | 0.01143768998 | 0.01052019789 | 0.01052019871 |
| Width RMSE | 0.000364241208 | 0.000340424510 | 0.000340424456 |
| Physical offset RMSE | 0.00200402519 | 0.00189953656 | 0.00189953718 |

兩種候選各有 3,192 次 shape calls；單次候選最大成本分別為 33／34 次方程評估，均低於凍結預算 128。
真值只在事後 scorer 使用。Failed-only 在此固定資料的 A/B/C RMSE 約下降 8.02%／6.54%／5.21%；
全 residual 政策增加大量精化與參考工作，但未帶來實質不同的結果。

Failed-only 的 recovery residual mean square 軌跡為：

| Attempt | 起始 residual | 採用 factor | trial residual |
| --- | ---: | ---: | ---: |
| 10 | 211.330615 | 1/8 | 172.423594 |
| 11 | 172.423594 | 1/32 | 164.434824 |
| 12 | 164.434824 | 1/64 | 160.626935 |
| 13 | 160.626935 | 1/128 | 158.771268 |
| 14 | 158.771268 | 無 | 全部八個 trial 被 best-objective-bound 拒絕 |

Fresh-residual 的對應軌跡為 `211.329177 → 172.422367 → 164.433674 → 160.625836 → 158.770182`。
兩者下降皆約 **24.87%**。各自合計 33 個 trials：4 次接受、29 次因 best-objective-bound 拒絕。
所有 recovery-current／已評估 trial 的 nominal offset solves 都合格；沒有新的 offset IRLS blocker。
最後八個 trial 在 objective 檢查即遭拒絕，其 residual 為 unavailable，不能推論它們沒有 residual 改善。

這個 24.87% 是 **recovery 軌跡上的 residual 降幅**，不是保存 state 的降幅。
最後仍依既有政策保存 best iteration 6，沒有套用 final polish。
其 final p99 約為 `[0.006222757, 0.008322470, 0.001002178]`，仍未通過 `1e-4`。
保存 state 的 audit objective 約為 `0.79594`，也不能因 RMSE 降低就稱該 objective 比原基準更好。

因此，本輪證明了可靠內層資格讓 recovery 可以工作，但**沒有證明完整固定點收斂**。
14 attempts 並不等於通過原本的 25 輪收斂驗收；品質 gate 仍是 `uncalibrated`。
下一步應沿 failed-only 的較小改動研究既有方向與固定 best-objective 上界的相容性；
沒有證據支持直接升級 production 預設或全面精化原成功結果。

## 驗證與證據

- 15 項候選／共用數值實驗／唯讀 probe 測試通過，涵蓋已知第二根、alpha=0、exact-fit、
  近零噪聲、rank／denominator／floor、預算耗盡、終點 covariance 及政策觸發規則。
- 197 項既有相關測試中，196 項通過，原外部 replay 測試因未設定路徑先 skip；
  隨後提供實際 337 shape＋1 offset capture 重跑該測試並通過，沒有把 skip 當作驗收。
- Audit OFF 的 8 項候選與 numerical probe 測試、34 項 Python integration tests、文件連結檢查及 `lint_repo` 通過。
- 新 legacy／原基準／audit OFF j1 的 168 筆完整持久化 records、品質、summary 與 certificate 完全一致。
- 兩種候選各自的 audit ON j4/v4、平行 j4/v3、audit OFF j1/v3，完整持久化 records、品質、summary、
  certificate 以及候選／參考工作量完全一致；初始 snapshot 的數值內容相同。
- `BUILD_TESTING=OFF` 建置成功，symbol inspection 未發現新增實驗、capture 或候選函式。

版本化精簡證據為 [完整階段結果與 recovery trials](figures/endpoint-refinement/results.json)
及 [341 例逐分支核對表](figures/endpoint-refinement/branches.csv)。
它們保存結果與完整本機產物雜湊，不取代原始 DB、snapshot、逐 solve 證據與 provenance。

## 重跑方式與產物

建立 Debug／SYSTEM／OpenMP 測試版，關閉 Python bindings 與 UMAP，啟用 second-stage audit。
建置 `rhbm_gem_cli`、`rhbm_tests`、`mdpde_experiment`。另建 audit OFF 的 testing 版，
以及 `BUILD_TESTING=OFF` 的 production 隔離檢查版。

先用 [既有 capture runner](../../tests/integration/mdpde_experiment.py) 重建原始基準。
以下 `CAPTURES` 指向該 run 的 `solver-failures`，`MODEL`／`MAP` 指向前述固定輸入；
每次使用全新 output 目錄，不覆寫歷史資料。

```sh
python3 tests/integration/endpoint_refinement.py calibrate \
  --executable build/endpoint-refinement/on/bin/mdpde_experiment \
  --captures "$CAPTURES" --output build/endpoint-refinement/stage-a

python3 tests/integration/endpoint_refinement.py run --mode compare \
  --stage-a build/endpoint-refinement/stage-a/stage-a.json \
  --reference build/endpoint-refinement/baseline-run \
  --executable build/endpoint-refinement/on/bin/RHBM-GEM \
  --model "$MODEL" --map "$MAP" --output build/endpoint-refinement/compare

# 對 failed-only 與 fresh-residual 各自執行；v4 保存完整 recovery audit。
python3 tests/integration/endpoint_refinement.py run --mode failed-only \
  --stage-a build/endpoint-refinement/stage-a/stage-a.json \
  --stage-b build/endpoint-refinement/compare/stage-b.json \
  --reference build/endpoint-refinement/compare --verbosity 4 \
  --executable build/endpoint-refinement/on/bin/RHBM-GEM \
  --model "$MODEL" --map "$MAP" --output build/endpoint-refinement/failed-only-audit
```

使用 v3／j4 重跑驗證平行路徑；audit OFF／j1 檢查隔離性。
Debug logging 會依既有行為停用部分外層平行工作，因此 v4 的 j4 本身不算完整平行驗證。
Runner 要求 A 通過才可 compare，指定政策通過 B 才可 closed loop。
科學門檻未通過會保存結果並停止升階；缺案例、重播不符、資料或程式執行中改變則回傳失敗。

完整本機產物保存在 `build/endpoint-refinement/`，不納入版本控制：
獨立基準來源／建置、各 run 的 DB／log／score／certificate、source／binary／input hashes、
A 的逐例 root／reference／branch 結果、凍結政策、B 的 before／after snapshots 與逐 atom 比較、
C 的逐 solve JSONL（拒絕時附精確 dataset）及 audit 軌跡。應另行備份此目錄。
