# Fixed-target Halo FullABC P4：Production Promotion Gate

日期：2026-10-21
分支：`develop`  期別：P4

## 決定

**Do Not Promote。** Production 保留 `geometry-local-target-median-v1`；FullABC 與 Recovery 繼續作為明確 opt-in research path。此決定由已觀測到的數值/統計 gate failure 支持，不依賴完整 6Z6U Chain A run 是否完成；未測的 resource/integration gates也不能視為通過。

## Promotion Gate

| 領域 | 判定 | 證據 |
|---|---|---|
| 第一輪 Target endpoint凍結 | Pass（實驗範圍） | 每個 paired case只估計一次 Target；P2 Target A/B/C drift = 0 |
| `n_h<3`唯一事前固定 | Pass（程式/實驗） | 343/519 Halos預先固定；176/519均先嘗試FullABC；沒有rank/chemistry pre-pruning |
| FullABC qualification | Fail | P2五個 local cases aggregate FullABC qualified 0/5；P2 raw FullABC不可用 |
| Recovery independently qualified | Partial / Fail gate | 5/5 reconstruction available且assembly replay通過，但boundary case仍保留巨大振幅端點；qualification不足以表示可信參數或科學精度 |
| KKT、width、replay、endpoint | Partial | P0/P1受控小型回歸通過；boundary Recovery replay/objective精確一致仍出現 A≈8.35e122、B≈0.0934 Å，證明replay本身不足以排除端點病態 |
| 真實幾何 | Pass（幾何驗證） | 保留 1.179–1.191 Å bonded contacts；無 exact collisions或nonbonded <1.0 Å；support overlaps已量測 |
| Halo signal accuracy | Partial | 三短鍵seed component-local Recovery median signal RMSE 0.00652，優於kNN-4 0.07363；但單一boundary worst case 0.04251且其他策略更差 |
| Peeling / worst case | Fail | 三短鍵seed recovery median Peeling 約0.6284–0.6289；boundary recovery Peeling RMSE約1.91e57、Tail Ratio MAE約8.83e56 |
| Local regression | Partial | 小型案例median Local A/B約0.273–0.280/0.0215–0.0217；boundary recovery Local A約0.36–0.40，且Peeling失控 |
| Group Alpha / Group posterior | Not measured for training | 五個 local scenes均0 group members，alpha及posterior為null；完整Chain A run結果尚待 |
| 4 GB大型Halo RSS / Scale B | Fail / Not measured | 完整 Chain A runner運行約2小時未輸出、以Ctrl-C中止(exit 130)；完成run的peak RSS未量測。Scale D既有1.471 GB數字是zero-Halo prepare/preflight，不能認證Halo FullABC |
| Scale C大型Halo workload | Not measured | 無更大含Halo selection的完整測量 |
| SQLite / UMAP / Two Stage integration | Not measured this round | 本輪沒有接入Production policy或schema；未跑完整 save/reload、UMAP與Two Stage promotion integration |
| Production regression | Partial / Fail gate | targeted Halo tests 53/53通過；full CTest 31項中28項通過，`joint_production_bounded_synthetic`重建並隔離重跑後通過，`joint_component_regression`與`joint_component_physical_smoke`各觸及600秒timeout。這兩項無法作為完整integration pass |

## 生產契約與未來條件

本輪沒有改動Production Halo policy、SQLite schema或Target fit。Pre-fixed/FullABC/Recovery provenance仍區分；recovered endpoint不標成FullABC solver success。即使未來只升格某一 Recovery policy，需先增加能拒絕窄width/巨大振幅的參數穩定性qualification與測試，重新跑 all-case accuracy/worst-case/Group Alpha，並完成含 Halo的 Scale B/Scale C RSS、runtime、SQLite save/reload、UMAP、Two Stage及回歸。

完整數值證據與P2逐案來源見 [P2 physical validation report](fixed-target-halo-fullabc-p2-physical-validation-20261021.md) 及 `docs/research/results/halo_fullabc_p2/`。

## 本輪回歸執行結果

- Build：`cmake --build build/debug --target joint_halo_fullabc_validation --parallel 4`成功；重建`joint_production_bounded_synthetic`也成功。所有建構最多4個parallel jobs。
- Targeted：`RHBM-GEM-TEST --gtest_filter='HaloMedianReconstructionTest.*:HaloFullABCJointEstimationTest.*:PureHaloSupportExclusionTest.*'`，53/53通過，19.488秒。
- Full CTest：`ctest --test-dir build/debug --output-on-failure -j1`完成31項，28 pass、3 fail，總計約3,586秒。`joint_production_bounded_synthetic`第一次出現bus error；重建後直接執行及單項CTest均通過（3.22秒）。`joint_component_regression`與`joint_component_physical_smoke`觸及各自600秒timeout，未產生完整summary。
- 額外回歸重試：無CTest timeout的`joint_component_regression`執行約11分鐘仍停在第一個C++ fixture且沒有result JSON，Ctrl-C結束，exit 130。`joint_component_physical_smoke`直接重試逾10分鐘後仍只寫出double precision結果，沒有float32或Model/Map reload結果；Ctrl-C結束，exit 130。
- CTest timeout清理輸出`/bin/ps: Operation not permitted`屬執行環境權限限制；它不代表測試本身通過或造成 solver 故障。未能完成的integration gate仍列為未通過。
