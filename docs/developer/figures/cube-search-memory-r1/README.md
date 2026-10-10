# Cube search memory frontier r1

This closed campaign compared SPQR ordering choices for the historical
OperatorPcg search route on cube workloads, including diagonal and Schwarz
preconditioners. BEST and METIS reduced fill and peak memory on completed
512/1024 cases, but no ordering completed cube-2048 within the 4 GiB gate; the
multi-block result therefore remains insufficient evidence and production
ordering was unchanged. The 6 GiB diagnostic is excluded from the formal gate.

`cube-memory-summary.csv` and `cube-memory-summary.json` retain the per-case
results; `cube-memory-analysis.json` retains the comparison and decision.
The 15 individual execution JSONs (784,282 bytes) were pruned. Output paths in
the manifest and summary remain historical provenance, not retained files.
