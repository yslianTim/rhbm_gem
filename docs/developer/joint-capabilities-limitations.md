# Joint capability and ownership guide

A saved result and a qualified result are different states. The command can
save candidate diagnostics and a stop reason, but only a qualified endpoint
publishes formal Target estimates and starts downstream fitting.

## Current production behavior

| Area | Current behavior | Limit |
| --- | --- | --- |
| Joint fit | Target-only Pure Halo-Support Exclusion with FixedNeighbor and bounded Active-set LM | Halo parameters are not optimized by Joint fitting |
| Qualification | Complete-sweep stationarity, endpoint verification, Independent Replay, and component/assembly checks | A matching replay alone does not establish convergence |
| Point quality | v2 combines numerical convergence, valid parameters and Local Conditioning | Conditioning warnings do not establish global identifiability or reject a qualified endpoint |
| Halo values | Geometry-local median from the nearest four qualified Target estimates | Reconstruction is unavailable when there are no qualified donors |
| Peeling | Uses the Original observation domain and qualified Halo reconstruction | Retained Halo rows are excluded from Joint fitting |
| Local and Group fitting | Operate on selected Targets with their own qualification rules | Joint does not provide global parameter covariance or certify FirstStage Halo accuracy |
| Historical results | Only current SQLite/JSON schemas are accepted | Historical payloads require the corresponding Git version |

The default CTest groups cover the current contract, workflow, numerics,
persistence, and CLI smoke behavior. The
[runtime guide](joint-component-runtime.md) owns the detailed production
contract.
