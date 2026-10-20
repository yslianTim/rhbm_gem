# Joint capability and ownership guide

Joint result handling separates command execution, global numerical convergence,
local point quality, and scientific interpretation. A zero command exit code
means the result was saved; it does not by itself mean convergence passed.
Poor local conditioning or parameter sensitivity is reported as point-quality
evidence and does not invalidate a numerically converged Joint state.

## Current owners

| Stage | Question | Owner | Build category | Output meaning |
| --- | --- | --- | --- | --- |
| A | Did the accepted Joint endpoint converge and replay consistently? | `JointComponentEstimator` and `joint_component_runtime.py` | CORE | Global stationarity, endpoint verification, component and assembly consistency |
| B | Are the current numerical and production-workflow invariants covered? | Joint contract/workflow/numerical tests, bound-active synthetic, persistence, and CLI regressions | CORE | Current Active-set LM behavior and workflow contracts |
| C | Did the full command path run and what resources did it use? | [`joint_benchmark.py`](../../tests/integration/joint_benchmark.py), `--profile command` | BENCHMARK | Analysis, save/reload, JSON/CSV export, and process/resource classification |

There is no combined Stage A/B/C runner. The fixed-seed noise/mismatch study
and its 450 saved outcomes are historical evidence from an earlier solver
path, not scientific validation of the current Active-set LM. Its campaign
runner is retired; a deterministic generator contract preserves the noise
stream and observation hash. Process completion does not replace numerical or
persistence tests.

## Point-quality scope

After convergence, local conditioning describes the accepted state within each
FixedNeighbor block while external neighbors are held fixed. Parameter
stability combines recent sweep changes with bounded, sensitivity-ranked local
perturbations of selected blocks. Conditioning and stability can identify
conditional numerical weakness or parameter sensitivity. Target point
validation maps block evidence through compact atom-to-block indexes. These
diagnostics do not certify global parameter identifiability, and local
perturbations do not replace global multi-start or block-order studies.

Retired global-rank and weak-direction tools are not current entry points.
Historical rank records remain readable with their original meaning. They are
not converted into local conditioning or stability evidence.

## Downstream boundaries

Production selected target points use the qualified target-only Joint state
from Pure Halo-Support Exclusion. Retained non-target atoms exclude fitting
rows and do not provide Joint parameters or FirstStage halo values. Point-
quality warnings do not silently remove target contributors. After a qualified
endpoint, Halo median reconstruction is a separate availability-checked step
for original-domain peeling; local MDPDE covariance and Group Fitting remain
target-only and retain their own qualification and resource controls.

An unqualified candidate keeps diagnostics and stop reason but never becomes a
formal target estimate and does not start downstream inference. A qualified
Pure Exclusion endpoint can still have unavailable Halo reconstruction when no
qualified target donor exists. The explicit Fixed-Halo research/compatibility
APIs retain their separate FirstStage `(A, B, 0)` conditional contract. Joint
results do not supply a global parameter covariance or certify FirstStage halo
statistical accuracy.
