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
| B | How do fixed-seed noise and position mismatch affect estimates? | [`joint_statistical_experiment.py`](../../tests/integration/joint_statistical_experiment.py) and `joint_statistical_experiment` | RESEARCH | Fixed-seed outcomes, with execution, convergence, and estimate error reported separately |
| C | Did the full command path run and what resources did it use? | [`joint_benchmark.py`](../../tests/integration/joint_benchmark.py), `--profile command` | BENCHMARK | Analysis, save/reload, JSON/CSV export, and process/resource classification |

There is no combined Stage A/B/C runner. Statistical research is not a
correctness test, and process completion does not replace numerical or
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

Retained halo and selected target points use the same converged Joint state.
Point-quality warnings do not silently remove retained contributors. Post-fit
peeling, local MDPDE covariance, and group fitting retain their own local
qualification and resource controls. Joint results do not supply a global
parameter covariance.
