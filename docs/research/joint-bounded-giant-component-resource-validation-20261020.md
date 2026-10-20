# Bounded Joint giant-component resource validation — 2026-10-20

> Historical comparison note: the prior FixedHalo-specific synthetic test was
> migrated into a target-domain memory regression with 10,000 non-target
> contributors. The migrated test no longer models fixed predictions.

## Decision

The exact requested shape now has a synthetic, single-component preparation
and endpoint-certification preflight with all 37,406 contributors selected as
targets. Its preflight and measured isolated-process RSS are below the 4 GB
target. The complete FirstStage → bounded search → local diagnostics → assembly
workflow was not run at this scale. The measured preflight is not evidence that
the full 37,406-target fit stays below 4 GB; large-workload promotion remains
blocked pending a complete run or equivalent stage-by-stage measurement.

## Memory units and exact-shape result

This report uses decimal GB (`1 GB = 1,000,000,000 bytes`) for the requested
4 GB ceiling. The FixedNeighbor memory planner's 3,200,000,000-byte budget is
3.2 GB decimal. GiB values are also shown for comparison.

The isolated GTest workload used 5,938,087 observation rows, 37,406 targets and
contributors, one connected target component, and 3,118 blocks at core size
12. It held 74,812 support entries, 74,812 A/C design columns, and 112,218
target parameters. The synthetic support joins adjacent targets through shared
rows. It does not represent measured physical density or real FirstStage donor
fits.

| Measure | Bytes | Decimal GB | GiB |
|---|---:|---:|---:|
| FixedNeighbor preparation estimate | 912,690,261 | 0.9127 | 0.8500 |
| Global state evaluation estimate | 975,208,628 | 0.9752 | 0.9082 |
| Endpoint verification estimate | 627,146,451 | 0.6271 | 0.5841 |
| Memory planning budget | 3,200,000,000 | 3.2000 | 2.9802 |
| Test-process baseline RSS | 24,887,296 | 0.0249 | 0.0232 |
| Measured process peak RSS | 1,471,479,808 | 1.4715 | 1.3700 |
| Requested RSS ceiling | 4,000,000,000 | 4.0000 | 3.7253 |

This isolated test run passed both preflight checks and endpoint certification.
It performs zero search sweeps and zero Active-set LM subproblem solves. The
accepted state has zero A/C coefficients on a zero-valued synthetic observation
vector; independent endpoint replay therefore validates the certification
path, not a fitted high-scale estimate. Four isolated all-target preflight
invocations measured process peaks from 1,286,668,288 to 1,471,479,808 bytes;
no complete-fit RSS was measured. At the time of this report,
`ProductionScaleSingleComponentPreparationStaysBelowBudget` tested 27,406
targets plus 10,000 Fixed Halos separately. That test was later migrated into
a target-domain memory regression retaining 10,000 non-target contributors;
it no longer constructs a fixed prediction.

The run allocates a 47,504,696-byte effective observation vector and retains
row identities, sparse support, target state, compact block mappings, and
endpoint replay workspaces. No global dense observation-by-target matrix is
allocated. Such a matrix for 5,938,087 rows by 74,812 A/C columns would require
3,553,921,317,152 bytes before other storage.

The historical FixedHalo comparison setup also passed independently at 27,406
targets plus 10,000 halos: its preparation estimate was 969,975,622 bytes,
endpoint estimate 624,026,451 bytes, and measured process peak RSS
1,458,044,928 bytes. Those measurements describe that saved test revision;
the current migrated regression runs no fixed prediction or search sweep.

## What was and was not measured

| Stage | Evidence at 37,406 targets |
|---|---|
| Problem/data construction | Synthetic input creation ran; isolated process RSS includes it. No Production `BuildPureHaloSupportExclusionProblem` memory estimate was captured. |
| FirstStage initialization | Not run. |
| Qualified-donor calibration | Not run. |
| FixedNeighbor structural preparation | Ran through zero-sweep preparation; estimated peak 912,690,261 B. |
| Global state evaluation | Preflight estimate 975,208,628 B. |
| Active-set LM | Not run; zero subproblem solves. |
| Local A/C block solve | Not run. |
| Local conditioning | Not run. |
| Conditional local reprofiling | Not run. |
| Parameter Stability | Not run. |
| Endpoint verification / Independent Replay | Ran on the zero state; estimated peak 627,146,451 B and certification passed. |
| Component assembly | Not run. |
| Persistence preparation and downstream | Not run. |

The whole-process peak RSS is below 4 GB for this preflight test only. It must
not be interpreted as a complete fitting RSS result.

## Medium-scale completed candidate runs

The Release non-instrumented bounded candidate completed both 1,418-target
6Z6U runs after the local diagnostic changes:

| Workload | Wall time | Process peak RSS | Selected local probes | Local solves / evaluations | Maximum local workspace estimate |
|---|---:|---:|---:|---:|---:|
| Simulation, 180,478 observations | 303.694 s | 708,755,456 B | 5 of 138 | 30 / 120 | 106,629,766 B |
| Unsharpened, 189,377 observations | 435.770 s | 1,389,789,184 B | 6 of 140 | 19 / 159 | 110,724,617 B |

Both runs passed RuntimeConvergence, bound feasibility, projected Width KKT,
and independent replay. The simulation had 2 Stable, 3 Unstable, and 133
NotAssessed blocks. The unsharpened run had 1 Unstable, 5 Unavailable, and 134
NotAssessed blocks; both bound-active blocks were selected. These are complete
medium-scale fits, not extrapolated estimates for 37,406 targets.

## Allocation and cost review

The bounded preparation uses sparse support entries and compact row mappings;
it does not introduce a global dense design, covariance, or per-block
full-length observation-row mapping. Endpoint verification is streaming. The
new local-conditioning values and probe-selection flags are fixed-size fields
stored with per-block evidence, so retained diagnostic metadata grows with the
number of blocks. The local stability perturbation budget remains eight blocks
per fit, with bounded per-block CPU, profile evaluations, and a local memory
preflight. On the two completed 6Z6U runs, the maximum per-block local workspace
estimate was about 107–111 MB.

No high-scale LocalConditioning or Parameter Stability allocation was made in
the 37,406-target test. Their total retained block metadata and worst local
workspace at 3,118 blocks are therefore not claimed as measured. The 4 GB
full-fit gate remains open until these stages and Active-set LM run at the
specified scale under a process-level RSS measurement.

## Validation

- Build: `cmake --build build --target rhbm_tests --parallel 4` passed.
- Direct exact-shape preflight:
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 build/bin/RHBM-GEM-TEST --gtest_filter='JointFixedNeighborBlockCoordinateTest.BoundedAllTargetLargeShapePreflightStaysBelowBudget'`
  passed.
- Existing Fixed-Halo scale test remains in place and was run independently.
- Focused CTest: `rhbm_tests_joint_component_numerical` and
  `joint_production_bounded_synthetic` passed 2/2.
- Full 37,406-target Active-set LM and downstream workflow: not run.

No Notion documents or public CLI options were changed.
