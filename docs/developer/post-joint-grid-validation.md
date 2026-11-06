# Post-Joint grid-domain validation

## P0 baseline (2026-11-06)

Initial local and remote `develop`: `99de29ed0ba8e64bdba7d72317376f128f9426f6`.
Clean working tree; no branch or worktree created. Notion algorithm page was
read only. Its schema descriptions predate the repository (current Joint JSON
24, Stage JSON 4, SQLite 20); repository contracts take precedence.

Release, Apple Clang, C++20, Eigen, OpenMP, ROOT, UMAP, existing system dependencies.
Build parallelism 4. CTest parallelism 4 with OpenMP thread limit 4.
19/21 CTest groups passed (14.95 s). Existing failures:

- DocumentationLayoutTest: two missing runtime-guide links
  (`bcd-refactor-validation.md`, `tau-naming-validation.md`).
- Eigen thread restoration test with OMP_NUM_THREADS=1; passes unchanged when
  rerun with OMP_NUM_THREADS=2 and thread limit 2.

All Joint contract, numerical, workflow, persistence and CLI groups passed.
No 6Z6U/EMD-11103 data found under the project directory or nearby project
folders searched; no real-map benchmark is claimed.

The existing 20-qualified-target, two-Halo synthetic workflow fixture has
truth A=2, tau=0.65, C=-0.2 for isolated donor targets. FirstStage stores 200
Fibonacci samples per target, 100 fall within the Local 1 Å range. The baseline
has 20/20 available locals, alpha=0, OLS/MDPDE A=2.003270672,
tau=0.651122504, information condition=48.178684. Peeling: 4000/4000 available.
All 20 groups have one eligible member and remain unavailable (`single-member`);
there is no group posterior or WEB timing to compare for this fixture.

One isolated Release measurement (not a stable performance estimate):

| Phase | seconds |
|---|---:|
| Context | 0.006252 |
| Peeling | 0.034073 |
| Own charge | 0.006065 |
| Local alpha | 0.034147 |
| Local fit | 0.000191 |
| Group alpha | 0.000000291 |
| WEB | 0 |
| Total postprocessing | 0.081158 |

RSS snapshots: before 50,413,568 B; context 53,706,752 B; peeling 58,556,416 B;
local 58,736,640 B; group 58,753,024 B. Largest snapshot delta 8,339,456 B;
this is not a measured process peak. Estimated payload 8,236,168 B.
Machine-readable evidence and test logs are in the ignored
`build/release/grid-domain-validation/` directory. The regression fixture emits
compact local, group, phase-time and memory evidence through GoogleTest XML.

## Data flow and ownership before the change

| Stage | Source / identity | Ownership / lifetime | Failure and coverage |
|---|---|---|---|
| Joint construction | Map voxels; decimal global grid IDs; atom support rows | JointProblem owns immutable OriginalInput and compact fitting input | Builder rejects malformed/nonfinite observations; dual-domain tests |
| Fitting mask | Excludes Halo-supported original rows | Original-to-fitting mappings retained in problem and snapshot | Pure Exclusion tests |
| Qualification | Certified target endpoint and point evidence | Joint snapshot and atom Second-stage source/run | Unqualified runs skip downstream; runtime/target tests |
| Halo reconstruction | Geometry kNN-4 qualified target medians | Snapshot provenance; borrowed donor workspace during postprocessing | Missing donor rejection; Halo tests |
| FirstStage sampling | Fibonacci positions, map interpolation, selection mask | LocalPotentialEntry raw buffer | Geometry/selection remain FirstStage-specific; sampler tests |
| Peeling | Original prediction plus reconstruction-only tricubic nodes at raw positions | Context borrows map/model/problem/snapshot; results own raw-index-paired responses | Missing contributor and stencil diagnostics; partial-selection tests |
| Charge correction | Tricubic own C/tau field at same raw positions | Temporary PostJointLocalSamples | raw_sample_index/count identity; interpolation reference tests |
| Local alpha/OLS/MDPDE | Positive corrected samples within 1 Å; log Gaussian design | Transient design/fit; compact result owns diagnostics and information | QR, solve, working covariance, SPD information qualification |
| Group alpha/WEB | Qualified Local points/information; fixed Joint charge separate | Model group summaries and member posterior | Source/run and local evidence validation; unavailable/singular states |
| Persistence/invalidation | Stage JSON 4 / SQLite 20; Joint JSON 24 | ModelAnalysisData, editor atomic batches | Raw edits clear dependent results; endpoint edits invalidate run; atom geometry edits clear derived model state |

Existing raw-index coupling also affects the Peeling ratio/UI views, Stage JSON,
LocalPotentialEntry reconstruction, editor batch validation and Group evidence
validation. MapObject is independent of ModelObject; map edits cannot notify an
unrelated model. This existing ownership boundary needs explicit handling when
reusing a model against changed map observations.
