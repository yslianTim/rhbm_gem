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

## P1 sample contract

Added PostJointGridObservation (Original row, global grid index, squared
radius, raw voxel response, physical position) and a support-based collector.
The collector rejects duplicate rows, invalid grid/mapping/geometry and missing
or changed original observations. It returns stable Original-row order and
includes the 2 Å boundary independently of the fitting mask. Production still
uses the baseline downstream path in this phase. Release build and all 30
selected grid-contract, Pure Exclusion and sampler tests pass.

## P2–P3 grid subtraction and charge correction

Production Peeling now owns only original-grid samples within 2 Å. Each sparse
Original support contribution is evaluated once; each Target then sums cached
non-self contributions in extended precision. This avoids repeated kernel
work and avoids subtracting a large own signal from a total prediction. No
dense observation-by-atom matrix, interpolation, stencil fallback or node cache
is used. The same voxel distance evaluates the Target's own C/tau charge.

Target points require Passed global/target/component convergence and parameter
validity. Local conditioning remains diagnostic for Target endpoints; Halo
donors retain their stricter Adequate-conditioning requirement. Halo points
remain geometry-local median reconstructions, not optimizer parameters.

The switch necessarily wires the editor, Local preparation, UI ratios, Group
validation and persistence together: a mixed raw-index/grid-index contract
would be invalid even for one commit. Stage JSON is now 5 and SQLite 21;
Joint JSON remains 24. Old Stage/SQLite contracts are rejected. Failures carry
Original row and global grid identities. Corrected, positive and retained
counts are separate; nonpositive finite responses are not reconstruction
failures. Production fitting remains 1 Å pending the radius study.

Release build succeeds. 20/21 CTest groups pass (36.64 s, serial CTest,
OMP_NUM_THREADS=2 / OMP_THREAD_LIMIT=4). The only failure is the two pre-existing
runtime-guide links recorded above. Independent Gaussian/erf references cover
single and overlapping Targets, signed nonzero charge, Halo medians, narrow
widths, off-grid atoms, missing contributor state and the total-prediction
identity. FirstStage sampler tests, Joint numerical tests, Stage/SQLite
round-trips, Group integration and CLI workflow tests pass.

## P4 radius decision

768 paired synthetic trials per path: spacing 0.35/0.7 Å, center alignment
0/0.17 Å (with unequal y/z shifts), base tau 0.4/0.7/0.9 Å, additive Gaussian
noise SD 0.5%/3% of the Gaussian peak, independent or spatially correlated
noise. Each of the 48 configurations has 16 trials with independent innovations
(seed 1729+trial); A and tau vary by deterministic independent sinusoidal
patterns. Signed charge is known and corrected at each voxel. The correlated
field averages a normalized 3x3x3 innovation neighborhood (adjacent-voxel
correlation approximately 2/3 away from map boundaries). This does not model
all cryo-EM noise. Neighbor/Halo subtraction is independently covered by the
numerical tests; this statistical experiment conditions on the true endpoint
and does not estimate Halo uncertainty.

Legacy comparison recreates the single-Target baseline Fibonacci map sampling
and tricubic own-charge field, then uses the same Gaussian estimator and 1 Å
range. Grid variants use the production Peeling/correction and Local fitting
qualification with a test-only radius argument. FirstStage defaults are
unchanged. Training, fitting and information now consume one selected dataset;
the common PreparedLocalGaussianDesign API defines the explicit radial range.

| Metric | Legacy Fibonacci 1 Å | Grid 1 Å | Grid 1.5 Å | Grid 2 Å |
|---|---:|---:|---:|---:|
|Qualified Local / 768|768/768|767/768|766/768|768/768|
|A bias|-0.015803|0.012441|-0.030434|-0.063734|
|A RMSE|0.121849|0.098797|0.152310|0.501155|
|tau bias (Å)|0.002373|0.000476|0.023456|0.115998|
|tau RMSE (Å)|0.016132|0.012598|0.057009|0.214714|
|Positive retained / radius candidates|99.98%|99.77%|93.42%|86.49%|
|Available synthetic groups / 48|48/48|48/48|48/48|48/48|
|Posterior nominal 95% ellipse coverage|6.25%|44.59%|27.81%|13.67%|
|Local beta MSE / mean working covariance trace|64.82|2.20|58.14|316.33|
|Mean information condition|59.62|68.35|24.60|19.59|

Groups contain the 16 qualified synthetic members per configuration and use
the current Group alpha, MDPDE and WEB kernels; production source/run eligibility
is tested separately in the workflow tests. Group alpha training retains its
existing random split policy, so individual Group summaries are exploratory,
not a deterministic covariance calibration certificate. Raw per-configuration
OLS/MDPDE bias, RMSE, width error SD, alpha, conditioning, sample counts,
posterior MSE/covariance traces, coverage, failures and timings are retained in
[the compact CSV](post-joint-grid-statistics.csv).

**Decision: retain 1 Å grid fitting; do not promote 2 Å.** For base tau=0.4 Å,
mean tau bias at 2 Å is +0.3416 Å with independent noise and +0.3343 Å with
correlated noise, versus +0.00083/+0.00407 Å at 1 Å. Broader tau=0.9 Å cases
benefit from the larger design leverage (2 Å bias -0.00016/-0.00115 Å), but this
does not offset the narrow-peak tail selection bias. Improved design condition
and 768/768 numerical qualifications at 2 Å do not establish statistical
correctness. No estimator qualification was relaxed. The three unavailable
grid outcomes (one at 1 Å, two at 1.5 Å) are `local-mdpde-unqualified`.

Working covariance is conditional and diagonal by observation. The experiment
shows substantial overconfidence even at 1 Å, especially with correlated noise;
the legacy Fibonacci interpolation has still worse dependence/calibration.
These are measured limitations, not a claim that voxels are independent or
that any nominal 95% coverage is calibrated. Joint charge and Halo uncertainty
remain conditioned out. Changing the Group statistical model is out of scope.

Across the 768 grid fits, alpha training totals are 1.35/2.53/3.73 s and fitting
0.0132/0.0239/0.0461 s for 1/1.5/2 Å in the first run. This isolates the extra
fitting-domain cost from removal of interpolation. The expanded four-path study
ran in 12.30 s with measured process peak RSS 32,030,720 B. This is a small
synthetic study peak, not the real-map pipeline peak.

N_grid is all candidate voxels at r²<=4. N_corrected counts finite corrected
responses across that candidate domain. N_positive counts positive corrected
responses across the same domain; N_retained is the positive finite subset
inside the chosen fitting radius actually used for alpha, fit and information.
The CSV names N_positive explicitly as `positive_all_2A`; retention fractions
use N_retained divided by the radius-specific candidate count. Nonpositive
responses are not Peeling failures.

## P5 integrity, invalidation and persistence

Release build and all 21 CTest groups pass (39.25 s, at most four threads).
The complete-source Group fixture varies A/tau across 20 separated Targets,
uses the actual workflow, and round-trips Gaussian posterior covariance and
fixed Joint charge through SQLite and Stage JSON. It replaces abstract Group
persistence fixtures without a Joint snapshot. Ordered grid IDs, Original row
mapping, map positions, response/reason pairing, charge identity and all four
Local counts are checked. Available Local evidence independently replays the
fixed charge to verify correction/positive/retained counts. Local Gaussian
offsets must be zero. Empty positive datasets retain zero actual solver rows;
the historical FirstStage/Two-Stage sentinel behavior remains unchanged.

A live map revision binds its analyses and model copies. Value replacement,
clearing and normalization invalidate them at next access; copying or destroying
the map does not invalidate independent/historical evidence. Position changes
clear analysis and spatial caches. Replacing/removing a Joint snapshot clears
downstream grid records. Grid geometry or known map/model content identity
cannot be changed by metadata-only updates. FirstStage raw edits retain grid
results. These dependencies have regression coverage.

The Peeling contributor universe now contains only Original contributors.
Every 2 Å candidate lies inside a Target's Original 2.5 Å support, so all atoms
with an actual contribution there are already included by the Original builder.
Unindexed contributors are rejected; no reconstruction-only interpolation
lookup or fallback remains. Sparse cached values require no dense matrix.

Current documentation describes Stage JSON 5 / SQLite 21 / Joint JSON 24.
The two baseline broken links are retained as plain historical report names
rather than links to nonexistent files; all documentation link checks pass.
The real files were subsequently located through `local/run_example.sh` in
`/Users/yslian/data`; the earlier negative search applied only to the project
folders searched. A paired real-map baseline was rebuilt from an ordinary
`git archive` source directory (not a Git worktree) at the initial commit.
