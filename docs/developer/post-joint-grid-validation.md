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

Trial t=0..15 uses A=2+0.2 sin(1.7t), tau=w[1+0.08 cos(1.1t)],
C=+0.2 for even t and -0.2 for odd t, and center (a,0.37a,-0.23a).
The 0.35 Å map is 19³ with origin (-3.15,-3.15,-3.15); the 0.7 Å map is
11³ with origin (-3.5,-3.5,-3.5). The same innovation field is paired across
radii. Width variability in the CSV is the SD of width error, not the
unconditional spread of the intentionally varied truth widths. Posterior
coverage uses the 2D beta-coordinate 95% chi-square threshold 5.991464547.

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
|Within-configuration tau error SD, pooled RMS (Å)|0.012357|0.012237|0.023061|0.049525|
|OLS A RMSE on qualified cases|0.100552|0.090896|0.195406|0.467861|
|OLS tau RMSE on qualified cases (Å)|0.013382|0.012490|0.057429|0.211591|
|Mean Local alpha|0.23841|0.34811|0.35783|0.38854|
|Positive retained / radius candidates|99.98%|99.77%|93.42%|86.49%|
|Available synthetic groups / 48|48/48|48/48|48/48|48/48|
|Posterior nominal 95% ellipse coverage|6.25%|44.59%|27.81%|13.67%|
|Local beta MSE / mean working covariance trace|64.82|2.20|58.14|316.33|
|Posterior beta MSE / mean posterior covariance trace|67.99|4.56|119.49|673.52|
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
At grid 1 Å, posterior ellipse coverage is 239/384 (62.24%) with independent
noise and 103/383 (26.89%) with correlated noise; Local beta MSE divided by
mean working covariance trace is 1.65 and 3.29 respectively. Every variant
has trials selecting alpha at the upper boundary 1.0; means alone do not
establish stable alpha selection. Per-configuration values are in the CSV.
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

## P6 real-map paired benchmark

6Z6U / EMD-11103 unsharpened map, all 70,190 model atoms retained; selected
non-H chain A sequence IDs 30–50: 181 Targets. Map 480³, spacing 0.4920000076 Å.
The workflow normalizes the map and supplies reference Gaussian parameters;
Q-score calculation is omitted because it does not enter Joint/Local fitting.
The archived initial source and final source use Release, Apple Clang, Eigen,
ROOT and OpenMP. The temporary baseline build disables optional UMAP/experimental
commands, which are absent from this benchmark and do not affect fitting.
Three paired independent processes per implementation, one numerical thread,
no concurrent builds/tests during measurement. Times include phase bookkeeping,
not external file import. They are wall times on this machine.

All six runs qualify. The 181 Target A/tau/C floating-point bit patterns match
across all runs (SHA-256 of big-endian doubles ordered by Atom ID:
`a5a57ea3ac5f6f20b1958c1a6e37dc4196d2e23f15c11ae8498b21e703acefdf`).
Joint objective is exactly 0.01566223279243328; Original rows 32,933, Fitting rows
19,757, reconstructed Halos 314. No optimizer, BCD, KKT or replay implementation
was changed.

| Median phase | Baseline seconds | Grid 1 Å seconds |
|---|---:|---:|
| Context | 0.005714291 | 0.004159625 |
| Peeling | 0.605913083 | 0.030749875 |
| Own charge / preparation | 0.030867829 | 0.007354877 |
| Local alpha | 0.571445703 | 0.346706747 |
| Local fit | 0.005666874 | 0.002600571 |
| Group alpha | 0.000000208 | 0.000000292 |
| Group WEB | 0.000129289 | 0.000132541 |
| Total postprocessing | 1.222804708 | 0.410749625 |

Median total postprocessing falls 66.4% (2.98x faster); Peeling is
19.7x and own-charge preparation 4.20x faster.
This comparison retains a 1 Å fitting range: it measures the observation-domain
change and interpolation removal. Separate P4 radius timings quantify the cost
of expanding to 1.5/2 Å.
Group WEB changes by 0.003 ms in absolute median time;
sub-millisecond timings do not establish a general WEB speed change.

Peeling availability: baseline 36,200/36,200 Fibonacci samples; final
50,936/50,936 voxel observations. Retained Local solver rows total 13,505 vs
6,344. Both have 178/181 qualified Local results and 25/85 available Groups.
Both exclude two invalid Local parameter outcomes and one unqualified MDPDE.
The identities change: 401 and 523 gain eligibility; 671 and 735 lose it;
737 remains unavailable. Failures never enter Group inference. These correspond
to TYR C, SER OG, ASP CB, LYS CE and LYS NZ in the input atom-site records;
the latter two have alternate occupancy 0.5 and high temperature factors.
Those attributes are observations, not demonstrated causes of failure.
Median available tau is 0.65357 vs 0.63716 Å; no real-map truth is known, so
these differences are not accuracy claims. On members qualified in both paths,
Target 400 changes from tau 0.14173 to 0.60954 Å. Some targets and Groups change
substantially despite unchanged aggregate coverage; see the member data.

Median measured process peak RSS: 1,498,513,408 B vs 1,498,923,008 B
(+409,600 B, +0.03%). This includes map import and is not isolated
to postprocessing. Median largest sampled postprocessing RSS increase is
17,596,416 B vs 9,601,024 B (baseline range
17,383,424–20,938,752 B; final range 9,437,184–9,732,096 B).
Allocation reuse makes snapshot deltas variable. Estimated postprocessing payload
falls 19,451,926 B to 8,658,344 B
(55.5%). No material RSS growth or dense observation-by-atom allocation was observed.

Detailed evidence: [phase/RSS repeats](post-joint-grid-real-performance.csv),
[Local and Joint points](post-joint-grid-real-local.csv),
[Group means, MDPDE, member IDs, posterior covariance and fixed charge](post-joint-grid-real-groups.json).
Real-map inference conditions on fitted Joint and reconstructed Halo profiles;
these data do not calibrate charge/Halo uncertainty or Group coverage.

### Overlapping Targets and Halo misspecification

An additional paired study uses spacing 0.3 Å, a 31x21x21 map at origin
(-3,-3,-3), two off-grid/overlapping Targets at (0.13,0.21,-0.07) and
(1.5,0.1,0), and an optional Halo at (3,0,0). Target charges are -0.3/+0.2.
Trial t uses A1=2+0.2 sin(1.7t), tau1=0.4[1+0.05 cos(1.1t)],
A2=3+0.15 cos(1.3t), tau2=0.8[1+0.05 sin(0.7t)]. Independent additive noise
has SD 0.5% of Target 1's Gaussian peak, seed 8128+t. The source endpoints
are controlled true Target points, not newly optimized noisy endpoints.
All 16 trials/radius/Target remain numerically qualified.

The no-Halo and exact-median Halo cases produce identical error summaries,
confirming that overlap, signed charge and correctly reconstructed Halo signals
are removed consistently. With exact subtraction, Target 1 width RMSE is
0.00143/0.02388/0.23306 Å for 1/1.5/2 Å; Target 2 is
0.00764/0.00393/0.01098 Å. The broad Target benefits at 1.5 Å but not uniformly
at 2 Å, while narrow-peak tail bias remains substantial.

In the misspecified case the actual Halo has A=1.2 times the Target median,
tau=1.15 times the median and C=median+0.1; Peeling deliberately retains the
required qualified-Target median rule. Target 1 width bias becomes
+0.00240/+0.05007/+0.34183 Å, and Target 2 +0.07199/+0.07931/+0.19921 Å.
These controlled offsets measure sensitivity, not a probabilistic Halo
uncertainty model. They reinforce rejection of 2 Å and show that 1 Å fitting
does not remove all Halo-model risk. Full data are in the
[overlap/Halo CSV](post-joint-grid-overlap-statistics.csv).

### Final regression scope

Post-Joint contexts reject uncertified Joint runs. Halo reconstruction is
revalidated even when the build workspace has already transferred its snapshot;
a stale workspace cannot bypass parameter/donor/run provenance checks. Coarse
equal-radius voxels fail Local rank qualification independently of upstream
flags. Retired unclamped-node/index/weight helpers are removed from the shared
interpolation header; FirstStage interpolation itself is unchanged and its
clamp/weight regression now uses an independent closed-form reference in tests.
No temporary statistical driver or benchmark instrumentation remains in
production or test sources.

Final Release build and all 21 CTest groups pass (37.72 s); final Debug build
and all 21 groups pass (446.57 s). Both use BUILD_TESTING=ON, C++20 and the
existing dependency set. Builds use `--parallel 4`; final CTest uses
`--parallel 1`, OMP_NUM_THREADS=2, OMP_THREAD_LIMIT=4 and
VECLIB_MAXIMUM_THREADS=1. Tests with explicit OpenMP settings remain capped
at four threads. Debug's unoptimized 300-atom numerical/replay test initially
exceeded its 180 s timeout; its final measured time is 233.60 s. The numerical
CTest group now has a Debug-only 900 s timeout, with Release still 180 s;
solver budgets, assertions and qualification thresholds are unchanged.
The initial documentation links are repaired, and the Eigen restoration test
passes with the documented two-thread environment. Final logs reside under
the ignored `build/release/grid-domain-validation/` directory.

### Production decision and limits

Grid Peeling, direct own-charge evaluation, unified Local datasets and the
new grid evidence/persistence contract are adopted. Post-Joint interpolation
is fully removed. Candidate observations remain inclusive at 2 Å; the
production fit uses the validated interim 1 Å range. The 2 Å candidate fails
truth recovery, tail bias and covariance checks despite numerical solvability.
Group eligibility and its statistical kernels retain their existing rules.
No CLI experiment switch, legacy fallback, temporary study driver or new
dependency remains. FirstStage and Two-Stage retain their 1 Å semantics.

Measured remaining risks are coarse-grid rank failure, low-SNR positive-tail
selection bias, overconfident covariance under correlated observations, and
sensitivity to Halo median misspecification. Rank/SPD qualification ensures
legal numerical evidence; it does not certify width identifiability or nominal
coverage. The real-map case establishes a performance improvement for this
selection, not truth accuracy or scaling to every whole-model selection.
Joint endpoint and charge/Halo reconstruction uncertainty remain conditional
inputs. A broader covariance model or a different tail estimator requires
separate statistical work; this change does not claim to resolve either.
