# Operator Guarded LM and overlapping Schwarz (PR2 / PR3)

Production still defaults to `SearchMethod::LegacyCompact`. Operator search is
available through the internal `EvaluationContext::search` and
`FitWithSearchPolicy`; no installed API or serialized schema was extended.
Promotion requires both backend regression and benchmark gates below. Search
completion is independent of the returned state's convergence evidence.

## Search and ownership

The new path solves `(J'J + mu Gamma'Gamma) step = -J'(r/s)` using the full PR1
operator, including residual correction. It never constructs a global normal
matrix, reduced Jacobian or dense LM QR. PR1's transient p-by-p rank SVD remains;
A/C, reference, assessment, assembly and uncertainty retain their existing paths.

Each accepted state owns a dedicated operator factor. Trial evaluation uses a
separate mutable workspace; rejection retains the accepted linearization.
Accepting a trial rebuilds state-dependent numerical data. Only structural block
topology survives accepted updates. There is no state/face cache.

Raw-width column norms `d_j=||D_j||/s` provide a deliberately approximate search
metric. Positive floors are `sqrt(epsilon)*max(d)` (all-zero norms use ones), and
accepted updates take coordinatewise maxima. All three preconditioner controls
use the same metric. This approximation changes trajectory, not objective,
KKT, rank or endpoint stationarity thresholds.

PCG starts from zero and verifies the true metric-scaled residual at relative
1e-10, every 32 iterations and before success. A residual replacement restarts
the recurrence. The iteration limit is min(4m,1000); internal test overrides may
only reduce it. A fixed damping, metric, operator and inverse action are used
throughout each solve. Curvature, nonfinite action, stale context and exhaustion
have explicit failure reasons, never a dense or diagonal fallback.
The fixed SPD inverse follows the [PCG requirements](https://www.netlib.org/linalg/old_html_templates/subsubsection2.6.3.1.1.html).

The initial radius is 0.1*||Gamma eta||, or 0.1 if zero. Initial damping is 1e-3,
floored at 1e-12. Bracket search returns a step inside the radius, with at most 20
damping solves per accepted linearization; internal overrides can only reduce
this limit. Invalid/untrusted trials quarter the
radius and quadruple damping. Acceptance uses the full trial objective and
`-g.dot(step)-0.5*||J step||^2`, requires ratio >=1e-4, and passes scalar replay.
Radius adaptation uses 0.25/0.75 ratio boundaries. Reduction and step tolerances
remain 1e-14 and 1e-12; gradient stopping uses the full gradient infinity norm
1e-12. These are search stops, not endpoint evidence. Profile/update budgets
remain 200/100. Search performs zero reference solves. Endpoint certification
and newest-first fallback still check the saved coefficients independently.

## Structural and numerical preconditioners

Partitioning uses informative frozen-support incidence, not numerical nonzeros.
Stable atom-identity ordering defines BFS core groups of at most 128 atoms.
Every core receives a complete one-ring structural overlap. Parent atom/row
indices are retained; width and canonical free-A/C mappings are separate.
The partition owns the immutable snapshot. Legacy hand-built numerical fixtures
without snapshot identities get a bounded copied profile snapshot for testing.
The production indexed-domain path retains the original snapshot and mappings.

No clique adjacency is materialized. Row-to-atom CSR incidence and visited marks
bound topology workspace. These blocks never enter the exact `SolveLinear`
block-decomposition path.

Schwarz approximates local Gauss-Newton curvature after local A/C elimination.
It uses the global state's beta, active face and parent scale, without local
nonlinear fits. For normalized local free design Z and
`E=D*Gamma_local^-1/s`, it computes sparse products
`A=Z'Z`, `B=Z'E`, `C=E'E`, then `S=C-B'(A+lambda I)^-1 B`.
`lambda=sqrt(epsilon)*max(1,||A||_infinity)` affects only this approximation.
It intentionally does not replace the full global residual-corrected Jacobian.

The symmetrized S is factored with damping and a local numerical shift:
`S + (mu+tau) I`. Tau starts at
`64*epsilon*max(1,||S||_infinity,mu)` and may grow tenfold for at most six LLT
attempts. Local rank deficiency therefore does not become a global rank verdict.
Maximum lambda/tau and build counts are recorded. Resource-audit mode also retains
each local lambda and every shifted-factor tau attempt, tagged by local build,
factor build and stable block ordinal. Damping changes refactor only
shifted matrices; metric/state/face changes require new local numerical data.

Both restriction and scatter use weights `1/sqrt(coordinate membership)`.
Local solves are additionally scaled on both sides by `Gamma_local^-1`.
The resulting additive inverse action is fixed and SPD during a PCG solve.
This uses symmetric restriction/scatter, consistent with the basic additive
form distinguished from restricted Schwarz in [PETSc's PCASM documentation](https://petsc.org/main/manualpages/PC/PCASMType/).
Identity and raw-width diagonal inverses provide controls. Local inverses are
never reused as reference factors or covariance.

Limits are 512 atoms per block, 512 MiB held topology/model/factor storage and
256 MiB conservative construction scratch. Bounds exclude the shared problem,
global evaluation and PR1 factor, which remain subject to the process watchdog.
Limit failures report unavailable without truncating overlap or changing the
solver. Search retains its last accepted state for normal endpoint certification.
Shape probes are known allocations, not an exhaustive allocator trace. No coarse
level, iterative A/C solver or LSMR backend is implemented here.

## Reproduce validation

See the [closed validation report](joint-operator-search-validation.md) for
completed measurements, user-cancelled work and the decision to retain the
legacy production default.

Build EIGEN and SPQR Release variants with BUILD_TESTING=ON,
RHBM_GEM_DEP_PROVIDER=SYSTEM, RHBM_GEM_ENABLE_UMAP=OFF,
RHBM_GEM_ROOT_MODE=OFF and RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON.
Targets: rhbm_tests, joint_sparse_benchmark, joint_component_runtime,
joint_validation and rhbm_gem_cli.

The baseline is a git archive of
`c6c869cd4f4939104dfde96984ae17f736a8ce4a`. Its production sources must remain
unchanged. Copy the candidate `tests/experiments/joint_sparse_benchmark.cpp`
into that archive and add `PR23_BASELINE_DRIVER` to this target's compile
definitions in tests/CMakeLists.txt. This compiles the same search/assessment
measurement wrapper against pristine numerical code. Baseline mode only accepts
legacy search; its production source hash is verified against git archive.
Driver, generator and runner hashes are separately retained.

Use `tests/integration/joint_search_validation.py` with a shared `--work-dir`
and `--input-dir`. Run stages baseline (pristine --eigen/--spqr build directories),
eigen and spqr (candidate build directories), large-local (candidate), then
compare. The first four invocations are the four campaigns; each receipt refuses
replacement. Baseline generates missing single-128/512 inputs with the unchanged
existing generator and initializer; fixtures are unpacked with catalogue hashes.
Optional --model/--map select the 6Z6U control. Reusing prior frozen inputs is
permitted, but baseline/candidate input fingerprints must match exactly.

All measurements use one numerical thread, a 600-second pipeline limit, 4 GiB
sampled process-tree RSS, and an 80-minute campaign deadline. Sampling may
overshoot; process and OS peak values are separate. There is no automatic restart
or budget increase. Three independent completed runs produce small-case medians;
a resource stop terminates that case's repetitions. Remaining campaign entries
are explicitly not-run-budget after deadline.

`joint_sparse_benchmark ... --search legacy|identity|diagonal|schwarz --resources`
measures search followed by unchanged endpoint certification. Search costs are
saved before assessment, along with PCG/operator and preconditioner counters.
`synthetic chain|cube N workflow ... --search schwarz` additionally exercises the
internal full runtime path, capture and uncertainty. Public runtime/CLI regression
continues to exercise the production default.

`synthetic chain|cube N local ... --resources` is the only new large-case mode.
It generates supplied beta/eta, basis and residual, then builds/applies local
preconditioners. It never calls global profile solve, rank, operator, search,
reference, assessment or uncertainty. Requests for large fixed/workflow runs are
rejected. Input hashes, block counts, SPD/symmetry checks and zero forbidden-work
counters are retained for 128,512,2k,5k,10k in both arrangements and backends.

## Acceptance and promotion

Core tests cover the damped dense oracle, inverse symmetry/SPD, permutation,
coverage, parent/nuisance mappings, local rank regularization, resource failure,
stale contexts, budgets, absence of global derivative reduction, replay-only
search, and full runtime uncertainty/persistence. PR1 tests retain cancellation,
active-face and factor lifetime checks. Both backend regressions must pass.

Search trajectories may differ. Returned rank, availability and necessary evidence
must match; qualified converged A/C/B endpoints use scaled 1e-10 and normalized
objectives absolute 1e-12. Uncertified states are not passing endpoints. Explicit
matching initial rank failures are recorded as equivalent unavailability.

For each backend, single-512 Schwarz search median must be <=1.10 times pristine
legacy. Schwarz must also beat diagonal, including every build/rebuild cost, in
at least one common overlap case on both backends. A reduction in Krylov iterations
alone is insufficient. Incomplete comparisons never promote. The runner reports
benchmark eligibility but never edits the production default. Evidence checks
and the final default decision must be recorded explicitly.

6Z6U is an unchanged public-workflow resource control, including in candidate
campaigns while legacy remains default. Prior RSS stops provide no numerical
verdict and do not establish operator search parity on this dataset. This
increment does not claim 10k search or full-workflow scalability.
