# Fixed-state profile operator and PR0 resource campaign

This increment adds an internal operator and PR3 contracts. Production search
still uses `PrepareDerivative` / `ReduceDerivative` and the existing Guarded LM.
Public results, active-set decisions, reference independence, endpoint checks,
uncertainty and serialized schemas are unchanged.

## Representation and ownership

`ProfileJacobianOperator` freezes the canonical free face (every C, A > 0),
column-normalized free design Z, sparse raw width derivative D, observation
scale, and the contractions T. Each free column has exactly one width owner,
so T is an owner index plus scalar per free column, not a p-by-m matrix.

The actions are `(ProjectComplement(D*v) - PseudoInverseTranspose(T*v))/s`
and `(D.transpose()*ProjectComplement(w) - T.transpose()*LeastSquares(w))/s`.
They include the nonzero-residual correction. The factor implements both
pseudoinverse actions with the same Q, triangular factor and permutation.
Complement projection zeroes the leading free coordinates of Q-transformed
vectors; no RHS-dependent approximation or cancellation fallback is selected.

The operator creates a dedicated factor, never borrowing the mutable primary
workspace. Destroying the input evaluation or factoring a trial cannot invalidate
it. Copies refer to the same immutable linearization identity. New preparations
receive different identities even at numerically equal states. No persistent
state/face cache is introduced. Factors are serial-use (SPQR owns mutable scratch).

Both configure-time backends are supported. Eigen's new factor is used only by
the operator; its production active-set row-reduction route remains intact.
Eigen extracts the existing R and restores its column permutation for the
compact rank check rather than repeating Q-transpose actions over every column.
Construction still performs the existing compact SVD rank check using original
context rows, not informative/reduced row counts. Rank failure is unavailable,
not an approximate derivative. The p-by-p compact is transient; this is **not**
a scalable rank implementation. Persistent operator storage is sparse matrices
and factors plus O(p+m) vectors; each action uses O(N+p+m) vector workspace.
Sparse fill can still be large. No p-by-m derivative coefficients/correction,
N-by-m dense Jacobian or m-by-m reduced Jacobian is stored by the operator.

Contribution-only halo remains analytically eliminated. Callers use the existing
profile domain/context; no nuisance width coordinates are introduced and parent
scale/original rank rows remain in force.

## PR3 contract (no automatic partitioner or Schwarz solver yet)

`PreconditionerPartition` owns immutable problem/layout identities, stable block
IDs, disjoint core atoms, overlap atoms, informative rows and reverse membership.
Indices in structural blocks refer to the parent snapshot. Solver coordinates
are distinct: `SolverBlockMapping` supplies local-to-global mappings for Width
or FreeAC. FreeAC mappings must be rebuilt for a new canonical active face.
Numerical factors are never shared across coordinate spaces or with reference.

Block coordinate construction rejects missing coverage, duplicates and invalid
indices. Restriction and scatter both apply weight `1/sqrt(membership)`;
therefore an identity local action assembles to global identity. Future SPD local
factors give `sum R_k^T W_k H_k^-1 W_k R_k`. Local regularization affects only
this preconditioner, never objective, rank evidence, or marginal covariance.

`PreconditionerContext` binds linearization identity, coordinate space, positive
metric diagonal and nonnegative damping. `IdentityPreconditioner` takes an owned
copy and rejects stale identities, changed metric/damping/space and bad RHSs.
Future implementations must provide the same fixed, linear, symmetric positive
inverse action throughout one Krylov solve. Operator and preconditioner remain
separate arguments of the future solver; no preconditioner lives in J.

PR3 starts with deterministic 128-core-atom blocks and one structural overlap
ring. Numerical zeros do not remove edges. Resource exhaustion must return
partition/preconditioner unavailable rather than truncate global coupling.
Automatic partitioning, local numerical factors, rebuild policy integration,
coarse spaces and PCG/LM are intentionally outside this increment.

## Reproducible measurements

Build pristine `15e71e3fc230205e6efd7f60c8a915eb1989733b` and the candidate
separately for EIGEN and SPQR, Release, `BUILD_TESTING=ON`,
`RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=ON`. Build `rhbm_tests` (candidate),
`rhbm_gem_cli`, `joint_component_runtime`, `joint_sparse_benchmark`, and
`joint_validation`. Runtime ROOT/UMAP are unnecessary; use identical settings.

Run `tests/integration/joint_operator_validation.py` with `--stage baseline`,
then `candidate`, then `large`, finally `compare`. Each invocation takes
`--work-dir` and (except compare) `--eigen` / `--spqr` build directories.
It refuses to overwrite the stage receipt. Optional `--model` and `--map`
select the real-data control. Receipts contain source/binary/input fingerprints,
configuration and linked libraries, single-worker environment, commands,
process-tree samples, OS peak RSS and numerical statuses. Run campaigns without
concurrent builds or unrelated numerical jobs when comparing elapsed times.

Each campaign is capped at 4800 seconds; each subprocess at 600 seconds and
4 GiB sampled tree RSS. These are watchdog limits, not hard OS allocation limits;
post-exit RSS and sampling gaps are retained. Three completed fixed-state runs
produce medians. A stopped run is not a numerical verdict or a successful gate.
Comparison exits 1 on numerical mismatch and 2 for incomplete/unavailable
measurements; unavailable comparisons have `passed: null`, never true.
Inclusive phase times overlap and must not be summed. Dense shape probes report
known matrix sizes/roles and maximum bytes per matrix, not a complete allocator
trace or simultaneously live memory. Factor nnz is Eigen R nnz / the existing
SPQR bound, not a byte count for the entire factor.
Candidate `--operator` measurements retain the primary state and legacy
derivative oracle while preparing/testing the new operator. Their process RSS
therefore measures the combined audit, not isolated operator memory or a future
production search. Setup/factor cost is included; no production speedup is claimed.

The synthetic `frozen-lattice-v1` workload builds support directly, avoiding the
old fixed map box. Integer half-angstrom coordinates, 7-step atom spacing and
radius-5 spheres give 515 memberships per atom. Unique voxels are sorted z/y/x.
Chain centers advance along x; cube side is ceil(cuberoot(atom count)), filled
x, then y, then z. Truth is A=2, C=.2, B=.5; start B=.55. Noise is
`1e-5*sin(.13*z+.17*y+.19*x)` in integer coordinates. The generator source hash
plus canonical manifest hash identifies the recipe. Preparation also streams
SHA-256 over all identities, observations and memberships using length-prefixed
bytes, little-endian uint64 indices and exact double bits (`frozen-lattice-input-v1`).
Real/fixture input files have their own SHA-256 hashes.

Large runs (2k/5k/10k) only construct input, structural partition, parameter
layout, sparse basis and supplied-state residual. Their state is not a profiled
or certified fit. The driver rejects large fixed/workflow requests before work.
Preparation receipts explicitly list all unexecuted solve/evidence stages.
128/512 operator parity uses bounded oracle tiles, never a full dense Jacobian.

## Acceptance

The [measured acceptance report](figures/joint-profile-operator/acceptance.md)
contains the completed campaign results, retained receipts and regression logs.
All operator audits and preparation cases passed; the 6Z6U numerical comparison
remains unverified because both baseline and candidate reached the RSS watchdog.

Both backends must pass same-state derivative relative 1e-8, adjoint/linearity
scaled 1e-12, gradient 1e-13 + 2e-9*abs(reference), and cancellation absolute
1e-12. Rank/availability decisions use the existing thresholds. Tests include
nonzero correction, active A/signed C, permutation, original parent rank rows,
noncontiguous nuisance-row elimination, expired trial factors, snapshot lifetime,
SPD weighted overlap, identity/stale-context rejection and deterministic support.
Existing joint unit/runtime/CLI regressions retain their numerical contracts.

Passing these checks establishes the operator foundation, not 10k search or
full-workflow scalability. PR2/PR3 will integrate the solver/preconditioner;
PR4/PR5/PR6 still own factor, rank/endpoint and uncertainty scalability.
