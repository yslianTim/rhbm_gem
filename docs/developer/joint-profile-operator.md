# Fixed-state profile operator

This guide owns the current tiled derivative, profile-Jacobian operator, and
structural partition contracts. Historical operator measurements, incomplete
campaigns, and retrieval data belong in the
[canonical historical evidence](joint-component-evidence.md).
Production search selection and solver gates are documented in
[joint-operator-search.md](joint-operator-search.md).

The profile operator is internal. Production continues to use the current
Guarded LM path unless the internal search policy explicitly selects the
operator path. Public results, active-set decisions, reference independence,
endpoint checks, uncertainty, and serialized schemas are shared.

## Tiled derivative and assessment

The derivative and assessment paths process observations in 8,192-row tiles.
They form temporary tile-by-width work, reduce with compact QR factors, and do
not materialize the full observation-by-width dense Jacobian. The full residual
norm, including its orthogonal tail, remains in the objective, actual-reduction,
and stopping calculations. Tile boundaries change storage scheduling, not the
mathematical Jacobian or residual norm.

The immutable observation/support snapshot and original parent rank context
remain attached through component views. Informative-row profiling does not
replace the parent row count used by rank policy. Atom and row identities remain
in parent coordinates, and no component edge is cut. Dense derivatives exist
only in test support as a parity reference.

Derivative/LM temporary dense workspace is O(tile * atoms + atoms^2), in
addition to sparse designs, sparse-factor storage and O(rows) residual vectors.
This does not bound sparse fill-in or parent-global assembled assessment as the
number of atoms grows. No dense runtime fallback is used.

## Profile Jacobian operator

ProfileJacobianOperator freezes the canonical free face (every C and A > 0),
column-normalized free design Z, sparse raw width derivative D, observation
scale, and contractions T. Each free column has exactly one width owner, so T
is an owner index plus scalar per free column, not a p-by-m matrix.

The actions are

    Jv = (ProjectComplement(D*v) - PseudoInverseTranspose(T*v))/s
    J'w = (D.transpose()*ProjectComplement(w) - T.transpose()*LeastSquares(w))/s

They include the nonzero-residual correction. The factor implements both
pseudoinverse actions with the same Q, triangular factor and permutation.
Complement projection zeroes the leading free coordinates of Q-transformed
vectors. No RHS-dependent approximation or cancellation fallback is selected
by the operator.

The operator creates a dedicated factor and never borrows the mutable primary
workspace. Destroying the input evaluation or factoring a trial cannot
invalidate it. Copies refer to the same immutable linearization identity; new
preparations receive different identities even at numerically equal states.
There is no persistent state/face cache. Factors are serial-use because SPQR
owns mutable scratch.

Both configure-time backends are supported. Eigen's operator factor is used
only by the operator; its production active-set row-reduction route remains
intact. Eigen extracts the existing R and restores its column permutation for
the compact rank check rather than repeating Q-transpose actions over every
column. Construction uses original context rows, not informative/reduced row
counts. Rank failure is unavailable, not an approximate derivative. The
p-by-p compact is transient; this is not a scalable rank implementation.
Persistent operator storage is sparse matrices and factors plus O(p+m)
vectors; each action uses O(N+p+m) vector workspace. Sparse fill can still be
large. The operator stores no p-by-m derivative coefficients/correction,
N-by-m dense Jacobian, or m-by-m reduced Jacobian.

Contribution-only halos remain analytically eliminated. The operator uses the
existing profile domain/context, introduces no nuisance width coordinates, and
preserves parent scale and original rank rows.

## Structural partition and preconditioner coordinates

PreconditionerPartition owns immutable problem/layout identities, stable block
IDs, disjoint core atoms, overlap atoms, informative rows, and reverse
membership. Indices in structural blocks refer to the parent snapshot. Solver
coordinates are separate: SolverBlockMapping supplies local-to-global mappings
for Width or FreeAC. FreeAC mappings are rebuilt for a new canonical active
face. Numerical factors are never shared across coordinate spaces or with the
reference.

Block construction rejects missing coverage, duplicates, and invalid indices.
Restriction and scatter both apply weight 1/sqrt(membership), so an identity
local action assembles to global identity. The current additive SPD local action
is sum R_k' W_k H_k^-1 W_k R_k. Local regularization affects only this
preconditioner, never the objective, rank evidence, or marginal covariance.

PreconditionerContext binds linearization identity, coordinate space, positive
metric diagonal, and nonnegative damping. IdentityPreconditioner owns a copy
and rejects stale identities, changed metric/damping/space, and invalid right
hand sides. Any inverse action held through one Krylov solve must be fixed,
linear, symmetric and positive definite. Operator and preconditioner are
separate solver arguments.

The current partition uses deterministic breadth-first core groups of at most
128 atoms and a complete one-ring structural overlap. Numerical zeros do not
remove edges. It does not materialize clique adjacency; row-to-atom CSR
incidence and visited marks bound topology workspace. Resource exhaustion
returns unavailable rather than truncating overlap or global coupling.

Current limits are 512 atoms per block, 512 MiB held topology/model/factor
storage, and 256 MiB conservative construction scratch. These bounds exclude
the shared problem, global evaluation, and profile factor; process RSS remains
a separate resource boundary. Shape probes are known allocations, not an
allocator trace.

## Current validation and measurement

Permanent operator tests cover Jv, J-transpose-w, adjoint and linearity,
normal-action and dense parity, nonzero correction, active A/signed C,
permutation, parent rank rows, nuisance-row elimination, factor lifetime,
weighted overlap, stale-context rejection, and deterministic support.
Same-state numerical controls retain the current derivative, adjoint, gradient,
cancellation, rank, and availability tolerances described by the
[runtime guide](joint-component-runtime.md).

Use the unified benchmark profiles for current measurements:

    python3 tests/integration/joint_benchmark.py \
      --profile prepare --case chain-8 --build-dir build/joint-eigen \
      --output build/operator-prepare.json
    python3 tests/integration/joint_benchmark.py \
      --profile fixed --case chain-8 --build-dir build/joint-eigen \
      --output build/operator-fixed.json

These controls establish the operator and partition contracts. They do not
establish 10,000-atom search, full-workflow scalability, or a complete 6Z6U
numerical outcome; those historical boundaries are summarized in the
[evidence index](joint-component-evidence.md).
