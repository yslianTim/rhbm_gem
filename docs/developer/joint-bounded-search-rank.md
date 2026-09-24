# Bounded search revalidation and free-design rank prototype

This increment follows `5f61bb6e`. Production still uses `LegacyCompact` and
operator search still uses the dense rank backend. No A/C, damping, PCG stopping,
independent reference, endpoint, assembly, or uncertainty policy is changed.
See [the validation results](joint-bounded-search-rank-validation.md) for outcomes
and retained limitations.

## Fixed-step and search evidence

The fixed-state comparator and campaign now return
0 for complete numerical comparisons, 1 for an observed numerical failure,
3 for incomplete evidence, and 2 for malformed input or invocation. Failure has
precedence over incompleteness; both flags remain in the JSON. These codes do
not promote search or turn historical phase subtotals into wall-clock timings.
The fixed comparator accepts a requested subset of preconditioners for the new
bounded runner; its historical default still requires all three.

Each new fixed-step measurement owns a fresh operator and independently times
preparation, gradient, metric, applicable partition/local factors, PCG, and
prediction. `fixed_step_wall_seconds` is one continuous interval through the
prediction, excluding supplied-state reconstruction, separate parity actions,
serialization and teardown. The separate audit factor is released before these
intervals. `total_seconds` remains the sum of the listed phases. Inverse/action
subtimers are nested in the solve and must not be added again.

The bounded experiment compared legacy and operator methods built from the same
source, with serial execution, one numerical thread and a 4 GiB process-tree
RSS limit. Its total budget was 60 minutes, split into fixed-step, search and
rank stages without transfers or extensions. Resource stops prevented further
repetitions or larger cases in the affected group.

Qualified endpoints must retain the existing rank, availability and evidence
contracts, scaled A/C/B tolerance 1e-10 and normalized objective tolerance 1e-12.
Equivalent unavailability or unconverged non-regression is reported separately
and cannot qualify a performance comparison. The 512 gate requires three
qualified paired repetitions and Schwarz search median <= 1.10 times legacy.
Search, assessment and full process time are separate. The report has independent
fixed/search/rank statuses, and `promote` is always false. This bounded schedule
does not establish Schwarz superiority over diagonal or satisfy all historical
promotion requirements.

## SPQR rank prototype

`EvaluateFreeDesignRank` takes the normalized sparse design, an optional immutable
factor, the existing `RankRequest`, and `RankBudget`. The factor may be absent
for structural checks. Its result distinguishes `FullRank`, `Deficient` and
`Unavailable`, with rank bounds, spectral/threshold bounds, the verification
error, weak-direction bound, work count, workspace bound and elapsed time.
Exact rank is reported only when established. Original problem dimensions,
absolute overrides and `SvdNative`/`StrictGreater` policies remain at the caller.

The internal operator constructor can explicitly select `SpqrBounds`; its
default is `Dense`. Only offline/fixed-state tests select the prototype. Eigen
has no prototype factor view and returns unavailable for nonstructural cases.
The dense `EvaluateRank` remains the oracle and is never a hidden fallback.

For the public SPQR representation `Z P = Q [R; 0] + E`:

1. Sparse column, row and Frobenius norms enclose the largest singular value,
   hence the original rank threshold. The absolute-override enclosure includes
   the dense backend's division/multiplication rounding; unresolved subnormal
   cutoffs are unavailable.
2. Exact zero/duplicate columns and deficient dimensions provide structural
   bounds. Nonzero problems at a zero threshold remain unavailable rather than
   promising agreement with roundoff-sensitive SVD equality. The zero matrix
   and an insufficient row dimension have direct rank bounds.
3. Up to three smallest-pivot triangular directions are verified using the
   original sparse design. Only a direct-action upper bound strictly below the
   threshold lower bound establishes numerical deficiency. Pivot size alone
   never decides rank.
4. Let `M(R)` have diagonal `abs(R_ii)` and off-diagonal `-abs(R_ij)`. Nonnegative
   triangular solves with ones bound `||R^-1||_infinity` and `||R^-1||_1`.
   Their geometric mean bounds `||R^-1||_2`. An already insufficient comparison
   bound returns `rank-bound-too-wide` without an expensive futile verification.
5. For each stored reflector `I - tau h h'`, enclose its exceptional eigenvalue
   `1 - tau ||h||^2`. The product of the smaller of its absolute lower bound and
   one bounds `sigma_min(Q)` below; permutations preserve it.
6. Apply the stored reflectors to each R column using outward-rounded intervals,
   compare to the corresponding Z column, and accumulate a Frobenius upper bound
   `delta` on E. Only one observation column is retained. Full rank requires
   `q_lower / inverse_norm_upper - delta > threshold_upper`, including the native
   SVD positive-minimum guard. Equality, overlap, overflow and budget exhaustion
   remain unavailable.

The inverse bound follows [Higham's triangular comparison-matrix bounds](https://nhigham.com/2021/03/30/bounds-for-the-norm-of-the-inverse-of-a-triangular-matrix/).
Each elementary bound is expanded with `nextafter`; IEEE binary arithmetic and
round-to-nearest are required, and fast-math is rejected. No guessed safety
factor or stochastic estimate is used as a certificate.

Defaults are 120 seconds, 100 million charged sparse-entry work units and a conservative
`48*n + 192*p` byte workspace bound capped at 256 MiB. The latter excludes the
input design and immutable factor; the process watchdog includes their storage,
construction scratch and allocator overhead. Sparse exported-array bytes are
reported separately. Shape probes are known allocations, not an allocator trace.

The prototype allocates no global dense compact, normal matrix, inverse or
singular-vector matrix. It does retain sparse QR and its fill-in. Streaming
verification can cost O(p * nnz(H) + n*p), and comparison bounds can be very
conservative. This is a memory-representation/certificate prototype, not a claim
of linear runtime, complete rank recovery, or 2k/5k/10k workflow scalability.
