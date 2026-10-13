# Joint Component Assessment Memory Census

This report records the current correction-free assessment preflight and the
one-way structural rank obstruction for the production-scale synthetic
single-component shape. The structural check does not change numerical
thresholds or produce an exact numerical rank.

## Input and structural shape

| Measure | Value |
| --- | ---: |
| Observations | 5,938,087 |
| Contributors / full atoms | 37,406 |
| Targets | 37,406 |
| Halo atoms | 0 |
| Connected components | 1 |
| Support memberships | 74,811 |
| Prepared blocks | 3,118 |
| Total affected rows | 40,523 |
| Total block support entries | 74,811 |
| Maximum block rows | 13 |
| Maximum block support entries | 24 |
| Design columns | 74,812 |
| Free A/C columns at the tested state | 37,406 |
| Actual full-design structural nonzeros | 149,622 |
| Unique rows supported by the full design | 37,406 |
| Full-design structural rank upper bound | 37,406 |
| Full-design structural matching | NotRun; row-count obstruction is sufficient |
| Free-face structural matching | 37,406 of 37,406 columns |
| Width columns | 37,406 |
| Target columns | 112,218 |
| Full Jacobian shape | 5,938,087 × 112,218 |

The fixture is deliberately sparse: each atom contributes one or two support
memberships. It is not a physical cryo-EM support graph. The tested state has
zero amplitudes, so only the 37,406 charge columns are free in the rank-shape
estimate; target geometry contains 112,218 columns for the complete A/C/log-B
parameter vector.
At the tested widths, each support membership contributes a nonzero Gaussian
and charge design entry. The unique supported-row count is therefore 37,406,
not the 40,523 sum of affected-row entries. Since 37,406 is below the 74,812
full-design columns, `rank(X) <= 37,406 < 74,812` proves full-design rank
deficiency. The free face has 37,406 columns and 37,406 supported rows, so
that row-count inequality does not prove free-face deficiency; no full-design
conclusion is transferred to the free face. A separate sparse free-face graph
has a perfect structural matching, which still does not certify its numerical
rank or conditioning.

## Estimated memory

The planning budget is 3,200,000,000 bytes (3.2 GB); the hard limit is
4,000,000,000 bytes (4 GB). The stage estimates retain their own resident data
and temporary allocations. The parent rank and target estimates use the
maximum sequential stage peak, rather than adding every stage's temporary
storage.

| Stage | Factor / matrix shape | Resident bytes | Temporary bytes | Estimated peak bytes | Outcome |
| --- | ---: | ---: | ---: | ---: | --- |
| Input snapshot | — | 330,746,835 | — | 330,746,835 | measured estimate |
| FixedNeighbor preparation | — | — | — | 910,595,445 | Pass |
| Global state evaluation | — | — | — | 972,989,100 | Pass |
| Prepared row mappings | — | — | — | 2,892,736 mapping bytes | prepared |
| Endpoint verification | streaming | — | — | 625,051,683 | Pass |
| Structural design check | 5,938,087 rows × 74,812 columns | 574,853,787 | 3,665,740 | 578,519,527 | Pass; full-design deficiency proven |
| DesignSpectrum | 74,812 × 74,812 | 574,853,787 | 188,909,277,376 | 189,484,131,163 | ResourceLimit; NotRun |
| Projected-width assessment | 37,406 × 37,406 | 574,853,787 | 38,491,073,152 | 39,065,926,939 | ResourceLimit; NotRun |
| Target geometry | 112,218 × 112,218 | 477,151,131 | 316,944,928,128 | 317,422,079,259 | ResourceLimit; NotRun |

The structural check uses support traversal and compact row/column IDs; its
preflight is O(P + E + U), where U is the unique supported-row count, and includes pair storage,
supported-row IDs and per-column marks. Matching is skipped when the
supported-row count already proves deficiency. The bound is an upper bound,
not an exact numerical rank; no singular values or condition number are
reported for this result. Target geometry is
the largest estimated stage. Its full Jacobian factor is estimated at
100,743,036,192 bytes (100.743 GB); the right-vector/null-space storage is
included in the stage estimate. The `right_null_space_bytes` field describes a
subset of `right_singular_vector_bytes`; do not add both fields together.

DesignSpectrum uses a dense `ReferenceQR` rank factorization. Its 74,812-square
factor is 44,774,682,752 bytes (44.775 GB). The dense QR workspace is
144,129,806,720 bytes (144.130 GB), and sparse design copies are 4,787,904
bytes. The resulting 189.484 GB stage estimate is rejected before the factor
or workspace is allocated. Production global reference-solve calls remain
zero; this does not mean DesignSpectrum would allocate no dense QR when its
preflight passes.

Projected-width assessment includes one free-column × width coefficient
array, square width factors/right vectors, sparse derivative copies, and
8,192-row tiled workspace. For this tested shape it estimates
11,193,670,688 bytes for that array, 22,387,341,376 bytes for the square
factor/right-vector pair, 4,902,879,232 bytes for tiled workspace, and
7,181,856 bytes for sparse copies. The full stage peak including resident data
is 39.066 GB. Qualification no longer allocates the duplicated projected
Jacobian or correction workspaces.

The prior LocalCorrection stage estimate was 61,453,268,315 bytes and the
prior TargetLocalCorrection estimate was 95,087,470,723 bytes; both production
stages and their factor/step workspaces are removed. The projected-width stage
remains necessary and drops from 61,453,268,315 to 39,065,926,939 estimated
bytes. The rank aggregate still peaks at DesignSpectrum (189,484,131,163
bytes), while target geometry remains the largest target stage (317,422,079,259
bytes). These are planner estimates, not equivalent Peak RSS reductions:
preflight rejected the dense stages before they allocated their modeled
workspaces.

Target geometry is the only target-identifiability stage. Its full Jacobian
factor is estimated at 100,743,036,192 bytes (100.743 GB); right-vector and
null-space storage are included in the 317.422 GB stage peak. The endpoint
evaluation remains resident. The `right_null_space_bytes` field is a subset of
`right_singular_vector_bytes`; do not add both fields together.

The original aggregate planner regression values are unchanged and asserted
exactly: 189,288,489,472 bytes for the frozen rank-planner control and
317,385,746,272 bytes for the frozen target-planner control. Arithmetic
overflow is checked; overflow remains distinct from a budget rejection.

## Timing, RSS, and runtime status

The following data came from an isolated serial Release test process on macOS.
`Current RSS` is sampled before and after the stage. `Process peak RSS` is the
process-wide high-water mark, not a stage-local allocation.

| Stage | Elapsed | Current RSS before → after | Process peak RSS | Status |
| --- | ---: | ---: | ---: | --- |
| Primary endpoint evaluation | 0.016791 s | 474,726,400 → 479,019,008 bytes | 479,019,008 bytes | Entered; completed |
| Independent streaming replay | 0.049350 s | 479,019,008 → 622,542,848 bytes | 622,542,848 bytes | Entered; completed |
| Rank assessment and structural check | 0.003648 s | 622,542,848 → 624,967,680 bytes | 624,967,680 bytes | Structural deficiency proven; dense preflight rejected |
| Repeated structural check on retained fixture | 0.003758 s | 675,872,768 → 678,363,136 bytes | 678,363,136 bytes | Upper bound confirmed |
| Target memory preflight calculation | 0.006283 s | 678,395,904 → 678,412,288 bytes | 678,412,288 bytes | Preflight rejected |
| DesignSpectrum factorization | NotRun | Not sampled | Not allocated | ResourceLimit |
| Projected-width factorization | NotRun | Not sampled | Not allocated | ResourceLimit |
| Target geometry factorization | NotRun | Not sampled | Not allocated | ResourceLimit |

The final process peak was 678,412,288 bytes (0.678 GB), about
316,743,666,971 bytes below the largest estimate. That gap is expected here:
the 317.422 GB target estimate and the 189.484 GB rank estimate were rejected,
so those allocations never happened. The observed RSS is evidence for the
input, preparation, streaming verification and preflight path only; it cannot
validate the predicted RSS of an unexecuted dense stage.

Runtime statuses were: endpoint verification `Passed`; full-design rank
`Failed` by the structural upper bound; the rank assessment stopped at the
structural obstruction; projected-width rank and numerical identifiability
were `Unavailable`; target-geometry preflight was `ResourceLimit`; target
geometry was `NotRun`; and `RuntimeConvergence` was `Unavailable`. The rank
preflight recorded zero production
reference solves and zero `ReferenceQR` factorization calls because it stopped
before DesignSpectrum. A separate small `chain-8` benchmark records zero
production reference solves and one `ReferenceQR` call, demonstrating that the
counters describe different work.

## Scope and limitations

The planner models explicit sparse storage, dense factors, right vectors and
tiled workspaces with checked arithmetic. It cannot measure allocator, Eigen,
BLAS or operating-system overhead for a stage that safety preflight forbids
from running. The estimates are resource guards, not measured allocations or
proof that a physical input fits under 4 GB. No user-provided physical model or
map was found in the workspace; only repository test fixtures are present.
**Physical input not measured.**

RSS functions return unavailable on platforms without a supported reader
rather than returning zero. macOS and Linux readers are centralized in
`ResourceWork`; Windows uses the unavailable path. The Windows/MSVC build was
not run in this environment.

## Experimental SparseQR measurements

Natural and COLAMD ran in isolated serial processes on the synthetic
chain/cube designs. The conservative SparseQR preflight includes a dense
worst-case `N × P` workspace term and 5% for Eigen and allocator overhead.
The completed 512 cases stayed under the 3.2 GB planning budget, but their
observed process peaks were higher than `DesignSpectrum`'s estimated peak on
the same shape. That comparison is an observed SparseQR peak against an
unexecuted dense-rank estimate, not a direct factorization-to-factorization
measurement.

| Case | Planned peak, Natural / COLAMD | Observed peak RSS, Natural / COLAMD | Outcome |
| --- | ---: | ---: | --- |
| chain-128 | 177,972,838 / 177,923,686 bytes | 152,043,520 / 151,994,368 bytes | Completed |
| cube-128 | 154,999,603 / 154,934,067 bytes | 143,278,080 / 143,278,080 bytes | Completed |
| chain-256 | 590,659,584 / 590,675,968 bytes | 491,241,472 / 435,994,624 bytes | Completed |
| cube-256 | 490,420,428 / 490,453,196 bytes | 470,007,808 / 444,055,552 bytes | Completed |
| chain-512 | 2,199,014,604 / 2,199,047,372 bytes | 1,392,148,480 / 1,308,246,016 bytes | Completed |
| cube-512 | 1,767,156,940 / 1,767,124,172 bytes | 1,166,622,720 / 1,164,247,040 bytes | Completed |
| chain-1024 | 8,539,393,228 / 8,539,458,764 bytes | 99,139,584 / 99,205,120 bytes | ResourceLimit; factorization NotRun |
| cube-1024 | 6,744,236,851 / 6,744,269,619 bytes | 95,059,968 / 95,092,736 bytes | ResourceLimit; factorization NotRun |

The 1,024-atom shapes were rejected before `analyzePattern()` or
`factorize()`. Those RSS readings cover synthetic input construction and the
preflight only; they are not SparseQR factor peaks. The fixed 37,406-atom
shape has 5,938,087 rows, 74,812 columns and 149,622 actual nonzeros. The
production-fixture SparseQR preflight estimated 4,154,773,326,379 temporary
bytes and a 4,155,449,133,611-byte peak with its measured resident state; its
process peak RSS remained 675,823,616 bytes because factorization was `NotRun`.
A separate 600,000,000-byte resident-control preflight estimates a
4,155,373,326,379-byte peak after the 5% margin. The production fixture
separately proves full-design deficiency from its 37,406 supported rows
before any rank factorization.

Full per-ordering rank, fill, factor storage and timing data are in the
[bounded global rank research report](joint-bounded-rank-research.md). The
next stage remains memory-bounded projected-width assessment and target
identifiability; those algorithms are unchanged here.
