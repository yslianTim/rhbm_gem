# Joint Component Assessment Memory Census

This report records the current reference-free assessment preflight and the
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
estimate; target correction conservatively budgets all 74,812 A/C columns.
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
| Projected-width assessment | 37,406 × 37,406 | 574,853,787 | 60,878,414,528 | 61,453,268,315 | ResourceLimit; NotRun |
| Local correction | 37,406 × 37,406 | 574,853,787 | 60,878,414,528 | 61,453,268,315 | ResourceLimit; NotRun |
| Target geometry | 112,218 × 112,218 | 477,151,131 | 316,944,928,128 | 317,422,079,259 | ResourceLimit; NotRun |
| Target local correction | 37,406 × 37,406 | 628,044,131 | 94,459,426,592 | 95,087,470,723 | ResourceLimit; NotRun |

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

Projected-width and local-correction estimates include the free-column by
width arrays, square width factors/right vectors, sparse derivative copies,
and 8,192-row tiled workspace. For this tested state, each stage estimates
33,581,012,064 bytes for the three free-column × width arrays,
22,387,341,376 bytes for the square factor and right-vector pair,
4,902,879,232 bytes for tiled workspace, and 7,181,856 bytes for sparse
copies. The full stage peak including resident data is 61.453 GB.

After weak directions are copied, the projected-width right-vector matrix is
released before local correction. Target geometry similarly releases its
full right-vector matrix and column scales after target leakage evidence is
copied, before target-local correction begins. The endpoint evaluation remains
resident for target-local correction because that stage consumes it. Its
estimate therefore includes endpoint sparse matrices and vectors, but does not
double count the released geometry null-space matrix.

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
| Primary endpoint evaluation | 0.023681 s | 474,431,488 → 478,773,248 bytes | 478,773,248 bytes | Entered; completed |
| Independent streaming replay | 0.054572 s | 478,773,248 → 622,592,000 bytes | 622,592,000 bytes | Entered; completed |
| Rank assessment | 0.004313 s | 622,592,000 → 622,592,000 bytes | 622,592,000 bytes | Preflight rejected |
| Target memory preflight calculation | 0.005581 s | 622,608,384 → 622,624,768 bytes | 622,624,768 bytes | Preflight rejected |
| DesignSpectrum factorization | NotRun | Not sampled | Not allocated | ResourceLimit |
| Projected-width factorization | NotRun | Not sampled | Not allocated | ResourceLimit |
| Local-correction factorization | NotRun | Not sampled | Not allocated | ResourceLimit |
| Target geometry / target-local correction | NotRun | Not sampled | Not allocated | ResourceLimit |

The final process peak was 622,624,768 bytes (0.623 GB), about
316,799,454,491 bytes below the largest estimate. That gap is expected here:
the 317.422 GB target estimate and the 189.484 GB rank estimate were rejected,
so those allocations never happened. The observed RSS is evidence for the
input, preparation, streaming verification and preflight path only; it cannot
validate the predicted RSS of an unexecuted dense stage.

Runtime statuses were: endpoint verification `Passed`; full-design rank
`Failed` by the structural upper bound; dense design-rank preflight
`ResourceLimit`; projected-width and local-correction stages `NotRun`;
numerical identifiability `Unavailable`; target-identifiability preflight
`ResourceLimit`; target geometry and target-local correction `NotRun`; and
`RuntimeConvergence` `Unavailable`. The rank preflight recorded zero production
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

The next research step is to evaluate ways to bound the dense rank and target
geometry stages while preserving the current numerical contracts. No bounded
rank method is introduced here.
