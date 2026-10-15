# Joint component memory census

## Current production planning

`ResourceWork` keeps stage-specific preflights for Joint problem construction,
structural construction, FixedNeighbor preparation, global streaming
state/objective evaluation, endpoint verification, reference solves, local
conditioning, and parameter stability. Estimates use checked arithmetic,
resident bytes, support entries, affected rows, block dimensions, and the
largest local workspace. The planner budget is 3,200,000,000 estimated peak
bytes; the hard limit is 4,000,000,000 bytes. Both values are limits for
estimates, not observed RSS limits.

Production no longer creates a global free-design factor, a global projection
coefficient matrix, or a global projected-width factor. The accepted search and
endpoint checks use streaming support traversal. Local conditioning and
parameter perturbation allocate workspace from block dimensions. Construction
and FixedNeighbor preparation have separate preflights and can still reject a
large input before fitting.

## Synthetic large-component measurement

The retained deterministic shape has 5,938,087 observation rows, 37,406
contributors, and 3,118 blocks, with 74,811 support memberships. A prior
no-sweep preparation run reported 2,892,736 prepared-mapping bytes and an
876,119,666-byte global-evaluation estimated peak. The whole test process
recorded 474,546,176 bytes peak RSS (0.475 GB, 0.442 GiB).

That run disabled sweeps and endpoint verification. It is a synthetic
preparation and streaming-memory control. It does not measure a physical-data
search, local conditioning on an actual component, perturbation stability, or
the physical support graph. Its estimated peak and observed RSS are different
measurements and should not be compared as if they described the same
allocation.

The current `ProductionScaleSingleComponentPreparationStaysBelowBudget`
regression was rerun with 5,938,087 rows, 37,406 contributors, 74,811 support
entries, and 3,118 blocks. It reported a 972,989,100-byte preparation peak
estimate and a 625,051,683-byte endpoint-verification peak estimate. The whole
test process measured 622,460,928 bytes peak RSS and completed in 263 ms. The
synthetic test performs a no-sweep preparation check and a zero-coefficient
streaming endpoint check; it does not fit the component or execute local
conditioning and stability over all blocks. No physical input was used.

The preparation test retains actual support and affected-row counts and
checks the planner before expensive allocation. A synthetic passing preflight
does not claim that a physical-data fit completed.

The point-diagnostics regression separately exercises 37,406 point identities
and 3,118 blocks. It compares component-local with assembled-global evidence
and verifies exact linear work counts: `4P` membership visits while building
both indexes and `3P` atom lookups for conditioning, stability, and selected
targets. That fixture has one observation row; it isolates the point-mapping
cost and is not a 5,938,087-row fit or a peak-RSS measurement. The test prints
the two construction times for reference, while operation counts remain the
deterministic complexity check.

## Retired global-rank estimates

The estimates below describe the deleted production rank path and are retained
only to explain the architectural change. They are not current preflights or
allocations.

- The old free-design and projected-width shape at 37,406 columns estimated
  67,162,024,128 bytes of BDC-SVD workspace and 78,355,694,816 bytes at the
  isolated SVD peak before retained sparse inputs and projection coefficients.
- The retired dense `DesignSpectrum` and combined projected-width planner
  estimated 189,484,131,163 bytes and 39,065,926,939 bytes, respectively.
- An earlier target-geometry planner estimated 317,422,079,259 bytes and the
  global covariance preflight estimated 418,399,065,443 bytes. Those stages
  were not part of the current production estimator and were removed.

No current production result depends on these rank estimates, singular values,
threshold decisions, or global weak directions. Local SparseQR rank checks and
CompactSvd utilities remain where the block solver and profile derivative use
them.

## Interpretation

Memory estimates are planning evidence; they are not process RSS or guarantees
about operating-system behavior. The prior RSS values are whole-process
high-water marks, not allocations attributed to one stage. No real physical
large-component fit was run for this census. Construction, preparation,
endpoint verification, local conditioning, and stability each require their
own measured result before making physical-data resource claims.

The historical global-rank experiments and their derivations are archived in
[bounded global rank research](joint-bounded-rank-research.md).
