# Joint component memory planning

`ResourceWork` preflights Joint problem construction, structural construction,
FixedNeighbor preparation, global streaming evaluation, endpoint verification,
local conditioning, and parameter stability. It uses checked arithmetic and
accounts for resident bytes, support entries, affected rows, block dimensions,
and the largest local workspace.

The estimated peak budget is 3.2 GB and the hard estimate limit is 4.0 GB.
These are planner limits, not measured RSS guarantees. Construction and
FixedNeighbor preparation have separate checks and may reject a large input
before fitting.

Production evaluates the accepted state and endpoint through streaming support
traversal. It does not allocate a global free-design factor, projection
coefficient matrix, or projected-width factor. Local conditioning and stability
workspaces are bounded by block dimensions.

The large synthetic measurements did not complete a full search, local
conditioning and stability over all blocks, or downstream workflow. The full
physical-data resource gate remains open. The
[Joint evidence record](joint-component-evidence.md) retains the measured
figures and their scope.
