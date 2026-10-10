# FixedNeighbor workspace attribution

The repeated workspace-residency study is closed. It compared persistent
prepared-block workspaces with fresh block-local workspaces on chain and cube
cases from 256 through 1024 atoms. Both paths used the same historical
core-64, Forward, 30-sweep search policy.

Persistent residency produced no material wall-time improvement. It increased
sampled peak RSS by 15.6%–34.8% in the repeated cases. Exact numeric-factor
reuse opportunities were limited to initial-profile evaluations and covered
1.17%–15.00% of the measured inputs; the treatment still performed fresh
numeric factorization.

Derivative preparation and reduction accounted for roughly 82%–86% of search
time in the measured cases. The result did not justify a workspace-policy
switch, numeric-factor cache, or candidate-replay optimization. Production
semantics remain the frozen core-12 FixedNeighbor route with local Profile LM.

This concise conclusion is the retained record. Closed campaign drivers and
machine-readable outputs are recoverable from Git history and are intentionally
not retained in the current tree. Current measurements are documented in the
[Joint benchmark guide](joint-benchmark.md).
