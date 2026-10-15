# FixedNeighbor local profile recovery

## Policy

A block with no trusted accepted local profile gets one retry on the next
complete sweep. The first failed visit leaves the component parameters,
prediction, residual, and objective unchanged for that block; the remaining
blocks still run in forward Gauss-Seidel order. At the retry visit the local
effective response is rebuilt from the current component state, including any
neighbor updates accepted later in the previous sweep.

The retry candidate must pass the existing local constrained A/C feasibility
and KKT checks, local trust replay, global objective replay, and monotonicity
gate before it can update the component state. A successful retry is recorded
only after those checks pass. If the retry again has no trusted local profile,
the search stops with `block-recovery-exhausted`. If the sweep budget ends
while a retry is pending, the reason is `block-recovery-sweep-budget`.

A sweep containing a failed block is not eligible for first-order or confirmed
stationarity. Retaining the accepted state or exhausting the retry budget does
not qualify a component. Any final component state still needs the unchanged
global KKT `1e-10`, width-gradient `1e-12`, confirmed eta-change `1e-10`, and
independent endpoint verification checks.

## Relationship to precision recovery

The retry is a bounded response to later Gauss-Seidel neighbor updates. It does
not substitute old A/C values for a failed local solution and does not retry
the same local response with a relaxed KKT threshold. P2's Block 21 evidence
does not establish a multiprecision candidate that can be converted to double
and independently certified, so no precision fallback is wired into this
path. An unresolved local numerical failure terminates after the single
retry, preserving diagnostics and the last accepted component state.

The tests inject a transient no-profile result to verify that later blocks can
change the state before the fresh-response retry, and inject a persistent
failure to verify bounded exhaustion, no stationarity promotion, and state
preservation.
