# Retired fixed-action investigation

Status: **Closed**

The fixed-state factor, normal-action and operator-action experiments were
part of the retired OperatorPcg research surface. They are historical evidence,
not a current implementation or benchmark profile.

The current production Pure Exclusion estimator is bounded FixedNeighbor with
Active-set LM. Fixed-Halo compatibility retains the Profile LM
`SearchProfile()` path. The existing trust, global replay, stationarity,
endpoint-certification and `RuntimeConvergence` contracts remain in force.
Shared numerical components remain where the bounded or compatibility paths
still use them.

The historical fixed-state comparison used two Q/Q′ calls for normal actions,
compared with six in the composed control. That result did not establish a
different full-search default or large-workflow scalability. The conclusion
is summarized in [joint-component-evidence.md](joint-component-evidence.md);
the experiment-specific reports and figures are retired.
