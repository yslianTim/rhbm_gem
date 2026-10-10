# Retired fixed-action investigation

Status: **Closed**

The fixed-state factor, normal-action and operator-action experiments were
part of the retired OperatorPcg research surface. They are historical evidence,
not a current implementation or benchmark profile.

The current Joint estimator is FixedNeighbor. Each structural block invokes
the existing local Profile LM `SearchProfile()` operation, then
uses the existing trust, global replay, stationarity, endpoint-certification
and `RuntimeConvergence` contracts. Shared numerical components remain when
FixedNeighbor or its local profile search still calls them.

The historical fixed-state comparison used two Q/Q′ calls for normal actions,
compared with six in the composed control. That result did not establish a
different full-search default or large-workflow scalability. The conclusion
is summarized in [joint-component-evidence.md](joint-component-evidence.md);
the experiment-specific reports and figures are retired.
