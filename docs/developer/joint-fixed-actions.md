# Retired fixed-action investigation

Status: **Closed / Removed**

The fixed-state factor, normal-action and operator-action experiments were
part of the retired OperatorPcg research surface. They are historical evidence,
not a current implementation or benchmark profile.

The current Joint estimator is FixedNeighbor. Each structural block invokes
the existing local `SearchProfile()` / `LegacyCompact` profile solver, then
uses the existing trust, global replay, stationarity, endpoint-certification
and `RuntimeConvergence` contracts. Shared numerical components remain when
FixedNeighbor or its local profile search still calls them.

Historical fixed-action parity, factor-residency and preparation results remain
available through the compact reports and figures indexed by
[`joint-component-evidence.md`](joint-component-evidence.md). They are not
recreated by the current tree.
