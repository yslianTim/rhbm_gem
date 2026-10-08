# Optimized FixedNeighbor matched route frontier

This P6 campaign reuses the P1 LegacyCompact, OperatorPcg, and original Full/128-core FixedNeighbor reports and runs the qualified OneAccepted/64-core FixedNeighbor policy on the same chain and cube cases. All measurements are search-only with one Eigen thread, the SPQR backend, and a 4 GiB RSS limit. The P1 records are linked by the manifest; they are not rerun in this campaign.

Formal envelope: 600.0 seconds. Diagnostic envelope: 7200.0 seconds.

Route position: **OperatorPcg-dominant**.

Promotion gate: **deferred**. OperatorPcg remained fastest and passed the formal envelope at every matched point; the historical OneAccepted/core64 FixedNeighbor route had lower RSS, but this campaign did not show a formal resource-survival case that warrants an internal-route promotion. The compact evidence remains; raw per-case telemetry was pruned after the route was superseded by the qualified core12 policy.
