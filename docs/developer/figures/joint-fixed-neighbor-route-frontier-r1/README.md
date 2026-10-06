# Matched Joint route frontier

This campaign measures LegacyCompact, OperatorPcg, and FixedNeighbor in separate processes with the same frozen synthetic workload, initial widths, SPQR backend, one Eigen thread, resource monitor, and 4 GiB RSS limit. Search-only is the primary comparison. Full endpoint runs are limited to 512 and 768 atoms. Formal runs use a 600.0-second cap; diagnostic reruns use a 7200.0-second wall cap with the same RSS limit; retry policy: deferred.

P1 route position: **OperatorPcg-dominant**.
