# FixedNeighbor bounded local-work qualification

This P4 campaign runs the P3-selected OneAcceptedLocalUpdate policy at chain-1024 and cube-1024 in search-only mode. FullLocalSearch values are the existing matched search-only frontier from `joint-fixed-neighbor-scaling-r1`; the candidate uses the same SPQR backend, 128-atom core, forward order, one Eigen thread, and frozen stationarity checks.

`local_work_seconds` is the comparable legacy local profile-work field. The candidate also reports the newer `profile_factor_seconds` telemetry; the pre-P3 Full baseline does not contain that field. No full endpoint claim is made for these 1024 points; endpoint correctness was gated separately in P3.

Qualification gate: **passed**. Selected policy: **OneAccepted**.
