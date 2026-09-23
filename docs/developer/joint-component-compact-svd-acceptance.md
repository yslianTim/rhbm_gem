# Joint compact SVD acceptance

Baseline: `ce58c89747d4091f91967e7f4bd9c7890d51c203` (production source
fingerprint `9718531067e72cd3ecf180d12ae8033b8e39360b7fe1d49d24ae847de4240fa0`).

Acceptance was cancelled at the user's request on 2026-09-23. Implementation
and the existing results are retained; no further acceptance runs are planned
as part of this delivery. The performance gates and Single 512 end-to-end gate
are not certified by this delivery.

Implemented: values-only derivative spectra, private automatic Jacobi/BDCSVD
selection with threshold-sensitive Jacobi retry, independent reference solves,
cost counters, compact capture/replay and a bounded measurement runner.

Completed verification: EIGEN and SPQR Release numerical/regression tests,
12 runner tests, repository guards, and a testing-disabled SPQR build with
Single 128 analysis, persistence and export. The full CTest logs include a
documentation-link failure subsequently resolved by the separate documentation
rechecks. See the [verification record](figures/joint-compact-acceptance/verification.json)
and its associated logs for the exact scope and results.

The [implementation and reproduction guide](joint-component-compact-svd.md)
defines the numerical controls, unchanged policies and separate compact and
complete-command gates. Partial measurement outputs remain locally under
`build/joint-compact-measurements-final`; these are not a completed acceptance
package and must not be used to claim that the required performance gates passed.
