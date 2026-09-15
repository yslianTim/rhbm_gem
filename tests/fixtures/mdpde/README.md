# Captured MDPDE experiment inputs

These version-1 fixtures were captured from actual fitting calls using the fixed
fold-168 map at baseline `d169e070`, with testing-only instrumentation. They are
not reconstructed from terminal diagnostics. `index.json` records identities,
roles, hashes and the hashes of the full state/background contexts retained with
the local experiment evidence.

The three shape fixtures contain X, log response, alpha, execution options and
the exact expected beta/variance/status. The offset fixture contains the original
sparse system (including explicit zeros), ridge, anchor and exact expected result.
Replay uses no model, map or simulation manifest.

Run the `mdpde_experiment solve` executable with this directory, or set
`RHBM_TEST_REPLAY_DIR` to it when running
`ProductionFittingTest.ReplaysCapturedSolverFailureWithoutSimulationInputs`.
These inputs record failures of the unchanged production solver; improved
experimental roots do not replace the expected production outputs.
