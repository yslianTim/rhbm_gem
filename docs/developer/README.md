# Developer Documentation

Use these guides when changing the codebase, validating a build, or working
with project internals.

## Build and repository workflow

- [Build and configuration](build-and-configuration.md) describes dependencies,
  active CMake options, build modes, and validation commands.
- [Development guidelines](development-guidelines.md) covers engineering rules,
  test labels, generated documentation, and repository checks.
- [Command architecture](architecture/command-architecture.md),
  [object architecture](architecture/object-architecture.md), and
  [data-object I/O architecture](architecture/dataobject-io-architecture.md)
  describe the main implementation boundaries.
- [Adding data-object operations](adding-dataobject-operations.md) and
  [adding a command](adding-a-command.md) provide implementation checklists.
- [Command developer notes](commands/README.md) indexes command-specific
  documentation.

## Production fitting

- [First-stage fitting](estimate-local-gaussian-with-offset.md) describes local
  Gaussian initialization.
- [Production fitting](production-fitting.md) owns per-atom MDPDE solves,
  failed-only refinement, covariance, numerical failure handling, and solver
  replay.
- [Second-stage local fitting](second-stage-local-fitting.md) owns outer
  iteration, candidate acceptance, recovery, convergence, and finalization.

## Joint Component estimator

- [Runtime guide](joint-component-runtime.md) is the production API and
  regression authority.
- [Stage integration](joint-stage-integration.md) describes integration with
  the fitting workflow.
- [Capabilities and limitations](joint-capabilities-limitations.md) identifies
  supported behavior and current evidence limits.
- [Width recovery closeout](joint-width-closeout.md) records recovery integration,
  full-model Width sensitivity, and downstream eligibility evidence.
- [Target point validation](joint-target-estimability.md) defines v2 quality and
  schema 22 persistence; the [cleanup record](joint-diagnostics-cleanup.md)
  records numerical equivalence and the deliberate quality migration.
- [Memory planning](joint-component-memory-census.md) documents production
  preflight estimates and their limits.
- [Historical evidence](joint-component-evidence.md) summarizes retired
  decisions and scientific limits. Closed implementations and machine-readable
  results are recoverable from Git history and are not retained in the tree.
- [Certification contract](joint_abc_certification_contract.md),
  [component contract](joint_abc_components_contract.md), and the
  [compact-SVD implementation contract](joint-component-compact-svd.md) own
  the maintained mathematical and numerical rules.

## Testing and release

- [Test organization](../../tests/README.md) describes default correctness
  tests, Joint validation and production benchmarks.
- Read [release compliance](release-compliance.md) before preparing a source
  or binary release.
