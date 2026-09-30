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
- [Benchmark guide](joint-benchmark.md) documents current benchmark profiles
  and measurement semantics.
- [Capabilities and limitations](joint-capabilities-limitations.md) identifies
  supported behavior and current research boundaries.
- [Historical evidence](joint-component-evidence.md) summarizes retained
  comparisons, counterexamples, fixture ownership, and artifact retrieval.
- [Certification contract](joint_abc_certification_contract.md),
  [component contract](joint_abc_components_contract.md), and
  [compact-SVD acceptance](joint-component-compact-svd-acceptance.md) retain
  independent mathematical and numerical contracts.
- [Endpoint-reference acceptance](joint-component-reference-acceptance.md)
  preserves its distinct historical comparison results and provenance.
- [Joint v1 acceptance](joint-component-v1-acceptance.md) records the
  integration and metadata contract.

## Testing and release

- [Test organization](../../tests/README.md) describes default correctness
  tests, optional Joint validation, benchmarks, and research tools.
- Read [release compliance](release-compliance.md) before preparing a source
  or binary release.
