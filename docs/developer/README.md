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

- [Production fitting](production-fitting.md) owns per-atom MDPDE solves,
  failed-only refinement, covariance, numerical failure handling, and solver
  replay.
- [Two-Stage local fitting](second-stage-local-fitting.md) owns outer
  iteration, candidate acceptance, recovery, convergence, and finalization.

## Joint Component Estimator

- [Runtime](joint-component-runtime.md) defines the Target model, observation
  domains, BCD updates and global qualification.
- [Numerics](joint-component-numerics.md) explains constrained A/C, Variable
  Projection and the bounded LM solve.
- [Stage Integration](joint-stage-integration.md) follows FirstStage
  initialization through Target publication, Halo reconstruction, peeling,
  Local and Group fitting, and persistence.

## Testing and release

- [Test organization](../../tests/README.md) describes default correctness
  tests and Joint validation.
- Read [release compliance](release-compliance.md) before preparing a source
  or binary release.
