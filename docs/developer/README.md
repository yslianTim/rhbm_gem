# Developer Documentation

Start here if you are changing the codebase, validating build configurations, or working with project internals.

## Recommended Reading Order

1. Read [`docs/developer/build-and-configuration.md`](/docs/developer/build-and-configuration.md) for CMake parameters, dependency selection, coverage, and feature-mode validation commands.
2. Read [`docs/developer/development-guidelines.md`](/docs/developer/development-guidelines.md) for repository-wide engineering rules, test/label expectations, command-manifest sync requirements, and quality-check alignment (`lint_repo`, formatter/tidy checks).
3. Read [`docs/developer/architecture/command-architecture.md`](/docs/developer/architecture/command-architecture.md) and [`docs/developer/architecture/dataobject-io-architecture.md`](/docs/developer/architecture/dataobject-io-architecture.md) when you need architecture context for commands, data object traversal policies, typed dispatch, I/O, or persistence.
4. Read [`docs/developer/adding-dataobject-operations-and-iteration.md`](/docs/developer/adding-dataobject-operations-and-iteration.md) when you need implementation checklists for extending operations/iteration on existing `DataObject` types.
5. Read [`docs/developer/adding-a-command.md`](/docs/developer/adding-a-command.md) when you need a concrete implementation template for a new command.
6. Browse [`docs/developer/commands/README.md`](/docs/developer/commands/README.md) for command-specific developer notes generated alongside individual commands.
7. Read [`docs/developer/release-compliance.md`](/docs/developer/release-compliance.md) before preparing source or binary releases.
