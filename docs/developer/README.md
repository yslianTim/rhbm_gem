# Developer Documentation

Start here if you are changing the codebase, validating build configurations, or working with project internals.

## Recommended Reading Order

1. Read [`build-and-configuration.md`](build-and-configuration.md) for CMake parameters, dependency selection, coverage, and feature-mode validation commands.
2. Read [`development-guidelines.md`](development-guidelines.md) for repository-wide engineering rules, test/label expectations, command-manifest sync requirements, and quality-check alignment (`lint_repo`, formatter/tidy checks).
3. Read [`architecture/command-architecture.md`](architecture/command-architecture.md) and [`architecture/dataobject-io-architecture.md`](architecture/dataobject-io-architecture.md) when you need architecture context for commands, data object traversal policies, typed dispatch, I/O, or persistence.
4. Read [`adding-dataobject-operations-and-iteration.md`](adding-dataobject-operations-and-iteration.md) when you need implementation checklists for extending operations/iteration on existing `DataObject` types.
5. Read [`adding-a-command.md`](adding-a-command.md) when you need a concrete implementation template for a new command.
6. Read [`release-compliance.md`](release-compliance.md) before preparing source or binary releases.
