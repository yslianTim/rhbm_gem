# Documentation

Use this directory as the entry point for the project documentation.
The content is split by audience so you can go straight to the guide that matches your task.

## Directory Ownership

- Keep narrative documentation, onboarding, architecture notes, and release/process guidance under `docs/`.
- Keep runnable examples, reusable templates, and helper scripts under [`../resources/`](../resources/README.md).
- Do not move repository-reading guides into `resources/`; that tree is reserved for assets that users or tooling execute directly.

## User Documentation

- Set up your environment and install RHBM-GEM.
- Verify the Python bindings.
- Run the example workflows.
- Start with [`user/README.md`](user/README.md) and [`user/getting-started.md`](user/getting-started.md).
- Browse shared examples and helper resources in [`../resources/README.md`](../resources/README.md).

## Developer Documentation

- Change the codebase and follow repository conventions.
- Validate build configurations and dependency choices.
- Read architecture notes and release guidance.
- Start with [`developer/README.md`](developer/README.md) and [`developer/build-and-configuration.md`](developer/build-and-configuration.md).
- [`developer/development-guidelines.md`](developer/development-guidelines.md)
- [`developer/architecture/command-architecture.md`](developer/architecture/command-architecture.md)
- [`developer/architecture/dataobject-io-architecture.md`](developer/architecture/dataobject-io-architecture.md)
- [`developer/adding-dataobject-operations-and-iteration.md`](developer/adding-dataobject-operations-and-iteration.md)
- [`developer/adding-a-command.md`](developer/adding-a-command.md)
- [`developer/commands/README.md`](developer/commands/README.md)
