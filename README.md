# Robust Hierarchical Bayesian Model for Gaussian Estimation of cryo-EM Maps (RHBM-GEM)
RHBM-GEM is a robust hierarchical Bayesian framework for Gaussian modeling and parameter estimation of atom-associated cryo-EM map distributions.



## Getting Started

User-facing setup, build, install, Python, coverage, and CMake configuration details are maintained in [`docs/getting-started.md`](docs/getting-started.md).

- [`docs/getting-started.md`](docs/getting-started.md) for platform-specific installation and build instructions.
- [`examples/python/00_quickstart.py`](examples/python/00_quickstart.py) for a minimal Python workflow once the project is installed.
- [`examples/python/01_end_to_end_from_test_data.py`](examples/python/01_end_to_end_from_test_data.py) for an end-to-end example using test data.

## Verbosity Levels

The command line option `-v, --verbose` controls how much information is printed.
The logger supports five levels:

| Level |  Name   | Description                       |
|-------|---------|-----------------------------------|
| 0     | Error   | Only error messages are shown.    |
| 1     | Warning | Warnings and errors are shown.    |
| 2     | Notice  | Notice messages are shown.        |
| 3     | Info    | Informational messages are shown. |
| 4     | Debug   | Detailed debug output.            |

Passing a larger value to `--verbose` enables all lower levels as well.



## License and Third-Party Notices

- Project license: Apache License 2.0 (see `LICENSE`).
- Third-party licenses and attribution summary: see `THIRD_PARTY_NOTICES.md`.
- This repository can use either system packages or bundled third-party dependencies under `third_party/` depending on CMake configuration.

**Distribution policy (default release profile)**
1. Release builds define `EIGEN_MPL2_ONLY` to constrain Eigen header usage to MPL-2.0-only subsets.
2. ROOT and Boost are optional and not required for non-plot workflows; release artifacts may be shipped without these linked binaries.
3. If binaries are distributed, include `LICENSE` and `THIRD_PARTY_NOTICES.md` in the package.
4. If large `.sqlite` datasets are distributed, verify their source-data license terms separately.

**Compliance checklist (source and binary release)**
1. Include `LICENSE`.
2. Include `THIRD_PARTY_NOTICES.md`.
3. Ensure bundled third-party license file paths referenced in `THIRD_PARTY_NOTICES.md` remain valid.
4. Validate build variants used for release:
   - `USE_SYSTEM_LIBS=OFF`
   - `USE_SYSTEM_LIBS=ON`
   - `CMAKE_DISABLE_FIND_PACKAGE_ROOT=TRUE`
   - `CMAKE_DISABLE_FIND_PACKAGE_Boost=TRUE`


# Developer Documentation

Repository-internal design, extension rules, and developer-facing architecture notes are maintained under [`docs/`](docs/).

- [`docs/development-guidelines.md`](docs/development-guidelines.md) for repository-wide engineering rules and change expectations.
- [`docs/architecture/command-architecture.md`](docs/architecture/command-architecture.md) for command registration, execution flow, and extension patterns.
- [`docs/architecture/dataobject-io-architecture.md`](docs/architecture/dataobject-io-architecture.md) for DataObject file I/O, SQLite persistence, and related constraints.
