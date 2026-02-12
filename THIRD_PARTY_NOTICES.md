# Third-Party Notices

This project is licensed under Apache License 2.0. This file summarizes third-party components that may be used when building or testing the project.

## Scope and Dependency Selection

- The build option `USE_SYSTEM_LIBS=ON` prefers system packages for Eigen3, CLI11, pybind11, and SQLite3.
- If a system package is unavailable, the build falls back to bundled copies in `third_party/`.
- Test dependencies (for example GoogleTest) are only used when `BUILD_TESTING=ON`.

## Bundled Components

| Component | Usage in this project | License | Local path | Upstream license file(s) |
|---|---|---|---|---|
| CLI11 | Command-line parsing | BSD-3-Clause | `third_party/CLI11` | `third_party/CLI11/LICENSE` |
| Eigen | Linear algebra templates | MPL-2.0 (primary), with additional notices in upstream package | `third_party/eigen` | `third_party/eigen/COPYING.README`, `third_party/eigen/COPYING.MPL2`, `third_party/eigen/COPYING.BSD`, `third_party/eigen/COPYING.LGPL`, `third_party/eigen/COPYING.GPL`, `third_party/eigen/COPYING.APACHE`, `third_party/eigen/COPYING.MINPACK` |
| pybind11 | Python bindings | BSD-3-Clause | `third_party/pybind11` | `third_party/pybind11/LICENSE` |
| SQLite (amalgamation) | Embedded database backend fallback | Public Domain (upstream declaration) | `third_party/sqlite` | Declaration in `third_party/sqlite/sqlite3.h` and `third_party/sqlite/sqlite3ext.h` |
| GoogleTest | Unit tests only | BSD-3-Clause | `third_party/googletest` | `third_party/googletest/LICENSE` |

## Optional System Dependencies (Not Bundled)

- ROOT (optional plotting features): discovered via `find_package(ROOT ...)` and not bundled in this repository.
- OpenMP runtime (optional parallelism): provided by compiler toolchain or system packages and not bundled in this repository.

## Redistribution Notes

- Keep this file and the main `LICENSE` file in source and binary distributions.
- Preserve copyright and license notices required by third-party licenses.
- If shipping binaries, include attribution/notice documents in the package and/or accompanying documentation.
- If distributing large data artifacts (for example `.sqlite` datasets), verify data-source-specific license terms separately from software license compliance.

