# Third-Party Notices

This project is licensed under Apache License 2.0. This file summarizes third-party components that may be used when building or testing the project.

## Scope and Dependency Selection

- The build option `RHBM_GEM_DEP_PROVIDER=SYSTEM` requires system packages for Eigen3, CLI11, pybind11, SQLite3, and (when tests are enabled) GoogleTest.
- The build option `RHBM_GEM_DEP_PROVIDER=FETCH` uses CMake FetchContent for Eigen3, CLI11, pybind11, SQLite3, and GoogleTest from pinned release sources with SHA256 checksums.
- Test dependencies (for example GoogleTest) are only used when `BUILD_TESTING=ON`.

## Components

| Component | Usage in this project | License | Source mode | Upstream license file(s) |
|---|---|---|---|---|
| CLI11 | Command-line parsing | BSD-3-Clause | `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent | `https://github.com/CLIUtils/CLI11/blob/v2.5.0/LICENSE` |
| Eigen | Linear algebra templates | MPL-2.0 (primary), with additional notices in upstream package | `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent | `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.README`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.MPL2`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.BSD`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.LGPL`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.GPL`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.APACHE`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.MINPACK` |
| pybind11 | Python bindings | BSD-3-Clause | `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent | `https://github.com/pybind/pybind11/blob/v3.0.2/LICENSE` |
| SQLite (amalgamation) | Embedded database backend fallback | Public Domain (upstream declaration) | `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent (`sqlite-amalgamation-3490100`) | `https://www.sqlite.org/copyright.html` |
| GoogleTest | Unit tests only | BSD-3-Clause | `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent | `https://github.com/google/googletest/blob/v1.17.0/LICENSE` |

## Optional System Dependencies (Not Bundled)

- ROOT (optional plotting features): discovered via `find_package(ROOT ...)` and not bundled in this repository.
- OpenMP runtime (optional parallelism): provided by compiler toolchain or system packages and not bundled in this repository.

## Redistribution Notes

- Keep this file and the main `LICENSE` file in source and binary distributions.
- Preserve copyright and license notices required by third-party licenses.
- If shipping binaries, include attribution/notice documents in the package and/or accompanying documentation.
- If distributing large data artifacts (for example `.sqlite` datasets), verify data-source-specific license terms separately from software license compliance.
