# Third-Party Notices

RHBM-GEM
Copyright 2024-2026 Yu-Hsiang Lien

Unless otherwise noted, the original source code, build scripts, example
programs, and documentation in this repository are licensed under the MIT
License. See `LICENSE`.

This file summarizes third-party components that may be used when building,
testing, installing, or redistributing the project. It is provided for
attribution and compliance guidance and does not modify the license terms of
this project or of any third-party component.

## Scope and Dependency Selection

- The build option `RHBM_GEM_DEP_PROVIDER=SYSTEM` requires system packages for Eigen3, CLI11, SQLite3, and Boost. It also requires `pybind11` plus a Python 3 interpreter/development module when `BUILD_PYTHON_BINDINGS=ON`, and GoogleTest when `BUILD_TESTING=ON`.
- The build option `RHBM_GEM_DEP_PROVIDER=FETCH` uses CMake `FetchContent` with pinned release sources and SHA256 checksums for Eigen3, CLI11, SQLite3, and Boost. It additionally fetches `pybind11` when `BUILD_PYTHON_BINDINGS=ON` and GoogleTest when `BUILD_TESTING=ON`.
- Release builds define `EIGEN_MPL2_ONLY` to constrain Eigen header usage to MPL-2.0-only subsets.
- In `FETCH` mode, the install/export flow redistributes fetched Eigen, CLI11, Boost, and SQLite artifacts as part of the installed package.

## Components

| Component | Usage in this project | License | Source mode | Upstream license file(s) |
|---|---|---|---|---|
| Boost | Math/statistics helper headers (`boost::math`) | Boost Software License 1.0 | Required in both modes: `SYSTEM` mode uses system package discovery; `FETCH` mode provides pinned header fallback | `https://github.com/boostorg/boost/blob/boost-1.90.0/LICENSE_1_0.txt` |
| CLI11 | Command-line parsing | BSD-3-Clause | `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent | `https://github.com/CLIUtils/CLI11/blob/v2.5.0/LICENSE` |
| Eigen | Linear algebra templates | MPL-2.0 (primary), with additional notices in upstream package | `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent; release builds define `EIGEN_MPL2_ONLY` | `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.README`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.MPL2`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.BSD`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.LGPL`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.GPL`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.APACHE`, `https://gitlab.com/libeigen/eigen/-/blob/5.0.1/COPYING.MINPACK` |
| pybind11 | Python bindings | BSD-3-Clause | Used only when `BUILD_PYTHON_BINDINGS=ON`: `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent | `https://github.com/pybind/pybind11/blob/v3.0.2/LICENSE` |
| SQLite (amalgamation) | Database backend: system SQLite library in `SYSTEM` mode; embedded amalgamation (`sqlite3.c`/`sqlite3.h`) in `FETCH` mode | Public Domain (upstream declaration) | `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent (`sqlite-amalgamation-3490100`) | `https://www.sqlite.org/copyright.html` |
| GoogleTest | Unit tests only | BSD-3-Clause | Used only when `BUILD_TESTING=ON`: `SYSTEM` provider uses system package; `FETCH` provider uses pinned FetchContent | `https://github.com/google/googletest/blob/v1.17.0/LICENSE` |

## Optional System Dependencies (Not Bundled)

- ROOT (optional plotting features): discovered via `find_package(ROOT ...)` and not bundled in this repository.
- OpenMP runtime (optional parallelism): provided by compiler toolchain or system packages and not bundled in this repository.
- Python 3 interpreter and development module (required only when `BUILD_PYTHON_BINDINGS=ON`): discovered via `find_package(Python ...)` and not bundled in this repository.

## Redistribution Notes

- Keep `LICENSE` and this file in source and binary distributions.
- In `FETCH` installs or binary packages, account for the fetched Eigen, CLI11, Boost, and SQLite artifacts redistributed by the install/export rules.
- Preserve copyright and license notices required by third-party licenses.
- If shipping binaries, include attribution/notice documents and the applicable upstream third-party license texts in the package and/or accompanying documentation.
- If enabling optional system integrations such as ROOT, OpenMP, or Python, verify the selected local packages' license terms separately because they are not bundled by this repository.
- If distributing large data artifacts (for example `.sqlite` datasets), verify data-source-specific license terms separately from software license compliance.
