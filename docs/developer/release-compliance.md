# Release and Compliance

This guide is for maintainers preparing source or binary releases.

## License and Third-Party Notices

- Project license: Apache License 2.0 (see `../../LICENSE`).
- Third-party licenses and attribution summary: see `../../THIRD_PARTY_NOTICES.md`.
- This repository can use system packages or FetchContent fallbacks (Eigen3/CLI11/pybind11/SQLite3/GoogleTest), depending on CMake configuration.

## Distribution Policy (Default Release Profile)

1. Release builds define `EIGEN_MPL2_ONLY` to constrain Eigen header usage to MPL-2.0-only subsets.
2. ROOT and Boost are optional and not required for non-plot workflows; release artifacts may be shipped without these linked binaries.
3. If binaries are distributed, include `LICENSE` and `THIRD_PARTY_NOTICES.md` in the package.
4. If large `.sqlite` datasets are distributed, verify their source-data license terms separately.

## Compliance Checklist

1. Include `LICENSE`.
2. Include `THIRD_PARTY_NOTICES.md`.
3. Ensure version pins and notices in `THIRD_PARTY_NOTICES.md` match the configured dependency sources.
4. Validate build variants used for release:
   - `USE_SYSTEM_LIBS=OFF`
   - `USE_SYSTEM_LIBS=ON`
   - `CMAKE_DISABLE_FIND_PACKAGE_ROOT=TRUE`
   - `CMAKE_DISABLE_FIND_PACKAGE_Boost=TRUE`
