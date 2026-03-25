# Release and Compliance

This guide is for maintainers preparing source or binary releases.

## License and Third-Party Notices

- Project license: Apache License 2.0 (see `../../LICENSE`).
- Project notice file and copyright statement: see `../../NOTICE`.
- Third-party licenses and attribution summary: see `../../THIRD_PARTY_NOTICES.md`.
- This repository supports two dependency providers selected by `RHBM_GEM_DEP_PROVIDER`: `SYSTEM` or `FETCH`.

## Distribution Policy (Default Release Profile)

1. Release builds define `EIGEN_MPL2_ONLY` to constrain Eigen header usage to MPL-2.0-only subsets.
2. ROOT remains optional, but Boost support is required in all release builds (via system Boost or bundled FETCH fallback headers).
3. If binaries are distributed, include `LICENSE`, `NOTICE`, and `THIRD_PARTY_NOTICES.md` in the package.
4. If `FETCH` dependencies are redistributed, include the applicable upstream third-party license texts and notices with the package.
5. If large `.sqlite` datasets are distributed, verify their source-data license terms separately.

## Compliance Checklist

1. Include `LICENSE`.
2. Include `NOTICE`.
3. Include `THIRD_PARTY_NOTICES.md`.
4. Ensure version pins and notices in `THIRD_PARTY_NOTICES.md` match the configured dependency sources.
5. For `RHBM_GEM_DEP_PROVIDER=FETCH` releases, include the applicable upstream third-party license texts for redistributed artifacts.
6. Validate build variants used for release:
   - `RHBM_GEM_DEP_PROVIDER=FETCH`
   - `RHBM_GEM_DEP_PROVIDER=SYSTEM`
   - `CMAKE_DISABLE_FIND_PACKAGE_ROOT=TRUE`
