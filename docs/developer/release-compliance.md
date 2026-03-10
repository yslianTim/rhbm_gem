# Release and Compliance

This guide is for maintainers preparing source or binary releases.

## License and Third-Party Notices

- Project license: Apache License 2.0 (see `../../LICENSE`).
- Third-party licenses and attribution summary: see `../../THIRD_PARTY_NOTICES.md`.
- This repository supports two dependency providers selected by `RHBM_GEM_DEP_PROVIDER`: `SYSTEM` or `FETCH`.

## Distribution Policy (Default Release Profile)

1. Release builds define `EIGEN_MPL2_ONLY` to constrain Eigen header usage to MPL-2.0-only subsets.
2. ROOT and Boost are optional and not required for non-plot workflows; release artifacts may be shipped without these linked binaries.
3. If binaries are distributed, include `LICENSE` and `THIRD_PARTY_NOTICES.md` in the package.
4. If large `.sqlite` datasets are distributed, verify their source-data license terms separately.

## Legacy v1 Migration Lifecycle

- `RHBM_GEM_LEGACY_V1_SUPPORT` is `ON` by default to preserve migration compatibility.
- Releases may set `RHBM_GEM_LEGACY_V1_SUPPORT=OFF` only after explicit deprecation notice.
- Recommended retirement gate: observe no production need for legacy migration across at least 2 consecutive release cycles, then remove migration code and fixtures in a major update.

## Compliance Checklist

1. Include `LICENSE`.
2. Include `THIRD_PARTY_NOTICES.md`.
3. Ensure version pins and notices in `THIRD_PARTY_NOTICES.md` match the configured dependency sources.
4. Validate build variants used for release:
   - `RHBM_GEM_DEP_PROVIDER=FETCH`
   - `RHBM_GEM_DEP_PROVIDER=SYSTEM`
   - `CMAKE_DISABLE_FIND_PACKAGE_ROOT=TRUE`
   - `CMAKE_DISABLE_FIND_PACKAGE_Boost=TRUE`
