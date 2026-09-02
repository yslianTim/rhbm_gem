# Release and Compliance

This guide is for maintainers preparing source or binary releases.

## License and Third-Party Notices

- Project-owned source code, build scripts, examples, and documentation are licensed under the MIT License (see [`LICENSE`](/LICENSE)).
- Third-party licenses and attribution summary: see [`THIRD_PARTY_NOTICES.md`](/THIRD_PARTY_NOTICES.md).
- This repository supports two dependency providers selected by `RHBM_GEM_DEP_PROVIDER`: `SYSTEM` or `FETCH`.

## Distribution Policy (Default Release Profile)

1. Release builds define `EIGEN_MPL2_ONLY` to constrain Eigen header usage to MPL-2.0-only subsets.
2. ROOT remains optional, but Boost support is required in all release builds (via system Boost or bundled FETCH fallback headers).
3. The `umap_embedding` command is enabled by default but remains optional. A fresh `SYSTEM` build strictly requires installed core packages, prefers the installed UMAP package stack, and fetches only missing UMAP components at fixed versions. The standard presets use `FETCH`; set `RHBM_GEM_ENABLE_UMAP=OFF` for a build without UMAP lookup or downloads. Fetched umappp targets and headers are build-only and excluded from the installed RHBM-GEM package.
4. If binaries are distributed, include [`LICENSE`](/LICENSE) and [`THIRD_PARTY_NOTICES.md`](/THIRD_PARTY_NOTICES.md) in the package.
5. If `FETCH` dependencies are redistributed, include the applicable upstream third-party license texts and notices with the package.
6. If large `.sqlite` datasets are distributed, verify their source-data license terms separately.

## Compliance Checklist

1. Include [`LICENSE`](/LICENSE).
2. Include [`THIRD_PARTY_NOTICES.md`](/THIRD_PARTY_NOTICES.md).
3. Ensure version pins and notices in [`THIRD_PARTY_NOTICES.md`](/THIRD_PARTY_NOTICES.md) match the configured dependency sources.
4. For `RHBM_GEM_DEP_PROVIDER=FETCH` releases, include the applicable upstream third-party license texts for redistributed artifacts.
5. Validate build variants used for release:
   - `RHBM_GEM_DEP_PROVIDER=FETCH`
   - `RHBM_GEM_DEP_PROVIDER=SYSTEM`
   - default `RHBM_GEM_ENABLE_UMAP=ON` with each dependency provider
   - explicit `RHBM_GEM_ENABLE_UMAP=OFF` with each dependency provider
   - `CMAKE_DISABLE_FIND_PACKAGE_ROOT=TRUE`
