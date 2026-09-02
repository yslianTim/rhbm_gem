# Build and Configuration

This guide is for contributors and maintainers who need the full build surface: dependency selection, CMake parameters, feature-mode validation, and coverage workflows.

If you only want to install and run the project, use [`/docs/user/getting-started.md`](/docs/user/getting-started.md). For the end-user workflow, start with [`/docs/user/getting-started.md#environment-setup`](/docs/user/getting-started.md#environment-setup), then continue to [`/docs/user/getting-started.md#installation`](/docs/user/getting-started.md#installation), [`/docs/user/getting-started.md#python-bindings`](/docs/user/getting-started.md#python-bindings), and [`/docs/user/getting-started.md#python-examples`](/docs/user/getting-started.md#python-examples).

Top-level CMake logic is split into modular files under `/cmake/` (`RHBMGemOptions`, `RHBMGemDependencies`, `RHBMGemInstall`, `RHBMGemDevTools`) to keep maintenance localized while keeping test and lint entrypoints co-located.
All executable runtime targets from a build tree (CLI and C++ test executables) are emitted under `<build-dir>/bin/`.

## Dependency Strategy

This project requires CMake 3.24 or newer and uses C++20 with GNU extensions enabled by default. A single dependency-provider switch controls third-party resolution:

- `RHBM_GEM_DEP_PROVIDER=SYSTEM`: strictly require system packages for Eigen3 `>=5.0.0,<6.0.0`, CLI11, SQLite3, and Boost; also require `pybind11` plus Python development headers when bindings are enabled, and GTest when tests are enabled. UMAP is the only dependency group with the system-preferred fallback described below.
- `RHBM_GEM_DEP_PROVIDER=FETCH`: use pinned `FetchContent` sources for Eigen3 5.0.0, CLI11, SQLite3, and Boost; additionally fetch `pybind11` when bindings are enabled and GTest when tests are enabled.

The stable `umap_embedding` command is enabled by default. Set `RHBM_GEM_ENABLE_UMAP=OFF` to omit the command and avoid resolving its dependencies. When enabled, the build attaches the header-only `libscran::umappp` target to the internal build interface of `rhbm_gem`:

- With the `FETCH` provider, CMake downloads pinned releases of umappp 3.3.2, aarand 1.1.0, irlba 3.1.0, subpar 0.5.0, sanisizer 0.2.0, and knncolle 3.1.0. These build-only dependencies are excluded from the RHBM-GEM installation.
- With the `SYSTEM` provider, CMake first resolves and validates the real system Eigen3 package, then resolves the UMAP stack in dependency order while preferring each installed package. This ordering avoids unrecoverable nested `REQUIRED` lookups in an incomplete upstream package stack. If umappp or one of its package dependencies (`ltla_aarand`, `ltla_irlba`, `ltla_subpar`, `ltla_sanisizer`, and `knncolle_knncolle`) is unavailable, only the missing UMAP components are fetched at the pinned versions above; a complete installation causes no UMAP download. A build-tree-only package redirect lets fetched targets satisfy later package lookups. The fallback never downloads Eigen, CLI11, SQLite3, Boost, pybind11, or GTest; a qualifying system Eigen3 remains mandatory and is reused by umappp and irlba. Network access is required when a missing UMAP component is not already cached. No patched umappp or irlba package is required. Eigen3 6.x is outside the supported range. Upstream umappp 3.3.2 does not install its generated `ConfigVersion.cmake`, so its system lookup is unversioned and the documented system version remains v3.3.2.
- UMAP remains an internal dependency: installed consumers of `RHBM_GEM::rhbm_gem` do not need umappp. UMAP-enabled installations export the `RHBM_GEM_ENABLE_UMAP` compile definition so their public `UmapEmbeddingRequest` type matches the linked library.

To force fetched dependency sources:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_DEP_PROVIDER=FETCH
```

To force system packages:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_DEP_PROVIDER=SYSTEM
```

## Coverage (`gcov`)

This project supports text-based coverage reports using compiler-provided `gcov` only.

1. Configure a coverage build:

```bash
cmake -S . -B build-cov -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON
```

2. Build tests:

```bash
cmake --build build-cov -j
```

3. Run the coverage target:

```bash
cmake --build build-cov --target coverage
```

4. Read reports:

```bash
cat build-cov/coverage/coverage_summary.txt
cat build-cov/coverage/coverage_detail.txt
```

Notes:

1. Default summary includes only core `/src/` files and excludes `/src/python/`.
2. To include `/tests/`, configure with `-DCOVERAGE_INCLUDE_TESTS=ON`.
3. Coverage artifacts are generated under the build directory.
4. `ENABLE_COVERAGE=ON` requires `BUILD_TESTING=ON`.

## Test Execution and Labels

Build all C++ test targets:

```bash
cmake --build build --target tests_all -j
```

The C++ unit tests are linked into a single executable at `build/bin/RHBM-GEM-TEST`. The `ctest`
entries remain split by domain/intent grouping and invoke filtered subsets from that one binary.

Run all tests:

```bash
ctest --test-dir build --output-on-failure
```

Run tests by domain label:

```bash
ctest --test-dir build -L domain:data --output-on-failure
```

Run tests by intent label:

```bash
ctest --test-dir build -L intent:migration --output-on-failure
```

Supported labels:

- domain: `core`, `data`, `utils`, `integration`
- intent: `contract`, `command`, `validation`, `io`, `schema`, `migration`, `algorithm`, `bindings`

## Repository Lint Checks

Repository guard checks (style/structure/hygiene/fixture tracking/absolute-path/install consumer smoke) are executed through `lint_repo`:

```bash
cmake --build build --target lint_repo
```

`lint_repo` also exercises the command catalog indirectly through compile-time typed visitors and
contract tests.

## Static Quality Checks (Targeted)

Clang-tidy check for painter/parser directories:

```bash
bash resources/tools/developer/run_clang_tidy_check.sh build
```

Clang-tidy baseline guard (no-regression policy):

```bash
bash resources/tools/developer/run_clang_tidy_check.sh build --baseline resources/tools/developer/clang_tidy_baseline.json
```

CTest with failure classification output:

```bash
bash resources/tools/developer/run_ctest_with_classification.sh build -j8
```

This repository does not use GitHub-hosted CI. Run quality checks locally, or wire these commands into an external CI system if you need automation.

Suggested local validation bundles:

1. Fast validation: `lint_repo` + focused `ctest` labels for the area you changed, plus targeted format/tidy checks when touching covered painter/parser paths
2. Full validation: full `ctest` run + feature-mode configuration checks (`OPENMP` / `ROOT` / `LEGACY`) + failure classification output when investigating broader regressions

## Feature Mode Checks (`AUTO` / `OFF` / `ON`)

1. Force OpenMP OFF:

```bash
cmake -S . -B build-openmp-off -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_OPENMP_MODE=OFF
cmake --build build-openmp-off -j
./build-openmp-off/bin/RHBM-GEM --help
```

2. Force OpenMP ON on macOS with Homebrew `libomp`:

```bash
cmake -S . -B build-openmp-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_OPENMP_MODE=ON -DOpenMP_ROOT=/opt/homebrew/opt/libomp -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/libomp
cmake --build build-openmp-on -j
./build-openmp-on/bin/RHBM-GEM --help
```

Notes:

- On macOS with AppleClang, the project auto-probes Homebrew `libomp` in `/opt/homebrew/opt/libomp` and `/usr/local/opt/libomp` when `OpenMP_ROOT` is not set.
- If VSCode CMake Tools still reports missing OpenMP after installing `libomp`, clear the old configure cache and configure again (for example: `CMake: Delete Cache and Reconfigure`).

3. Force OpenMP ON on Linux:

```bash
cmake -S . -B build-openmp-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_OPENMP_MODE=ON
cmake --build build-openmp-on -j
./build-openmp-on/bin/RHBM-GEM --help
```

4. Force ROOT OFF even if ROOT is installed:

```bash
cmake -S . -B build-root-off -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ROOT_MODE=OFF
cmake --build build-root-off -j
./build-root-off/bin/RHBM-GEM --help
```

5. Force ROOT ON:

```bash
cmake -S . -B build-root-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ROOT_MODE=ON
```

6. If OpenMP is missing on Linux, install the package set that matches your compiler toolchain (Ubuntu/Debian + Clang example):

```bash
sudo apt install -y libomp-dev
```

For GCC-based builds, ensure `g++`/`build-essential` is installed so `libgomp` is available.

## CMake Parameters

Beginner / common:

| Parameter | Default | Description |
| --- | --- | --- |
| `CMAKE_BUILD_TYPE` | generator-dependent | Build type for single-config generators such as `Release` or `Debug`. |
| `CMAKE_INSTALL_PREFIX` | platform-dependent | Base install path for `cmake --install`. |
| `BUILD_TESTING` | `ON` | Build test targets (aggregate target: `tests_all`). |
| `ENABLE_COVERAGE` | `OFF` | Enable `gcov` coverage instrumentation and the `coverage` target. |
| `COVERAGE_INCLUDE_TESTS` | `OFF` | Include `/tests/` files in the coverage summary when coverage is enabled. |
| `RHBM_GEM_DEP_PROVIDER` | `SYSTEM` | Dependency provider mode: `SYSTEM` or `FETCH`. |
| `BUILD_SHARED_LIBS` | `ON` | Build shared libraries instead of static libraries. |
| `BUILD_PYTHON_BINDINGS` | `ON` | Build the pybind11 module in `/src/python/`. |
| `RHBM_GEM_ENABLE_UMAP` | `ON` | Resolve umappp and compile the stable `umap_embedding` command. |
| `RHBM_GEM_OPENMP_MODE` | `AUTO` | OpenMP mode control: `AUTO`, `ON`, or `OFF`. |
| `RHBM_GEM_ROOT_MODE` | `AUTO` | ROOT mode control: `AUTO`, `ON`, or `OFF`. |
| `RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` | `OFF` | Enable experimental features across the project. |
| `RHBM_GEM_ENABLE_FOLD_168_REGRESSION` | `OFF` | Enable the opt-in external-data 168-atom simulation regression. |
| `RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT` | `OFF` | Enable developer-only frozen-IRLS trust-model instrumentation; requires tests. |
| `RHBM_GEM_FOLD_168_MODEL` | empty | Path to the hash-verified fold-168 CIF input. |
| `RHBM_GEM_FOLD_168_MAP` | empty | Path to the hash-verified fold-168 map input. |
| `RHBM_GEM_PYTHON_INSTALL_LAYOUT` | `SITE_PREFIX` | Python module install layout: `SITE_PREFIX` or `LIBDIR`. |
| `RHBM_GEM_PYTHON_INSTALL_DIR` | empty | Explicit install directory for the Python extension module. |

Advanced / environment control:

| Parameter | Example | Description |
| --- | --- | --- |
| `CMAKE_PREFIX_PATH` | `/opt/homebrew;/opt/homebrew/opt/libomp` | Extra package search prefixes for `find_package(...)`. |
| `OpenMP_ROOT` | `/opt/homebrew/opt/libomp` | Help CMake find OpenMP on macOS/Homebrew. |
| `Boost_ROOT` | `/opt/homebrew/opt/boost` | Help CMake find system Boost when automatic detection fails. |
| `Python_EXECUTABLE` | `/opt/homebrew/bin/python3` | Select the interpreter used to build bindings and derive the install layout. |

Notes:

1. `RHBM_GEM_DEP_PROVIDER=FETCH` is recommended when system dependencies are unavailable.
2. The project-specific mode flags (`RHBM_GEM_OPENMP_MODE`, `RHBM_GEM_ROOT_MODE`) are preferred over `CMAKE_DISABLE_FIND_PACKAGE_*`.
3. `FETCH` mode requires network access unless archives are already cached.
4. A fresh non-preset configure defaults to `SYSTEM + UMAP ON`. Core dependencies, including Eigen3, must be installed; the build prefers the installed UMAP stack and fetches only missing UMAP components at fixed versions. The Debug and Release presets explicitly use `FETCH + UMAP ON` instead.
5. Changing the source default does not overwrite an existing CMake cache. A manually configured build cached with UMAP `OFF` must be reconfigured with `-DRHBM_GEM_ENABLE_UMAP=ON` or recreated; the standard presets explicitly apply `ON`.
6. Set `RHBM_GEM_ENABLE_UMAP=OFF` to skip both the UMAP system probe and its network fallback.

## CMake Preset Workflow

The project provides Debug and Release configure/build presets.
The presets require CMake 3.24 or newer and the Ninja generator. They use
separate build directories, explicitly select `FETCH + UMAP ON`, disable tests
and Python bindings, and build the CLI target. The first configure downloads
the pinned dependencies unless they are already cached:

| Preset | Build type | Executable |
| --- | --- | --- |
| `debug` | `Debug` | `build/debug/bin/RHBM-GEM` |
| `release` | `Release` | `build/release/bin/RHBM-GEM` |

To configure and incrementally build both variants in order:

```bash
cmake -P cmake/BuildAllPresets.cmake
```

The script stops at the first configure or build failure and identifies the
failed preset. Running it again reuses the two existing build directories.

To configure or build one variant independently:

```bash
cmake --preset debug
cmake --build --preset debug

cmake --preset release
cmake --build --preset release
```

In VS Code, `CMake: Build Debug + Release` is the default build
task, so `Cmd/Ctrl + Shift + B` runs the same all-variants build script.

`CMakePresets.json` is the single source of truth for the standard configure and
build settings. The VS Code default build task is the recommended project
workflow; workspace-level CMake Tools settings should not override the presets.
The shared Ninja progress format is defined by the build presets so the same
output is used from VS Code and the command line.

VS Code IntelliSense is provided by the CMake Tools configuration provider in
`.vscode/c_cpp_properties.json`. Select a Debug or Release
configure preset before editing. Do not add compiler paths or third-party
include paths there; CMake Tools supplies them from the active CMake target.
The common configure preset also exports `compile_commands.json` for each
standard build directory, which keeps auxiliary tooling aligned with CMake.

When the LLVM `clangd` extension is enabled, `.clangd` points it to the Debug
core build database. The standard presets share the same compiler, C++20 mode,
defines, and dependency paths; only optimization/debug flags differ. Keep the
Debug preset configured when using clangd for the default core source tree.

For coverage, external regression data, or other option-specific validation,
continue to use explicit `cmake -S/-B` build directories as shown in the
sections above. This keeps those specialized configurations independent from
the two standard preset directories.

## Validation Examples

These examples are for validating configuration behavior and CMake options. For platform-specific end-user setup and install flows, use [`/docs/user/getting-started.md#environment-setup`](/docs/user/getting-started.md#environment-setup), [`/docs/user/getting-started.md#installation`](/docs/user/getting-started.md#installation), [`/docs/user/getting-started.md#python-bindings`](/docs/user/getting-started.md#python-bindings), and [`/docs/user/getting-started.md#python-examples`](/docs/user/getting-started.md#python-examples).

```bash
# Release build, no tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF

# Force fetched dependencies
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_DEP_PROVIDER=FETCH

# Force system dependencies
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_DEP_PROVIDER=SYSTEM

# Disable the default UMAP embedding command and its dependency stack
cmake -S . -B build-no-umap \
  -DCMAKE_BUILD_TYPE=Release \
  -DRHBM_GEM_ENABLE_UMAP=OFF

# Pure C++ build (skip Python bindings)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=OFF

# Force ROOT/OpenMP requirements
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ROOT_MODE=ON -DRHBM_GEM_OPENMP_MODE=ON

# Enable project-wide experimental features
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE=ON

# Enable the external 168-atom regression benchmark
cmake -S . -B build-fold-168 \
  -DRHBM_GEM_ENABLE_FOLD_168_REGRESSION=ON \
  -DRHBM_GEM_FOLD_168_MODEL=/path/to/fold_test_model_0.cif \
  -DRHBM_GEM_FOLD_168_MAP=/path/to/sim_map_gaus_grid0.10_charge1_bw0.50.map
cmake --build build-fold-168 --target tests_all -j
ctest --test-dir build-fold-168 -R fold_168_simulation_regression --output-on-failure

# Enable the developer-only frozen-IRLS trust-model experiment
cmake -S . -B build-trust-model \
  -DBUILD_TESTING=ON \
  -DBUILD_PYTHON_BINDINGS=OFF \
  -DRHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT=ON
cmake --build build-trust-model --target tests_all -j

# Run the production-only convergence corpus (no experiment flag required)
cmake --build build --target convergence_exposure_corpus

# Install Python module into <prefix>/<CMAKE_INSTALL_LIBDIR>/pythonX.Y/site-packages (default layout)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
cmake --install build --prefix "$HOME/.local"
PYVER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PYTHONPATH="$HOME/.local/lib/python${PYVER}/site-packages:$HOME/.local/lib64/python${PYVER}/site-packages" python3 -c "import rhbm_gem_module"

# Keep the old libdir install style for the Python module
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_PYTHON_INSTALL_LAYOUT=LIBDIR
```

Note: The Python examples here demonstrate layout validation only. For the user-facing install and example flow, follow [`/docs/user/getting-started.md#python-bindings`](/docs/user/getting-started.md#python-bindings) and [`/docs/user/getting-started.md#python-examples`](/docs/user/getting-started.md#python-examples).

The fold-168 runner verifies the two external SHA-256 identities before it
starts the command. It passes a nonexistent database path in a temporary
directory to `potential_analysis`; the normal persistence layer creates the
current schema and stores the benchmark output. The database is therefore an
ephemeral output, not a fitting input. The blocking gate requires the complete
168 serial-ID set, finite valid atom parameters, and truth-based
amplitude/width/offset RMSE plus maximum absolute offset no worse than 105% of
the schema-4 reference metrics. The parsed second-stage final summary must
report no more than ten accepted iterations. The single wall-time measurement
remains diagnostic/report-only.
The required fixture hashes are:

- CIF: `156d35aa326f0d4408d726a999329d2ffede775489aeaa5d99a2cc9b9f663cab`
- map: `5e0dbb13fc3a76f8a944e6e2b18393d1896fafc2ec9020457cca8e8a421f120e`

Run artifacts are written under
`<build>/benchmark-results/fold_168/{run.log,actual.json,report.json}`. A valid
algorithm change must update the checked-in reference metrics only after
manually reviewing `actual.json`; the runner never overwrites the baseline.
`run.log` is retained for diagnosis, and only the stable second-stage final
summary is parsed. The temporary
database is deleted after its 168 atom results have been read and is not
retained as an artifact.

After installation, downstream CMake projects can consume this project with:

```cmake
find_package(RHBM_GEM CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE RHBM_GEM::rhbm_gem)
```

## Migration Notes

If you are upgrading from the previous CMake interface:

- `USE_SYSTEM_LIBS=ON` -> `RHBM_GEM_DEP_PROVIDER=SYSTEM`
- `USE_SYSTEM_LIBS=OFF` -> `RHBM_GEM_DEP_PROVIDER=FETCH`
- `RHBM_GEM::core` (and split library targets) -> `RHBM_GEM::rhbm_gem`
- repository/style guard tests in `ctest` -> `cmake --build <build> --target lint_repo`

If you consume the C++ algorithm headers directly:

- `NormalizedChange.hpp`, `ParameterChange.hpp`, and `ParameterChangeStats.hpp`
  -> `Convergence.hpp`
- `WeightedRidgeSystem.hpp` -> `WeightedRidgeSolver.hpp`; `WeightedRidgeSystem`
  keeps the same name and fields
- `WeightedRidgeNormalEquation.hpp`, `WeightedRidgeNormalEquation`, and
  `BuildWeightedRidgeNormalEquation(...)` were internalized into
  `WeightedRidgeSolver`
- `CalculateRobustLoss(RobustLossKind::Cauchy, residual, cutoff)`
  -> `CalculateCauchyLoss(residual, cutoff)`
- `CalculateRobustWeight(RobustLossKind::Cauchy, residual, scale, multiplier)`
  -> `CalculateCauchyWeight(residual, scale, multiplier)`; the unused Huber
  variant was removed

To inspect the full cache after configure:

```bash
cmake -S . -B build -LAH
```
