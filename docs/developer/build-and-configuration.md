# Build and Configuration

This guide is for contributors and maintainers who need the full build surface: dependency selection, CMake parameters, feature-mode validation, and coverage workflows.

If you only want to install and run the project, use [`../user/getting-started.md`](../user/getting-started.md). For the end-user workflow, start with [`../user/getting-started.md#environment-setup`](../user/getting-started.md#environment-setup), then continue to [`../user/getting-started.md#installation`](../user/getting-started.md#installation), [`../user/getting-started.md#python-bindings`](../user/getting-started.md#python-bindings), and [`../user/getting-started.md#python-examples`](../user/getting-started.md#python-examples).

Top-level CMake logic is split into modular files under `cmake/` (`RHBMGemOptions`, `RHBMGemDependencies`, `RHBMGemInstall`, `RHBMGemDevTools`) to keep maintenance localized while keeping test and lint entrypoints co-located.
All runtime targets from a build tree (CLI, GUI, and C++ test executables) are emitted under `<build-dir>/bin/`.

## Dependency Strategy

This project uses CMake + C++17 with a single dependency-provider switch:

- `RHBM_GEM_DEP_PROVIDER=SYSTEM`: require system packages for Eigen3, CLI11, SQLite3, pybind11 (if bindings enabled), and GTest (if tests enabled).
- `RHBM_GEM_DEP_PROVIDER=FETCH`: use pinned `FetchContent` sources for Eigen3, CLI11, SQLite3, Boost, pybind11, and GTest.

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

1. Default summary includes only `src/` files.
2. To include `tests/`, configure with `-DCOVERAGE_INCLUDE_TESTS=ON`.
3. Coverage artifacts are generated under the build directory.
4. `ENABLE_COVERAGE=ON` requires `BUILD_TESTING=ON`.

## Test Execution and Labels

Build all C++ test targets:

```bash
cmake --build build --target tests_all -j
```

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

`lint_repo` also exercises the command manifest indirectly through compile-time expansion and
contract tests.

## Static Quality Checks (Targeted)

Format check for painter/parser directories:

```bash
bash scripts/developer/run_clang_format_check.sh
```

Clang-tidy check for painter/parser directories:

```bash
bash scripts/developer/run_clang_tidy_check.sh build
```

Clang-tidy baseline guard (no-regression policy):

```bash
bash scripts/developer/run_clang_tidy_check.sh build --baseline scripts/developer/clang_tidy_baseline.json
```

CTest with failure classification output:

```bash
bash scripts/developer/run_ctest_with_classification.sh build -j8
```

CI layout:

1. Pull request / push: fast gate (`lint_repo` + `core/data/integration` tests + format/tidy checks)
2. Nightly schedule: full test run with feature-mode matrix (`OPENMP` / `ROOT` / `LEGACY`) + failure classification (`build` / `contract` / `runtime`)

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
| `COVERAGE_INCLUDE_TESTS` | `OFF` | Include `tests/` files in the coverage summary when coverage is enabled. |
| `RHBM_GEM_DEP_PROVIDER` | `SYSTEM` | Dependency provider mode: `SYSTEM` or `FETCH`. |
| `BUILD_SHARED_LIBS` | `ON` | Build shared libraries instead of static libraries. |
| `BUILD_PYTHON_BINDINGS` | `ON` | Build the pybind11 module in `bindings/`. |
| `RHBM_GEM_BUILD_GUI` | `ON` | Build the Qt6 GUI executable (`RHBM-GEM-GUI`) when Qt6 Core/Widgets are available. |
| `RHBM_GEM_OPENMP_MODE` | `AUTO` | OpenMP mode control: `AUTO`, `ON`, or `OFF`. |
| `RHBM_GEM_ROOT_MODE` | `AUTO` | ROOT mode control: `AUTO`, `ON`, or `OFF`. |
| `RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS` | `OFF` | Enable the experimental bond-analysis workflow hook inside `PotentialAnalysisCommand`. |
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

## Standard Build Directory Workflow

Use explicit source/build directory commands instead of CMake presets:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build
```

## Validation Examples

These examples are for validating configuration behavior and CMake options. For platform-specific end-user setup and install flows, use [`../user/getting-started.md#environment-setup`](../user/getting-started.md#environment-setup), [`../user/getting-started.md#installation`](../user/getting-started.md#installation), [`../user/getting-started.md#python-bindings`](../user/getting-started.md#python-bindings), and [`../user/getting-started.md#python-examples`](../user/getting-started.md#python-examples).

```bash
# Release build, no tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF

# Force fetched dependencies
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_DEP_PROVIDER=FETCH

# Force system dependencies
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_DEP_PROVIDER=SYSTEM

# Pure C++ build (skip Python bindings)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=OFF

# Build the Qt6 GUI executable (optional; target exists only when Qt6 Core/Widgets is found)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BUILD_GUI=ON
cmake --build build --target RHBM_GEM_gui -j

# Force ROOT/OpenMP requirements
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ROOT_MODE=ON -DRHBM_GEM_OPENMP_MODE=ON

# Enable the experimental PotentialAnalysis bond workflow
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS=ON

# Install Python module into <prefix>/<CMAKE_INSTALL_LIBDIR>/pythonX.Y/site-packages (default layout)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
cmake --install build --prefix "$HOME/.local"
PYVER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PYTHONPATH="$HOME/.local/lib/python${PYVER}/site-packages:$HOME/.local/lib64/python${PYVER}/site-packages" python3 -c "import rhbm_gem_module"

# Keep the old libdir install style for the Python module
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_PYTHON_INSTALL_LAYOUT=LIBDIR
```

Note: The Python examples here demonstrate layout validation only. For the user-facing install and example flow, follow [`../user/getting-started.md#python-bindings`](../user/getting-started.md#python-bindings) and [`../user/getting-started.md#python-examples`](../user/getting-started.md#python-examples).

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

To inspect the full cache after configure:

```bash
cmake -S . -B build -LAH
```
