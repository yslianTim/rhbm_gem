# Build and Configuration

This guide is for contributors and maintainers who need the full build surface: dependency selection, CMake parameters, feature-mode validation, and coverage workflows.

If you only want to install and run the project, use [`../user/getting-started.md`](../user/getting-started.md).

## Dependency Strategy

This project uses CMake + C++17. By default it prefers system-installed Eigen3, CLI11, pybind11, and SQLite3. If any are missing, CMake automatically falls back to the bundled copies in `third_party/`. Boost support is controlled independently with `RHBM_GEM_BOOST_MODE` (`AUTO` by default).

To force bundled dependencies where supported:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_LIBS=OFF
```

Boost has no bundled fallback. Keep `RHBM_GEM_BOOST_MODE=AUTO` or set `RHBM_GEM_BOOST_MODE=OFF` if Boost is unavailable.

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

## Feature Mode Checks (`AUTO` / `OFF` / `ON`)

1. Force OpenMP OFF:

```bash
cmake -S . -B build-openmp-off -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_OPENMP_MODE=OFF
cmake --build build-openmp-off -j
./build-openmp-off/RHBM-GEM --help
```

2. Force OpenMP ON on macOS with Homebrew `libomp`:

```bash
cmake -S . -B build-openmp-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_OPENMP_MODE=ON -DOpenMP_ROOT=/opt/homebrew/opt/libomp -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/libomp
cmake --build build-openmp-on -j
./build-openmp-on/RHBM-GEM --help
```

3. Force OpenMP ON on Linux:

```bash
cmake -S . -B build-openmp-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_OPENMP_MODE=ON
cmake --build build-openmp-on -j
./build-openmp-on/RHBM-GEM --help
```

4. Force ROOT OFF even if ROOT is installed:

```bash
cmake -S . -B build-root-off -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ROOT_MODE=OFF
cmake --build build-root-off -j
./build-root-off/RHBM-GEM --help
```

5. Force ROOT ON:

```bash
cmake -S . -B build-root-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ROOT_MODE=ON
```

6. Force Boost OFF:

```bash
cmake -S . -B build-boost-off -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=OFF
```

7. Force Boost ON:

```bash
cmake -S . -B build-boost-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=ON
```

8. Force Boost ON with required components:

```bash
cmake -S . -B build-boost-components -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=ON -DRHBM_GEM_BOOST_COMPONENTS="filesystem;system"
```

9. If OpenMP is missing on Linux:

```bash
sudo apt install -y libomp-dev
```

## CMake Parameters

Beginner / common:

| Parameter | Default | Description |
| --- | --- | --- |
| `CMAKE_BUILD_TYPE` | generator-dependent | Build type for single-config generators such as `Release` or `Debug`. |
| `CMAKE_INSTALL_PREFIX` | platform-dependent | Base install path for `cmake --install`. |
| `BUILD_TESTING` | `ON` | Build unit tests (`unit_tests` target). |
| `ENABLE_COVERAGE` | `OFF` | Enable `gcov` coverage instrumentation and the `coverage` target. |
| `COVERAGE_INCLUDE_TESTS` | `OFF` | Include `tests/` files in the coverage summary when coverage is enabled. |
| `USE_SYSTEM_LIBS` | `ON` | Prefer system packages for Eigen3, CLI11, pybind11, and SQLite3. |
| `BUILD_SHARED_LIBS` | `ON` | Build shared libraries instead of static libraries. |
| `BUILD_PYTHON_BINDINGS` | `ON` | Build the pybind11 module in `bindings/`. |
| `RHBM_GEM_OPENMP_MODE` | `AUTO` | OpenMP mode control: `AUTO`, `ON`, or `OFF`. |
| `RHBM_GEM_BOOST_MODE` | `AUTO` | Boost mode control: `AUTO`, `ON`, or `OFF`. |
| `RHBM_GEM_BOOST_COMPONENTS` | empty | Semicolon-separated Boost components required when Boost is enabled. |
| `RHBM_GEM_ROOT_MODE` | `AUTO` | ROOT mode control: `AUTO`, `ON`, or `OFF`. |
| `RHBM_GEM_PYTHON_INSTALL_LAYOUT` | `SITE_PREFIX` | Python module install layout: `SITE_PREFIX` or `LIBDIR`. |
| `RHBM_GEM_PYTHON_INSTALL_DIR` | empty | Explicit install directory for the Python extension module. |

Advanced / environment control:

| Parameter | Example | Description |
| --- | --- | --- |
| `CMAKE_PREFIX_PATH` | `/opt/homebrew;/opt/homebrew/opt/libomp` | Extra package search prefixes for `find_package(...)`. |
| `OpenMP_ROOT` | `/opt/homebrew/opt/libomp` | Help CMake find OpenMP on macOS/Homebrew. |
| `Boost_ROOT` | `/opt/homebrew/opt/boost` | Help CMake find Boost when automatic detection fails. |
| `Python_EXECUTABLE` | `/opt/homebrew/bin/python3` | Select the interpreter used to build bindings and derive the install layout. |

Notes:

1. For Eigen3, CLI11, pybind11, and SQLite3, `-DUSE_SYSTEM_LIBS=OFF` is usually simpler than setting multiple `CMAKE_DISABLE_FIND_PACKAGE_*` flags.
2. The project-specific mode flags (`RHBM_GEM_OPENMP_MODE`, `RHBM_GEM_BOOST_MODE`, `RHBM_GEM_ROOT_MODE`) are preferred over `CMAKE_DISABLE_FIND_PACKAGE_*`.

## Validation Examples

```bash
# Release build, no tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF

# Force bundled third-party dependencies
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_LIBS=OFF

# Pure C++ build (skip Python bindings)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=OFF

# Force ROOT/OpenMP/Boost requirements
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ROOT_MODE=ON -DRHBM_GEM_OPENMP_MODE=ON -DRHBM_GEM_BOOST_MODE=ON

# Force Boost with specific components
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=ON -DRHBM_GEM_BOOST_COMPONENTS="filesystem;system"

# Install Python module into <prefix>/lib/pythonX.Y/site-packages (default)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
cmake --install build --prefix "$HOME/.local"
PYVER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PYTHONPATH="$HOME/.local/lib/python${PYVER}/site-packages" python3 -c "import rhbm_gem_module"

# Keep the old libdir install style for the Python module
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_PYTHON_INSTALL_LAYOUT=LIBDIR
```

After installation, downstream CMake projects can consume this project with:

```cmake
find_package(RHBM_GEM CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE RHBM_GEM::core)
```

To inspect the full cache after configure:

```bash
cmake -S . -B build -LAH
```
