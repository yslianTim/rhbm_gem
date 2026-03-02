# Robust Hierarchical Bayesian Model for Gaussian Estimation of cryo-EM Maps (RHBM-GEM)
RHBM-GEM is a robust hierarchical Bayesian framework for Gaussian modeling and parameter estimation of atom-associated cryo-EM map distributions.



## Getting Started (macOS / Linux / Windows)

This project uses CMake + C++17. By default it prefers system-installed Eigen3 / CLI11 / pybind11 / SQLite3; if any are missing, CMake automatically falls back to the bundled copies in `third_party/`. Boost support is controlled independently with `RHBM_GEM_BOOST_MODE` (`AUTO` by default). To force bundled deps for the bundled set, pass `-DUSE_SYSTEM_LIBS=OFF` during configure.

**macOS (Homebrew)**
1. Install Xcode Command Line Tools:
```bash
xcode-select --install
```
2. Install Homebrew if `brew` is not available:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
3. Install dependencies (adjust as needed):
```bash
brew install cmake eigen sqlite3 python pybind11 cli11 libomp root boost
```

**Linux (Ubuntu/Debian example)**
1. Install dependencies (adjust as needed):
```bash
sudo apt update
sudo apt install -y build-essential cmake pkg-config python3 python3-dev libsqlite3-dev libeigen3-dev libboost-dev
```
2. If you want to build Python bindings, also install:
```bash
sudo apt install -y pybind11-dev
```

**Windows (PowerShell + Visual Studio 2022)**
1. Install prerequisites:
   - Visual Studio 2022 (or Build Tools) with the `Desktop development with C++` workload
   - CMake (>= 3.18), Python 3, and Git
2. Quick start using bundled dependencies (recommended on Windows):
```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DUSE_SYSTEM_LIBS=OFF -DRHBM_GEM_ROOT_MODE=OFF
cmake --build build --config Release
```
3. Optional: use vcpkg system dependencies instead of bundled copies:
```powershell
git clone https://github.com/microsoft/vcpkg $env:USERPROFILE\vcpkg
& $env:USERPROFILE\vcpkg\bootstrap-vcpkg.bat
& $env:USERPROFILE\vcpkg\vcpkg.exe install eigen3 sqlite3 pybind11 cli11 boost:x64-windows
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="$env:USERPROFILE\vcpkg\scripts\buildsystems\vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows
cmake --build build --config Release
```

**ROOT (recommended for plotting features)**
This project can build without ROOT, but ROOT is strongly recommended if you need figure/plot output.

1. Install ROOT on macOS (Homebrew):
```bash
brew install root
```
2. Install ROOT on Linux (choose one):
```bash
# Option A: package manager (if available on your distro)
sudo apt update
sudo apt install -y root-system

# Option B: conda-forge (works on most Linux setups)
conda create -n emmap-root -c conda-forge root
conda activate emmap-root
```
3. Install ROOT on Windows (conda-forge):
```powershell
conda create -n emmap-root -c conda-forge root
conda activate emmap-root
```
4. Verify ROOT installation:
```bash
root-config --version
root-config --prefix
```
On Windows + conda, verify with PowerShell:
```powershell
where.exe root.exe
root.exe -b -q
```
5. If CMake cannot find ROOT automatically, pass ROOT prefix explicitly:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$(root-config --prefix)"
```
On Windows + conda, use PowerShell:
```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_PREFIX_PATH="$env:CONDA_PREFIX;$env:CONDA_PREFIX\Library"
```

Without ROOT, build still succeeds, but plotting code paths are compiled out (`HAVE_ROOT` disabled). In practice:
1. `potential_display` (all painter modes: `gaus`, `model`, `comparison`, `demo`, `atom`) will not generate ROOT-based figures.
2. `map_visualization` will not generate the 2D histogram plot output from `LocalPainter::PaintHistogram2D`.
3. `potential_analysis` study plots from `LocalPainter::PaintTemplate1` (for example alpha-scan figure outputs) will not be generated.
4. Non-plotting workflows (data processing / analysis / simulation / dumping) remain available.

**Boost (optional)**
Boost is optional. In `AUTO` mode, CMake enables Boost features only when Boost is found.

1. Build with default Boost behavior (`AUTO`):
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=AUTO
```
2. Force Boost ON (configure fails if missing):
```bash
cmake -S . -B build-boost-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=ON
```
3. Force Boost OFF:
```bash
cmake -S . -B build-boost-off -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=OFF
```
4. Require specific Boost components:
```bash
cmake -S . -B build-boost-components -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=ON -DRHBM_GEM_BOOST_COMPONENTS="filesystem;system"
```

**Installation**
1. From the project root:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```
2. Run the main executable:
```bash
./build/RHBM-GEM --help
```
3. Install to the default prefix (`CMAKE_INSTALL_PREFIX`, usually `/usr/local`):
```bash
cmake --install build
```
4. Install to a user-specified path (choose one):
```bash
# Option A: set install prefix at configure time
cmake -S . -B build-local -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build-local -j
cmake --install build-local

# Option B: override prefix at install time
cmake --install build --prefix "$HOME/.local"
```
5. Run the installed executable:
```bash
"$HOME/.local/bin/RHBM-GEM" --help
```
Windows PowerShell equivalent:
```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DUSE_SYSTEM_LIBS=OFF -DRHBM_GEM_ROOT_MODE=OFF
cmake --build build --config Release
.\build\Release\RHBM-GEM.exe --help
cmake --install build --config Release --prefix "$env:USERPROFILE\AppData\Local\RHBM-GEM"
& "$env:USERPROFILE\AppData\Local\RHBM-GEM\bin\RHBM-GEM.exe" --help
```

**Troubleshooting**
1. If you do not have Eigen / SQLite3 / pybind11 / CLI11 installed, you can skip installing them; CMake will use the bundled versions in `third_party/`.
2. Boost has no bundled fallback. Keep `RHBM_GEM_BOOST_MODE=AUTO` (default) or set `RHBM_GEM_BOOST_MODE=OFF` if Boost is not available.
3. To force bundled deps (where available), configure with:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_LIBS=OFF
```

**Feature mode checks (`AUTO` / `OFF` / `ON`)**
1. Force OpenMP OFF (macOS / Linux):
```bash
cmake -S . -B build-openmp-off -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_OPENMP_MODE=OFF
cmake --build build-openmp-off -j
./build-openmp-off/RHBM-GEM --help
```
2. Force OpenMP ON on macOS (Homebrew `libomp`):
```bash
cmake -S . -B build-openmp-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_OPENMP_MODE=ON -DOpenMP_ROOT=/opt/homebrew/opt/libomp -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/libomp
cmake --build build-openmp-on -j
./build-openmp-on/RHBM-GEM --help
```
3. Force OpenMP ON on Linux (no `OpenMP_ROOT` needed):
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
5. Force ROOT ON (configure fails if ROOT is unavailable):
```bash
cmake -S . -B build-root-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_ROOT_MODE=ON
```
6. Force Boost OFF:
```bash
cmake -S . -B build-boost-off -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=OFF
```
7. Force Boost ON (configure fails if Boost is unavailable):
```bash
cmake -S . -B build-boost-on -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=ON
```
8. Force Boost ON with required components:
```bash
cmake -S . -B build-boost-components -DCMAKE_BUILD_TYPE=Release -DRHBM_GEM_BOOST_MODE=ON -DRHBM_GEM_BOOST_COMPONENTS="filesystem;system"
```
9. If OpenMP is missing on Linux, install:
```bash
sudo apt install -y libomp-dev
```

**Python bindings module name**
1. The canonical Python extension module name is `rhbm_gem_module`.

**Python examples**
Recommended single-prefix quickstart (build + install + verify + run):
```bash
INSTALL_PREFIX="$HOME/.local"

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_PYTHON_BINDINGS=ON \
  -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}"
cmake --build build -j
cmake --install build

PYVER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PY_SITE="${INSTALL_PREFIX}/lib/python${PYVER}/site-packages"

# Verify module import before running examples.
PYTHONPATH="${PY_SITE}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 -c "import rhbm_gem_module as rgm; print(rgm.__file__)"

# Run examples from source tree.
PYTHONPATH="${PY_SITE}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 examples/python/00_quickstart.py
PYTHONPATH="${PY_SITE}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 examples/python/01_end_to_end_from_test_data.py --workdir /tmp/rhbm_py_demo

# Run installed example.
PYTHONPATH="${PY_SITE}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 "${INSTALL_PREFIX}/share/RHBM_GEM/examples/python/01_end_to_end_from_test_data.py" \
  --workdir /tmp/rhbm_py_demo_installed
```

Expected outputs for the pipeline example:
1. SQLite database: `<workdir>/demo.sqlite`
2. Simulated map files: `<workdir>/maps/sim_map_*.map`
3. Dump files: `<workdir>/dump/*`

Troubleshooting:
1. `ModuleNotFoundError: No module named 'rhbm_gem_module'`  
   Verify `BUILD_PYTHON_BINDINGS=ON`, run `cmake --install`, and ensure `PYTHONPATH` points to your install site-packages.
2. Install prefix mismatch (`/usr/local` vs `$HOME/.local`)  
   If you install under `/usr/local`, do not use `$HOME/.local/.../site-packages` in `PYTHONPATH`.
   Reinstall with a single `INSTALL_PREFIX` (recommended), or update `PYTHONPATH` to match your actual install prefix.
3. `Could not find sample model ...`  
   Pass `--model /path/to/your_model.cif` explicitly.

**Coverage test without additional third-party libraries (`gcov`)**
This project supports text-based coverage reports using compiler-provided `gcov` only.

1. Configure a coverage build:
```bash
cmake -S . -B build-cov -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON
```
2. Build tests:
```bash
cmake --build build-cov -j
```
3. Run coverage target:
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

**CMake parameters (quick reference)**
Beginner (most common):

| Parameter | Default | Description |
|---|---|---|
| `CMAKE_BUILD_TYPE` | generator-dependent | Build type for single-config generators (for example `Release` / `Debug`). |
| `CMAKE_INSTALL_PREFIX` | platform-dependent (usually `/usr/local`) | Base install path for `cmake --install`; change this to install into a custom location. |
| `BUILD_TESTING` | `ON` | Build unit tests (`unit_tests` target). |
| `ENABLE_COVERAGE` | `OFF` | Enable `gcov` coverage instrumentation and the `coverage` target. Requires GNU/Clang/AppleClang with `BUILD_TESTING=ON`. |
| `COVERAGE_INCLUDE_TESTS` | `OFF` | Include `tests/` files in coverage summary when `ENABLE_COVERAGE=ON`. |
| `USE_SYSTEM_LIBS` | `ON` | Prefer system packages for Eigen3 / CLI11 / pybind11 / SQLite3. If `OFF`, use bundled `third_party` copies. |
| `BUILD_SHARED_LIBS` | `ON` | Build shared libraries instead of static libraries. |
| `BUILD_PYTHON_BINDINGS` | `ON` | Build the pybind11 module in `bindings/`. Set `OFF` for pure C++ builds without Python dependency. |
| `RHBM_GEM_OPENMP_MODE` | `AUTO` | OpenMP mode control: `AUTO`, `ON`, or `OFF`. `ON` fails configure if unavailable. |
| `RHBM_GEM_BOOST_MODE` | `AUTO` | Boost mode control: `AUTO`, `ON`, or `OFF`. `ON` fails configure if unavailable. |
| `RHBM_GEM_BOOST_COMPONENTS` | empty | Semicolon-separated Boost components required when Boost is enabled (for example `filesystem;system`). |
| `RHBM_GEM_ROOT_MODE` | `AUTO` | ROOT mode control: `AUTO`, `ON`, or `OFF`. `ON` fails configure if unavailable. |
| `RHBM_GEM_PYTHON_INSTALL_LAYOUT` | `SITE_PREFIX` | Python module install layout: `SITE_PREFIX` or `LIBDIR`. |
| `RHBM_GEM_PYTHON_INSTALL_DIR` | empty | Explicit install directory for Python extension module (overrides layout). |

Advanced (debugging / environment control):

| Parameter | Example | Description |
|---|---|---|
| `CMAKE_PREFIX_PATH` | `/opt/homebrew;/opt/homebrew/opt/libomp` | Extra package search prefixes for `find_package(...)`. |
| `OpenMP_ROOT` | `/opt/homebrew/opt/libomp` | Help CMake find OpenMP on macOS/Homebrew. |
| `Boost_ROOT` | `/opt/homebrew/opt/boost` | Help CMake find Boost when automatic detection fails. |
| `Python_EXECUTABLE` | `/opt/homebrew/bin/python3` | Select Python interpreter used to build bindings and derive major/minor version for install layout. |

Notes:
1. For Eigen3 / CLI11 / pybind11 / SQLite3, `-DUSE_SYSTEM_LIBS=OFF` is usually simpler than setting multiple `CMAKE_DISABLE_FIND_PACKAGE_*` flags.
2. The project-specific mode flags (`RHBM_GEM_OPENMP_MODE`, `RHBM_GEM_BOOST_MODE`, `RHBM_GEM_ROOT_MODE`) are preferred over `CMAKE_DISABLE_FIND_PACKAGE_*`.

Examples:

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

# Keep old libdir install style for Python module
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

## Verbosity Levels

The command line option `-v, --verbose` controls how much information is printed.
The logger supports five levels:

| Level |  Name   | Description                       |
|-------|---------|-----------------------------------|
| 0     | Error   | Only error messages are shown.    |
| 1     | Warning | Warnings and errors are shown.    |
| 2     | Notice  | Notice messages are shown.        |
| 3     | Info    | Informational messages are shown. |
| 4     | Debug   | Detailed debug output.            |

Passing a larger value to `--verbose` enables all lower levels as well.



## License and Third-Party Notices

- Project license: Apache License 2.0 (see `LICENSE`).
- Third-party licenses and attribution summary: see `THIRD_PARTY_NOTICES.md`.
- This repository can use either system packages or bundled third-party dependencies under `third_party/` depending on CMake configuration.

**Distribution policy (default release profile)**
1. Release builds define `EIGEN_MPL2_ONLY` to constrain Eigen header usage to MPL-2.0-only subsets.
2. ROOT and Boost are optional and not required for non-plot workflows; release artifacts may be shipped without these linked binaries.
3. If binaries are distributed, include `LICENSE` and `THIRD_PARTY_NOTICES.md` in the package.
4. If large `.sqlite` datasets are distributed, verify their source-data license terms separately.

**Compliance checklist (source and binary release)**
1. Include `LICENSE`.
2. Include `THIRD_PARTY_NOTICES.md`.
3. Ensure bundled third-party license file paths referenced in `THIRD_PARTY_NOTICES.md` remain valid.
4. Validate build variants used for release:
   - `USE_SYSTEM_LIBS=OFF`
   - `USE_SYSTEM_LIBS=ON`
   - `CMAKE_DISABLE_FIND_PACKAGE_ROOT=TRUE`
   - `CMAKE_DISABLE_FIND_PACKAGE_Boost=TRUE`


# Overall Information for Developer

For the repository-specific development conventions and code habits that are already in use, see [docs/development-guidelines.md](docs/development-guidelines.md).
It summarizes project structure, naming/style patterns, command/data-layer design, logging, testing, and documentation update expectations.


## Data Object Manager

Commands rely on `DataObjectManager` to load and store data objects. The manager instance is owned by `CommandBase` and shared with derived commands. A command lazily creates a `DataObjectManager` and initializes it with a database path via `SetDatabaseManager`.

`DataObjectManager` is not accessed concurrently. The database pointer and the internal data map are guarded by simple `std::mutex` members to keep the implementation straightforward.

`DatabaseManager` is safe for concurrent use. Operations like `SaveDataObject` and `LoadDataObject` lock an internal mutex so multiple threads can access the database simultaneously without interference.

### File IO and persistence work through cooperating managers:

- **DataObjectManager** selects a reader or writer based on file extension. It uses `FileProcessFactoryRegistry` to instantiate the appropriate factory. Commands load and save objects through it, delegating file operations and database access to `DatabaseManager`.
- **DatabaseManager** wraps the SQLite database. DAO objects are created via `DataObjectDAOFactoryRegistry` and cached for reuse.

Together these managers keep data loading and persistence separate from command logic.



## Unit Tests (via GoogleTest v1.17.0)

The test suite automatically includes any immediate `.cpp` file in `tests/` directory.
