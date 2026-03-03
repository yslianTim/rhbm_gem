# Getting Started

This guide is for end users who want to install RHBM-GEM, build it on their platform, verify the Python bindings, and run the supported examples.

If you are validating build matrices, collecting coverage, or tuning advanced CMake settings, use [`../developer/build-and-configuration.md`](../developer/build-and-configuration.md) instead.

## Platforms

This project uses CMake + C++17. By default it prefers system-installed Eigen3, CLI11, pybind11, and SQLite3. If any of those are missing, CMake automatically falls back to the bundled copies in `third_party/`. Boost support is controlled independently with `RHBM_GEM_BOOST_MODE` (`AUTO` by default).

## macOS (Homebrew)

1. Install Xcode Command Line Tools:

```bash
xcode-select --install
```

2. Install Homebrew if `brew` is not available:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

3. Install dependencies as needed:

```bash
brew install cmake eigen sqlite3 python pybind11 cli11 libomp root boost
```

## Linux (Ubuntu/Debian example)

1. Install dependencies:

```bash
sudo apt update
sudo apt install -y build-essential cmake pkg-config python3 python3-dev libsqlite3-dev libeigen3-dev libboost-dev
```

2. If you want Python bindings, also install:

```bash
sudo apt install -y pybind11-dev
```

## Windows (PowerShell + Visual Studio 2022)

1. Install prerequisites:
   - Visual Studio 2022 (or Build Tools) with the `Desktop development with C++` workload
   - CMake (>= 3.18), Python 3, and Git
2. Quick start using bundled dependencies:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DUSE_SYSTEM_LIBS=OFF -DRHBM_GEM_ROOT_MODE=OFF
cmake --build build --config Release
```

3. Optional: use vcpkg dependencies instead of bundled copies:

```powershell
git clone https://github.com/microsoft/vcpkg $env:USERPROFILE\vcpkg
& $env:USERPROFILE\vcpkg\bootstrap-vcpkg.bat
& $env:USERPROFILE\vcpkg\vcpkg.exe install eigen3 sqlite3 pybind11 cli11 boost:x64-windows
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="$env:USERPROFILE\vcpkg\scripts\buildsystems\vcpkg.cmake" -DVCPKG_TARGET_TRIPLET=x64-windows
cmake --build build --config Release
```

## ROOT

RHBM-GEM can build without ROOT, but ROOT is recommended if you need plot or figure output.

1. Install ROOT on macOS:

```bash
brew install root
```

2. Install ROOT on Linux:

```bash
# Option A: package manager (if available on your distro)
sudo apt update
sudo apt install -y root-system

# Option B: conda-forge
conda create -n emmap-root -c conda-forge root
conda activate emmap-root
```

3. Install ROOT on Windows:

```powershell
conda create -n emmap-root -c conda-forge root
conda activate emmap-root
```

4. Verify ROOT:

```bash
root-config --version
root-config --prefix
```

On Windows + conda:

```powershell
where.exe root.exe
root.exe -b -q
```

Without ROOT, build still succeeds, but ROOT-based plotting paths are compiled out.

## Boost

Boost is optional. In `AUTO` mode, CMake enables Boost-backed features only when Boost is found. For a typical user install, the default `AUTO` mode is the safest choice.

## Installation

Choose the workflow that matches your environment. Complete the platform-specific prerequisite section above first, then follow only one of the installation flows below.

### macOS and Linux end users

For most macOS and Linux users, installing under `~/.local` is the simplest option because it avoids administrator privileges and keeps the installed files in your home directory.

1. Configure and build:

```bash
cmake -S . -B build-local \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build-local -j
```

2. Verify the executable from the build tree:

```bash
./build-local/RHBM-GEM --help
```

3. Install to your user-local prefix:

```bash
cmake --install build-local
```

4. Verify the installed executable:

```bash
"$HOME/.local/bin/RHBM-GEM" --help
```

### macOS and Linux administrators or shared systems

If you are installing for a shared machine and want CMake's default install prefix instead of a user-local path, use a separate build directory for that workflow.

1. Configure and build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

2. Verify the executable from the build tree:

```bash
./build/RHBM-GEM --help
```

3. Install to the default prefix:

```bash
cmake --install build
```

Depending on your platform and install prefix, this step may require administrator privileges.

### Windows users (PowerShell + Visual Studio 2022)

Use this flow if you are building from a Visual Studio 2022 developer environment or a PowerShell session with the required MSVC tools available.

1. Choose a per-user install location:

```powershell
$InstallPrefix = "$env:USERPROFILE\AppData\Local\RHBM-GEM"
```

2. Configure and build:

```powershell
cmake -S . -B build `
  -G "Visual Studio 17 2022" `
  -A x64 `
  -DUSE_SYSTEM_LIBS=OFF `
  -DRHBM_GEM_ROOT_MODE=OFF `
  -DCMAKE_INSTALL_PREFIX="$InstallPrefix"
cmake --build build --config Release
```

3. Verify the executable from the build tree:

```powershell
.\build\Release\RHBM-GEM.exe --help
```

4. Install and verify the installed executable:

```powershell
cmake --install build --config Release
& "$InstallPrefix\bin\RHBM-GEM.exe" --help
```

If you configured Windows builds with `vcpkg` in the platform setup section above, add the same toolchain arguments to the configure command in this workflow.

## Python Bindings

The canonical Python extension module name is `rhbm_gem_module`.

Recommended quickstart:

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

PYTHONPATH="${PY_SITE}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 -c "import rhbm_gem_module as rgm; print(rgm.__file__)"
```

## Python Examples

Run the examples from the source tree:

```bash
PYTHONPATH="${PY_SITE}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 examples/python/00_quickstart.py
PYTHONPATH="${PY_SITE}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 examples/python/01_end_to_end_from_test_data.py --workdir /tmp/rhbm_py_demo
```

Run the installed example:

```bash
PYTHONPATH="${PY_SITE}${PYTHONPATH:+:${PYTHONPATH}}" \
python3 "${INSTALL_PREFIX}/share/RHBM_GEM/examples/python/01_end_to_end_from_test_data.py" \
  --workdir /tmp/rhbm_py_demo_installed
```

Expected outputs for the pipeline example:

1. SQLite database: `<workdir>/demo.sqlite`
2. Simulated map files: `<workdir>/maps/sim_map_*.map`
3. Dump files: `<workdir>/dump/*`

## Troubleshooting

1. If you do not have Eigen, SQLite3, pybind11, or CLI11 installed, you can skip installing them; CMake will use the bundled versions in `third_party/`.
2. Boost has no bundled fallback. Keep `RHBM_GEM_BOOST_MODE=AUTO` or set `RHBM_GEM_BOOST_MODE=OFF` if Boost is unavailable.
3. `ModuleNotFoundError: No module named 'rhbm_gem_module'`
   Verify `BUILD_PYTHON_BINDINGS=ON`, run `cmake --install`, and ensure `PYTHONPATH` points to your actual install site-packages path.
4. Install prefix mismatch (`/usr/local` vs `$HOME/.local`)
   If you install under `/usr/local`, do not use `$HOME/.local/.../site-packages` in `PYTHONPATH`.
5. `Could not find sample model ...`
   Pass `--model /path/to/your_model.cif` explicitly.

## Verbosity Levels

The command line option `-v, --verbose` controls how much information is printed.

| Level | Name | Description |
| --- | --- | --- |
| 0 | Error | Only error messages are shown. |
| 1 | Warning | Warnings and errors are shown. |
| 2 | Notice | Notice messages are shown. |
| 3 | Info | Informational messages are shown. |
| 4 | Debug | Detailed debug output. |

Passing a larger value to `--verbose` enables all lower levels as well.
