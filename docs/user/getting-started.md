# Getting Started

Use this guide to set up your environment, install RHBM-GEM, verify the Python bindings, and run the example workflows.

Follow the sections in this order:

1. Complete **Environment Setup** for your platform.
2. Choose one **Installation** workflow.
3. If needed, continue with **Python Bindings** and **Python Examples**.

If you need coverage, build-matrix validation, or advanced CMake settings, use [`../developer/build-and-configuration.md`](../developer/build-and-configuration.md) instead. For advanced dependency selection, feature-mode control, and custom Python install layouts, see [`../developer/build-and-configuration.md#dependency-strategy`](../developer/build-and-configuration.md#dependency-strategy), [`../developer/build-and-configuration.md#cmake-parameters`](../developer/build-and-configuration.md#cmake-parameters), and [`../developer/build-and-configuration.md#feature-mode-checks-auto--off--on`](../developer/build-and-configuration.md#feature-mode-checks-auto--off--on).

## Environment Setup

RHBM-GEM uses CMake + C++17. Choose your platform first, then install any optional dependencies you need.

| If you need... | Prepare... |
| --- | --- |
| A standard C++ build | A compiler toolchain, CMake, and the platform prerequisites for your OS |
| Python interface or Python examples | Python 3; on Linux also install `python3-dev` |
| Plots or figure output | ROOT on the platform where you will build |

Some dependencies are provider-dependent or optional:

- Core C++ dependencies (`Eigen3`, `SQLite3`, `CLI11`, and Boost) are selected by `RHBM_GEM_DEP_PROVIDER`.
- `pybind11` and Python development headers are required only when `BUILD_PYTHON_BINDINGS=ON`.
- `GTest` is required only when `BUILD_TESTING=ON`.
- Use `RHBM_GEM_DEP_PROVIDER=SYSTEM` to require system packages.
- Use `RHBM_GEM_DEP_PROVIDER=FETCH` to use pinned `FetchContent` sources (network access required during configure).
- `ROOT` is optional. If it is not available, the build still succeeds, but ROOT-based plotting paths are compiled out.
- Runtime database default path is `${HOME}/.rhbmgem/data/database.sqlite`. Set `RHBM_GEM_DATA_DIR` to change the default root.

For the full dependency policy and override flags such as `RHBM_GEM_DEP_PROVIDER`, `OpenMP_ROOT`, `Boost_ROOT`, and `Python_EXECUTABLE`, see [`../developer/build-and-configuration.md#dependency-strategy`](../developer/build-and-configuration.md#dependency-strategy) and [`../developer/build-and-configuration.md#cmake-parameters`](../developer/build-and-configuration.md#cmake-parameters).

### macOS (Homebrew)

1. Install Apple command line tools:

```bash
xcode-select --install
```

2. Install Homebrew if `brew` is not available:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

3. Install the base toolchain:

```bash
brew install cmake python libomp
```

4. Optional (Recommended when using `RHBM_GEM_DEP_PROVIDER=SYSTEM`): install system packages:

```bash
brew install eigen sqlite3 pybind11 cli11 boost
```

If you plan to keep `BUILD_TESTING=ON`, install GTest too:

```bash
brew install googletest
```

5. Optional: install these only if you need the related features:

```bash
brew install root
```

6. If you installed ROOT, verify it before running CMake:

```bash
root-config --version
root-config --prefix
```

Notes:

- Install `root` only if you need plotting or figure output.

### Linux (Ubuntu/Debian example)

1. Install the base toolchain:

```bash
sudo apt update
sudo apt install -y build-essential cmake pkg-config python3
```

2. Optional: if you want distro packages with `RHBM_GEM_DEP_PROVIDER=SYSTEM`, install the core development packages:

```bash
sudo apt install -y libsqlite3-dev libeigen3-dev libboost-dev
```

If you plan to use Python bindings, also install:

```bash
sudo apt install -y python3-dev pybind11-dev
```

If you plan to keep `BUILD_TESTING=ON`, install GTest too:

```bash
sudo apt install -y libgtest-dev
```

3. Optional: install these only if you need the related features:

If you need ROOT, install it with either your distro package manager or conda:

```bash
# Option A: package manager (if available on your distro)
sudo apt install -y root-system

# Option B: conda-forge
conda create -n emmap-root -c conda-forge root
conda activate emmap-root
```

4. If you installed ROOT, verify it before running CMake:

```bash
root-config --version
root-config --prefix
```

Notes:

- If your distro does not package `CLI11`, use `-DRHBM_GEM_DEP_PROVIDER=FETCH`.
- If you installed ROOT through conda, keep that environment active while you configure, build, and install the project.

### Windows (PowerShell + Visual Studio 2022)

1. Install the required tools:
   - Visual Studio 2022 (or Build Tools) with the `Desktop development with C++` workload
   - CMake 3.18 or newer
   - Python 3 if you plan to use Python bindings or examples
   - Git if you plan to install optional packages with `vcpkg`

2. For a basic build without system packages, configure with `-DRHBM_GEM_DEP_PROVIDER=FETCH`.

3. Optional: if you prefer `vcpkg` packages, or you want a `SYSTEM`-mode dependency set on Windows, prepare them now:

```powershell
git clone https://github.com/microsoft/vcpkg $env:USERPROFILE\vcpkg
& $env:USERPROFILE\vcpkg\bootstrap-vcpkg.bat
& $env:USERPROFILE\vcpkg\vcpkg.exe install eigen3 sqlite3 pybind11 cli11 boost:x64-windows
```

4. Optional: if you need ROOT-based plots or figure output, install ROOT with conda-forge:

```powershell
conda create -n emmap-root -c conda-forge root
conda activate emmap-root
```

5. If you installed ROOT, verify it before running CMake:

```powershell
where.exe root.exe
root.exe -b -q
```

Notes:

- The installation workflow later in this guide does not require `vcpkg` for a basic build if you set `-DRHBM_GEM_DEP_PROVIDER=FETCH`.
- If you use `vcpkg`, reuse the same toolchain arguments when you configure the project later.
- If you installed ROOT through conda, keep that environment active while you configure, build, and install the project.
- The project defaults are still documented in [`../developer/build-and-configuration.md#dependency-strategy`](../developer/build-and-configuration.md#dependency-strategy). This guide chooses the simplest Windows path instead of listing every supported configuration.

## Installation

Finish **Environment Setup** first, then choose one installation workflow and follow only that workflow.

The commands below intentionally use `RHBM_GEM_DEP_PROVIDER=FETCH`, `BUILD_TESTING=OFF`, and `BUILD_PYTHON_BINDINGS=OFF` for a predictable first install. This avoids requiring system `Eigen3`/`SQLite3`/`CLI11`/`Boost`/`GTest` packages during onboarding.
Runtime executables from any build directory are placed under `<build-dir>/bin/`.

### macOS and Linux: user-local install

For most people on macOS and Linux, `~/.local` is the easiest choice. It avoids administrator privileges and keeps the installed files in your home directory.

1. Configure and build:

```bash
cmake -S . -B build-local \
  -DCMAKE_BUILD_TYPE=Release \
  -DRHBM_GEM_DEP_PROVIDER=FETCH \
  -DBUILD_TESTING=OFF \
  -DBUILD_PYTHON_BINDINGS=OFF \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build-local -j
```

2. Verify the executable from the build tree:

```bash
./build-local/bin/RHBM-GEM --help
```

3. Install to your user-local prefix:

```bash
cmake --install build-local
```

4. Verify the installed executable:

```bash
"$HOME/.local/bin/RHBM-GEM" --help
```

### macOS and Linux: shared install

Use this workflow on a shared machine, or when you want CMake's default install prefix instead of a user-local path.

1. Configure and build:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DRHBM_GEM_DEP_PROVIDER=FETCH \
  -DBUILD_TESTING=OFF \
  -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build -j
```

2. Verify the executable from the build tree:

```bash
./build/bin/RHBM-GEM --help
```

3. Install to the default prefix:

```bash
cmake --install build
```

Note: Depending on your platform and install prefix, this step may require administrator privileges.

### Windows (PowerShell + Visual Studio 2022)

Use this workflow in a Visual Studio 2022 developer environment, or in PowerShell after the required MSVC tools are available.

1. Choose a per-user install location:

```powershell
$InstallPrefix = "$env:USERPROFILE\AppData\Local\RHBM-GEM"
```

2. Configure and build:

```powershell
cmake -S . -B build `
  -G "Visual Studio 17 2022" `
  -A x64 `
  -DRHBM_GEM_DEP_PROVIDER=FETCH `
  -DRHBM_GEM_ROOT_MODE=OFF `
  -DBUILD_TESTING=OFF `
  -DBUILD_PYTHON_BINDINGS=OFF `
  -DCMAKE_INSTALL_PREFIX="$InstallPrefix"
cmake --build build --config Release
```

3. Verify the executable from the build tree:

```powershell
.\build\RHBM-GEM.exe --help
```

4. Install and verify the installed executable:

```powershell
cmake --install build --config Release
& "$InstallPrefix\bin\RHBM-GEM.exe" --help
```

Notes:

- If you prepared Windows builds with `vcpkg` in **Environment Setup**, add the same toolchain arguments to the configure command in this workflow.
- If you installed ROOT, remove `-DRHBM_GEM_ROOT_MODE=OFF` or replace it with `-DRHBM_GEM_ROOT_MODE=AUTO`.
- This workflow intentionally uses `-DRHBM_GEM_DEP_PROVIDER=FETCH`, `-DBUILD_TESTING=OFF`, `-DBUILD_PYTHON_BINDINGS=OFF`, and `-DRHBM_GEM_ROOT_MODE=OFF` to keep the first setup simple. For the project defaults and advanced alternatives, see [`../developer/build-and-configuration.md#dependency-strategy`](../developer/build-and-configuration.md#dependency-strategy) and [`../developer/build-and-configuration.md#feature-mode-checks-auto--off--on`](../developer/build-and-configuration.md#feature-mode-checks-auto--off--on).

## CLI Quickstart (macOS and Linux)

Use [`../../resources/examples/cli/00_quickstart.sh`](../../resources/examples/cli/00_quickstart.sh) when you want a Bash example that downloads one model/map pair and runs the CLI end to end.

This script is supported on macOS and Linux only. Windows users should keep using the manual `RHBM-GEM.exe` commands from the **Installation** section above.

The script resolves the executable in this order:

1. `--executable /path/to/RHBM-GEM`
2. `RHBM_GEM_EXECUTABLE`
3. `RHBM-GEM` on `PATH`
4. Installed prefix inferred from the script location
5. Common source-tree build directories such as `build-local/bin/` and `build/bin/`

If you built from the source tree, run the quickstart from the repository root after `cmake --build`:

```bash
bash resources/examples/cli/00_quickstart.sh --workdir /tmp/rhbm_cli_demo
```

If you installed to `~/.local`, you can run the installed copy directly:

```bash
bash "$HOME/.local/share/RHBM_GEM/resources/examples/cli/00_quickstart.sh" \
  --workdir /tmp/rhbm_cli_demo_installed
```

If `~/.local/bin` is not on `PATH`, pass the executable explicitly:

```bash
bash "$HOME/.local/share/RHBM_GEM/resources/examples/cli/00_quickstart.sh" \
  --workdir /tmp/rhbm_cli_demo_installed \
  --executable "$HOME/.local/bin/RHBM-GEM"
```

For offline reruns with pre-staged input files under `<workdir>/data`, add `--skip-download`.

## Python Bindings

The Python module is named `rhbm_gem_module`.

Complete **Installation** first. Then follow the Python binding workflow for your platform. Each workflow builds the module, installs it under your chosen prefix, and confirms that Python can import it.

### macOS and Linux

For most users, reuse the same `~/.local` install prefix from **Installation**.

1. Configure, build, and install with Python bindings enabled:

```bash
INSTALL_PREFIX="$HOME/.local"

cmake -S . -B build-py \
  -DCMAKE_BUILD_TYPE=Release \
  -DRHBM_GEM_DEP_PROVIDER=FETCH \
  -DBUILD_TESTING=OFF \
  -DBUILD_PYTHON_BINDINGS=ON \
  -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}"
cmake --build build-py -j
cmake --install build-py
```

2. Find the installed `site-packages` directory, export `PYTHONPATH` for the current shell, and verify the import:

```bash
PYVER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PY_SITE_LIB="${INSTALL_PREFIX}/lib/python${PYVER}/site-packages"
PY_SITE_LIB64="${INSTALL_PREFIX}/lib64/python${PYVER}/site-packages"
PY_SITE_PATHS="${PY_SITE_LIB}:${PY_SITE_LIB64}"
export PYTHONPATH="${PY_SITE_PATHS}${PYTHONPATH:+:${PYTHONPATH}}"

python3 -c "import rhbm_gem_module as rgm; print(rgm.__file__)"
```

Note: If you installed the project to a different prefix, update `INSTALL_PREFIX` before computing `PY_SITE_PATHS`.
Keep this shell session open for **Python Examples**. If you start a new shell later, rerun the `export PYTHONPATH=...` command first.

For custom Python layouts such as `RHBM_GEM_PYTHON_INSTALL_LAYOUT=LIBDIR` or an explicit `RHBM_GEM_PYTHON_INSTALL_DIR`, see [`../developer/build-and-configuration.md#cmake-parameters`](../developer/build-and-configuration.md#cmake-parameters) and [`../developer/build-and-configuration.md#validation-examples`](../developer/build-and-configuration.md#validation-examples).

### Windows (PowerShell)

Reuse the same install prefix from **Installation**, unless you intentionally chose a different location.

1. Configure, build, and install with Python bindings enabled:

```powershell
$InstallPrefix = "$env:USERPROFILE\AppData\Local\RHBM-GEM"

cmake -S . -B build-py `
  -G "Visual Studio 17 2022" `
  -A x64 `
  -DRHBM_GEM_DEP_PROVIDER=FETCH `
  -DRHBM_GEM_ROOT_MODE=OFF `
  -DBUILD_TESTING=OFF `
  -DBUILD_PYTHON_BINDINGS=ON `
  -DCMAKE_INSTALL_PREFIX="$InstallPrefix"
cmake --build build-py --config Release
cmake --install build-py --config Release
```

2. Find the installed `site-packages` directory and verify the import:

```powershell
$PyVer = python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
$PySite = "$InstallPrefix\lib\python$PyVer\site-packages"
$env:PYTHONPATH = if ($env:PYTHONPATH) { "$PySite;$env:PYTHONPATH" } else { $PySite }

python -c "import rhbm_gem_module as rgm; print(rgm.__file__)"
```

Notes:

- If you prepared Windows builds with `vcpkg` in **Environment Setup**, add the same toolchain arguments to the configure command here.
- If you installed ROOT, remove `-DRHBM_GEM_ROOT_MODE=OFF` or replace it with `-DRHBM_GEM_ROOT_MODE=AUTO`.
- For custom Python install layouts, see [`../developer/build-and-configuration.md#cmake-parameters`](../developer/build-and-configuration.md#cmake-parameters).

## Python Examples

Once the import check succeeds, continue with the example workflow for your platform. Each workflow runs the quickstart script, then the end-to-end pipeline from the source tree, and finally the installed copy of the same pipeline example.

### macOS and Linux

These commands assume `INSTALL_PREFIX`, `PY_SITE_PATHS`, and `PYTHONPATH` are still set in the same shell from **Python Bindings**.

1. Run the quickstart example from the source tree:

```bash
python3 resources/examples/python/00_quickstart.py
```

2. Run the end-to-end pipeline example from the source tree:

```bash
python3 resources/examples/python/01_end_to_end_from_test_data.py --workdir /tmp/rhbm_py_demo
```

3. Run the installed copy of the pipeline example:

```bash
python3 "${INSTALL_PREFIX}/share/RHBM_GEM/resources/examples/python/01_end_to_end_from_test_data.py" \
  --workdir /tmp/rhbm_py_demo_installed
```

### Windows (PowerShell)

These commands assume `$InstallPrefix`, `$PySite`, and `$env:PYTHONPATH` are still set in the same PowerShell session from **Python Bindings**.

1. Run the quickstart example from the source tree:

```powershell
python .\resources\examples\python\00_quickstart.py
```

2. Run the end-to-end pipeline example from the source tree:

```powershell
python .\resources\examples\python\01_end_to_end_from_test_data.py --workdir "$env:TEMP\rhbm_py_demo"
```

3. Run the installed copy of the pipeline example:

```powershell
python "$InstallPrefix\share\RHBM_GEM\resources\examples\python\01_end_to_end_from_test_data.py" `
  --workdir "$env:TEMP\rhbm_py_demo_installed"
```

The pipeline example should create:

1. SQLite database: `<workdir>/demo.sqlite`
2. Simulated map outputs under `<workdir>/maps/`
3. Dump outputs under `<workdir>/dump/`

Compatibility note:

- The installed legacy path `share/RHBM_GEM/examples/python/...` is still shipped temporarily for compatibility, but `share/RHBM_GEM/resources/examples/python/...` is now the canonical location.

## Troubleshooting

1. Missing `Eigen3`, `SQLite3`, `pybind11`, `CLI11`, or `Boost`
   Fallback sources are used only with `-DRHBM_GEM_DEP_PROVIDER=FETCH`. In `SYSTEM` mode, install the required packages first.
2. Missing `GTest` during configure
   If you do not need tests, set `-DBUILD_TESTING=OFF` (the installation workflows in this guide already do this).
3. `ModuleNotFoundError: No module named 'rhbm_gem_module'`
   Verify `BUILD_PYTHON_BINDINGS=ON`, run `cmake --install`, and make sure `PYTHONPATH` points to the actual install `site-packages` path.
4. Install prefix mismatch
   If you install under a different prefix, recompute the `site-packages` path from that prefix instead of reusing the examples for `~/.local` or `%USERPROFILE%\AppData\Local\RHBM-GEM`.
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
