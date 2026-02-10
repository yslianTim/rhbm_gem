# em_map_analysis

## Getting Started (macOS / Linux)

This project uses CMake + C++17. By default it prefers system-installed Eigen3 / CLI11 / pybind11 / SQLite3; if any are missing, CMake automatically falls back to the bundled copies in `third_party/`. To force bundled deps, pass `-DUSE_SYSTEM_LIBS=OFF` during configure.

**macOS (Homebrew)**
1. Install Xcode Command Line Tools:
```bash
xcode-select --install
```
1. Install dependencies (adjust as needed):
```bash
brew install cmake eigen sqlite3 python pybind11 cli11 libomp
```

**Linux (Ubuntu/Debian example)**
1. Install dependencies (adjust as needed):
```bash
sudo apt update
sudo apt install -y build-essential cmake pkg-config python3 python3-dev libsqlite3-dev libeigen3-dev
```
1. If you want to build Python bindings, also install:
```bash
sudo apt install -y pybind11-dev
```

**Build and run from scratch**
1. From the project root:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```
1. Run the main executable:
```bash
./build/PotentialMap_run --help
```

**Troubleshooting**
1. If you do not have Eigen / SQLite3 / pybind11 / CLI11 installed, you can skip installing them; CMake will use the bundled versions in `third_party/`.
1. To force bundled deps, configure with:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_LIBS=OFF
```

**OpenMP mode checks (both OFF and ON)**
1. Force OpenMP OFF (macOS / Linux):
```bash
cmake -S . -B build-openmp-off -DCMAKE_BUILD_TYPE=Release -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE
cmake --build build-openmp-off -j
./build-openmp-off/PotentialMap_run --help
```
1. Force OpenMP ON on macOS (Homebrew `libomp`):
```bash
cmake -S . -B build-openmp-on -DCMAKE_BUILD_TYPE=Release -DOpenMP_ROOT=/opt/homebrew/opt/libomp -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/libomp
cmake --build build-openmp-on -j
./build-openmp-on/PotentialMap_run --help
```
1. Force OpenMP ON on Linux (no `OpenMP_ROOT` needed):
```bash
cmake -S . -B build-openmp-on -DCMAKE_BUILD_TYPE=Release
cmake --build build-openmp-on -j
./build-openmp-on/PotentialMap_run --help
```
1. If OpenMP is missing on Linux, install:
```bash
sudo apt install -y libomp-dev
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




# Overall Information for Developer


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
