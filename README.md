# em_map_analysis

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



## Tests

The test suite automatically includes any immediate subdirectory of `tests/` containing a `CMakeLists.txt`. To add a new set of tests, create a new folder under `tests` with its own `CMakeLists.txt`.
