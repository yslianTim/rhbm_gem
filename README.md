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



## Data Object Manager

Commands rely on `DataObjectManager` to load and store data objects. The manager instance is owned by `CommandBase` and shared with derived commands. A command lazily creates a `DataObjectManager` and assigns or injects a `DatabaseManager` using `SetDatabaseManager`.

`DataObjectManager` is thread-safe. The database pointer and the internal data map are protected by separate `std::shared_mutex` instances. Read-only operations acquire a `std::shared_lock` so multiple threads can query concurrently. When an operation needs both locks, the database mutex is taken before the map mutex to avoid deadlock.

The manager stores objects as `std::shared_ptr`. Functions such as `GetTypedDataObject` return shared pointers, ensuring the objects remain valid as long as callers keep a reference. File loading and writing are delegated to `FileIOManager`, which simplifies the class and improves testability. `DatabaseManager` caches DataObject DAO objects using `std::shared_ptr`; callers should simply hold the returned pointers without managing lifetime themselves.

`DatabaseManager` is safe for concurrent use. Operations like `SaveDataObject` and `LoadDataObject` lock an internal mutex so multiple threads can access the database simultaneously without interference.

### File IO and persistence work through cooperating managers:

- **FileIOManager** selects a reader or writer based on file extension. It uses `FileProcessFactoryRegistry` to instantiate the appropriate factory and owns no persistent state.
- **DatabaseManager** wraps the SQLite database. DAO objects are created via `DataObjectDAOFactoryRegistry` and cached for reuse.
- **DataObjectManager** coordinates the above helpers. Commands load and save objects through it, delegating file operations to `FileIOManager` and database access to `DatabaseManager`.

Together these managers keep data loading and persistence separate from command logic.