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

Commands rely on `DataObjectManager` to load and store data objects.
The manager instance is now owned by `CommandBase` and shared with derived commands.
A command lazily creates a `DataObjectManager` and assigns a `DatabaseManager` using the configured database path.
If the option changes later on, the same `DataObjectManager` updates its `DatabaseManager` through a setter.
`DataObjectManager` is thread-safe. Access to the internal data map is guarded by a `std::shared_mutex`. Read-only operations acquire a `std::shared_lock` so multiple threads can query concurrently.
Changing the database path through `SetDatabaseManager` now also locks this mutex so updates are safe to perform from multiple threads.
The manager stores objects as `std::shared_ptr`. Functions such as `GetTypedDataObject` return shared pointers, ensuring the objects remain valid as long as callers keep a reference.
`DatabaseManager` caches DataObject DAO objects using `std::shared_ptr`. The `CreateDataObjectDAO` method returns a shared pointer that shares ownership with the cache so callers should simply hold the returned `shared_ptr` without attempting to manage lifetime themselves.