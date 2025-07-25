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
`DataObjectManager` is thread-safe. All operations that access the internal data object map are synchronized using a mutex so commands may share the manager safely across threads.
`DatabaseManager` caches DataObject DAO objects using `std::shared_ptr`. The `CreateDataObjectDAO` method returns a shared pointer that shares ownership with the cache so callers should simply hold the returned `shared_ptr` without attempting to manage lifetime themselves.