# em_map_analysis

## Verbosity Levels

The command line option `-v, --verbose` controls how much information is printed.
The new logger supports four levels:

| Level |  Name   | Description                       |
|-------|---------|-----------------------------------|
| 0     | Error   | Only error messages are shown.    |
| 1     | Warning | Warnings and errors are shown.    |
| 2     | Info    | Informational messages are shown. |
| 3     | Debug   | Detailed debug output.            |

Passing a larger value to `--verbose` enables all lower levels as well.