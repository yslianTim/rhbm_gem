# Command Documentation

Use this directory for command-specific developer notes.

- Keep workflow-wide guidance in [`docs/developer/adding-a-command.md`](/docs/developer/adding-a-command.md).
- Add one page per command under `docs/developer/commands/` when a command needs design notes, caveats, or maintenance guidance beyond the shared checklist.
- [`docs/developer/commands/potential_analysis.md`](/docs/developer/commands/potential_analysis.md) is the current deep-dive example for a command that combines file-backed input, command-local preprocessing, and database persistence.
- The scaffold at [`resources/tools/developer/command_scaffold.py`](/resources/tools/developer/command_scaffold.py) writes new command-note stubs into this directory.
- Keep runnable examples under [`resources/README.md`](/resources/README.md); this directory is for narrative developer documentation only.
