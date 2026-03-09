# Data Directory Policy

- `data/local/` is reserved for local-only datasets and large generated databases.
- Files under `data/local/` must not be committed.
- Shared test fixtures belong under `tests/data/` and must be tracked.
