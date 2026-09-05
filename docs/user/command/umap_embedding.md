# `umap_embedding`

`umap_embedding` loads one saved model from SQLite, reconstructs its local-fitting
features, and creates a two-dimensional UMAP embedding. The stable command is
enabled by default and is omitted only when the project is configured with
`RHBM_GEM_ENABLE_UMAP=OFF`.

## Build Configuration

The simplest self-contained configuration uses the pinned FETCH dependencies. No
UMAP option is needed because it defaults to `ON`:

```bash
cmake -S . -B build-umap \
  -DRHBM_GEM_DEP_PROVIDER=FETCH
cmake --build build-umap --target rhbm_gem_cli
```

SYSTEM builds prefer an installed umappp 3.3.2 stack and fetch only missing UMAP
components at fixed versions; Eigen3 and the other core dependencies must still
be installed. A configure needs network access if a missing UMAP archive is not
cached. See the [build and configuration guide](/docs/developer/build-and-configuration.md#dependency-strategy),
or set `RHBM_GEM_ENABLE_UMAP=OFF` to build without the command, UMAP lookup, or
fallback downloads.

## SQLite Model Contract

The database must already exist, be a regular file, and use the currently
supported SQLite schema (schema v15). `--model-key` must identify exactly one
model saved by `potential_analysis`; the command neither enumerates nor combines
models.

The saved model must contain at least three selected atoms and complete persisted
local-fitting analysis for every selected atom:

- raw and peeling samples sufficient to calculate both `[0, 1)` Å and
  `[1, 2]` Å peeling ratios;
- the peeling neighbor count;
- finite first-, second-, and third-stage MDPDE amplitude, width, and offset.

Missing samples, an undefined ratio, or any non-finite reconstructed feature
makes the model incomplete and fails the command. No atom is silently skipped.

For each selected atom, the command reconstructs the same 22 features used by
the local-fitting CSV:

```text
neighbor count for peeling,neighbor count in 2A,signal peeling ratio,tail peeling ratio,amplitude 1st,amplitude 2nd,amplitude 3rd,width 1st,width 2nd,width 3rd,offset 1st,offset 2nd,offset 3rd,amplitude rank 1st,amplitude rank 2nd,amplitude rank 3rd,width rank 1st,width rank 2nd,width rank 3rd,offset rank 1st,offset rank 2nd,offset rank 3rd
```

`neighbor count in 2A` counts all atoms in the owning model within an
inclusive 2 Å radius and excludes the current atom. Each peeling ratio is
`(raw sum - peeling sum) / raw sum` in its distance range. Each rank compares
the current atom with up to its three nearest selected atoms; the largest value
has rank 1, equal values share a rank, and models with fewer than four selected
atoms use all available atoms.

Rows are ordered by atom serial ID. The current build passes these 12 features
to UMAP:

- peeling neighbor count, 2 Å neighbor count, and both peeling ratios;
- first- and second-stage amplitude;
- first- and second-stage width;
- second-stage offset;
- second-stage amplitude, width, and offset ranks.

The other ten reconstructed features remain in the result CSV but are not passed
to UMAP. With the default internal setting `kFilterUmapInputBySpot = false`,
all selected atoms are retained. Developers can enable the filter and rebuild to
retain only `C`, `CA`, `N`, `O`, and residue `HOH` with spot `O`.

## CLI Usage

```bash
RHBM-GEM umap_embedding \
  --database analysis.sqlite \
  --model-key model \
  --folder embedding_output \
  --neighbors 15 \
  --min-dist 0.1 \
  --epochs 0 \
  --seed 42 \
  --jobs 1
```

| Option | Default | Meaning |
| --- | --- | --- |
| `-d, --database` | platform data directory `database.sqlite` | Existing SQLite database containing the saved model. |
| `-k, --model-key` | required | Key of the single saved model to embed. |
| `--neighbors` | `15` | Nearest neighbors; must be at least 2 and is limited to `retained row count - 1` when necessary. |
| `--min-dist` | `0.1` | Minimum embedding distance in the inclusive range `[0, 1]`. |
| `--epochs` | `0` | Optimization epochs; `0` selects the umappp automatic value. |
| `--seed` | `42` | Seed used for initialization and layout optimization. |
| `-j, --jobs` | `1` | Threads for neighbor search, fuzzy-set work, and spectral initialization. Layout optimization remains single-threaded. |
| `-o, --folder` | current directory | Output directory. |
| `-v, --verbose` | `3` | Standard RHBM-GEM verbosity level. |

The removed `-i, --input` option is not accepted and is reported as an unknown
CLI option.

## Processing and Output

Each selected feature is independently standardized with a population Z-score
before UMAP runs. A constant selected feature is replaced with zeros and
reported as a warning. If all 12 selected features are constant, the command
fails without writing an output file.

Unlike the legacy CSV input path, features are restored directly from SQLite as
full `double` values and are not quantized to two decimal places. The output
uses round-trip-safe floating-point precision for continuous features and UMAP
coordinates.

The CSV contains the identifiers `serial id,residue,spot`, all 22 reconstructed
features, and `umap x,umap y`: 27 columns in total. Its name is
`<folder>/umap_embedding_<sanitized-model-key>.csv`. ASCII letters, digits,
`.`, `_`, and `-` are retained in the key; other bytes become `_`.

When ROOT support is available, the command also writes a PDF with the same
stem. Categories `C`, `CA`, `N`, and `O` retain their existing colors;
residue `HOH` with spot `O` is plotted separately as `HOH`; unmatched atoms
are plotted in the gray `Other` category. Empty categories are omitted from
the legend. Builds without ROOT skip the PDF and report a warning.

The database and model are fully validated, features are reconstructed and
standardized, and UMAP completes before the CSV is opened. Database path,
opening, and schema failures are attributed to `--database`; missing keys and
incomplete model features to `--model-key`; output failures to `--folder`.
UMAP coordinates are meaningful relative to other rows in the same run; their
absolute rotation, reflection, and scale have no fixed interpretation.

## C++ and Python

The enabled build exposes the same request through C++:

```cpp
rhbm_gem::core::UmapEmbeddingRequest request;
request.database_path = "analysis.sqlite";
request.model_key_tag = "model";
request.output_dir = "embedding_output";
const auto result = rhbm_gem::core::RunCommand(request);
```

The Python binding uses the same executor:

```python
import rhbm_gem_module as rgm

request = rgm.UmapEmbeddingRequest()
request.database_path = "analysis.sqlite"
request.model_key_tag = "model"
request.output_dir = "embedding_output"
result = rgm.RunCommand(request)
```

For reproducible reruns, keep the saved model, UMAP parameters, seed, and job
count unchanged. Changing the spectral-initialization thread count can change
floating-point round-off and therefore the final coordinates.

## Migrating from CSV Input

The local-fitting CSV produced by `potential_analysis` remains available for
reporting, but UMAP no longer reads it. Replace:

```bash
RHBM-GEM umap_embedding --input local_fitting_result_model.csv
```

with the database and key used when the analysis was saved:

```bash
RHBM-GEM umap_embedding \
  --database analysis.sqlite \
  --model-key model
```

There is no compatibility alias or automatic CSV import. If only the old CSV
remains, rerun `potential_analysis` to persist a complete model in the current
SQLite schema before embedding it.
