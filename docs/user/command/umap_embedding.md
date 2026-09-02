# `umap_embedding`

`umap_embedding` converts one `local_fitting_result_*.csv` file into a two-dimensional UMAP embedding. The stable command is enabled by default and is omitted only when the project is configured with `RHBM_GEM_ENABLE_UMAP=OFF`.

## Build Configuration

The simplest self-contained configuration uses the pinned FETCH dependencies. No UMAP option is needed because it defaults to `ON`:

```bash
cmake -S . -B build-umap \
  -DRHBM_GEM_DEP_PROVIDER=FETCH
cmake --build build-umap --target rhbm_gem_cli
```

SYSTEM builds prefer an installed umappp 3.3.2 stack and fetch only missing UMAP components at fixed versions; Eigen3 and the other core dependencies must still be installed. A configure needs network access if a missing UMAP archive is not cached. See the [build and configuration guide](/docs/developer/build-and-configuration.md#dependency-strategy), or set `RHBM_GEM_ENABLE_UMAP=OFF` to build without the command, UMAP lookup, or fallback downloads.

## Input Contract

The input must use this exact 23-column header and order:

```text
serial id,residue,spot,neighbor count,peeling ratio,amplitude 1st,amplitude 2nd,amplitude 3rd,width 1st,width 2nd,width 3rd,offset 1st,offset 2nd,offset 3rd,amplitude rank 1st,amplitude rank 2nd,amplitude rank 3rd,width rank 1st,width rank 2nd,width rank 3rd,offset rank 1st,offset rank 2nd,offset rank 3rd
```

`serial id`, `residue`, and `spot` are preserved as identifiers. All remaining 20 columns are parsed as features. The current build passes these 16 features to UMAP:

- neighbor count and peeling ratio
- first- and second-stage amplitude
- first- and second-stage width
- second-stage offset
- first-, second-, and third-stage amplitude, width, and offset ranks

Each rank compares the current atom with up to its three nearest selected atoms. The largest parameter value has rank 1, equal values share a rank, and models with fewer than four selected atoms use all available atoms.

The third-stage amplitude and width plus the first- and third-stage offset are preserved but are not currently passed to UMAP. Every data row must contain exactly 23 comma-separated fields. The 20 feature fields must parse completely as finite numbers; empty values, `nan`, and infinity are errors. Quoted fields and commas embedded inside fields are not supported. Both LF and CRLF line endings are accepted, and at least three data rows are required.

## CLI Usage

```bash
RHBM-GEM umap_embedding \
  --input local_fitting_result_model.csv \
  --folder embedding_output \
  --neighbors 15 \
  --min-dist 0.1 \
  --epochs 0 \
  --seed 42 \
  --jobs 1
```

| Option | Default | Meaning |
| --- | --- | --- |
| `-i, --input` | required | One local fitting result CSV file. |
| `--neighbors` | `15` | Nearest neighbors; must be at least 2 and is limited to `row count - 1` when necessary. |
| `--min-dist` | `0.1` | Minimum embedding distance in the inclusive range `[0, 1]`. |
| `--epochs` | `0` | Optimization epochs; `0` selects the umappp automatic value. |
| `--seed` | `42` | Seed used for both initialization and layout optimization. |
| `-j, --jobs` | `1` | Threads for neighbor search, fuzzy-set work, and spectral initialization. Layout optimization remains single-threaded. |
| `-o, --folder` | current directory | Output directory. |
| `-v, --verbose` | `3` | Standard RHBM-GEM verbosity level. |

## Processing and Output

Each selected feature is independently standardized with a population Z-score before UMAP runs. A constant selected feature is replaced with zeros and reported as a warning. If all 16 selected features are constant, the command fails without writing an output file.

The Euclidean VP-tree UMAP result has two dimensions. Input row order and all 23 original values are preserved, and the output appends `umap x,umap y` using round-trip-safe floating-point precision, producing 25 columns per row.

For `local_fitting_result_model.csv`, the output is `<folder>/umap_embedding_model.csv`. Other input names use their complete stem, for example `custom.csv` becomes `umap_embedding_custom.csv`.

The output is opened only after parsing, standardization, and UMAP all succeed. UMAP coordinates are meaningful relative to other rows in the same run; their absolute rotation, reflection, and scale have no fixed interpretation.

## C++ and Python

The enabled build exposes the same request through C++:

```cpp
rhbm_gem::core::UmapEmbeddingRequest request;
request.input_csv_path = "local_fitting_result_model.csv";
request.output_dir = "embedding_output";
const auto result = rhbm_gem::core::RunCommand(request);
```

The Python binding uses the same executor:

```python
import rhbm_gem_module as rgm

request = rgm.UmapEmbeddingRequest()
request.input_csv_path = "local_fitting_result_model.csv"
request.output_dir = "embedding_output"
result = rgm.RunCommand(request)
```

For reproducible reruns, keep the input row order, feature values, UMAP parameters, seed, and job count unchanged. Changing the spectral-initialization thread count can change floating-point round-off and therefore the final coordinates.
