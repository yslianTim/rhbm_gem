# `umap_embedding`

`umap_embedding` converts one `local_fitting_result_*.csv` file into a two-dimensional UMAP embedding. The command is stable, but is available only in builds configured with `RHBM_GEM_ENABLE_UMAP=ON`.

## Enable the Command

The simplest self-contained configuration uses the pinned FETCH dependencies:

```bash
cmake -S . -B build-umap \
  -DRHBM_GEM_DEP_PROVIDER=FETCH \
  -DRHBM_GEM_ENABLE_UMAP=ON
cmake --build build-umap --target rhbm_gem_cli
```

SYSTEM builds require umappp 3.3.2 and its package dependencies as described in the [build and configuration guide](/docs/developer/build-and-configuration.md#dependency-strategy).

## Input Contract

The input must use this exact 14-column header and order:

```text
serial id,residue,spot,neighbor count,peeling ratio,amplitude 1st,amplitude 2nd,amplitude 3rd,width 1st,width 2nd,width 3rd,offset 1st,offset 2nd,offset 3rd
```

`serial id`, `residue`, and `spot` are preserved as identifiers. UMAP uses the remaining 11 columns:

- neighbor count and peeling ratio
- first-, second-, and third-stage amplitude
- first-, second-, and third-stage width
- first-, second-, and third-stage offset

Every data row must contain exactly 14 comma-separated fields. The 11 feature fields must parse completely as finite numbers; empty values, `nan`, and infinity are errors. Quoted fields and commas embedded inside fields are not supported. Both LF and CRLF line endings are accepted, and at least three data rows are required.

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

Each feature is independently standardized with a population Z-score before UMAP runs. A constant feature is replaced with zeros and reported as a warning. If all 11 features are constant, the command fails without writing an output file.

The Euclidean VP-tree UMAP result has two dimensions. Input row order and all original values are preserved, and the output appends `umap x,umap y` using round-trip-safe floating-point precision.

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
