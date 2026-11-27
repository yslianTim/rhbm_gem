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

`RHBM_GEM_ROOT_MODE=AUTO` (the default) enables ROOT plots when ROOT is available.
Use `ON` to require ROOT at configure time or `OFF` to disable plots.
Both modes still support CSV, metadata JSON, and bin summaries when UMAP is enabled.

## SQLite Model Contract

The database must already exist, be a regular file, and use a supported SQLite
schema (v28 only). Existing v28 results can be used directly. For earlier
schemas, rerun `potential_analysis` into a new database. `--model-key` selects
one model saved by `potential_analysis`.

UMAP uses the selected Output Targets. At least three rows must have every
enabled First UMAP feature. By default these are neighbor distance sum,
closest-neighbor distance, Signal peeling ratio, and Tail peeling ratio.
Missing required values exclude the row with a reported reason;
fewer than three remaining rows fails the command. Missing ancillary values are
written as empty CSV fields. Empty regions, incomplete Original-grid coverage,
and zero total component magnitude yield an unavailable ratio.

The adjacent `.csv.metadata.json` records estimator, feature names, peeling mode,
contributor reconstruction contract, donor policy, SQLite schema, map normalization (null when unknown), standardization, included/excluded row
counts and exclusions by reason. Its `column_semantics`
field identifies the two raw ratio columns and describes `amplitude joint`
as Gaussian amplitude A, `width joint` as the
shared Gaussian/charge-kernel width tau, and `offset joint` as the charge-kernel
coefficient C. It also describes the three Joint model response columns.
`local_fitting_features` lists the ten raw columns in CSV order, and
`local_fitting_feature_values` records that their exported values are raw.
`features` and `second_umap_features` list the selected inputs in matrix order.
`model_response_features` records the Joint point source, Signal radius in
angstroms, and spherical volume measure. The `standardization` field reports
`mad-robust-scaling` for both embeddings. `standardization_parameters` records the center, scale,
normal-consistency factor, and policies for zero MAD and constant columns.
`first_umap_zero_mad_fallback_features` and
`second_umap_zero_mad_fallback_features` list the columns that used the
nonzero-deviation median in each run.

For each selected atom, the command reconstructs ten local-fitting features.
All ten are First UMAP candidates in the following order.

| Index | Atom feature | Enabled by default |
| --- | --- | --- |
| 0 | `neighbor distance sum in 2A` | Yes |
| 1 | `distance to closest neighbor` | Yes |
| 2 | `signal peeling ratio` | Yes |
| 3 | `tail peeling ratio` | Yes |
| 4 | `amplitude joint` | No |
| 5 | `width joint` | No |
| 6 | `offset joint` | No |
| 7 | `joint model center height` | No |
| 8 | `joint signal gaussian integral` | No |
| 9 | `joint signal charge integral` | No |

`neighbor distance sum in 2A` sums the Euclidean distances to non-hydrogen
atoms in the owning model within an inclusive 2 Å radius, including unselected
atoms and excluding the current atom by identity. The result is in Å and is zero
when there are no neighbors. For neighbor distances of 1 Å and 1.5 Å, the sum is
2.5 Å. It shares the full-model non-hydrogen KD-tree used for the closest-neighbor
distance. A selected hydrogen atom also queries non-hydrogen neighbors.

`distance to closest neighbor` uses KNN to find the nearest non-hydrogen atom
in the entire owning model, including unselected atoms and excluding the current
atom by identity, without a radius limit. It is the Euclidean distance in Å;
distinct atoms at the same position have distance zero. A selected hydrogen atom
also queries non-hydrogen neighbors. This feature always excludes hydrogen
candidates; no eligible neighbor makes the feature unavailable.

The **L1 Component Occupancy Peeling Ratio** uses the saved Original-grid
voxel values `y` and their paired peeled responses `p` before own-charge
correction. Signal uses `0 ≤ r² ≤ 1 Å²` (0 through 1 Å). Tail uses
`1 < r² ≤ 4 Å²` (above 1 through 2 Å). The 1 Å boundary belongs to Signal
and the 2 Å boundary belongs to Tail.

```text
neighbor_magnitude = sum_v |y_v - p_v|
peeled_magnitude   = sum_v |p_v|
peeling_ratio      = neighbor_magnitude / (neighbor_magnitude + peeled_magnitude)
```

The numerator uses the modeled net neighbor contribution at each voxel.
Contributions from multiple neighbors can cancel within a voxel. The ratio
describes the relative magnitudes of the two modeled components and
naturally lies in [0, 1]. Signed responses and a zero signed raw sum are
valid. Cancellation can make the denominator differ from `sum_v |y_v|`.
The peeled response may contain target signal, model residuals, map noise,
and neighbor reconstruction error. Physical electron density and occupancy
probabilities require separate evidence. A zero ratio indicates zero
observable net neighbor magnitude in the region. A one indicates zero
peeled magnitude. Over-peeling severity cannot be inferred from the ratio alone.

Each region requires every Original-grid observation to have a valid,
finite raw/peeled pair from the same Joint run. Empty regions, missing
responses, or two zero component magnitudes yield an unavailable value.
The command reports missing required Signal or Tail values in metadata exclusions.
The saved pairs reflect the contributor reconstruction applied by
`potential_analysis`, which uses accepted Fitting Target profiles and
reconstructs External Halos from qualified donors. The embedding command
derives ratios from those pairs and takes Joint A, width and charge features
from the saved Joint endpoints.
Atom labels remain available for CSV and plot annotations.

Both commands compute from the same binary64 pairs; a neighbor contribution
that rounded away when forming the peeled value cannot be recovered from
those pairs. The metadata field `peeling_ratio_definition` identifies
`l1-component-occupancy-v1`, with the bounds in `peeling_ratio_range` and
`peeling_ratio_regions`. `peeling_mode` identifies how the samples were
produced. The same CSV column names now carry the raw occupancy ratios
identified by this metadata. The MAD scaling recorded in `standardization`
applies to selected UMAP inputs. The [0, 1] bound applies to the raw ratios.

`amplitude joint`, `width joint`, and `offset joint` read the original
accepted `joint.point` from Stage JSON. A is the Gaussian amplitude, τ is the
shared width in angstroms, and C is the signed charge-kernel coefficient.
The three response features use that same point. An accepted point at a width
boundary or with a Local Conditioning Warning remains eligible.

The model is $f(r)=A G(r;\tau)+C K(r;\tau)$ with a normalized 3D Gaussian $G$
and $K(r;\tau)=\operatorname{erf}(r/(\sqrt{2}\tau))/r$. Both integrals use
the Signal radius $R_s=1.0$ Å and the spherical volume measure $4\pi r^2\,dr$.

- `joint model center height` is the signed response at $r=0$.

  $$ H_0 = \frac{A}{(2\pi)^{3/2}\tau^3} + \sqrt{\frac{2}{\pi}}\frac{C}{\tau}. $$

- `joint signal gaussian integral` includes A and lies in $[0,A]$.

  $$ I_G = A\left[\operatorname{erf}\!\left(\frac{R_s}{\sqrt{2}\tau}\right)
  -\sqrt{\frac{2}{\pi}}\frac{R_s}{\tau}e^{-R_s^2/(2\tau^2)}\right]. $$

- `joint signal charge integral` includes C and preserves its sign. It measures
  integrated charge-kernel response within the Signal sphere. Interpreting it as
  physical charge, electron count, or occupancy would require a separate model.

  $$ I_C = 2\pi C\left[(R_s^2-\tau^2)\operatorname{erf}\!\left(\frac{R_s}{\sqrt{2}\tau}\right)
  +\sqrt{\frac{2}{\pi}}R_s\tau e^{-R_s^2/(2\tau^2)}\right]. $$

A missing or invalid Joint point leaves all three response features unavailable.
A response outside the finite double range leaves its own field unavailable.
An underflowing response can round to zero. The three quantities are derived
from A, τ, and C and affect distances when selected for a standardized matrix.

To change the first embedding's inputs, edit `include_in_umap` in
`kFeatureDefinitions` in `src/core/classifier/UmapEmbeddingFeatures.cpp` and
rebuild with at least one feature enabled. Enabled columns are required for
row retention. Feature reconstruction, bin statistics, and CSV output include
all ten features regardless of the switches. Rows are ordered by atom serial ID.
Input rows are not filtered by spot or residue name. Selected Joint Targets with
all required features are retained; rows with missing required features are
excluded as described above. Spot and residue names remain available for plot
classification.

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
| `-j, --jobs` | `1` | Threads for neighbor search, fuzzy-set work, spectral initialization, and mean shift when OpenMP is available. Layout optimization remains single-threaded. |
| `-o, --folder` | current directory | Output directory. |
| `-v, --verbose` | `3` | Standard RHBM-GEM verbosity level. |

The removed `-i, --input` option is not accepted and is reported as an unknown
CLI option.

## Processing and Output

Both UMAP runs scale each enabled column using its atom rows. The first
uses all retained rows after required-feature validation, and the second
uses the full column after median imputation. The median `m` is the middle
value for odd counts and the midpoint of the two middle values for even
counts. The median absolute deviation (MAD) is `median(abs(x - m))`. When
MAD is positive, scaling uses:

```text
z = (x - m) / (1.482602218505602 * MAD)
```

When MAD is zero and the column varies, the denominator uses the same
factor times the median of strictly positive absolute deviations.
This `NonzeroDeviationMedian` policy emits a warning for the column and
has different robustness properties from standard MAD. Constant columns
become zero. The first UMAP warns for constant columns and fails if every
enabled column is constant.

The command validates that standardization preserves the required floating-point
range and produces finite values before either UMAP run. It checks a
bounding-box upper bound for the Euclidean sum of squared differences and
checks for gaps whose squares would be below the normal double range. These
checks can reject extreme inputs even when the bounding-box diagonal is
unattainable by any pair of rows or other columns would separate a close pair.
The command reports numerical failures with the UMAP stage, model key,
column, and reason.

MAD scaling retains finite outliers without clipping. Those rows can have
large standardized distances. CSV feature values come directly from the
original reconstructed rows and are written with round-trip-safe precision,
as are the UMAP coordinates.

The per-atom CSV has 20 columns in order: `serial id,residue,spot`, the ten
reconstructed features, `umap x,umap y`, `mean shift x,mean shift y`, then
`bin number,second umap x,second umap y`. Its name is
`<folder>/umap_embedding_<sanitized-model-key>.csv`. ASCII letters, digits,
`.`, `_`, and `-` are retained in the key; other bytes become `_`.

Mean shift starts from the original two-dimensional UMAP coordinates and applies
100 updates per point. A KD-tree built once per embedding indexes the fixed
original points. Each update uses every original point within the two-dimensional
Euclidean radius, including points on the boundary. For the original coordinates'
bounding-box diagonal `D`, the search radius is `R = 0.1D` and the Gaussian
bandwidth is `sigma = 0.5R`. The first and second UMAP derive both values
separately from their own original coordinates.

Each query moves to the weighted mean of its neighbors' original coordinates.
With `alpha = 2.0`, each iteration computes the weight as
`exp(-4 * (d / R)^2)`, where `d` is the Euclidean distance from the current
query position to the neighbor's original position. The kernel bandwidth is
determined solely by the embedding geometry. Joint width is a separate CSV feature.
If a neighborhood is empty or its weights sum to zero, the query keeps its
current position. The command reports these events in a single warning per
embedding. Non-finite updates fail the command before CSV output. At zero radius,
the query stays fixed without division, including when the radius rounds to zero.
All retained rows participate regardless of plot visibility, including rows with
duplicate coordinates. OpenMP builds use `--jobs` threads for the independent
point trajectories. Builds without OpenMP process the trajectories serially.

### Mean shift bin summary

The command reports `UMAP mean shift bin count: N` at Info level and writes
`umap_embedding_<sanitized-model-key>_mean_shift_bins.csv` alongside the per-atom
CSV. Each final mean shift position is assigned to a square of width
`h = gamma * R = 0.05D`, where `gamma = 0.5` controls cell width independently
of Gaussian bandwidth. The origin `(ox, oy)` is the lower-left corner of the
original UMAP bounding box. The `bin x` and `bin y` columns are relative grid
indices `floor((x - ox) / h)` and `floor((y - oy) / h)`, where `(x, y)`
is the final mean shift position. Each round derives its own origin and width.
Each occupied square gets one row, sorted by `bin x` and then `bin y`. The first
column, `bin number`, numbers the rows consecutively from 1 in that order. The
per-atom CSV uses the first summary's bin numbers. Overlapping points count
individually.

If `D = 0`, all points share grid cell `(0, 0)`. At tiny nonzero scales,
binning uses a rescaled diagonal to retain a usable grid when the radius or
width would underflow. For large finite axis spans, the radius can still be
computed by scaling the spans before taking their hypotenuse. Non-finite
coordinates or axis spans that exceed the finite range fail the command.

The first summary has 20 columns. After `point count` it records the three Joint
parameter medians, then `joint model center height median of bin`,
`joint signal gaussian integral median of bin`, and
`joint signal charge integral median of bin`. The remaining columns start with
`spot class count`. Each median uses the finite values for its column. An even
sample count uses the average of the two middle values, and a column with no
finite values is blank. The spot classes match the plot legend and include all
points in the bin. The three most frequent classes are listed with their names,
counts, and count divided by `point count`. Ties are ordered by class name.
Absent ranks have a blank name and zero count and ratio.

### Second UMAP

After the first binning, the command computes a median for each of the ten
local-fitting features once per occupied bin. Each feature uses its own finite
values, excluding NaN and both infinities. An even sample count uses the midpoint
of the two middle values. A feature with no finite values has a missing median.
Each atom inherits the medians of its first bin. The first bin summary uses the
same statistics for the three Joint parameters and the three model responses.
`bin number` identifies the summary row. The geometric bin center and bin
number are excluded from the second feature matrix.

The second embedding defines sixteen candidate columns in this order. The ten
enabled columns appear in the feature matrix and `second_umap_features`
metadata in the same order.

| Candidate | Source | Local feature index | Enabled |
| --- | --- | --- | --- |
| `neighbor distance sum in 2A median of bin` | First bin median | 0 | Yes |
| `distance to closest neighbor median of bin` | First bin median | 1 | Yes |
| `signal peeling ratio median of bin` | First bin median | 2 | Yes |
| `tail peeling ratio median of bin` | First bin median | 3 | Yes |
| `amplitude joint median of bin` | First bin median | 4 | No |
| `width joint median of bin` | First bin median | 5 | No |
| `offset joint median of bin` | First bin median | 6 | No |
| `joint model center height median of bin` | First bin median | 7 | No |
| `joint signal gaussian integral median of bin` | First bin median | 8 | No |
| `joint signal charge integral median of bin` | First bin median | 9 | No |
| `amplitude joint` | Atom | 4 | Yes |
| `width joint` | Atom | 5 | Yes |
| `offset joint` | Atom | 6 | Yes |
| `joint model center height` | Atom | 7 | Yes |
| `joint signal gaussian integral` | Atom | 8 | Yes |
| `joint signal charge integral` | Atom | 9 | Yes |

To select a subset, edit `include_in_umap` in `kSecondUmapFeatureDefinitions` in
`src/core/classifier/UmapEmbeddingFeatures.cpp` and rebuild with at least one
feature enabled.

For each enabled column, the command fills missing or non-finite values with
the median of the finite atom rows, then applies the shared MAD scaler to
the full column. A bin contributes its median once per atom to both
imputation and scaling. An entirely missing column is filled with zeros
and reported in a warning. The per-atom CSV and bin summaries use the
raw feature values and the bin medians before imputation.

Both UMAP runs use the same request settings. The command then applies
mean shift to the second embedding. The updated coordinates are used for the
fourth scatter plot and binned by the same rules. Feature medians are computed
once for each of the second embedding's bins. The command writes
`umap_embedding_<sanitized-model-key>_second_umap_mean_shift_bins.csv` with the
20 columns of the first summary followed by six ratio columns in this order,
for 26 columns total.

```text
local carbon ratio,local nitrogen ratio,local oxygen ratio,global carbon ratio,global nitrogen ratio,global oxygen ratio
```

Each ratio uses the number of points in the bin with the corresponding atom
`Element` as its numerator. The local denominator is `point count`, which
includes points with other or unknown elements and counts overlapping points
individually. The global denominator is the number of atoms with that same
element among the rows retained for UMAP. For example, 10 carbon points in a
bin and 100 carbon atoms retained for UMAP give a global carbon ratio of 0.1.
For each element present in the embedding, its global ratios sum to 1 across
all bins, within floating-point precision. All six ratios range from 0 to 1
and are 0 when the bin contains no points of that element. If UMAP receives
no atoms of an element, its global ratio is 0 in every bin.

The second summary has its own bin numbering starting at 1. The command reports
`Second UMAP mean shift bin count: N` at Info level. Both summaries are written
in builds without ROOT.

When all enabled columns are constant after imputation, the command reports a warning, leaves
`second umap x` and `second umap y` blank, and writes only the header to the
second bin summary, including all six ratio columns, with a reported bin count
of zero. This replaces any earlier second summary at the same path.

### Scatter plots

When ROOT support is available, the command writes the original scatter plot
as `umap_embedding_<sanitized-model-key>.pdf` and the updated scatter plot as
`umap_embedding_<sanitized-model-key>_mean_shift.pdf`. Both use the original
plot's axis ranges. When the second embedding is available, the command writes
`umap_embedding_<sanitized-model-key>_second_umap.pdf` and
`umap_embedding_<sanitized-model-key>_second_umap_mean_shift.pdf`. Both use axis
ranges based on the displayed points from the second embedding before mean shift.
All plots use the same 10% margins, minimum-range
handling, category colors, markers, legend, and layout. Categories are `C`, `CA`,
`N`, `O`, and `O (Water)` for residue `HOH` with spot `O`; unmatched atoms are
grouped by element as `Side-Chain <element>` or `Other Unknown`. Empty categories
are omitted from the legend. Builds without ROOT write the CSV files and
metadata JSON, and report that plotting is unavailable.

The database and model are fully validated, features are reconstructed and
standardized, and the embeddings, mean shift, and binning are computed before the
CSV is opened. Database path,
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

`potential_analysis` no longer produces a local-fitting CSV. The UMAP output
CSV includes the local-fitting features alongside the embedding coordinates.
UMAP reads the saved model from the database. Replace:

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
