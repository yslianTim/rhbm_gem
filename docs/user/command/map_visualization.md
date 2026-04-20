# Map Visualization Command

`map_visualization` generates a 2D map slice around one atom and writes the result as a PDF heatmap.
The command is registered only when `RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE` is enabled in the build.

Related implementation files:

- [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp)
- [`include/rhbm_gem/core/command/CommandManifest.def`](/include/rhbm_gem/core/command/CommandManifest.def)
- [`src/core/command/detail/CommandRequestSchema.hpp`](/src/core/command/detail/CommandRequestSchema.hpp)
- [`src/core/command/MapVisualizationCommand.hpp`](/src/core/command/MapVisualizationCommand.hpp)
- [`src/core/command/MapVisualizationCommand.cpp`](/src/core/command/MapVisualizationCommand.cpp)
- [`include/rhbm_gem/utils/math/GridSampler.hpp`](/include/rhbm_gem/utils/math/GridSampler.hpp)
- [`src/utils/math/GridSampler.cpp`](/src/utils/math/GridSampler.cpp)
- [`src/core/command/detail/MapSampling.hpp`](/src/core/command/detail/MapSampling.hpp)

## Request Surface

`MapVisualizationRequest` inherits the shared fields from `CommandRequestBase` and adds:

- `model_file_path`
- `map_file_path`
- `atom_serial_id`
- `sampling_size`
- `window_size`

CLI options are wired in [`src/core/command/detail/CommandRequestSchema.hpp`](/src/core/command/detail/CommandRequestSchema.hpp):

- `-a, --model`
- `-m, --map`
- `-i, --atom-id`
- `-s, --sampling`
- `--window-size`

Current defaults from [`include/rhbm_gem/core/command/CommandApi.hpp`](/include/rhbm_gem/core/command/CommandApi.hpp):

- `atom_serial_id = 1`
- `sampling_size = 100`
- `window_size = 5.0`

## What The Command Produces

The command loads one model and one map, normalizes the map values, selects one target atom, constructs a local 2D sampling plane near that atom, samples interpolated map values on a square grid, and writes the result as a PDF heatmap.

The output file name comes from [`MapVisualizationCommand::BuildOutputFilePath()`](/src/core/command/MapVisualizationCommand.cpp#L281):

```text
map_slice_<model_name>_atom<atom_serial_id>.pdf
```

## End-To-End Flow

`MapVisualizationCommand` currently follows this sequence:

1. `NormalizeRequest()` validates required paths, forces `atom_serial_id > 0`, and resets invalid `sampling_size` and `window_size` values to defaults.
2. `BuildDataObject()` loads the model and map objects.
3. `RunMapObjectPreprocessing()` normalizes the map array.
4. `RunModelObjectPreprocessing()` selects all atoms and bonds so the command can derive local bond context.
5. `RunAtomMapValueSampling()` finds the requested atom, finds one usable bond attached to it, builds a local sampling plane, samples interpolated values, and renders the 2D histogram PDF.

In practice, the geometry work happens entirely inside [`MapVisualizationCommand::RunAtomMapValueSampling()`](/src/core/command/MapVisualizationCommand.cpp#L170).

## Geometry Overview

The command does not sample on the global `x-y` plane.
Instead, it builds a local plane around the chosen atom using three directions:

- `u`: an in-plane direction aligned with one bond attached to the atom
- `n`: the plane normal
- `v`: the second in-plane direction, perpendicular to `u`

The final sampling plane is the plane spanned by `u` and `v`.
Every sampled point has the form:

```text
point = atom_position + a * u + b * v
```

where `a` and `b` move across a square window.

## How `u`, `n`, And `v` Are Built

### Step 1: Choose `u` From A Bond

Inside [`RunAtomMapValueSampling()`](/src/core/command/MapVisualizationCommand.cpp#L206), the command scans the bonds connected to the target atom and picks the first non-degenerate bond vector:

```text
reference_u_vector = bond_vector
u = normalize(reference_u_vector)
```

This means the horizontal direction of the local slice is anchored to a chemically meaningful direction: one bond leaving the atom.

### Step 2: Build A Stable Plane Normal `n`

The command next needs a vector perpendicular to `u` so it can define a plane.
It first tries the global `z` axis:

```text
n = u x z
```

This is implemented in [`MapVisualizationCommand.cpp`](/src/core/command/MapVisualizationCommand.cpp#L228).

If `u` happens to be parallel to `z`, then `u x z = 0`, so that choice cannot define a plane normal.
The command then falls back to the global `x` axis:

```text
n = u x x_axis
```

This fallback is implemented in [`MapVisualizationCommand.cpp`](/src/core/command/MapVisualizationCommand.cpp#L230).

### Step 3: Reconstruct `v` Inside `GridSampler`

`MapVisualizationCommand` passes:

- `u` through `SetReferenceUVector(...)`
- `n` through `SampleMapValues(..., plane_normal = n)`

Inside [`GridSampler::GenerateSamplingPoints()`](/src/utils/math/GridSampler.cpp#L31), the sampler then does:

```text
n_unit = normalize(plane_normal)
u_proj = reference_u_vector - dot(reference_u_vector, n_unit) * n_unit
u_unit = normalize(u_proj)
v_unit = n_unit x u_unit
```

Since `u` chosen by `MapVisualizationCommand` is already perpendicular to `n`, the projection step keeps the same direction in normal cases.

## Right-Hand Rule Version

The cross products in the code follow the right-hand rule.
For vectors `a x b`, point your right-hand index finger along `a`, your middle finger along `b`, and your thumb gives the direction of the result.

The command uses:

```text
n = u x reference_axis
v = n x u
```

This order matters.

### Why `u x z` And Not `z x u`

Cross products are anti-commutative:

```text
a x b = -(b x a)
```

So:

```text
u x z = -(z x u)
```

Swapping the order flips the normal direction.
The plane itself would still be the same geometric plane, but its orientation would reverse.

That reversal matters because the next step also uses a cross product to build `v`.
If `n` flips sign, then `v = n x u` also flips sign.
So changing `u x z` to `z x u` would mirror the local coordinate frame.

### Why `v = n x u`

Once `n` and `u` are fixed, there are only two possible perpendicular directions for `v`: one is the negative of the other.

The sampler chooses:

```text
v = n x u
```

That gives a consistent right-handed local frame:

```text
u x v = n
```

You can verify it with the vector triple-product identity.
If `u` is perpendicular to `n` and both are unit vectors, then:

```text
u x (n x u) = n
```

So the frame ordering is:

```text
u, v, n
```

and not:

```text
u, -v, n
```

### Sign Changes And What They Do

The sign of each axis affects orientation as follows:

- Flip `u`: the grid is mirrored left-to-right.
- Flip `v`: the grid is mirrored top-to-bottom.
- Flip `n`: the plane stays the same, but both handedness and the sign choice for `v` change if `v` is recomputed from `n x u`.

So even when the sampled geometric plane is unchanged, changing cross-product order can rotate or mirror the final heatmap.

## Precise ASCII Sketch

The local frame intended by the current implementation is:

```text
                  +n
                   ^
                   |
                   |
                   O -----> +u
                  /
                 /
              +v
```

with:

```text
u x v = n
```

and equivalently:

```text
v = n x u
```

The square sampling window is laid out in the `u-v` plane:

```text
           v+
            ^
            |
   (-,-)    |    (+,-)
            |
   -------- O --------> u+
            |
   (-,+)    |    (+,+)
            |
```

The code iterates `v` in the outer loop and `u` in the inner loop in
[`GridSampler.cpp`](/src/utils/math/GridSampler.cpp#L64), so the flattened output is row-major in local `u-v` coordinates.

## One Concrete Example

If the selected bond points along global `+x`:

```text
u = (1, 0, 0)
z = (0, 0, 1)
n = u x z = (0, -1, 0)
v = n x u = (0, 0, 1)
```

So:

- `u` points along global `+x`
- `v` points along global `+z`
- `n` points along global `-y`

The command therefore samples an `x-z` slice through the atom, not an `x-y` slice.

## How `GridSampler` Uses The Frame

After the frame is defined, [`GridSampler::GenerateSamplingPoints()`](/src/utils/math/GridSampler.cpp#L31) creates a square of evenly spaced points:

```text
u_coord in [-window_size/2, +window_size/2]
v_coord in [-window_size/2, +window_size/2]
```

with step size:

```text
step = window_size / (sampling_size - 1)
```

and each sampled point:

```text
shift = u_coord * u_unit + v_coord * v_unit
position = reference_position + shift
```

The sampler returns:

- the Euclidean in-plane distance `|shift|`
- the 3D sample position

`SampleMapValues()` then interpolates the map value at each 3D position in
[`src/core/command/detail/MapSampling.hpp`](/src/core/command/detail/MapSampling.hpp).

## Important Interpretation Notes

- The drawn axes are labeled "`u Position`" and "`v Position`", so they are local coordinates, not global Cartesian axes.
- The histogram data is filled using the sampler's flattened output order, so the local frame orientation directly controls image orientation.
- With the default `sampling_size = 100`, the grid has an even number of points per side, so the exact atom center is between four bins rather than exactly on one bin center.
- If a target atom has no usable non-zero bond vector, the command stops with an error because it cannot define a stable local frame.
- If the selected bond is parallel to global `z`, the command switches to global `x` as the helper axis to avoid a zero cross product.

## Summary

`map_visualization` should be read as:

- choose one atom
- choose one attached bond as local `u`
- use a global helper axis to derive a plane normal `n`
- reconstruct `v` so that `u x v = n`
- sample a square in that local `u-v` plane
- interpolate map values and render them as a PDF heatmap

This is why the command is best understood as a bond-oriented local slice visualizer rather than a generic Cartesian map slicer.
