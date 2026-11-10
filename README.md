# Robust Hierarchical Linear (Bayesian) Model for Gaussian Estimation of cryo-EM Maps (RHBM-GEM)

RHBM-GEM provides a robust hierarchical Bayesian framework to extract atom-specific signal features (amplitude and width) from cryo-EM maps, explicitly accounting for noise and contamination via Gaussian modeling.

Read the manuscript of method here: [https://doi.org/10.1101/2025.07.10.664269](https://doi.org/10.1101/2025.07.10.664269).

## Documentation

For the full documentation map, start at [`docs/README.md`](/docs/README.md).

### User Documentation

- Set up your environment and install RHBM-GEM.
- Verify the Python bindings.
- Run the supported example workflows.
- Start with [`docs/user/README.md`](/docs/user/README.md) and [`docs/user/getting-started.md`](/docs/user/getting-started.md).

### Developer Documentation

- Change the repository and follow project conventions.
- Extend commands, data I/O, or other internal components.
- Validate build configurations and prepare releases.
- Use the default-enabled `umap_embedding` command. In `SYSTEM` mode, installed UMAP packages are preferred and only missing UMAP components are fetched at fixed versions; disable it with `RHBM_GEM_ENABLE_UMAP=OFF` to avoid UMAP resolution and downloads.
- Start with [`docs/developer/README.md`](/docs/developer/README.md) and [`docs/developer/build-and-configuration.md`](/docs/developer/build-and-configuration.md).

## Examples

Browse runnable examples and helper assets in [`resources/README.md`](/resources/README.md).

## Joint component workflow

Potential Analysis uses the Joint Component Estimator (JCE). Results use
**SQLite v23**. The repository accepts only the current database schema.

```sh
RHBM-GEM potential_analysis -a model.cif -m map.mrc -d joint.sqlite -k example
RHBM-GEM result_dump --printer joint -d joint.sqlite -k example -o results
```

Joint output contains A/tau/C, component availability and runtime convergence.
The C++/Python `CommandResult.succeeded` field reports whether execution
completed and the result was saved. Its `outcome` distinguishes initialization
failure, an unqualified fit, unavailable Halo reconstruction, and a qualified
fit. CLI exit codes are `1` for execution or saving failure, `2` when Joint
initialization failed after saving diagnostics, `3` when Joint fitting
did not qualify after saving diagnostics, `4` when Joint qualified but Halo
reconstruction was unavailable, and `0` for a qualified Joint fit. See the
[joint workflow contract](docs/developer/commands/potential_analysis.md#joint-component-workflow)
for options and JSON/CSV output. JointStage stores the fitted endpoints. Plots, comparisons, Gaussian/outlier dumps and UMAP use the
available saved data; missing uncertainty or peeling coverage is never replaced
with zero. Qualified Target A/tau/C points provide geometry-local kNN-4 Halo
medians for downstream peeling. The current formats are Joint JSON schema 25
and Stage JSON version 7. See [Joint Stage Integration](docs/developer/joint-stage-integration.md).

## License, Citation, and Third-Party Notices

- Project-owned source code, build scripts, examples, and documentation in this repository are licensed under the MIT License (see [`LICENSE`](/LICENSE)).
- Third-party licenses and attribution summary: see [`THIRD_PARTY_NOTICES.md`](/THIRD_PARTY_NOTICES.md).
- If you use this work in research, cite the manuscript linked above; software license terms and scholarly citation expectations are separate.
- Release and compliance guidance for maintainers: [`docs/developer/release-compliance.md`](/docs/developer/release-compliance.md).
