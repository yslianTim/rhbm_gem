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
- Evaluate fold-168 parameter errors against recorded simulation charges with the [truth scorer](docs/developer/build-and-configuration.md#fold-168-parameter-truth-scoring); its new accuracy thresholds are not yet calibrated.
- Use the default-enabled `umap_embedding` command. In `SYSTEM` mode, installed UMAP packages are preferred and only missing UMAP components are fetched at fixed versions; disable it with `RHBM_GEM_ENABLE_UMAP=OFF` to avoid UMAP resolution and downloads.
- Start with [`docs/developer/README.md`](/docs/developer/README.md) and [`docs/developer/build-and-configuration.md`](/docs/developer/build-and-configuration.md).

## Examples

Browse runnable examples and helper assets in [`resources/README.md`](/resources/README.md).

## Joint component workflow (opt-in)

The default estimator remains the two-stage workflow. New results use **SQLite v18**. Valid v17 databases are read without modification
and upgraded transactionally on their first write; v16 and earlier are rejected.

```sh
RHBM-GEM potential_analysis --estimator joint-components -a model.cif -m map.mrc -d joint.sqlite -k example
RHBM-GEM result_dump --printer joint -d joint.sqlite -k example -o results
```

Joint output contains A/B/C, component availability, runtime convergence and
separate offline-certificate status. A successful command means the outcome was
saved, including incomplete or unconverged outcomes. See the
[joint workflow contract](docs/developer/commands/potential_analysis.md#joint-component-opt-in)
for options, C++/Python examples and JSON/CSV output. Joint endpoints populate the
common Second stage. Plots, comparisons, Gaussian/outlier dumps and UMAP use the
available saved data; missing uncertainty or peeling coverage is never replaced
with zero.

## License, Citation, and Third-Party Notices

- Project-owned source code, build scripts, examples, and documentation in this repository are licensed under the MIT License (see [`LICENSE`](/LICENSE)).
- Third-party licenses and attribution summary: see [`THIRD_PARTY_NOTICES.md`](/THIRD_PARTY_NOTICES.md).
- If you use this work in research, cite the manuscript linked above; software license terms and scholarly citation expectations are separate.
- Release and compliance guidance for maintainers: [`docs/developer/release-compliance.md`](/docs/developer/release-compliance.md).
