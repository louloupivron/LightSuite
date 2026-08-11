# Changelog

All notable changes to the LightSuite Python package are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `lightsuite brain|spinal|multires run` with `--from`, `--through`, `--resume`
- `lightsuite brain|spinal|multires stages` checkpoint status tables
- `lightsuite multires build-manifest` for mesoSPIM / SmartSPIM pair JSON
- `lightsuite config init|explain|schema` helpers
- `lightsuite workflow list` decision guide
- Combined multires → brain guide (`docs/usage_combined_multires_brain.md`)
- Actionable `LightsuiteConfigError` messages from config loader

### Changed

- Removed `lightsuite brain refine-auto-points` and AP auto-pair filtering step
- Removed legacy `lightsuite mesospim` CLI and direct-mesoSPIM YAML workflow (use `lightsuite multires` instead)
- `pyproject.toml` PyPI metadata: authors, URLs, classifiers

### Known gaps vs MATLAB

See [Python vs MATLAB](docs/python_vs_matlab.md): cell detection, slice module, CZI reader,
spinal cohort NNMF, and GPU-accelerated detection remain MATLAB-only.

## [0.1.0] - TBD

Initial public Python release: brain, spinal, and multiresolution CLI pipelines with
JSON checkpoints, Napari GUIs, and MATLAB interoperability helpers.
