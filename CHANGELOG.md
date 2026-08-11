# Changelog

All notable changes to the LightSuite Python package are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Canonical copy-paste templates: `examples/brain_lightsheet.yaml`, `examples/spinal_cord.yaml`
- `examples/README.md` workflow index mapping templates to real sample configs
- GitHub Actions CI: ruff, pytest, MkDocs strict build
- `devtools/parity/` for MATLAB/Python parity scratch scripts (moved from repo root)

### Changed

- README quick start uses canonical example paths (no feature-branch checkout required for templates)
- `pyproject.toml` PyPI metadata: authors, URLs, classifiers

### Known gaps vs MATLAB

See [Python vs MATLAB](docs/python_vs_matlab.md): cell detection, slice module, CZI reader,
spinal cohort NNMF, and GPU-accelerated detection remain MATLAB-only.

## [0.1.0] - TBD

Initial public Python release: brain, spinal, and multiresolution CLI pipelines with
JSON checkpoints, Napari GUIs, and MATLAB interoperability helpers.
