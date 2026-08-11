# Python vs MATLAB

The Python pipeline (`lightsuite` CLI) covers **brain lightsheet registration**, **multiresolution overview ↔ ROI alignment**, and **spinal cord registration** with export and per-sample region statistics. The original MATLAB workflows remain in the repository for features not yet ported.

Use this page when deciding which runtime to use for a given task.

---

## Not implemented in Python

These capabilities exist in the MATLAB LightSuite demos but **have no Python equivalent today**. Use the MATLAB entry scripts listed in the repository `README.md`.

| Feature | MATLAB entry point | Notes |
|---------|-------------------|--------|
| **Built-in 3D cell detection** | `demos/ls_analyze_lightsheet_volume.m` | Band-pass filtering + local maxima; maps cells to atlas regions. In Python, set `detection.enabled: false` and import external spot lists with [`brain import-annotations`](annotation_import.md). |
| **Widefield coronal slice module** | `demos/ls_analyze_slice_volume.m` | Per-slice 2D registration for widefield / AxioScan data. See [Slices (MATLAB)](usage_slice.md). |
| **CZI volume reader** | MATLAB slice / import paths | Python accepts **TIFF stacks** (`channelperfile`, `planeperfile`). CZI must be exported to TIFF first. |
| **Spinal cord cohort analysis** | `example_analysis_spinal_cord.m` | Cross-subject NNMF normalization and cohort plotting. Python provides per-sample `spinal region-stats` only. |
| **GPU-accelerated cell detection** | MATLAB `cell_counting/` | `compute.use_gpu` in YAML is ignored by the Python pipeline. |

---

## Python scope (what is supported)

| Workflow | CLI namespace | Status |
|----------|---------------|--------|
| Brain lightsheet | `lightsuite brain` | Preprocess → orientation → init → match-points → register → export |
| Multiresolution (mesoSPIM / SmartSPIM) | `lightsuite multires` | Manifest-driven overview ↔ ROI registration |
| Spinal cord | `lightsuite spinal` | Preprocess → straighten → register → export → `region-stats` |
| External annotations | `lightsuite brain import-annotations` | Native `points_csv` / `mask_tiff` at sample resolution |

Registered volumes are written as **multi-page TIFF** (`export.registered_volume_format: tiff`). OME-Zarr output is not supported.

---

## Intentional differences (brain B-spline)

When comparing Python and MATLAB registration numerically, see [Migrating from MATLAB](usage_lightsheet_brain.md#migrating-from-matlab) for documented Elastix parameter differences (bending energy, histogram bins, ASGD step estimation).

---

## See also

- [Home](index.md) — quick start and documentation map
- [Installation](installation.md)
- [Brain lightsheet usage](usage_lightsheet_brain.md)
- [Spinal cord usage](usage_spinal_cord.md)
