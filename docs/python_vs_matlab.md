# Python vs legacy MATLAB LightSuite

This repository ships the **Python pipeline** (`lightsuite` CLI) only. The original MATLAB toolbox is no longer included here.

Use this page when deciding whether LightSuite Python covers your workflow, or when you need features that were never ported.

---

## Not implemented in Python

| Feature | Notes |
|---------|--------|
| **Built-in 3D cell detection** | Band-pass filtering + local maxima in the legacy MATLAB toolbox. In Python, set `detection.enabled: false` and import external spot lists with [`brain import-annotations`](annotation_import.md). |
| **Widefield coronal slice module** | Per-slice 2D registration for widefield / AxioScan data. |
| **CZI volume reader** | Python accepts **TIFF stacks** (`channelperfile`, `planeperfile`). CZI must be exported to TIFF first. |
| **Spinal cord cohort analysis** | Cross-subject NNMF normalization and cohort plotting. Python provides per-sample `spinal region-stats` only. |
| **GPU-accelerated cell detection** | `compute.use_gpu` in YAML is ignored by the Python pipeline. |

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

## MATLAB interoperability (optional)

If you still run legacy MATLAB registration elsewhere, conversion helpers live under
[`devtools/matlab_interop/`](../devtools/matlab_interop/README.md) (control points and
`regopts` export). Parity scratch scripts are in [`devtools/parity/`](../devtools/parity/README.md).

---

## See also

- [Home](index.md) — quick start and documentation map
- [Installation](installation.md)
- [Brain lightsheet usage](usage_lightsheet_brain.md)
- [Spinal cord usage](usage_spinal_cord.md)
