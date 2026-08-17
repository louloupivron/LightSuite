# LightSuite

[![Documentation Status](https://readthedocs.org/projects/lightsuite/badge/?version=latest)](https://lightsuite.readthedocs.io/en/latest/)
[![Python](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**LightSuite** registers large microscopy volumes to standard anatomical atlases and exports atlas-space intensities, parcellation tables, and cell coordinates. The **Python pipeline** (`lightsuite` CLI) targets 100 GB+ lightsheet datasets and runs **without a MATLAB license**.

Full documentation: [lightsuite.readthedocs.io](https://lightsuite.readthedocs.io/en/latest/)

---

## What the Python tool does

LightSuite Python is a command-line pipeline with optional **Napari** GUIs for manual refinement. Configuration is YAML-based; stages are resumable via JSON checkpoints.

| Workflow | CLI | Summary |
|----------|-----|---------|
| **Brain lightsheet** | `lightsuite brain` | Preprocess → orientation check → coarse alignment → control-point matching → B-spline registration (Elastix) → export registered volumes and regional intensity tables. Supports Allen CCF v3 and Perens LSFM atlases. |
| **Spinal cord** | `lightsuite spinal` | Preprocess → straighten (central canal / A–P axis) → register to the Fiederling spinal atlas → export → per-sample region statistics and cord-specific plots. |
| **Multiresolution** | `lightsuite multires` | Align a high-resolution ROI stack to a low-magnification overview (mesoSPIM, SmartSPIM, etc.) using vendor geometry and optional landmarks — no atlas required. |
| **External annotations** | `lightsuite brain import-annotations` | Import spot lists (`points.csv`) or masks (TIFF) after registration and map them to atlas regions. |

**Inputs:** multi-page TIFF stacks (channel-per-file or plane-per-file layouts).  
**Outputs:** atlas-space registered volumes (multi-page TIFF), parcellation intensity CSVs, optional annotation transforms.

Interactive stages (orientation, match-points, registration review) are available through `lightsuite gui -c config.yaml` or per-stage Napari tools.

---

## Quick start

```bash
git clone https://github.com/dimokaramanlis/LightSuite.git
cd LightSuite
uv sync --extra dev --extra gui --extra registration

cp examples/brain_lightsheet.yaml my_sample.yaml
# edit paths and voxel size

uv run lightsuite doctor -c my_sample.yaml
uv run lightsuite brain run -c my_sample.yaml --through match-points
```

Example configs: [`examples/README.md`](examples/README.md) (brain, spinal cord, multiresolution).

**Requirements:** Python 3.11+, [Elastix 5.1.0](https://github.com/SuperElastix/elastix/releases/tag/5.1.0) on `PATH`, atlas NIfTI files. See [installation guide](https://lightsuite.readthedocs.io/en/latest/installation/).

---

## Python vs original MATLAB LightSuite

The [original LightSuite](https://lightsuite.readthedocs.io/en/latest/) is a MATLAB toolbox covering whole-brain lightsheet, spinal cord, widefield coronal slices, built-in cell detection, and probe/implant tracing. The Python release is a **focused port** of the registration core — not a full rewrite.

| | **Python (`lightsuite`)** | **MATLAB (legacy)** |
|---|---------------------------|---------------------|
| **Runtime** | Python 3.11+, CLI + Napari GUIs | MATLAB R2022b+, MATLAB GUIs |
| **Configuration** | YAML files | Script-based `opts` structs |
| **Brain registration** | Yes | Yes |
| **Spinal cord registration** | Yes (per-sample stats) | Yes (+ cohort NNMF analysis) |
| **Multiresolution overview ↔ ROI** | Yes (`lightsuite multires`) | Limited / vendor-specific paths |
| **Built-in 3D cell detection** | No — import external spots | Yes (SNR band-pass + local maxima) |
| **Widefield coronal slice module** | No | Yes (2D per-slice registration, CZI) |
| **Probe & implant tracing** | No | Yes (`probe_ccf` export) |
| **GPU-accelerated detection** | No | Yes |
| **CZI reader** | No — export to TIFF first | Yes |

For features still MATLAB-only, use the legacy demos in `demos/` (`ls_analyze_lightsheet_volume.m`, `ls_analyze_spinal_cord.m`, `ls_analyze_slice_volume.m`). See the detailed parity table in [Python vs MATLAB](https://lightsuite.readthedocs.io/en/latest/python_vs_matlab/).

---

## Documentation

- [Installation](https://lightsuite.readthedocs.io/en/latest/installation/)
- [Brain lightsheet usage](https://lightsuite.readthedocs.io/en/latest/usage_lightsheet_brain/)
- [Spinal cord usage](https://lightsuite.readthedocs.io/en/latest/usage_spinal_cord/)
- [Multiresolution registration](https://lightsuite.readthedocs.io/en/latest/usage_multiresolution/)
- [Annotation import](https://lightsuite.readthedocs.io/en/latest/annotation_import/)

---

## Support & license

- **Issues:** [GitHub Issues](https://github.com/dimokaramanlis/LightSuite/issues)
- **License:** [GPL-3.0](LICENSE)

If you use LightSuite in your research, please cite this repository. A preprint is in preparation.
