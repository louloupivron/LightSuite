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

**Requirements:** Python 3.11+, [Elastix 5.1.0](https://github.com/SuperElastix/elastix/releases/tag/5.1.0) on `PATH`, atlas NIfTI files. See [installation guide](https://lightsuite.readthedocs.io/en/latest/installation/) and [Hardware requirements](#hardware-requirements) below.

---

## Hardware requirements

The Python pipeline is **CPU + RAM + fast local disk** bound. **GPU is not used** (`compute.use_gpu` is reserved and ignored). Built-in cell detection is not implemented in Python.

Run `uv run lightsuite doctor -c <config.yaml>` before long jobs — it checks Elastix, scratch free space (≥50 GB), atlas paths, and common config mistakes.

### Workstation profiles

| Profile | CPU | RAM | Scratch disk |
|---------|-----|-----|--------------|
| Spinal cord, 1–2 channels | 4+ cores | 16 GB | 50 GB SSD |
| Whole brain, 20 µm, 1–2 channels | 8+ cores | 32–64 GB | 100–200 GB SSD |
| Large brain / wide field / dual-channel MI | 8–16 cores | 64–128 GB | 200–500 GB SSD |
| Multiresolution (large overlap crop) | 8+ cores | 32–64 GB+ | 100 GB+ SSD |

Scratch should be a **fast local SSD/NVMe** (`sample.scratch`), not a network mount. Elastix/transformix temps also use scratch (spinal: `sample.scratch/<sample.name>/`).

### Sizing formulas (brain / spinal preprocess)

| Quantity | Formula |
|----------|---------|
| XY downsample factor | `vx / registration.resolution_um` (and `vy`) |
| Z downsample factor | `vz / registration.resolution_um` |
| Preprocess XY scratch (per channel) | `out_y × out_x × nz × 2` bytes (`uint16`; `nz` = native Z during load) |
| Registration TIFF (per channel) | `out_y × out_x × out_z × 2` bytes |
| Scratch disk (plane-per-file, rough) | `2 × ny × nx × nz × (vx/registres)²` bytes + Elastix temps |

In RAM if scratch fits below `compute.max_in_memory_scratch_gb` (default **24 GB**); otherwise a memmap on `sample.scratch`. `planeperfile` stacks always use **1 worker** (parallel reads thrash disk); `channelperfile` uses `compute.workers` (default 4).

### Brain pipeline (`lightsuite brain`)

| Step | Command | CPU | RAM | Disk | If requirements are not met |
|------|---------|-----|-----|------|----------------------------|
| Environment | `doctor` | 1 core | &lt;1 GB | ≥50 GB free on scratch | Point `sample.scratch` to fast SSD; install Elastix 5.1.0 on `PATH`. |
| Preprocess | `preprocess` | High (`workers` for channel-per-file; **1 worker** for plane-per-file) | XY scratch ≤24 GB default in RAM, else memmap | Heavy read/write of TIFF stacks | Fast NVMe scratch; correct `sample.voxel_um`; lower `max_in_memory_scratch_gb` to force memmap if OOM; raise `workers` only for channel-per-file. |
| Check orientation | `check-orientation` | Low | Napari previews capped (~256 MB/volume) | Minimal | Set `registration.orientation` in YAML and skip GUI; `uv sync --extra gui`. |
| Init registration | `init-registration` | Moderate (BCPD / Open3D ICP) | Full 20 µm volume + atlas (often 8–32 GB) | Small | Install BCPD; raise `registration.resolution_um` and re-preprocess; `sample_content_crop: auto`. |
| Align slices | `align-slices` | Low | Napari loads registration volume | Minimal | `--headless` or skip; register works without correspondence. |
| Match points | `match-points` | Low | Napari | Minimal | Optional — register can use init-registration auto points only. |
| Register | `register` | **Very high** (Elastix B-spline, CPU-only) | **Peak step** — sample + atlas + warped copies (often 32–128 GB for large/wide fields) | `elastix_temp/` under `save_path` | `sample_content_crop: auto` or manual box; larger `bspline_spatial_scale_mm`; `--single-step`; drop `channel_secondary`; raise `registration.resolution_um` and re-preprocess. |
| Export | `export` | High (transformix per channel) | Full volume per channel during warp | Large if `save_registered_volume: true` | `save_registered_volume: false`; `save_sample_space_volume: false`; CSV-only export; subset `analysis.intensity_channels`. |
| Import annotations | `import-annotations` | Low–moderate | Small (points); mask + transformix for TIFF masks | Small outputs | Prefer `points_csv`; requires `transformix` on `PATH`. |
| View registration | `view-registration` | Low | Napari + optional TIFF overlays | Reads existing exports | Run `export --space sample` first; needs `--extra gui`. |

See [brain lightsheet usage](https://lightsuite.readthedocs.io/en/latest/usage_lightsheet_brain/) for B-spline tuning (`bspline_bending_weight`, wide-field crops).

### Spinal cord pipeline (`lightsuite spinal`)

| Step | Command | CPU | RAM | Disk | If requirements are not met |
|------|---------|-----|-----|------|----------------------------|
| Preprocess | `preprocess` | High (streaming planes) | Caches all channels in `cache/` | Heavy plane TIFF reads | Correct `sample.voxel_um` (doctor warns if no downsampling on load); fast scratch. |
| Check orientation | `check-orientation` | Low | Small | Minimal | Set orientation in YAML / `cord_orientation.txt`. |
| Straighten | `straighten` | Moderate (per-slice 2D warps) | Full registration-grid volume | Straightened cache TIFF | Cord volumes are smaller than whole brain; trim Z in preprocess if needed. |
| Init / match / register | `init-registration`, `match-points`, `register` | High at register (Elastix) | Typically 4–16 GB at register | Elastix under `sample.scratch/<name>/` | Same Elastix mitigations as brain; match-points optional. |
| Export / region-stats | `export`, `region-stats` | High at export (transformix) | Upsampled Fiederling grid | `volume_registered/` TIFFs | Export subset of channels; defer full TIFF export. |

### Multiresolution pipeline (`lightsuite multires`)

| Step | Command | CPU | RAM | Disk | If requirements are not met |
|------|---------|-----|-----|------|----------------------------|
| Check geometry | `check-geometry` | Moderate | Warns if overlap peak &gt; available RAM | Minimal | Shrink overlap (`overlap_margin_um` negative); increase `registration_bin`. |
| Register | `register` | High (Elastix on overlap crop) | Peak ≈ `2.5 × crop_gb + max_slab_bytes` (default slab 500 MB); stacks streamed, not loaded whole | Elastix + outputs under `save_path` | Raise `registration_bin` (e.g. 2); lower `max_slab_bytes`; `write_full_overview_canvas: false`; fast scratch SSD. |
| Export / preview | `export-preview`, etc. | Moderate | Streaming slabs | Registered canvases | Lower `max_slab_bytes`. |

Full stacks are **not** held in RAM; overlap crops are materialized in slabs. See [multiresolution usage](https://lightsuite.readthedocs.io/en/latest/usage_multiresolution/).

### Cross-cutting remediation

| Problem | What to try |
|---------|-------------|
| Scratch full or slow | Move `sample.scratch` to NVMe; free ≥50 GB; delete old `elastix_temp/` and scratch memmaps (`chan_*_xy_*.dat`). |
| RAM OOM at register/export | Content crop; higher `registration.resolution_um` + re-preprocess; disable registered volume export; single channel; `brain register --single-step`. |
| RAM OOM at preprocess | Lower `compute.max_in_memory_scratch_gb` (disk memmap); reduce `compute.workers` on channel-per-file. |
| Process killed, no Python traceback | Likely Linux OOM killer during Elastix/transformix — same RAM mitigations. |
| GUI won't start | `uv sync --extra gui`; more RAM; set orientation/control points in YAML and skip Napari stages. |

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
