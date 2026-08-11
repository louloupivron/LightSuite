# LightSuite (Python)

**LightSuite** registers large microscopy volumes to standard brain atlases and exports atlas-space intensities. The Python version runs as a command-line tool (`lightsuite`) with optional Napari GUIs for manual registration refinement.

## What works today

| Workflow | Status |
|----------|--------|
| **Brain lightsheet** (3D whole-brain volumes) | Preprocess → check orientation → init registration → match points → register → export |
| **Multiresolution** (overview ↔ ROI) | Manifest → match points (optional) → check geometry → register ([usage guide](usage_multiresolution.md)) |
| **Spinal cord lightsheet** | Preprocess → straighten → register → export → per-sample `region-stats` ([usage guide](usage_spinal_cord.md)) |
| **Widefield coronal slices** | [MATLAB only](usage_slice.md) |

See **[Python vs MATLAB](python_vs_matlab.md)** for features that are **not implemented in Python** (cell detection, slice module, CZI reader, spinal cohort analysis, GPU-accelerated detection).

### Brain pipeline capabilities

- TIFF discovery (channel-per-file and plane-per-file layouts)
- Downsampling to registration resolution (default 20 µm)
- Interactive orientation checker (Napari 3D volume GUI)
- Coarse similarity alignment (Open3D ICP)
- Interactive control-point matching (Napari dual-pane GUI)
- Deformable B-spline registration (Elastix 5.1)
- Export of atlas-space registered volumes (TIFF) and parcellation intensity tables
- Allen and Perens brain atlas providers
- External annotation import (`points_csv`, `mask_tiff`)

## Quick start

```bash
# Install (see Installation guide)
uv sync --extra dev --extra gui

# Verify environment
uv run lightsuite doctor -c examples/brain_lightsheet.yaml

# Run the brain pipeline (one stage at a time)
uv run lightsuite brain preprocess           -c my_sample.yaml
uv run lightsuite brain check-orientation    -c my_sample.yaml
uv run lightsuite brain init-registration    -c my_sample.yaml
uv run lightsuite brain match-points         -c my_sample.yaml
uv run lightsuite brain register             -c my_sample.yaml
uv run lightsuite brain export               -c my_sample.yaml --save-volume --write-csv

# Spinal cord (see usage_spinal_cord.md)
uv run lightsuite spinal preprocess           -c my_spinal.yaml
uv run lightsuite spinal straighten           -c my_spinal.yaml
uv run lightsuite spinal init-registration    -c my_spinal.yaml
```

Copy [`examples/brain_lightsheet.yaml`](../examples/brain_lightsheet.yaml), edit paths and voxel size, then follow the [brain lightsheet guide](usage_lightsheet_brain.md).

## Hardware recommendations

| Task | Recommendation |
|------|----------------|
| Preprocess + registration | Workstation with fast SSD scratch space (≥50 GB free; 500 GB+ for very large brains) |
| Match-points GUI | Display with enough resolution for dual-pane Napari; `--extra gui` install |

## Supported input formats (brain)

| Format | Status |
|--------|--------|
| Multi-page TIFF, channel per file (`channelperfile`) | Supported |
| One TIFF per Z plane (`planeperfile`) | Supported |
| CZI | Not implemented — export to TIFF first |
| Imaris `.ims` | Spot conversion helper for spinal cord; full volume reader not implemented |

Registered outputs are written as **TIFF** (not OME-Zarr).

## Documentation map

1. **[Installation](installation.md)** — Python, `uv`, Elastix, atlas files
2. **[Python vs MATLAB](python_vs_matlab.md)** — feature parity and gaps
3. **[Brain lightsheet usage](usage_lightsheet_brain.md)** — YAML config, CLI stages, outputs, GUI
4. **[Multiresolution registration](usage_multiresolution.md)** — overview ↔ ROI alignment (mesoSPIM, SmartSPIM)
5. **[Spinal cord](usage_spinal_cord.md)** — Python spinal pipeline
6. **[Slice module](usage_slice.md)** — MATLAB workflow (not in Python)

## Getting help

- Run `uv run lightsuite doctor` to diagnose missing dependencies
- Open an [issue on GitHub](https://github.com/dimokaramanlis/LightSuite/issues) with OS, Python version, config (redacted paths), and the full error message
