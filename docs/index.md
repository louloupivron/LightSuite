# LightSuite (Python)

**LightSuite** registers large microscopy volumes to standard brain atlases and exports atlas-space intensities and (eventually) cell coordinates. The Python version runs as a command-line tool (`lightsuite`) with optional Napari GUIs for manual registration refinement.

This documentation covers the **Python pipeline** on branch `feature/python-migration`. The original MATLAB workflows remain in the repository for spinal cord and slice data until those modules are ported.

## What works today

| Workflow | Status |
|----------|--------|
| **Brain lightsheet** (3D whole-brain volumes) | Preprocess → init registration → match points → register → export |
| **Spinal cord lightsheet** | MATLAB only ([usage guide](usage_spinal_cord.md)) |
| **Widefield coronal slices** | MATLAB only ([usage guide](usage_slice.md)) |

### Brain pipeline capabilities

- TIFF discovery (channel-per-file and plane-per-file layouts)
- Downsampling to registration resolution (default 20 µm)
- Coarse similarity alignment (Open3D ICP)
- Interactive control-point matching (Napari dual-pane GUI)
- Deformable B-spline registration (Elastix 5.1)
- Export of atlas-space registered volumes and parcellation intensity tables
- Allen and Perens brain atlas providers

### Not yet ported from MATLAB

- 3D cell detection and atlas mapping of cell coordinates
- Spinal cord and slice pipelines
- OME-Zarr / Imaris readers (planned plugin layer)
- GPU-accelerated detection

## Quick start

```bash
# Install (see Installation guide)
uv sync --extra dev --extra gui

# Verify environment
uv run lightsuite doctor -c examples/brain_lightsheet.yaml

# Run the brain pipeline (one stage at a time)
uv run lightsuite brain preprocess           -c my_sample.yaml
uv run lightsuite brain init-registration    -c my_sample.yaml
uv run lightsuite brain match-points         -c my_sample.yaml
uv run lightsuite brain register             -c my_sample.yaml
uv run lightsuite brain export               -c my_sample.yaml --save-volume --write-csv
```

Copy [`examples/brain_lightsheet.yaml`](../examples/brain_lightsheet.yaml), edit paths and voxel size, then follow the [brain lightsheet guide](usage_lightsheet_brain.md).

## Hardware recommendations

| Task | Recommendation |
|------|----------------|
| Preprocess + registration | Workstation with fast SSD scratch space (≥50 GB free; 500 GB+ for very large samples) |
| Match-points GUI | Display with enough resolution for dual-pane Napari; `--extra gui` install |
| Cell detection (future) | GPU recommended for large volumes |

## Supported input formats (brain)

| Format | Status |
|--------|--------|
| Multi-page TIFF, channel per file (`channelperfile`) | Supported |
| One TIFF per Z plane (`planeperfile`) | Supported |
| CZI, OME-Zarr, Imaris | Planned |

## Documentation map

1. **[Installation](installation.md)** — Python, `uv`, Elastix, atlas files
2. **[Brain lightsheet usage](usage_lightsheet_brain.md)** — YAML config, CLI stages, outputs, GUI
3. **[Spinal cord](usage_spinal_cord.md)** — MATLAB workflow (not yet in Python)
4. **[Slice module](usage_slice.md)** — MATLAB workflow (not yet in Python)

## Getting help

- Run `uv run lightsuite doctor` to diagnose missing dependencies
- Open an [issue on GitHub](https://github.com/dimokaramanlis/LightSuite/issues) with OS, Python version, config (redacted paths), and the full error message
