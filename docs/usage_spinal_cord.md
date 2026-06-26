# Spinal cord lightsheet analysis

The **Python spinal cord pipeline** (MVP) covers preprocess, straightening, init registration, match-points, register, and export. Post-registration cohort analysis remains MATLAB-only (`example_analysis_spinal_cord.m`).

For the brain pipeline, see [Brain lightsheet analysis](usage_lightsheet_brain.md).

---

## Python workflow (MVP)

Spinal cord volumes register against the [Fiederling et al. (2021)](https://www.sciencedirect.com/science/article/pii/S2667237521001260) atlas.

### Getting started

1. Install LightSuite for Python ([Installation](installation.md))
2. Download the Fiederling atlas from [Mendeley Data](https://data.mendeley.com/datasets/4rrggzv5d5/1)
3. Copy [`examples/spinal_cord.yaml`](../examples/spinal_cord.yaml) and edit paths
4. Run:

```bash
uv run lightsuite doctor -c examples/spinal_cord.yaml
uv run lightsuite spinal validate-config -c my_spinal.yaml
uv run lightsuite spinal preprocess           -c my_spinal.yaml
uv run lightsuite spinal straighten           -c my_spinal.yaml
uv run lightsuite spinal init-registration    -c my_spinal.yaml
uv run lightsuite spinal match-points         -c my_spinal.yaml
uv run lightsuite spinal register             -c my_spinal.yaml
uv run lightsuite spinal export               -c my_spinal.yaml
uv run lightsuite spinal view                 -c my_spinal.yaml
```

### Key differences from brain

- **Straightening GUI** (`spinal straighten`) — trace central canal and anterior/posterior axis per slice before registration. Select **center**, **anterior**, or **posterior** in Napari’s layer list, then click on the slice to place a point. A **fit preview** (yellow axis, cyan center cross, green/red predicted ant/pos) appears once at least two slices are annotated; use **Toggle fit (p)** or **P** to show/hide it. **Clear last** undoes the most recent point; **Clear all** removes every annotation. Navigate with the **slice slider/spinbox**, **←/→**, **Page Up/Down**, or **scroll wheel**; save with **Save (s)** or **S** (works even when the dock has focus).
- Cord-specific Fiederling atlas (TIFF + CSV), not Allen/Perens NIfTI
- Export uses slice-style TIFF layout under `volume_registered/` (not brain OME-Zarr)

### Sample TIFF layout and resolution

`sample.source.tiff_type` controls how raw data are read:

| Value | Layout |
|-------|--------|
| `auto` | Detect plane-per-file (Terastitcher) vs channel-per-file |
| `planeperfile` | One 2D TIFF per Z slice (large datasets) |
| `channelperfile` | One multi-page TIFF per channel |
| `multichannel_single` | Single Bioformats-style stack |

**Important:** `sample.voxel_um` is the **native microscope voxel size**, not the registration grid. Registration always uses `registration.resolution_um` (default 20 µm). For Terastitcher exports (~1.8×1.8×4 µm), set e.g. `voxel_um: [1.8, 1.8, 4]` so plane-per-file stacks are downsampled while loading instead of loaded at full resolution.

During **preprocess**, all channels are cached as `cache/chan_N_sample_register_20um.tif` (paths recorded in `regopts.json` → `regvolpaths`). **Export** loads these cached volumes, upsamples registered output to the native Fiederling template grid (10×10×20 µm), and writes `annotation_registered.tiff` and `template_registered.tiff` alongside the channel TIFFs under `volume_registered/`. **View** opens Napari with the template, all exported channels, and warped annotation labels in a shared axis layout — do not manually load the raw atlas `Template.tif` or pre-registration cached TIFFs into the same viewer (different grid and orientation).

Optional `sample.source.skip_corrupt_slices: true` drops unreadable plane TIFFs instead of failing (use only when a few slices are damaged).

`lightsuite doctor -c …` warns when a large plane-per-file folder has `voxel_um` equal to `registration.resolution_um`.

### Checkpoints

Under `<sample.save_path>/`:

| File / folder | Stage |
|------|-------|
| `regopts.json` | preprocess (+ updated by init-registration) |
| `spinal_alignment_opt.json` | straighten |
| `corresponding_points.json` | match-points |
| `transform_params.json` | register |
| `cache/` | preprocess + init-registration intermediates (registration-grid TIFFs, straightened volume, resampled atlas) |
| `transforms/` | elastix affine + inverted B-spline parameter files |
| `qc/` | registration overlay PNGs |
| `volume_registered/` | export (`chan*_channel*.tiff`, `template_registered.tiff`, `annotation_registered.tiff`) |

Ephemeral elastix / transformix workspaces are written under `<sample.scratch>/<sample.name>/` (not under `save_path`).

Legacy flat layouts (artifacts directly under `save_path/`, transformix folders next to checkpoints) are still read when present; re-running a stage writes the organized layout above.

### Parity validation

```bash
uv run lightsuite spinal validate-parity --fixture-root tests/fixtures/spinal_cord
```

---

## MATLAB workflow (legacy)

The original MATLAB entry script `demos/ls_analyze_spinal_cord.m` remains available for comparison and for analysis scripts not yet ported.
