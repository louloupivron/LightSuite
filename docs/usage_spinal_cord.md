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
uv run lightsuite spinal check-orientation    -c my_spinal.yaml
uv run lightsuite spinal preprocess           -c my_spinal.yaml
uv run lightsuite spinal straighten           -c my_spinal.yaml
uv run lightsuite spinal align-longitudinal   -c my_spinal.yaml
uv run lightsuite spinal init-registration    -c my_spinal.yaml
uv run lightsuite spinal match-points         -c my_spinal.yaml
uv run lightsuite spinal register             -c my_spinal.yaml
uv run lightsuite spinal export               -c my_spinal.yaml
uv run lightsuite spinal import-annotations   -c my_spinal.yaml
uv run lightsuite spinal region-stats         -c my_spinal.yaml
uv run lightsuite spinal inspect-imports      -c my_spinal.yaml
uv run lightsuite spinal view                 -c my_spinal.yaml
```

Convert Imaris spot exports before import:

```bash
uv run lightsuite spinal convert-imaris-spots \
  -s /path/to/Spot_OnePageMultiComponent_Detailed.csv \
  -o /path/to/converted \
  --voxel-um 1,1,1.8 \
  --shape-yxz 10234,2006,837
```

`--voxel-um` is the size of one LightSuite native voxel **in the same units as the
Imaris Position columns** (Image Properties → Geometry → Voxel Size), **not**
necessarily `sample.voxel_um` from the YAML.

Check the `.ims` calibration before converting:

| Imaris voxel size | LightSuite `sample_reference.json` | Typical `--voxel-um` |
|-------------------|------------------------------------|----------------------|
| Same as microscope (e.g. 1.8³ µm) | Same grid and voxel size | `1.8,1.8,1.8` |
| 1.00³ µm (Position ≈ voxel indices) | Same XYZ shape as Imaris | `1,1,1` |
| 1.00³ µm, but **Z plane count differs** from LightSuite (common: Imaris keeps full acquisition Z, LightSuite has fewer planes at true 1.8 µm) | XY shape matches; `nz_imaris / nz_lightsuite ≈ sample.voxel_um[2]` | **`1,1,1.8`** (XY as Imaris indices, Z scaled to the LightSuite grid) |

Using the microscope voxel size (`1.8,1.8,1.8`) when the `.ims` is calibrated at 1 µm shrinks coordinates toward the origin and places most spots **outside the cord** after import. Pass `--shape-yxz` from `sample_reference.json` to get a warning when the units look wrong; confirm placement with `inspect-imports`.

### Cell counts from imported spots

After `import-annotations`, bin atlas-space points into Fiederling regions and
rostrocaudal segments (`Segments.csv`):

```bash
uv run lightsuite spinal region-stats -c my_spinal.yaml
```

Configure which import labels to count (optional):

```yaml
analysis:
  count_points: true
  point_labels:
    - imaris_TAyellow
    - imaris_MG_cyan
```

Outputs under `volume_registered/`:

- `region_stats.csv` — combined long-form table (`cell_count`, `cell_density` per region × segment)
- `{label}_region_counts.csv` — per-import-label table

Each row includes `segment` (e.g. `C5`, `L3`), `parcellation_index`, region name/acronym,
and `hemisphere` = `whole` (cord has no left/right split).

### Inspect imports (Napari QC)

After `import-annotations`, open a Napari viewer with registered channels, atlas
annotation, and imported point layers overlaid:

```bash
uv run lightsuite spinal inspect-imports -c my_spinal.yaml
```

Use `--headless` to validate inputs without opening the GUI.

### Key differences from brain

- **Orientation (`spinal check-orientation`)** — a small Napari window with side-by-side **longitudinal max projections**: atlas template (left, rostral at top) and sample (right). Click **Rostrocaudal** or **Caudorostral** to flip the sample until rostral anatomy aligns with the atlas, then **Save**. This writes `cord_orientation.txt` in `save_path`, which `preprocess` reads to decide whether to flip the atlas (replaces the old automatic detection). You can also edit `cord_orientation.txt` by hand or set `registration.longitudinal_direction: rostrocaudal|caudorostral` in the YAML to skip the GUI. **Preprocess fails if neither is set.** On the first run, check-orientation also caches downsampled per-channel TIFFs under `cache/` (`register_cache.json` + `chan_*_sample_register_*um.tif`); **preprocess reuses these** so you do not load the raw Terastitcher stack twice.
- **Straightening GUI** (`spinal straighten`) — trace central canal and anterior/posterior axis per slice before registration. Select **center**, **anterior**, or **posterior** in Napari’s layer list, then click on the slice to place a point. A **fit preview** (yellow axis, cyan center cross, green/red predicted ant/pos) appears once at least two slices are annotated; use **Toggle fit (p)** or **P** to show/hide it. **Clear last** undoes the most recent point; **Clear all** removes every annotation. Navigate with the **slice slider/spinbox**, **←/→**, **Page Up/Down**, or **scroll wheel**; save with **Save (s)** or **S** (works even when the dock has focus).
- **Longitudinal align GUI** (`spinal align-longitudinal`) — match straightened sample z-slices to atlas z-planes when the imaged cord is shorter than the template. Sample transverse slice on the **left** (fixed), atlas on the **right** (scroll with **PgUp/PgDn**, spinbox, or wheel over the atlas panel). Confirm ~20 anchors with **Enter** or **Confirm anchor**, then **Save && Close** before init-registration.
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

**Multi-channel Terastitcher (like brain):** when each channel is its own folder of plane TIFFs, list them under `sample.source.channels` (order = channel index). `path` is optional and defaults to the first folder. Requires `tiff_type: planeperfile` or `auto`. When channels have different slice counts, LightSuite intersects planes by the numeric index embedded in each TIFF filename (Terastitcher trailing `_NNNNN.tif`, SmartSPIM `_NNNNN_ChN.tif`, etc.) and loads only planes present in every channel.

```yaml
sample:
  source:
    tiff_type: planeperfile
    channels:
      - /data/OP93M1/Ex_488_Em_488F_Ch1_stitched
      - /data/OP93M1/Ex_561_Em_561F_Ch2_stitched
  save_path: /data/OP93M1/registered

registration:
  channel_primary: 2   # e.g. structural / autofluorescence channel
```

All channel folders must share the same `(ny, nx, nz)` and matching slice ordering. Preprocess caches every channel; registration uses `channel_primary`; export warps all channels with the same transform.

**Important:** `sample.voxel_um` is the **native microscope voxel size**, not the registration grid. Registration always uses `registration.resolution_um` (default 20 µm). For Terastitcher exports (~1.8×1.8×4 µm), set e.g. `voxel_um: [1.8, 1.8, 4]` so plane-per-file stacks are downsampled while loading instead of loaded at full resolution.

During **preprocess**, all channels are cached as `cache/chan_N_sample_register_20um.tif` (paths recorded in `regopts.json` → `regvolpaths`). **Export** loads these cached volumes, upsamples registered output to the native Fiederling template grid (10×10×20 µm), and writes `annotation_registered.tiff` and `template_registered.tiff` alongside the channel TIFFs under `volume_registered/`. **View** opens Napari with the template, all exported channels, and warped annotation labels in a shared axis layout — do not manually load the raw atlas `Template.tif` or pre-registration cached TIFFs into the same viewer (different grid and orientation).

Optional `sample.source.skip_corrupt_slices: true` drops unreadable plane TIFFs instead of failing (use only when a few slices are damaged).

`lightsuite doctor -c …` warns when a large plane-per-file folder has `voxel_um` equal to `registration.resolution_um`.

### Checkpoints

Under `<sample.save_path>/`:

| File / folder | Stage |
|------|-------|
| `cord_orientation.txt` | check-orientation (rostrocaudal/caudorostral) |
| `regopts.json` | preprocess (+ updated by init-registration) |
| `spinal_alignment_opt.json` | straighten |
| `longitudinal_correspondence.json` | align-longitudinal |
| `corresponding_points.json` | match-points |
| `transform_params.json` | register |
| `sample_reference.json` | preprocess |
| `{label}_atlas_coords.npz` | import-annotations |
| `region_stats.csv` | region-stats |
| `{label}_region_counts.csv` | region-stats |
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
