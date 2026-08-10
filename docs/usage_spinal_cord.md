# Spinal cord lightsheet analysis

The **Python spinal cord pipeline** covers preprocess, straightening, registration, export, intensity parcellation, hierarchy rollups, point-based region stats, and cord-specific plots. Cohort NNMF normalization remains in `example_analysis_spinal_cord.m`.

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
uv run lightsuite spinal inspect-imports      -c my_spinal.yaml --space sample
uv run lightsuite spinal view                 -c my_spinal.yaml
uv run lightsuite spinal view                 -c my_spinal.yaml --space sample
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

### Intensity parcellation (registered channels)

After `export`, compute **median intensity**, **std**, and **volume_mm3** per
Fiederling region × segment from the exported `chan*_channel*.tiff` volumes
(same grid as `annotation_registered.tiff`):

```bash
uv run lightsuite spinal region-stats -c my_spinal.yaml
```

`region-stats` combines intensity parcellation and imported point counts into
one `region_stats.csv`. Disable either step:

```bash
uv run lightsuite spinal region-stats -c my_spinal.yaml --no-count-points
uv run lightsuite spinal region-stats -c my_spinal.yaml --no-parcellate-intensities
```

Configure in YAML:

```yaml
analysis:
  parcellate_intensities: true
  intensity_channels: [1, 2]   # omit to use all exported channels
  relative_intensity_to: none  # or background (per-segment id 0 reference)
  rollups: [division, structure]  # optional GM/WM and lamina/funiculus rollups
  count_points: true
  point_labels:
    - imaris_TAyellow
```

Outputs under `volume_registered/`:

- `region_stats.csv` — combined long-form table (intensities + cell counts)
- `chan{NN}_region_stats.csv` — per-channel intensity table
- `{label}_region_counts.csv` — per-import-label cell count table

Intensity metrics: `median_intensity`, `std`, `volume_mm3`, and optionally
`relative_median_intensity` when `relative_intensity_to: background`.

Optional hierarchy rollups aggregate finest-region stats to **division** (GM/WM)
or **structure** (combined laminas and funiculi) with volume-weighted means:

```yaml
analysis:
  rollups: [division, structure]
```

Rolled rows are appended to `region_stats.csv` with `rollup_level` =
`region`, `division`, or `structure`. Per-level sidecars:
`region_stats_division.csv`, `region_stats_structure.csv`.

Each row includes `segment` (e.g. `C5`, `L3`), `parcellation_index`, region name/acronym,
and `hemisphere` = `whole` (cord has no left/right split) unless
`analysis.split_hemispheres: true` is set (uses `Hemisphere_Annotation.tif` from the
atlas package; rows are labeled `left` / `right`, with optional `hemisphere_flip` and
`hemisphere_keep_whole`).

### Cord plots

After `region-stats`, generate matplotlib figures (requires the `gui` extra / matplotlib).

Write outputs under the sample **data** folder (`sample.save_path/plots/`), not inside the
LightSuite repository. Override with `analysis.plots_dir` in the YAML if needed.

```bash
PLOTS=/media/gbm/NVME2/ALICe-pipelines-data/spinal_cord/OP87F4/registered/plots
mkdir -p "$PLOTS"

uv run lightsuite analysis plot-cord-structure \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/structure.png" \
  --channel 1 --metric median_intensity

uv run lightsuite analysis plot-cord-division-profile \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/division_profile.png" \
  --channel 1 --metric median_intensity

uv run lightsuite analysis plot-cord-segment-bars \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/segment_counts.png" \
  --channel imaris_Coloc_pink_yellow --metric cell_count

uv run lightsuite analysis plot-cord-segment-grouped-bars \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/coloc_segment_counts.png"

uv run lightsuite analysis plot-cord-structure-panel \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/coloc_structure_panel.png"

uv run lightsuite analysis plot-cord-top-regions \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/top_regions.png" \
  --channel imaris_Coloc_pink_yellow --segment L5

uv run lightsuite analysis cord-coloc-overlap \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/coloc_overlap.png"

# Rexed laminae: % GM occupied by signal intensity vs. cell density
uv run lightsuite analysis plot-cord-laminae-pct-gm \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/laminae_pct_gm.png" \
  --intensity-channel 1 \
  --cell-channel imaris_Coloc_pink_yellow

# Optional: restrict to cervical segments
#   --segments C4,C5,C6,C7

# Rexed laminae × cord level (cervical / thoracic / lumbar means)
uv run lightsuite analysis plot-cord-laminae-level-bars \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/laminae_level_bars.png" \
  --channel 1 --metric median_intensity

# Dorsal funiculus subregions (dcs, cu, gr, psdc, df) × segment
uv run lightsuite analysis plot-cord-df-subregion-heatmap \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/df_subregion_heatmap.png" \
  --channel 1 --metric median_intensity

# Dorsal / ventral / central horn × segment
uv run lightsuite analysis plot-cord-horn-heatmap \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/horn_heatmap.png" \
  --channel imaris_Coloc_pink_yellow --metric cell_count --hemisphere left
```

| Command | Uses `rollup_level` | MATLAB equivalent |
|---------|---------------------|-------------------|
| `plot-cord-structure` | `structure` | Top-row `imagesc` heatmap |
| `plot-cord-structure-panel` | `structure` (multi-label) | Multi-channel `imagesc` row |
| `plot-cord-division-profile` | `division` | Bottom-row GM/WM line plot |
| `plot-cord-segment-bars` | `region` (summed per segment) | — |
| `plot-cord-segment-grouped-bars` | `region` (summed per segment, multiple labels) | — |
| `plot-cord-top-regions` | `region` | — |
| `plot-cord-top-regions-grouped` | `region` | Multi-label top region comparison |
| `plot-cord-laminae-pct-gm` | `structure` (laminae I–X) | % GM bar chart |
| `plot-cord-laminae-level-bars` | `structure` (laminae I–X) | Level-grouped bar chart |
| `plot-cord-laminae-grouped-bars` | `structure` (laminae I–X) | Multi-label lamina comparison |
| `plot-cord-df-subregion-heatmap` | `region` + `structure` (`df`) | Dorsal funiculus heatmap |
| `plot-cord-horn-heatmap` | `horn` (`DH`, `VH`, `C`) | Dorsal/ventral horn heatmap |
| `cord-coloc-overlap` | atlas-space spot coords | — |

`plot-cord-segment-grouped-bars` compares several imported spot labels on one chart
(e.g. Imaris coloc components). With `--config`, labels default to
`analysis.point_labels`; override with `--channels label_a,label_b`. Use
`--min-total 1` to hide segments with zero cells across all labels.

`cord-coloc-overlap` measures how many spots in one import label have a neighbor
within `--tolerance-voxels` (default 2) in another label, using
`*_atlas_coords.npz`. Writes a summary CSV and overlap bar chart.

`plot-cord-laminae-pct-gm` compares volume-weighted signal intensity (gray bars)
against cell-count share (red bars) across combined Rexed laminae I–X at
`structure` rollup. Use `--intensity-channel` for the imaging channel and
`--cell-channel` for an import label from `analysis.point_labels`. SEM error
bars appear automatically when multiple samples are present in the CSV.

`plot-cord-laminae-level-bars` averages a metric per lamina within each cord
level (default cervical / thoracic / lumbar). Override levels with
`--levels C,T,L,S` and filter segments with `--segments C4,C5,C6,C7`.

`plot-cord-laminae-grouped-bars` compares several import labels side-by-side at
each combined Rexed lamina (structure rollup). Labels default to
`analysis.point_labels` from `--config`; override with `--channels`. Use
`--segments L3,L4,L5,L6` to restrict rostrocaudal extent.

`plot-cord-df-subregion-heatmap` shows dorsal funiculus subregions
(`dcs`, `cu`, `gr`, `psdc`) at finest `region` rollup, plus the combined
`df` row from `structure` rollup. Disable the parent row with
`--no-include-parent-df`.

`plot-cord-horn-heatmap` shows dorsal horn (`DH`), ventral horn (`VH`), and
central (`C`) regions at `horn` rollup (requires `analysis.rollups` to include
`horn`). Horn rollup assigns each finest-level region to exactly one horn by
walking `parent_ID` ancestors (not by overlapping descendant sets). Use
`--hemisphere left|right` when `split_hemispheres: true`.

Enable left/right hemisegment stats in the YAML before running `region-stats`:

```yaml
analysis:
  split_hemispheres: true
  hemisphere_flip: false          # swap 0/255 assignment if needed
  hemisphere_keep_whole: false    # also emit whole-cord rows
```

Then plot one side or both:

```bash
PLOTS=/media/gbm/NVME2/ALICe-pipelines-data/spinal_cord/OP87F4/registered/plots

uv run lightsuite analysis plot-cord-structure \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/structure_right.png" \
  --channel 1 --metric median_intensity --hemisphere right

uv run lightsuite analysis plot-cord-structure-hemisphere-panel \
  -c examples/config/spinal_cord/OP87F4.yaml \
  -o "$PLOTS/structure_lr_panel.png" \
  --channel 1 --metric median_intensity
```

Re-run `lightsuite spinal region-stats` after upgrading to pick up `volume_mm3` on
point counts (enables structure/division `cell_density` rollups).

Use `--input region_stats.csv` instead of `--config` when plotting from a copied CSV.
Each command also writes a sidecar `.csv` next to the PNG.

### Inspect imports (Napari QC)

After `import-annotations`, open a Napari viewer with registered channels, atlas
annotation, and imported point layers overlaid:

```bash
uv run lightsuite spinal inspect-imports -c my_spinal.yaml
```

Sample space — straightened 20 µm grid with warped labels and `*_sample_coords.npz`
spot layers (requires `export --space sample`):

```bash
uv run lightsuite spinal inspect-imports -c my_spinal.yaml --space sample
```

Use `--headless` to validate inputs without opening the GUI.

When `cord_orientation.txt` sets `caudorostral` (`tofliprc: true`), sample-space inspect
flips the registration **Z** axis for display so rostrocaudal scrolling matches atlas-space
inspect. Left/right hemisegment placement uses the warped `hemisphere_in_sample_20um.tif`
(see [Registration spaces](registration_spaces.md)); that is separate from the longitudinal
orientation flip.

### View registration (Napari)

Atlas space (default) — sample warped onto the Fiederling export grid:

```bash
uv run lightsuite spinal view -c my_spinal.yaml
```

Sample space — straightened 20 µm registration grid with warped atlas labels
(requires `lightsuite spinal export --space sample` first):

```bash
uv run lightsuite spinal export -c my_spinal.yaml --space sample
uv run lightsuite spinal view -c my_spinal.yaml --space sample
```

Layers: straightened channel(s), warped atlas template, warped annotation labels.
Imported spot layers (`*_sample_coords.npz`) are overlaid when present.

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

During **preprocess**, all channels are cached as `cache/chan_N_sample_register_20um.tif` (paths recorded in `regopts.json` → `regvolpaths`). **Export** (atlas space, default) upsamples registered output to the native Fiederling template grid (10×10×20 µm) under `volume_registered/`. **Sample-space export** (`export.spaces: [sample]` or `--space both`) writes warped labels and straightened channels under `volume_registered/sample_space/` — see [Registration spaces](registration_spaces.md).

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
| `region_stats_division.csv` | region-stats (rollup) |
| `region_stats_structure.csv` | region-stats (rollup) |
| `chan{NN}_region_stats.csv` | region-stats (intensity) |
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
