# Brain lightsheet analysis (Python)

This guide walks through the Python brain pipeline: from stitched TIFF stacks to atlas-space registered volumes and regional intensity tables.

The workflow mirrors the MATLAB script `demos/ls_analyze_lightsheet_volume.m`, but uses a **YAML config file** and **`lightsuite brain`** CLI commands instead of editing an `opts` struct in MATLAB.

---

## Before you start

You will need:

- Stitched TIFF data (Terastitcher, BigStitcher, or similar)
- **Voxel size** in microns (`[x, y, z]`)
- Fast scratch space for intermediate files
- Elastix 5.1.0 and Allen (or Perens) atlas NIfTIs installed ([Installation](installation.md))
- Napari GUI extras for the match-points step: `uv sync --extra gui`

---

## Pipeline overview

| Step | CLI command | Type | MATLAB equivalent |
|:----:|-------------|------|-------------------|
| 0 | `lightsuite doctor` | Check | `check_lightsuite_installation.m` |
| 1 | `lightsuite brain preprocess` | Automated | `preprocessLightSheetVolume.m` |
| 2 | `lightsuite brain check-orientation` | **Manual (GUI)** | `getBrainOrientation.m` |
| 3 | `lightsuite brain init-registration` | Automated | `initializeRegistration.m` |
| 4 | `lightsuite brain align-slices` | **Manual (GUI)** | *(new — AP slice correspondence)* |
| 5 | `lightsuite brain refine-auto-points` | Automated | *(new — AP-filtered auto landmarks)* |
| 6 | `lightsuite brain match-points` | **Manual (GUI)** | `matchControlPoints_unified.m` |
| 7 | `lightsuite brain register` | Automated | `multiobjRegistration.m` |
| 8 | `lightsuite brain export` | Automated | `generateRegisteredBrainVolumes.m` |
| 9 | `lightsuite brain import-annotations` | Automated | `transformPointsToAtlas.m` |
| 10 | `lightsuite brain inspect-imports` | **Manual (GUI)** | *(Napari — atlas- or sample-space import QC)* |
| 11 | `lightsuite analysis region-stats` | Automated | *(new — tidy region table + cell counts)* |
| 12 | `lightsuite analysis group-stats` | Automated | *(new — cross-subject group summaries)* |
| 13 | `lightsuite analysis view-divisions` | **Manual (GUI)** | *(Napari — division-masked channel QC)* |
| 14 | `lightsuite analysis registration-qc` | Automated / **GUI** | *(unassigned-voxel registration score)* |

Built-in cell detection is **not yet ported**; set `detection.enabled: false` and use `import-annotations` with native `points.csv` / `mask.tif` exports (see [Annotation import](annotation_import.md)).

---

## Configuration

Create a YAML file for each sample. Start from the example:

```bash
cp examples/brain_lightsheet.yaml my_mouse.yaml
```

### Minimal example

```yaml
sample:
  name: M001
  source:
    format: tiff_stack
    path: /data/M001/stitched
    tiff_type: channelperfile    # or planeperfile
  scratch: /fastssd/lightsuite_scratch/M001
  save_path: /data/M001/lightsuite_results
  voxel_um: [5.26, 5.26, 5.0]    # required for preprocess

atlas:
  provider: allen                 # or perens
  resolution_um: 10
  atlas_dir: /atlases/allen_ccf_10um

registration:
  resolution_um: 20               # registration working resolution
  channel_primary: 1              # autofluorescence / structural channel
  channel_secondary: 2            # optional second channel for dual MI
  bspline_spatial_scale_mm: 0.64
  bspline_bending_weight: 2.0     # deformation smoothness; raise if contours look wavy
  control_point_weight: 0.2
  augment_points: false
  orientation: [1, 2, 3]          # or omit; uses brain_orientation.txt if present

detection:
  enabled: false                  # not yet implemented in Python

export:
  save_registered_volume: false
  write_cells_csv: true

import:
  write_csv: true
  annotations:
    - format: points_csv
      path: /data/annotations/cfos_cells.csv
      label: cfos_cells
    - format: mask_tiff
      path: /data/annotations/region_mask.tif
      label: hippocampus
```

See [Annotation import](annotation_import.md) for the native sample-space convention and tool conversion notes.

Validate before running:

```bash
uv run lightsuite brain validate-config -c my_mouse.yaml
uv run lightsuite doctor -c my_mouse.yaml
```

### Configuration reference

#### Sample

| Field | Description |
|-------|-------------|
| `sample.name` | Short sample identifier (used in output filenames) |
| `sample.source.path` | Folder containing TIFF files (optional when `source.channels` is set; defaults to first channel folder) |
| `sample.source.tiff_type` | `channelperfile` (BigStitcher-style) or `planeperfile` (Terastitcher-style) |
| `sample.source.channels` | Optional list of planeperfile roots (one folder per channel). Channel index follows list order. Requires `tiff_type: planeperfile`. |
| `sample.scratch` | Fast temp directory; **required disk space** for `planeperfile` XY-downsampled memmap (~`2 × ny × nx × nz × (vx/registres)²` bytes) plus detection scratch |
| `compute.workers` | Parallel slice workers for `channelperfile` stacks (default `4`). **`planeperfile` always uses 1 worker** — many TIFFs are faster read sequentially |
| `compute.max_in_memory_scratch_gb` | XY-downsampled scratch kept in RAM up to this size (default `24`); larger stacks use `sample.scratch` memmap |
| `sample.save_path` | All pipeline outputs and checkpoints |
| `sample.voxel_um` | Native voxel size `[x, y, z]` in µm — **required** for preprocess |

#### Atlas

| Field | Description |
|-------|-------------|
| `atlas.provider` | `allen` (10 µm) or `perens` (20 µm) |
| `atlas.resolution_um` | Atlas resolution in µm |
| `atlas.atlas_dir` | Directory with template + annotation NIfTIs; add `annotation_boundary_10.nii.gz` for Allen atlas boundary overlays (or use `LIGHTSUITE_ATLAS_PATH`) |

#### Registration

| Field | Description | Default |
|-------|-------------|---------|
| `registration.resolution_um` | Downsample target for registration volumes | `20` |
| `registration.channel_primary` | Primary channel for alignment | `1` |
| `registration.channel_secondary` | Optional second channel for dual-channel mutual information | `null` |
| `registration.bspline_spatial_scale_mm` | B-spline grid spacing in mm; smaller = finer warping, but also more free parameters to constrain | `0.64` |
| `registration.bspline_bending_weight` | Weight of the Elastix `TransformBendingEnergyPenalty` smoothness term; `0` disables it (MATLAB parity) | `2.0` |
| `registration.control_point_weight` | Landmark weight in Elastix (0–1) | `0.2` |
| `registration.augment_points` | Add thinned auto-landmarks to user control points | `false` |
| `registration.ap_pair_tolerance_vox` | AP residual tolerance for `refine-auto-points` (registration voxels) | `12.0` |
| `registration.use_slice_correspondence_affine` | Compose align-slices correspondence into the affine pre-warp before B-spline | `true` |
| `registration.use_slice_correspondence_landmarks` | Add align-slices anchors as extra B-spline landmarks in register | `true` |
| `registration.correspondence_landmark_weight` | Minimum landmark metric weight when correspondence landmarks are merged | `0.2` |
| `registration.correspondence_landmark_max_count` | Cap on total B-spline landmark pairs after merging correspondence anchors | `96` |
| `registration.ap_pair_min_kept` | Minimum auto pairs kept after AP filtering | `24` |
| `registration.orientation` | Axis permutation, e.g. `[1, 2, 3]`; flips use negative indices | auto |

#### Export

| Field | Description | Default |
|-------|-------------|---------|
| `export.save_registered_volume` | Write atlas-space TIFFs under `volume_registered/` | `false` |
| `export.write_cells_csv` | Write parcellation intensity CSVs | `true` |

---

## Step-by-step commands

Run all commands from the repository root with your config path.

### 0. Check environment

```bash
uv run lightsuite doctor -c my_mouse.yaml
```

Fix any failed checks (Elastix, atlas paths) before continuing.

### 1. Preprocess

Downsamples each channel to registration resolution and writes checkpoint files.

```bash
uv run lightsuite brain preprocess -c my_mouse.yaml
```

**Outputs** in `save_path`:

- `chan_{N}_sample_register_{20}um.tif` — one multi-page TIFF per channel
- `regopts.json` — volume metadata and paths (replaces `regopts.mat`)
- `sample_reference.json` — native sample-space grid for external segmentation exports

### 2. Check orientation (manual)

Compare your sample against the atlas before coarse registration. Opens a Napari viewer with the atlas volume on the left and the permuted sample on the right.

```bash
uv run lightsuite brain check-orientation -c my_mouse.yaml
```

Requires `uv sync --extra gui`.

**Workflow:**

1. Use the three dropdowns to map sample dimensions to atlas dimensions (with optional flips).
2. Click **Update preview** after each change.
3. Scroll through slices and use Napari's **change order of visible axes** control (**Ctrl+E**) to inspect different projections.
4. When anatomical axes align (AP, DV, LR), click **Save orientation && close**.

The viewer uses a **dual-panel layout** (atlas left, sample right) with one 3D volume per side. The sample is always shown at native resolution; the atlas is resampled to match for display when sizes differ.

**Output:** `brain_orientation.txt` — used by init-registration and downstream stages.

If you already know the permutation, set `registration.orientation` in YAML instead.

### 3. Initial registration

Coarse similarity alignment of sample to atlas using BCPD (or Open3D ICP fallback) on gradient point clouds.

```bash
uv run lightsuite brain init-registration -c my_mouse.yaml
```

Defaults match MATLAB (`cloudthres=5`, `pcdownsample=0.1` inside `extractSamplePoints.m`).
You normally do not set these in YAML. If `init_registration_diagnostics.json` reports a
**sparse sample cloud** and coarse alignment looks poor, see
[Advanced init-registration tuning](#advanced-init-registration-tuning) below.

**Outputs:**

- Updated `regopts.json` (`original_trans`, auto control point pairs)
- `brain_orientation.txt` (if orientation was set in config)
- `dim{1,2,3}_initial_registration.png` — eight sample slices per axis with warped atlas annotation edges overlaid (MATLAB `plotAnnotationComparison` style)

If orientation is wrong, set `registration.orientation` in YAML or edit `brain_orientation.txt`, then re-run init-registration.

#### Advanced init-registration tuning

These parameters affect only the **coarse BCPD step** (`extractSamplePoints.m` port). They are
omitted from example configs because MATLAB defaults work for most brains.

| Field | Default | When to change |
|-------|---------|----------------|
| `registration.cloud_threshold` | `5.0` (MATLAB `cloudthres`) | Lower (e.g. `3.0`) if diagnostics warn of a sparse sample cloud and `dim*_initial_registration.png` looks under-constrained |
| `registration.sample_cloud_subsample` | `0.1` (hardcoded in MATLAB) | Rarely change; fixed at MATLAB parity. Not listed in example YAMLs |

```yaml
registration:
  cloud_threshold: 3.0   # advanced only — more gradient points for sparse samples
```

### 4. Align slices (optional, recommended)

Opens a Napari dual-pane GUI to match **sample anatomy to the correct atlas plane** along **all three volume axes** (Y, X, Z), before placing control points.

```bash
uv run lightsuite brain align-slices -c my_mouse.yaml
```

Requires `uv sync --extra gui`.

**Workflow:**

1. Use the **Volume axis** control (1=Y, 2=X, 3=Z) to choose which axis you are aligning.
2. For each axis, confirm ~20 evenly spaced anchor slices:
   - Sample slice is shown on the **left** (fixed).
   - Scroll the atlas on the **right** with **PgUp** / **PgDn**, the atlas plane spinbox, or **mouse wheel over the atlas panel**.
3. When anatomy matches, press **Enter** or click **Confirm slice** — advances to the next anchor.
4. After the last anchor on an axis, you are prompted to continue on the next axis (or use **Next axis ▶**).
5. Use **←** / **→** to revisit anchors; switch axes any time with the dropdown.
6. Click **Save && Close** when finished (all three axes recommended).

**Output:**

- `slice_correspondence.json` (v2) — per-axis sample index ↔ atlas plane maps under `axes.{1,2,3}`

`match-points` and `refine-auto-points` read this file and pre-fill / filter using the curve for each cut axis. Older v1 files (single AP axis) are still loaded and migrated on save.

**Automated (all three axes, no GUI):**

```bash
uv run lightsuite brain align-slices -c my_mouse.yaml --headless
```

### 5. Refine auto points (optional)

Filters automatic landmark pairs from **init-registration** using the **multi-axis** slice maps from **align-slices**. A pair is removed if its atlas position disagrees with the correspondence curve on **any** confirmed axis. Useful for **auto-only** registration on stretched or sparse samples.

```bash
uv run lightsuite brain refine-auto-points -c my_mouse.yaml
```

Requires `slice_correspondence.json` with confirmed anchors (run **align-slices** first). To auto-estimate all three axes without the GUI:

```bash
uv run lightsuite brain refine-auto-points -c my_mouse.yaml --bootstrap-correspondence
```

Re-run with `--force` after changing correspondence or re-running init-registration.

**Outputs:**

- Updated `regopts.json` (`autocpsample` / `autocpatlas`, `auto_points_refined`, `auto_points_mode`)
- `auto_points_refine_stats.json` — pair counts, active axes, and median max-axis residual before/after

### 6. Match control points (optional)

Opens a Napari dual-pane GUI: sample on the left, atlas on the right.

**Optional:** you can skip this step and run **register** using only the auto control point pairs written by **init-registration** (same as MATLAB when `*tform.mat` is absent).

```bash
uv run lightsuite brain match-points -c my_mouse.yaml
```

Requires `uv sync --extra gui`.

**Workflow:**

1. Navigate slices with the **Navigation** panel on the right:
   - **←** / **→** arrow keys for previous/next chooselist slice (or use the panel buttons)
   - **Atlas plane** spinbox, **PgUp** / **PgDn**, or **mouse wheel over the atlas panel** to scroll along the atlas cut axis until anatomy matches the sample (like MATLAB)
   - Change **Slice #** and click **Show slice** to jump to a specific chooselist index
   - Toggle the atlas boundary overlay with **O** or the checkbox in Navigation
   - The spinbox also refreshes when you press Enter or change the value
   - Do **not** use Napari’s bottom dimension slider — it does not switch chooselist slices
2. Click corresponding landmarks on **sample** (left) and **atlas** (right). Markers are **numbered**; pair **1** on the sample matches pair **1** on the atlas, and so on.
3. When point counts match on a slice, the affine fit updates automatically and the **red atlas boundary overlay** on the sample refits (toggle with **O**).
4. **Backspace** removes the most recently placed point (on whichever side you clicked last).
5. Use **Clear slice points** to reset the current slice.
6. Click **Save && Close** when finished.

**Tips:**

- Use at least **16 matched pairs** spread along the anterior–posterior axis; more is better.
- Good landmarks: ventricle boundaries, corpus callosum, major nuclei outlines.
- Avoid damaged tissue; estimate where structures should be instead of bending the atlas into artifacts.

**Output:**

- `atlas2histology_tform.json` — control points and manual alignment (replaces `atlas2histology_tform.mat`)

For automated tests only:

```bash
uv run lightsuite brain match-points -c my_mouse.yaml --headless
```

### 7. Register (Elastix B-spline)

Runs affine + deformable registration using your control points and Elastix.

When `slice_correspondence.json` exists with confirmed anchors (from **align-slices**), **register**:

1. Composes a correspondence-informed affine correction before warping the atlas (step 1).
2. Adds anchor tissue centroids as extra B-spline landmarks, constraining through-plane alignment during Elastix (step 2).

This applies the axis-wise slice-index maps as geometry initialization rather than hard-filtering auto control points. Disable either step with `registration.use_slice_correspondence_affine: false` or `registration.use_slice_correspondence_landmarks: false`.

```bash
uv run lightsuite brain register -c my_mouse.yaml
```

Faster (lower quality) single-resolution schedule:

```bash
uv run lightsuite brain register -c my_mouse.yaml --single-step
```

**Outputs:**

- `transform_params.json` — full transform metadata (replaces `transform_params.mat`)
- `bspline_samp_to_atlas_20um.txt` — inverse B-spline (sample → atlas)
- `bspline_atlas_to_samp_20um.txt` — forward B-spline
- `registration_diagnostics.json` — registration checkpoint: affine/B-spline residuals, annotation overlap, GOOD/MODERATE/POOR status
- `affine_fit_stats.json` — affine landmark residuals (voxels): median/p95/max, auto vs manual, coarse baseline
- `correspondence_affine_stats.json` — slice-correspondence affine correction (anchor count, residual before/after)
- `correspondence_landmark_stats.json` — correspondence B-spline landmarks merged into Elastix
- `{name}_dim{1,2,3}_affine_registration.png` — eight sample slices per axis with warped atlas region outlines overlaid (same style as `dim*_initial_registration.png`)
- `{name}_dim{1,2,3}_bspline_registration.png` — same layout after B-spline
- `elastix_temp/` — Elastix working directory (keep until register finishes)

#### Wavy or over-warped B-spline contours

If the affine preview looks right but the B-spline preview shows tangled, wandering region
outlines, the deformation has too many degrees of freedom for the amount of image evidence
driving it — not a bad initialization.

Elastix draws ~5000 Mattes MI samples per iteration from a single small `SampleRegionSize`
cube. On a wide acquisition canvas (a 2048² mesoSPIM frame where the brain fills under a
third of the field, say) the B-spline grid has tens of thousands of parameters, most control
points sit over background, and each one receives only a handful of informative updates. The
rest of its trajectory is a random walk. Two things keep that in check:

- `registration.bspline_bending_weight` adds a smoothness prior. `2.0` is a good starting
  point; raise it if contours still wander. `0` reproduces MATLAB, which has no penalty term.
- `registration.bspline_spatial_scale_mm` controls how many parameters exist at all. Going
  from `0.64` to `1.0` mm cuts the grid roughly threefold. Counter-intuitively this usually
  *improves* landmark agreement on wide canvases, because the remaining control points are
  each much better constrained.

Check `elastix_temp/IterationInfo.0.R3.txt` to confirm the diagnosis. A `2:Metric0` column
that is frequently exactly `0.000000` means those iterations sampled a region with no atlas
overlap and contributed nothing; a `2:Metric1` (landmark) column that is flat across the
finest resolution means the last pyramid level is adding noise rather than alignment.

Cropping the registration volume to the brain with `registration.sample_content_crop: auto`
reduces the wasted grid too, though it is a smaller effect than the two settings above.

### 8. Export

Warps channels to atlas space (default) and optionally exports atlas labels warped onto the **registration grid** (sample space). See [Registration spaces](registration_spaces.md).

```bash
uv run lightsuite brain export -c my_mouse.yaml --save-volume --write-csv
uv run lightsuite brain export -c my_mouse.yaml --space both --save-volume --write-csv
```

Flags override YAML defaults:

- `--save-volume` / `--no-save-volume`
- `--write-csv` / `--no-write-csv`
- `--space atlas|sample|both` — output coordinate space(s)

**Outputs** in `volume_registered/`:

- `chan_{NN}_registered_atlas.tif` — if `--save-volume` (atlas space)
- `sample_space/` — warped annotation, template, division labels on registration grid
- `chan{NN}_intensities.csv` — regional median intensity, std, volume per hemisphere
- `chan{NN}_intensities.json` — same statistics in JSON form
- `chan{NN}_region_stats_sample.csv` / `region_stats_sample.csv` — sample-space stats when enabled

Allen parcellation CSV export requires `parcellation_to_parcellation_term_membership.csv` (see [Installation](installation.md)).

---

## Full command cheat sheet

```bash
export CONFIG=my_mouse.yaml

uv run lightsuite doctor -c $CONFIG
uv run lightsuite brain validate-config -c $CONFIG
uv run lightsuite brain preprocess -c $CONFIG
uv run lightsuite brain check-orientation -c $CONFIG
uv run lightsuite brain init-registration -c $CONFIG
uv run lightsuite brain align-slices -c $CONFIG
uv run lightsuite brain refine-auto-points -c $CONFIG
uv run lightsuite brain match-points -c $CONFIG
uv run lightsuite brain register -c $CONFIG
uv run lightsuite brain export -c $CONFIG --save-volume --write-csv
uv run lightsuite brain import-annotations -c $CONFIG
uv run lightsuite brain inspect-imports -c $CONFIG
uv run lightsuite brain inspect-imports -c $CONFIG --space sample
```

`inspect-imports` defaults to **atlas** space (Perens/Allen export grid). Use `--space sample` to overlay imported points/masks on the **20 µm registration grid** (`chan_*_sample_register_*um.tif` from preprocess). Warped atlas **template** and **annotation** labels are shown when present under `volume_registered/sample_space/` — run `lightsuite brain export --space sample` first (or set `export.spaces: [sample]`).

Atlas-space Napari QC applies the same canonical coronal orientation used in registration plots (Perens coronal is no longer upside-down vs Allen). Annotation volumes render as **label** layers rather than float images so CCF structure ids are fully visible.

---

## Output file layout

```
<save_path>/
├── regopts.json                          # Preprocess + init-registration state
├── sample_reference.json                 # Native grid for external segmentation
├── brain_orientation.txt                 # Axis permutation
├── slice_correspondence.json             # AP sample ↔ atlas plane map (align-slices)
├── auto_points_refine_stats.json         # AP filter stats (refine-auto-points)
├── correspondence_affine_stats.json      # Correspondence affine in register
├── correspondence_landmark_stats.json    # Correspondence B-spline landmarks
├── atlas2histology_tform.json            # Manual control points
├── transform_params.json                 # Final registration parameters
├── bspline_samp_to_atlas_20um.txt
├── bspline_atlas_to_samp_20um.txt
├── chan_1_sample_register_20um.tif       # Downsampled volumes (per channel)
├── dim1_initial_registration.png
├── {name}_dim1_affine_registration.png
├── {name}_dim1_bspline_registration.png
├── elastix_temp/                         # Elastix forward run
└── volume_registered/
    ├── chan_01_registered_atlas.tif
    ├── chan01_intensities.csv          # legacy wide intensity table
    ├── chan01_intensities.json
    ├── chan01_region_stats.csv         # tidy per-channel table (with region names)
    └── region_stats.csv                # combined tidy table (analysis region-stats)
```

Checkpoints use **JSON** instead of MATLAB `.mat` files. Legacy MATLAB outputs in the same folder are not read automatically — re-run the Python stages to produce JSON checkpoints.

---

## Resuming and re-running

Each stage reads checkpoints from `save_path`:

| Stage | Requires |
|-------|----------|
| preprocess | Valid YAML, TIFF source path |
| init-registration | `regopts.json` from preprocess |
| match-points | `regopts.json` with `original_trans` |
| register | `regopts.json` (with `original_trans` and auto pairs) + Elastix; `atlas2histology_tform.json` optional |
| export | `regopts.json` + `transform_params.json` + transformix |

Re-running a stage overwrites its outputs and downstream dependencies. After changing control points, re-run **register** and **export**.

---

## Dual-channel registration

If you have two structural contrasts at the same resolution (e.g. autofluorescence + fluorescent label), set:

```yaml
registration:
  channel_primary: 1
  channel_secondary: 2
  dual_channel_mi_weight_autofluor: 0.4
  dual_channel_mi_weight_signal: 0.4
```

Both channels are preprocessed in step 1; register uses dual fixed-image mutual information in Elastix.

### Multi-channel planeperfile (Terastitcher)

When each channel lives in its own folder of Z-plane TIFFs (e.g. `channel_488/RES(...)/`, `channel_561/RES(...)/`), list the roots under `source.channels` instead of running one config per channel:

```yaml
sample:
  source:
    tiff_type: planeperfile
    channels:
      - /data/Gilda/output/channel_488/RES(5690x5834x2052)
      - /data/Gilda/output/channel_561/RES(5690x5834x2052)
  save_path: /data/Gilda/output/registered_perens   # shared across channels

registration:
  channel_primary: 1
  channel_secondary: 2
```

All channel folders must have the same `(ny, nx, nz)` and matching slice ordering (files sorted by name). Export applies the same transform to every preprocessed channel.

---

## Import external annotations

After `register`, warp native sample-space coordinates or masks into atlas space.

Full specification: **[Annotation import](annotation_import.md)**.

```bash
uv run lightsuite brain import-annotations -c my_mouse.yaml
```

Supported formats (YAML `import.annotations`):

| `format` | Input | Notes |
|----------|-------|-------|
| `points_csv` | CSV with `x,y,z` columns | 1-based native voxel indices; see `sample_reference.json` |
| `mask_tiff` | 3D binary TIFF | Shape `(Y,X,Z)` at native resolution; matches `sample_reference.json` |

Outputs land in `volume_registered/` (`*_atlas_coords.npz`, optional CSV, mask TIFF).

---

## Region statistics and cell counts

After `export` (and optionally `import-annotations`), build a tidy region table:

```bash
uv run lightsuite analysis region-stats -c my_mouse.yaml
```

This reads the existing atlas-space outputs (no transformix run) and writes
`volume_registered/region_stats.csv`, a **long-form** table with one row per
`(sample, channel, parcellation_index, hemisphere, metric)` measurement:

| Column | Notes |
|--------|-------|
| `sample`, `channel`, `atlas` | provenance |
| `parcellation_index` | atlas label value (Allen ABC index, or Perens/Allen CCF id) |
| `acronym`, `name`, `structure`, `division` | region metadata in a common **Allen** vocabulary |
| `hemisphere` | `right` / `left` (same split as the intensity stats) |
| `metric` | `median_intensity`, `std`, `volume_mm3`, `cell_count`, `cell_density` |
| `value` | the measurement |

Intensity/std/volume rows come from each `chanXX_intensities.json`; `cell_count`
and `cell_density` come from binning each `*_atlas_coords.npz` point cloud into
atlas regions (gated by `analysis.count_points`). The Perens atlas is translated
into the Allen division grouping where CCF ids match.

`export` also writes a per-channel `chanXX_region_stats.csv` (same schema) unless
`analysis.write_tidy_csv: false`.

```yaml
analysis:
  write_tidy_csv: true       # emit chanXX_region_stats.csv during export
  count_points: true         # bin imported points into per-region counts/density
  point_labels: null         # restrict to specific import labels; null = all
```

---

## Cross-subject group analysis

After each subject has a `region_stats.csv` (step 11), define a **cohort YAML** listing
subjects, experimental groups, and comparison options. See `examples/cohort_region_stats.yaml`.

```bash
uv run lightsuite analysis validate-cohort -c examples/cohort_region_stats.yaml
uv run lightsuite analysis group-stats -c examples/cohort_region_stats.yaml
```

Each sample entry needs either:

| Field | Role |
|-------|------|
| `group` | Experimental group label (e.g. `control`, `treatment`) |
| `id` | Subject id (optional; defaults to brain config `sample.name`) |
| `region_stats` | Direct path to `volume_registered/region_stats.csv` |
| `config` | Brain pipeline YAML — `region_stats` resolved from `sample.save_path` |

### Cohort configuration

```yaml
name: amyloid_study
output_dir: ./cohort_results/amyloid_study

samples:
  - id: M001
    group: control
    config: examples/M001.yaml
  - id: M002
    group: treatment
    region_stats: /data/M002/volume_registered/region_stats.csv

group_analysis:
  channels: [1]                    # null = all channels
  metrics: [median_intensity]      # null = all metrics
  hemispheres: [right, left]       # null = both
  rollups: [division]              # optional division/structure summaries
  comparisons:
    - [control, treatment]
  test: mannwhitney
  fdr_alpha: 0.05
```

### Cohort outputs

Written under `output_dir`:

| File | Contents |
|------|----------|
| `cohort_region_stats_long.csv` | All subjects concatenated (filtered) with `group` column |
| `group_summary_by_region.csv` | Per-group mean, std, SEM, `n_subjects` at region level |
| `group_comparisons_by_region.csv` | Pairwise Mann–Whitney tests with BH-FDR `q_value` |
| `group_summary_by_division.csv` | Division rollup (if `rollups` includes `division`) |
| `group_comparisons_by_division.csv` | Division-level comparisons (optional) |

Division/structure rollups first average region values **within each subject**, then
summarize across subjects — so each subject contributes equally regardless of how many
fine regions fall in a division.

---

## Visualization

Simple matplotlib plots read `region_stats.csv`, legacy `chanXX_intensities.csv`, or cohort
`group_summary_by_division.csv`. No extra dependencies beyond the core package.

### Per-subject plots

Write figures under `sample.save_path/plots/` (not inside the LightSuite repository).

```bash
PLOTS=/media/gbm/NVME2/ALICe-pipelines-data/JulieBuron/registered_allen/plots
mkdir -p "$PLOTS"

# From brain config (resolves volume_registered/region_stats.csv)
uv run lightsuite analysis plot-division-bars \
  -c examples/config/mesoSPIM/JulieBuron.yaml -o "$PLOTS/division_bars.png" --channel 1

uv run lightsuite analysis plot-lr-scatter \
  -c examples/config/mesoSPIM/JulieBuron.yaml -o "$PLOTS/lr_scatter.png" \
  --keep-division Isocortex --keep-division Thalamus

uv run lightsuite analysis plot-top-regions \
  -c examples/config/mesoSPIM/JulieBuron.yaml -o "$PLOTS/top_regions.png" \
  --metric cell_count --top-n 15 --channel 1

# Or pass a CSV directly
uv run lightsuite analysis plot-division-bars \
  -i /data/M001/volume_registered/region_stats.csv \
  -o /data/M001/plots/division_bars.png --metric median_intensity --channel 1
```

| Command | Description |
|---------|-------------|
| `plot-division-bars` | Horizontal bars: left vs right **per division** (sum for `cell_count`, mean for intensity by default) |
| `plot-lr-scatter` | Left vs right **per region**, points coloured by division |
| `plot-top-regions` | Horizontal bars: left vs right for the **top N regions** ranked by left + right total (`--top-n`, default 10) |
| `plot-group-division` | Cohort **group means by division** (from `group_summary_by_division.csv`) |

Common options: `--metric` (`median_intensity`, `cell_count`, …), `--channel`, `--aggregate` (`auto`, `sum`, `mean`),
`--title`, `--exclude-division`, `--dpi`. Scatter also supports `--keep-division`, `--axis-min`, `--axis-max`.
Top-region plots support `--keep-division` and `--top-n` / `-n` (default 10).

`--aggregate auto` (default) sums fine regions within each division for `cell_count`, and averages them for
continuous metrics such as `median_intensity`.

Division bar plots also write a `.csv` aggregation next to the PNG.

### Cohort plot

After `group-stats` with `rollups: [division]`:

```bash
uv run lightsuite analysis plot-group-division \
  -i ./cohort_results/study/group_summary_by_division.csv \
  -o ./cohort_results/study/group_division_bars.png
```

Brain atlas heatmaps (brainglobe slice figures) are planned as an optional extra; these plots
cover the main QC views from the legacy notebook workflow.

### Napari division viewer

Interactive QC: toggle which Allen **divisions** contribute voxels and mask all registered
channel layers at once (ported from the legacy `ccf_division_napari_viewer` notebook).

Requires `uv sync --extra gui` and registered atlas-space volumes:

```bash
uv run lightsuite brain export -c my_mouse.yaml --save-volume
uv run lightsuite analysis build-division-map -c my_mouse.yaml   # optional pre-cache
uv run lightsuite analysis view-divisions -c my_mouse.yaml
```

| Option | Purpose |
|--------|---------|
| `--stride 2` | Downsample 2× per axis for faster interaction on large brains |
| `--force` | Rebuild cached division labels beside the atlas NIfTIs |
| `--headless` | Validate inputs only (for CI / scripting) |

**Cached atlas assets** (built once per atlas, shared across samples):

- `division_labels_10um.tif` / `division_labels_20um.tif`
- `division_id_legend_10um.csv` / `division_id_legend_20um.csv`

The dock panel lists divisions from the legend; unchecked divisions become transparent (NaN)
on every channel image layer. Enable the optional `division labels` layer for anatomical
reference.

### Registration QC (unassigned-voxel score)

Naive alignment proxy from the legacy `unassigned_registration_score` notebook: among
voxels above an intensity threshold, what fraction fall in **unassigned** division
labels (`division_id == 0`)? A higher percentage suggests signal spilling outside the
atlas — often a registration warning sign.

```bash
uv run lightsuite brain export -c my_mouse.yaml --save-volume

# Inspect threshold in Napari (magenta = id-0 signal above threshold), then score
uv run lightsuite analysis registration-qc -c my_mouse.yaml --inspect --threshold 500

# Headless score with auto-threshold + sensitivity sweep
uv run lightsuite analysis registration-qc -c my_mouse.yaml --channel 1 --sweep
```

| Option | Purpose |
|--------|---------|
| `--threshold` / `-t` | Global intensity cutoff for “signal” voxels; auto-estimated if omitted |
| `--inspect` | Open Napari to validate threshold before scoring (`uv sync --extra gui`) |
| `--sweep` | Save threshold sensitivity CSV + PNG (9 points, 0.5×–1.5× reference threshold) |
| `--headless` | Skip Napari; write score only |

**Outputs** under `<save_path>/registration_qc/`:

| File | Contents |
|------|----------|
| `chan01_unassigned_score.csv` | Score row + sample/channel metadata |
| `chan01_unassigned_score.json` | Same summary as JSON |
| `chan01_threshold_sweep.csv` | Optional multi-threshold sweep |
| `chan01_threshold_sweep.png` | Optional sensitivity plot |

**Caveats:** global threshold sensitivity; boundary biology and acquisition artifacts can
inflate the score. Treat as a quick QC flag, not a geometric error metric.

---

## Analysis roadmap (future improvements)

| Area | Next steps |
|------|------------|
| **Replicate combine** | Median-match LLS/RLS (or other replicates) before averaging intensities — deferred for now |
| **Brain heatmaps** | Optional `brainglobe-heatmap` slice figures on tidy/cohort tables (Allen atlas) |
| **Registration QC** | Adaptive/local thresholds; boundary-band restriction; morphology cleanup |
| **Statistics** | Effect sizes, volcano plots, cluster-corrected voxel-wise group tests |
| **Segmentation** | Wire cellpose/stardist imports → `cell_count` / `cell_density` in one `region-stats` pass |
| **Cross-atlas** | Full Perens ontology-tree rollup; explicit Allen↔Perens report views |
| **Outputs** | Parquet / AnnData export; dashboard notebook template for cohort studies |
| **Performance** | OME-Zarr registered volumes; lazy loading in Napari for full-resolution brains |
| **Automation** | `analysis run-all` meta-command: region-stats → plots → cohort summary from one config |

---

## Known limitations

- **Cell detection** — not implemented; preprocessing warns if `detection.enabled: true`
- **Orientation GUI** — `lightsuite brain check-orientation` (Napari); or set `registration.orientation` in YAML
- **OME-Zarr / Imaris** — planned; TIFF only today
- **Perens division names** — cross-atlas division grouping fills in only where Perens CCF ids match the Allen ontology; full ontology-tree rollup is a planned enhancement

---

## Migrating from MATLAB

| MATLAB | Python |
|--------|--------|
| `opts` struct in demo script | YAML config file |
| `regopts.mat` | `regopts.json` |
| `atlas2histology_tform.mat` | `atlas2histology_tform.json` |
| `transform_params.mat` | `transform_params.json` |
| `matchControlPoints_unified` | `lightsuite brain match-points` (Napari) |
| `check_lightsuite_installation.m` | `lightsuite doctor` |

You can keep using the same TIFF layouts, atlas files, and Elastix version as the MATLAB pipeline.

### Intentional differences in the B-spline step

`build_bspline_params` otherwise mirrors `performMultObjBsplineRegistration.m` — same metrics,
weights, transform, optimizer, pyramid schedule, iteration counts, spatial-sample count and
`SampleRegionSize` formula. It deviates in three places:

- **Bending energy.** `registration.bspline_bending_weight` appends a
  `TransformBendingEnergyPenalty` metric that MATLAB does not have. Set it to `0` for parity.
- **Histogram bins.** MATLAB inherits `NumberOfFixedHistogramBins`/`NumberOfMovingHistogramBins`
  = 32 from `matlab_elastix`'s `elastix_default.yml`, and Mattes MI prefers those over
  `NumberOfHistogramBins`. Python writes them explicitly at 32 rather than letting Elastix fall
  back to its own default.
- **ASGD step estimation.** Python sets `ASGDParameterEstimationMethod` to
  `DisplacementDistribution`; MATLAB leaves it at Elastix's built-in `Original`, whose larger
  steps can fold the grid.
