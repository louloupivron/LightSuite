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
| 4 | `lightsuite brain match-points` | **Manual (GUI)** | `matchControlPoints_unified.m` |
| 5 | `lightsuite brain register` | Automated | `multiobjRegistration.m` |
| 6 | `lightsuite brain export` | Automated | `generateRegisteredBrainVolumes.m` |

Cell detection and mapping cells to atlas coordinates are **not yet ported**; set `detection.enabled: false` for now.

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
  control_point_weight: 0.2
  augment_points: false
  orientation: [1, 2, 3]          # or omit; uses brain_orientation.txt if present

detection:
  enabled: false                  # not yet implemented in Python

export:
  save_registered_volume: false
  write_cells_csv: true
```

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
| `sample.source.path` | Folder containing TIFF files |
| `sample.source.tiff_type` | `channelperfile` (BigStitcher-style) or `planeperfile` (Terastitcher-style) |
| `sample.scratch` | Fast temp directory (large binary intermediates when detection is enabled) |
| `sample.save_path` | All pipeline outputs and checkpoints |
| `sample.voxel_um` | Native voxel size `[x, y, z]` in µm — **required** for preprocess |

#### Atlas

| Field | Description |
|-------|-------------|
| `atlas.provider` | `allen` (10 µm) or `perens` (20 µm) |
| `atlas.resolution_um` | Atlas resolution in µm |
| `atlas.atlas_dir` | Directory with template + annotation NIfTIs (or use `LIGHTSUITE_ATLAS_PATH`) |

#### Registration

| Field | Description | Default |
|-------|-------------|---------|
| `registration.resolution_um` | Downsample target for registration volumes | `20` |
| `registration.channel_primary` | Primary channel for alignment | `1` |
| `registration.channel_secondary` | Optional second channel for dual-channel mutual information | `null` |
| `registration.bspline_spatial_scale_mm` | B-spline grid spacing in mm; smaller = finer warping | `0.64` |
| `registration.control_point_weight` | Landmark weight in Elastix (0–1) | `0.2` |
| `registration.augment_points` | Add thinned auto-landmarks to user control points | `false` |
| `registration.orientation` | Axis permutation, e.g. `[1, 2, 3]`; flips use negative indices | auto |
| `registration.cloud_threshold` | Edge threshold for coarse ICP point extraction | `5.0` |

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

### 2. Check orientation (manual)

Compare mean projections of your sample against the atlas before coarse registration. Opens a Napari viewer with atlas projections (top) and permuted sample projections (bottom).

```bash
uv run lightsuite brain check-orientation -c my_mouse.yaml
```

Requires `uv sync --extra gui`.

**Workflow:**

1. Use the three dropdowns to map sample dimensions to atlas dimensions (with optional flips).
2. Pick a **Projection view** (axis 0/1/2) to compare atlas (left) and sample (right) at full resolution.
3. Click **Update preview** after each change.
4. When anatomical axes align (AP, DV, LR), click **Save orientation && close**.

The viewer uses a **dual-panel layout** (atlas left, sample right) so projections do not overlap. The sample is always shown at native pixel size; the atlas is upscaled to match when sizes differ.

**Output:** `brain_orientation.txt` — used by init-registration and downstream stages.

If you already know the permutation, set `registration.orientation` in YAML instead.

### 3. Initial registration

Coarse similarity alignment of sample to atlas using Open3D ICP on edge point clouds.

```bash
uv run lightsuite brain init-registration -c my_mouse.yaml
```

**Outputs:**

- Updated `regopts.json` (`original_trans`, auto control point pairs)
- `brain_orientation.txt` (if orientation was set in config)
- `dim{1,2,3}_initial_registration.png` — sanity-check images

If orientation is wrong, set `registration.orientation` in YAML or edit `brain_orientation.txt`, then re-run init-registration.

### 4. Match control points (manual)

Opens a Napari dual-pane GUI: sample on the left, atlas on the right.

```bash
uv run lightsuite brain match-points -c my_mouse.yaml
```

Requires `uv sync --extra gui`.

**Workflow:**

1. Navigate slices with the **Navigation** panel.
2. Click corresponding landmarks on **sample** (left) and **atlas** (right).
3. When point counts match on a slice, the affine fit updates automatically.
4. Use **Clear slice points** to reset the current slice.
5. Click **Save && Close** when finished.

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

### 5. Register (Elastix B-spline)

Runs affine + deformable registration using your control points and Elastix.

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
- `{name}_dim{1,2,3}_affine_registration.png`
- `{name}_dim{1,2,3}_bspline_registration.png`
- `elastix_temp/` — Elastix working directory (keep until register finishes)

### 6. Export

Warps all channels to atlas space and optionally writes parcellation statistics.

```bash
uv run lightsuite brain export -c my_mouse.yaml --save-volume --write-csv
```

Flags override YAML defaults:

- `--save-volume` / `--no-save-volume`
- `--write-csv` / `--no-write-csv`

**Outputs** in `volume_registered/`:

- `chan_{NN}_registered_atlas.tif` — if `--save-volume`
- `chan{NN}_intensities.csv` — regional median intensity, std, volume per hemisphere
- `chan{NN}_intensities.json` — same statistics in JSON form

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
uv run lightsuite brain match-points -c $CONFIG
uv run lightsuite brain register -c $CONFIG
uv run lightsuite brain export -c $CONFIG --save-volume --write-csv
```

---

## Output file layout

```
<save_path>/
├── regopts.json                          # Preprocess + init-registration state
├── brain_orientation.txt                 # Axis permutation
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
    ├── chan01_intensities.csv
    └── chan01_intensities.json
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
| register | `regopts.json` + `atlas2histology_tform.json` + Elastix |
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

---

## Known limitations

- **Cell detection** — not implemented; preprocessing warns if `detection.enabled: true`
- **Orientation GUI** — `lightsuite brain check-orientation` (Napari); or set `registration.orientation` in YAML
- **OME-Zarr / Imaris** — planned; TIFF only today
- **Parcellation names** — CSV includes numeric region IDs; Allen name/structure/division columns from MATLAB are not yet joined in Python export

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
