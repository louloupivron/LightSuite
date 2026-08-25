# Configuration reference

LightSuite configs are YAML files validated by Pydantic models in `src/lightsuite/config/models.py` (brain/spinal) and `src/lightsuite/multires/config_models.py` (multires). The Napari **Config** dock (`lightsuite gui`) edits a subset of keys; everything else can still be set in YAML and is **preserved on Save** when the GUI does not touch that key.

## GUI vs YAML legend

| Symbol | Meaning |
|--------|---------|
| **GUI** | Shown and editable in the Napari Config dock for that workflow |
| **YAML** | YAML only (not in the Config form; unchanged on Save unless you edit the file) |
| **GUI†** | GUI always writes this value on Save (no control in the form) |

**Workflows:** *brain* = whole-brain lightsheet, *spinal* = spinal cord, *multires* = overview↔ROI registration.

**Templates:** `examples/brain_lightsheet.yaml`, `examples/spinal_cord.yaml`, `examples/multiresolution.yaml`.

---

## `sample` — paths and imaging grid

| Key | Brain | Spinal | Multires | Default | Description |
|-----|:-----:|:------:|:--------:|---------|-------------|
| `sample.name` | GUI | GUI | GUI | *(required)* | Sample id (logs, exports, filenames) |
| `sample.scratch` | GUI | GUI | GUI | *(required)* | Fast disk for intermediates |
| `sample.save_path` | GUI | GUI | GUI | *(required)* | Results / checkpoints directory |
| `sample.voxel_um` | GUI | GUI | — | brain: optional; spinal: required | Native voxel size µm `[x, y, z]` |
| `sample.source.format` | YAML | YAML | — | `auto` / `tiff_stack` | Source reader (`auto`, `tiff_stack`, `czi`, `imaris`) |
| `sample.source.path` | GUI | GUI | — | — | Stitched TIFF folder (`channelperfile`) or first channel folder |
| `sample.source.tiff_type` | GUI | GUI | — | `channelperfile` | `channelperfile` or `planeperfile` (Terastitcher-style) |
| `sample.source.channels` | GUI | GUI | — | — | List of per-channel folders (`planeperfile` only) |
| `sample.source.skip_corrupt_slices` | — | YAML | — | `false` | Skip bad slice TIFFs instead of failing (spinal) |

---

## `atlas` — reference anatomy

### Brain (`atlas`)

| Key | GUI | Default | Description |
|-----|:---:|---------|-------------|
| `atlas.provider` | GUI | `allen` | `allen`, `perens`, `princeton`, `rat` (local files); BrainGlobe subset when `source: brainglobe` |
| `atlas.source` | GUI | `files` | `files` (local NIfTI) or `brainglobe` |
| `atlas.brainglobe_name` | GUI | inferred | Registry name when `source: brainglobe` |
| `atlas.resolution_um` | GUI | `10` | Native atlas voxel size µm (local files; auto from BrainGlobe catalog) |
| `atlas.atlas_dir` | GUI | — | Folder with template + annotation NIfTIs (local files) |
| `atlas.content_trim` | YAML | `off` | Trim atlas padding for registration: `off`, `auto`, `manual` |
| `atlas.content_margin_vox` | YAML | `8` | Margin around auto-detected atlas foreground |
| `atlas.content_box` | YAML | — | Manual atlas crop `[y0,y1,x0,x1,z0,z1]` when `content_trim: manual` |

### Spinal (`atlas`)

| Key | GUI | Default | Description |
|-----|:---:|---------|-------------|
| `atlas.atlas_dir` | GUI | *(required)* | Fiederling atlas folder (`Template.tif`, `Annotation.tif`, CSVs) |

---

## `registration` — alignment parameters

### Brain

| Key | GUI | Default | Description |
|-----|:---:|---------|-------------|
| `registration.resolution_um` | GUI | `20` | Isotropic working resolution µm |
| `registration.channel_primary` | GUI | `1` | 1-based registration channel |
| `registration.channel_secondary` | GUI | — | Optional second channel for dual MI (`none` in GUI = unset) |
| `registration.bspline_spatial_scale_mm` | GUI | `0.64` | B-spline grid spacing mm |
| `registration.control_point_weight` | GUI | `0.2` | Landmark weight in Elastix (0–1) |
| `registration.augment_points` | GUI | `false` | Add auto-landmarks to control points |
| `registration.dual_channel_mi_weight_primary` | GUI | `0.4` | Dual MI weight for primary channel (alias: `dual_channel_mi_weight_autofluor`) |
| `registration.dual_channel_mi_weight_secondary` | GUI | `0.4` | Dual MI weight for secondary channel (alias: `dual_channel_mi_weight_signal`) |
| `registration.orientation` | GUI | from `brain_orientation.txt` | Axis permutation, e.g. `[1, -3, 2]`; empty in GUI = omit (use check-orientation) |
| `registration.canvas_mode` | GUI | `off` | Working grid: `off`, `pad`, `crop`, `union` |
| `registration.bspline_bending_weight` | YAML | `2.0` | Elastix bending-energy penalty (0 = MATLAB parity) |
| `registration.cloud_threshold` | YAML | `5.0` | Init-registration cloud extraction gate |
| `registration.sample_cloud_subsample` | YAML | `0.1` | Fraction of gradient points kept |
| `registration.outlier_ratio` | YAML | `0.01` | Outlier rejection ratio |
| `registration.bcpd_path` | YAML | — | Path to BCPD binary |
| `registration.use_slice_correspondence_affine` | YAML | `true` | Correspondence-informed affine before B-spline |
| `registration.use_slice_correspondence_landmarks` | YAML | `true` | Add align-slices anchors as B-spline landmarks |
| `registration.correspondence_landmark_weight` | YAML | `0.2` | Min landmark weight for correspondence pairs |
| `registration.correspondence_landmark_max_count` | YAML | `96` | Max merged correspondence landmark pairs |
| `registration.sample_content_crop` | YAML | `off` | Crop registration TIFFs to sample foreground: `off`, `auto`, `manual` |
| `registration.sample_content_margin_vox` | YAML | `8` | Margin for auto sample crop |
| `registration.sample_content_box` | YAML | — | Manual sample crop box |
| `registration.sample_content_trim_z` | YAML | `true` | Drop sparse Z slices in auto crop |

### Spinal

| Key | GUI | Default | Description |
|-----|:---:|---------|-------------|
| `registration.resolution_um` | GUI | `20` | Working resolution µm after preprocess |
| `registration.channel_primary` | GUI | `1` | Registration channel |
| `registration.control_point_weight` | GUI | `0.2` | Landmark weight (0–1) |
| `registration.straightening_lambda_pos` | YAML | `5000` | Straightening position penalty |
| `registration.straightening_lambda_ang` | YAML | `5000` | Straightening angle penalty |
| `registration.target_orientation_deg` | YAML | `90` | Target orientation after straightening |
| `registration.longitudinal_direction` | YAML | from `cord_orientation.txt` | `rostrocaudal` or `caudorostral` override |

---

## `detection` — built-in cell detection (brain)

| Key | GUI | Default | Description |
|-----|:---:|---------|-------------|
| `detection.enabled` | GUI | `true` | Enable detection stage (Python implementation limited; prefer `import`) |
| `detection.backend` | YAML | `classical` | `classical`, `cellpose`, `stardist` |
| `detection.cell_diameter_um` | YAML | `14` | Expected cell diameter µm |
| `detection.thresholds` | YAML | `[0.5, 0.4]` | Detection thresholds (first ≥ second) |
| `detection.channel` | YAML | — | Channel index for detection |
| `detection.debug` | YAML | `false` | Debug outputs |
| `detection.save_cell_images` | YAML | `false` | Save cell image crops |
| `detection.write_to_csv` | YAML | `false` | Write detection CSV |

---

## `compute` — parallelism and memory

| Key | Brain | Spinal | Multires | Default | Description |
|-----|:-----:|:------:|:--------:|---------|-------------|
| `compute.workers` | GUI | GUI | GUI* | `4` | Parallel workers (*multires: optional `compute` block in YAML) |
| `compute.use_gpu` | YAML | YAML | — | `false` | Reserved (not used in Python) |
| `compute.max_in_memory_scratch_gb` | YAML | YAML | — | `24` | RAM threshold for in-memory scratch vs memmap |

---

## `export` — registered volumes and spaces

| Key | Brain | Spinal | Multires | Default | Description |
|-----|:-----:|:------:|:--------:|---------|-------------|
| `export.registered_volume_format` | YAML | YAML | — | `tiff` | Registered volume format |
| `export.write_pyramid` | YAML | YAML | — | `true` | Write image pyramids |
| `export.write_cells_csv` | YAML | YAML | — | `true` | Legacy parcellation CSV flag |
| `export.save_registered_volume` | YAML | YAML | — | `false` | Write full registered channel TIFFs |
| `export.spaces` | YAML | YAML | — | `[atlas]` | `atlas` and/or `sample` export spaces |
| `export.save_sample_space_volume` | YAML | YAML | — | `true` | Warp atlas template/labels to sample grid |

---

## `analysis` — region stats and parcellation

| Key | Brain | Spinal | Multires | Default | Description |
|-----|:-----:|:------:|:--------:|---------|-------------|
| `analysis.intensity_metrics` | GUI | GUI | — | `median_intensity`, `std`, `volume_mm3` | Metrics in `region_stats*.csv`: also `mean_intensity`, `variance` |
| `analysis.parcellate_intensities` | YAML† | GUI | — | `true` | Compute intensity stats from registered channels (brain: tied to export) |
| `analysis.count_points` | GUI† | GUI† | — | `true` | Bin imported points into `cell_count` / `cell_density` (GUI always saves `true`) |
| `analysis.write_tidy_csv` | YAML | YAML | — | `true` | Emit long-form `chanXX_region_stats.csv` (brain export) |
| `analysis.intensity_channels` | YAML | YAML | — | all | Subset of channel indices to parcellate |
| `analysis.relative_intensity_to` | YAML | YAML | — | `none` | Spinal: `none` or `background` → adds `relative_median_intensity` |
| `analysis.point_labels` | YAML | YAML | — | all | Filter import labels for counting (stem of `*_atlas_coords.npz`) |
| `analysis.stats_spaces` | GUI | GUI | — | `[atlas]` | Write stats in `atlas` and/or `sample` space |
| `analysis.hemisphere_flip` | YAML | YAML | — | `false` | Swap L/R for Fiederling hemisphere mask |
| `analysis.hemisphere_keep_whole` | YAML | YAML | — | `true` | Also emit whole-cord rows alongside left/right |
| `analysis.top_n_regions` | YAML | YAML | — | `10` | Write `region_stats_top{N}.csv` (0 disables) |
| `analysis.top_n_rank_by` | YAML | YAML | — | per channel: `cell_count` if present else first intensity metric | Preferred ranking metric for the top-N CSV |

**GUI intensity metric checkboxes:** Median, Mean, Std, Variance, Volume mm³.

**GUI stats space checkboxes:** Atlas and/or Sample (`analysis.stats_spaces`).

---

## `import` — external annotations

Each entry in `import.annotations[]`:

| Key | Brain | Spinal | Multires | Description |
|-----|:-----:|:------:|:--------:|-------------|
| `format` | GUI | GUI | GUI | `points_csv` or `mask_tiff` |
| `path` | GUI | GUI | GUI | Path to annotation file |
| `label` | GUI | GUI | GUI | Output filename stem (optional) |

| Key | Brain | Spinal | Multires | Default | Description |
|-----|:-----:|:------:|:--------:|---------|-------------|
| `import.write_csv` | YAML | YAML | YAML | `true` | Write coordinate CSV alongside NPZ |

**Notes:**

- Brain: `points_csv` and `mask_tiff` supported; point counts merged into `region_stats.csv` after export/import.
- Spinal: `points_csv` only (`mask_tiff` must be converted to points).
- Multires: warps ROI-native annotations to overview space (see [Annotation import](annotation_import.md)).

---

## `multires` — overview↔ROI registration

Top-level block in multires YAML only (not used by brain/spinal pipeline configs).

| Key | GUI | Default | Description |
|-----|:---:|---------|-------------|
| `multires.pair_label` | GUI | — | Short label for manifest and outputs |
| `multires.pair_manifest` | GUI | — | JSON pair manifest path |
| `multires.channels.<name>.overview` | GUI | — | Per-channel overview (low mag): TIFF or stitched folder |
| `multires.channels.<name>.roi` | GUI | — | Per-channel ROI (high mag): TIFF or stitched folder |
| `multires.channels.<name>.overview_meta_path` | GUI | — | Sidecar metadata (required for stitched overview folders) |
| `multires.channels.<name>.roi_meta_path` | GUI | — | Sidecar metadata (required for stitched ROI folders) |
| `multires.geometry_mode` | GUI | `metadata` | `metadata` or `hybrid` (landmark-refined geometry) |
| `multires.landmarks.fit_mode` | GUI | `similarity` | `similarity`, `affine`, or `rigid` |
| `multires.landmarks.session_path` | YAML | auto | Landmark session JSON path |
| `multires.landmarks.min_pairs` | YAML | `3` | Minimum landmark pairs |
| `multires.registration.reference_channel` | GUI | first channel | Channel for primary overview↔ROI transform |
| `multires.registration.overlap_margin_um` | GUI | `0` | Expand/shrink overlap canvas µm |
| `multires.registration.write_full_overview_canvas` | GUI | `true` | Write full-overview registered canvas |
| `multires.registration.registration_bin` | YAML | `1` | Binning for registration |
| `multires.registration.max_slab_bytes` | YAML | `500000000` | Memory limit per slab |
| `multires.registration.experiment_name` | YAML | `default` | Experiment id in outputs |
| `multires.registration.elastix_stages` | YAML | `[translation, rigid]` | Elastix stage list |
| `multires.registration.apply_transform_to` | YAML | — | Channels receiving reference transform |
| `multires.overview_meta_path` | YAML | — | Global overview metadata |
| `multires.mesospim_geometry.overview.*` | YAML | — | `stage_xy_is_center`, `lateral_flip`, `itk_lateral_dim0_motor` |
| `multires.mesospim_geometry.roi.*` | YAML | — | Same fields for ROI volume |

---

## `multires` link — brain view-registration overlay

Block alias: top-level `multires:` on **brain** configs (not the multires workflow file).

| Key | GUI | Default | Description |
|-----|:---:|---------|-------------|
| `multires.config` | YAML | — | Path to multires YAML (reads `multires_regopts.json`) |
| `multires.checkpoint` | YAML | — | Direct path to multires checkpoint JSON |
| `multires.use_full_overview` | YAML | `false` | Load full-overview registered ROIs in brain GUI |

---

## GUI behavior notes

1. **Round-trip:** Keys marked **YAML** are left unchanged when you Save from the Config dock, as long as they already exist in the file.
2. **Templates:** *New → From template* loads example YAML; edit paths, then **Save as…**.
3. **Validation:** Save runs config validation; invalid paths or missing files are reported in the dock status line.
4. **BrainGlobe:** When `atlas.source` is BrainGlobe, the GUI hides atlas directory and resolution (taken from the catalog).
5. **TIFF layout:** `planeperfile` shows channel folder rows; `channelperfile` shows a single source path.
6. **Orientation:** Brain orientation can be set in the GUI or left empty to use `brain_orientation.txt` from `check-orientation`.
7. **Stage-specific GUIs:** Match-points, align-slices, check-orientation, and view-registration use separate Napari panels; their outputs are checkpoint files, not YAML fields.

---

## Related docs

- [Brain lightsheet usage](usage_lightsheet_brain.md)
- [Spinal cord usage](usage_spinal_cord.md)
- [Multiresolution usage](usage_multiresolution.md)
- [Annotation import](annotation_import.md)
- [Registration spaces](registration_spaces.md)
