# Annotation import (LightSuite Sample Space v1)

Import cell coordinates or segmentation masks from **any** external tool into atlas space. All tools must export into a single native sample-space convention before running `import-annotations`.

---

## Workflow overview

```mermaid
flowchart LR
    A[preprocess] --> B[sample_reference.json]
    B --> C[Segment in any tool]
    C --> D[Convert to LSS v1]
    D --> E[register]
    E --> F[import-annotations]
```

1. **`lightsuite brain preprocess`** — writes `regopts.json` and **`sample_reference.json`** (native grid contract).
2. **Segment** at **native resolution** on the same grid as the stitched sample (not the 20 µm registration preview TIFF).
3. **Convert** tool output to `points.csv` or `mask.tif` (see below). Vendor-specific formats (Arivis, LCT, etc.) need a one-time local converter.
4. **`lightsuite brain register`** — produces `transform_params.json`.
5. **`lightsuite brain import-annotations`** — warps annotations into atlas space.

---

## Native sample-space reference

After preprocess, read `<save_path>/sample_reference.json`:

```json
{
  "format": "lightsuite_sample_space_v1",
  "sample_name": "Gilda",
  "shape_yxz": [5834, 5690, 2052],
  "voxel_um": [2.03, 2.03, 3.0],
  "index_base": 1,
  "axis_order": "xyz",
  "coordinate_units": "voxel_indices",
  "orientation_applied": false
}
```

| Field | Meaning |
|-------|---------|
| `shape_yxz` | Volume shape **(Y, X, Z)** = **(ny, nx, nz)** at native resolution |
| `voxel_um` | Voxel size **[x, y, z]** in µm (same as `sample.voxel_um` in YAML) |
| `index_base` | Always **1** — first voxel center is at `(1, 1, 1)` |
| `axis_order` | Always **`xyz`** for point coordinates |
| `orientation_applied` | Always **false** — do not apply `registration.orientation` yourself; LightSuite applies it during the warp |

**Do not** segment on `chan_*_sample_register_20um.tif` unless you rescale coordinates back to native indices first. The reference grid is the **full-resolution** stitched stack described in `regopts.json`.

---

## Point coordinates — `points.csv`

CSV with a header row. Required columns: **`x`**, **`y`**, **`z`** (case-insensitive).

```csv
x,y,z
1863.0,102.0,2.0
1864.5,103.2,2.0
```

- **Units:** 1-based voxel indices at native resolution (floats allowed for center-of-mass).
- **Bounds:** `1 ≤ x ≤ nx`, `1 ≤ y ≤ ny`, `1 ≤ z ≤ nz` (from `shape_yxz`).
- **Extra columns** (e.g. `intensity`, `volume_um3`) — numeric values are preserved as optional features; text columns (e.g. segment names) are skipped.

Example file: `examples/annotation_sample/points.csv`.

---

## Segmentation mask — `mask.tif`

3D binary TIFF at **native resolution**:

- **Shape:** `(Y, X, Z)` = `shape_yxz` from `sample_reference.json`
- **Values:** `0` / `255` or `0` / `1`
- **Layout:** Z-stack of 2D pages (one page per Z slice) or a single 3D TIFF page
- **Voxel size:** must match `voxel_um` in `sample_reference.json` (no resampling at import)

Example: segment on a native-resolution export of the structural channel, save as `mask.tif`.

---

## YAML configuration

```yaml
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

| Key | Description |
|-----|-------------|
| `format` | `points_csv` or `mask_tiff` |
| `path` | Path to the CSV or TIFF file |
| `label` | Output filename stem (defaults to file stem) |

Run:

```bash
uv run lightsuite brain import-annotations -c my_mouse.yaml
```

---

## Outputs

Written to `<save_path>/volume_registered/`:

| Kind | Files |
|------|-------|
| Points | `{label}_atlas_coords.npz`, optional `{label}_atlas_coords.csv` |
| Mask | `{label}_registered_atlas.tif` |
| Summary | `import_annotations_summary.json` |

NPZ arrays:

- `atlasptcoords` — atlas-space coordinates
- `sampleptcoords` — native sample-space coordinates (input)

---

## Converting from external tools

Each site maintains a small converter script. Common translations:

| Source | Conversion |
|--------|------------|
| **0-based indices** | Add 1 to each coordinate |
| **`[z, y, x]` order** | Reorder to `[x, y, z]` |
| **Arivis Blob Finder CSV** | Use COM columns as `x,y,z`; add 1 if 0-based |
| **Imaris Statistics CSV** | `Position X/Y/Z` in µm → `int(pos / voxel_um) + 1` per axis; see `examples/notebooks/convert_imaris_to_lightsuite.ipynb` |
| **Imaris mask TIFF series** | One label slice per Z (`*_Z####.tif`, 0-based in filename); binarize `(plane > 0)` and stack to multi-page TIFF |
| **LCT JSON** `[[z,y,x],…]` | Reorder to `x,y,z`; add 1 |
| **Downsampled segmentation** | Resample mask/coordinates to native `shape_yxz` before import |

Validate against `sample_reference.json` before import:

- Point count and coordinate ranges
- Mask shape `(ny, nx, nz)` and voxel size

---

## Validation errors

| Error | Fix |
|-------|-----|
| Missing `sample_reference.json` | Run `preprocess` |
| Mask shape mismatch | Resample mask to native `(Y, X, Z)` |
| All points out of bounds | Check axis order and index base (+1 for 0-based tools) |
| Missing `transform_params.json` | Run `register` before import |

---

## Limitations

- **Native resolution only** — registration-resolution (20 µm) masks/coordinates are not accepted.
- **Allen cell-count parcellation** — atlas coordinates are written; region counts are not yet joined in Python.
- **OME-Zarr masks** — use TIFF at native resolution; zarr pyramids are not supported in v1.

---

## See also

- [Brain lightsheet usage](usage_lightsheet_brain.md) — full pipeline and command cheat sheet
- `examples/annotation_sample/` — minimal example CSV
