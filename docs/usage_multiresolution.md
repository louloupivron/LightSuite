# Multiresolution registration (overview ↔ ROI)

This guide walks through registering a **high-resolution ROI stack** (e.g. 9× mosaic or single FOV) to a **low-resolution overview** (e.g. 1.6× or 1×) acquired on the same sample at different magnifications.

Unlike the brain atlas pipeline, multiresolution registration does **not** use an atlas. It aligns two volumes in a shared physical coordinate frame using stage metadata (and optionally manual landmarks), then runs Elastix on the shared overlap region.

---

## Before you start

You will need:

- An **overview** volume and a matching **ROI** volume from the same sample
- Vendor metadata that defines each stack’s physical frame (origin, spacing, orientation)
- Fast scratch space for Elastix intermediates (`sample.scratch`)
- A results directory (`sample.save_path`)
- Elastix on your `PATH` ([Installation](installation.md))
- Python extras for registration and (for hybrid mode) the landmark GUI:

```bash
uv sync --extra registration --extra gui
```

Supported input layouts today:

| Layout | Example |
|--------|---------|
| Single multi-page TIFF | mesoSPIM `1-561-1x.tif` or 1.25X ROI tile |
| Plane-per-file folder | SmartSPIM `All_Channels/` or mesoSPIM TeraStitcher `RES(...)` mosaic |

---

## Pipeline overview

| Step | CLI command | Type | Purpose |
|:----:|-------------|------|---------|
| 0 | *(conversion)* | Automated | Build a **pair manifest** JSON from vendor exports |
| 1 | `lightsuite multires validate-config` | Check | Load YAML + manifest; confirm paths |
| 2 | `lightsuite multires inspect-geometry` | **Manual (GUI)** | Trial `lateral_flip` on mesoSPIM pairs *(when metadata alignment looks mirrored)* |
| 3 | `lightsuite multires match-points` | **Manual (GUI)** | Place landmark pairs *(hybrid mode only)* |
| 4 | `lightsuite multires check-geometry` | Automated / QC | FOV overlap plots and geometry report |
| 5 | `lightsuite multires export-preview` | Optional | Small TIFF crops for visual alignment checks |
| 6 | `lightsuite multires register` | Automated | Elastix translation / rigid on overlap crops |
| 7 | `lightsuite multires inspect-registration` | **Manual (GUI)** | Overlay the registered ROI on the overview in Napari |
| 8 | `lightsuite multires import-annotations` | Optional | Warp ROI-native segmentation into overview-native space |

The pipeline streams planes from disk. Full overview and ROI stacks are **not** loaded into RAM at once.

Orchestrate automated steps with `lightsuite multires run -c my.yaml` (see `lightsuite multires stages` for checkpoint status). Interactive stages (inspect-geometry, match-points, inspect-registration) also appear in the unified GUI: `lightsuite gui -c my.yaml`.

---

## Step 1 — Publish a pair manifest

The multires pipeline reads a **LightSuite Multires Pair v1** JSON manifest. Vendor-specific parsing lives in conversion notebooks and helper scripts — the core `lightsuite multires` commands only consume the manifest.

Each manifest describes:

- `overview` — low-magnification stack path + physical geometry (`shape_zyx`, `spacing_um`, `origin_um`, `direction`)
- `roi` — high-magnification stack path + geometry
- `pair_label` — short identifier used in output folders and default landmark filenames

### mesoSPIM (TIFF + `*_meta.txt`)

Use the CLI, conversion notebook, the vendor helper, or the sample build script:

```bash
# CLI (preferred for new acquisitions):
lightsuite multires build-manifest \
  --vendor mesospim \
  --sample-name OP39M2 \
  --pair-label spinal_cord_488_561 \
  --channels-json '{"488":{"overview":"/path/2.5X/ch488","roi":"/path/1.25X/ch488.tif"},...}' \
  --reference-channel 488 \
  --overview-meta /path/anchor_tile_meta.txt \
  -o /path/OP39M2_pair.json

# Notebook: examples/notebooks/convert_mesospim_to_multires.ipynb
# OP39M2 spinal cord (stitched 2.5X overview + 1.25X ROI, 488 + 561):
uv run python scripts/build_op39m2_multires_manifest.py
```

Single-channel hyperstack pair:

```python
from pathlib import Path
from lightsuite.multires.vendor.mesospim import build_mesospim_pair_manifest

manifest = build_mesospim_pair_manifest(
    sample_name="JulieBuron",
    pair_label="baseline_561_TC",
    overview_path=Path("/data/JulieBuron/1X/1-561-1x.tif"),
    roi_path=Path("/data/JulieBuron/3.2X/1-561-3.2x_TC.tif"),
    output_manifest_path=Path("/data/JulieBuron/multiresolution_results/converted/baseline_561_TC_pair.json"),
)
```

Stitched mesoSPIM overview (TeraStitcher `RES(...)` folder) plus ROI tile:

```python
from lightsuite.multires.vendor.mesospim import build_mesospim_pair_manifest

manifest = build_mesospim_pair_manifest(
    sample_name="OP39M2",
    pair_label="spinal_cord_488",
    overview_path=Path("/data/OP39M2/2.5X/output/channel_488/RES(12486x2162x1163)"),
    roi_path=Path("/data/OP39M2/1.25X/..._Ch488_....tiff"),
    overview_meta_path=Path("/data/OP39M2/2.5X/..._Mag2.5x_ch488_Tile0.tiff_meta.txt"),
    output_manifest_path=Path("/data/OP39M2/multiresolution_results/converted/op39m2_488_pair.json"),
)
```

Multichannel paths in the YAML config (recommended for OP39M2-style data):

```yaml
multires:
  pair_label: spinal_cord_488_561
  pair_manifest: /data/OP39M2/multiresolution_results/converted/OP39M2_spinal_cord_pair.json
  overview_meta_path: /data/OP39M2/2.5X/..._Tile0.tiff_meta.txt
  channels:
    "488":
      overview: /data/OP39M2/2.5X/output/channel_488/RES(...)
      roi: /data/OP39M2/1.25X/..._Ch488_....tiff
    "561":
      overview: /data/OP39M2/2.5X/output/channel_561/RES(...)
      roi: /data/OP39M2/1.25X/..._Ch561_....tiff
  registration:
    reference_channel: "488"
    apply_transform_to: ["561"]
```

`validate-config` / `check-geometry` / `register` rebuild the pair manifest from `multires.channels` automatically. The standalone builder remains available:

```bash
uv run python scripts/build_op39m2_multires_manifest.py
```

Or call the helper directly:

```python
from lightsuite.multires.vendor.mesospim import build_mesospim_multichannel_pair_manifest

manifest = build_mesospim_multichannel_pair_manifest(
    sample_name="OP39M2",
    pair_label="spinal_cord_488_561",
    reference_channel="488",
    channels={
        "488": {
            "overview": Path("/data/OP39M2/2.5X/output/channel_488/RES(...)"),
            "roi": Path("/data/OP39M2/1.25X/..._Ch488_....tiff"),
        },
        "561": {
            "overview": Path("/data/OP39M2/2.5X/output/channel_561/RES(...)"),
            "roi": Path("/data/OP39M2/1.25X/..._Ch561_....tiff"),
        },
    },
    overview_meta_path=Path("/data/OP39M2/2.5X/..._Tile0.tiff_meta.txt"),
    output_manifest_path=Path("/data/OP39M2/multiresolution_results/converted/op39m2_pair.json"),
)
```

### SmartSPIM / ASI (plane-per-file + `metadata.txt` / `metadata.json`)

For the bundled **Multi_RES_SCANs** example dataset, rebuild all three pair manifests with:

```bash
uv run python scripts/rebuild_multi_res_manifests.py
```

This writes manifests such as `cortex_9x_561_pair.json`, `cerebellum_9x_561_pair.json`, and `single_fov_561_pair.json` under `registration_results/converted/`.

For new SmartSPIM acquisitions, adapt `scripts/rebuild_multi_res_manifests.py` or use `lightsuite.multires.vendor.smartspim.build_smartspim_pair_manifest` with your overview/ROI paths and stage metadata. Both legacy tab-separated `metadata.txt` and JSON `metadata.json` (`sample_metadata` + `tiles`) are accepted.

---

## Step 2 — Create a YAML config

Copy one of the examples under `examples/config/multiresolution/` and edit paths.

### Minimal example (metadata-only geometry)

Use this when stage metadata alone places the ROI correctly inside the overview FOV (typical for mesoSPIM pairs with reliable `*_meta.txt` files):

```yaml
sample:
  name: JulieBuron
  save_path: /data/JulieBuron/multiresolution_results
  scratch: /fastssd/lightsuite_scratch

multires:
  pair_manifest: /data/JulieBuron/multiresolution_results/converted/baseline_561_TC_pair.json

  geometry_mode: metadata

  registration:
    overlap_margin_um: 0.0
    registration_bin: 1          # use 2 for faster smoke tests
    experiment_name: baseline_561_TC
    elastix_stages: [translation, rigid]
    write_full_overview_canvas: true
```

### Hybrid geometry (metadata crop + manual landmarks)

Use this when metadata gives a reasonable starting crop but fine placement needs manual correction (typical for SmartSPIM mosaic ↔ overview pairs):

```yaml
sample:
  name: Multi_RES_SCANs
  save_path: /data/Multi_RES_SCANs/registration_results
  scratch: /fastssd/lightsuite_scratch

multires:
  pair_manifest: /data/Multi_RES_SCANs/registration_results/converted/cortex_9x_561_pair.json

  geometry_mode: hybrid
  landmarks:
    session_path: null           # default: save_path/multires_landmarks_<pair_label>.json
    fit_mode: similarity         # similarity | rigid | affine
    min_pairs: 3

  registration:
    overlap_margin_um: 0.0
    registration_bin: 2
    experiment_name: cortex_9x_561_hybrid
    elastix_stages: [translation, rigid]
    write_full_overview_canvas: false
```

### Configuration reference

| Key | Default | Description |
|-----|---------|-------------|
| `multires.pair_manifest` | auto under `save_path/converted/` | Pair JSON path. Required unless `multires.channels` is set; when channels are set, the JSON is rebuilt from those paths |
| `multires.channels` | — | Optional map of channel → `{overview, roi}` paths (mesoSPIM). Preferred way to declare multichannel inputs in YAML |
| `multires.overview_meta_path` | — | Anchor tile `*_meta.txt` for stitched overview folders |
| `multires.pair_label` | from experiment/sample | Short identifier used in output folders and default landmark filenames |
| `multires.geometry_mode` | `metadata` | `metadata` — overlap from stage geometry only; `hybrid` — metadata crop + landmark fit |
| `multires.mesospim_geometry.overview.lateral_flip` | `[1, -1]` | Per-volume axis sign when building manifests from `multires.channels` (mesoSPIM only). **Must match on overview and ROI** |
| `multires.mesospim_geometry.roi.lateral_flip` | `[1, -1]` | Same as overview; see [mesoSPIM lateral_flip](#mesospim-lateral_flip-axis-convention) below |
| `multires.landmarks.session_path` | auto | Path to landmark JSON; when `null`, uses `multires_landmarks_<pair_label>.json` under `save_path` |
| `multires.landmarks.fit_mode` | `similarity` | Transform fitted from landmark pairs: `similarity`, `rigid`, or `affine` |
| `multires.landmarks.min_pairs` | `3` | Minimum matched pairs required before geometry / registration |
| `multires.registration.overlap_margin_um` | `0.0` | Expand or shrink the physical overlap box (µm) |
| `multires.registration.registration_bin` | `1` | Bin factor for Elastix (full-res transforms are re-applied when > 1) |
| `multires.registration.experiment_name` | `default` | Slug for output subfolder and filenames |
| `multires.registration.elastix_stages` | `[translation, rigid]` | Elastix parameter maps to chain |
| `multires.registration.write_full_overview_canvas` | `true` | Write ROI embedded in full overview grid |
| `multires.registration.max_slab_bytes` | `500000000` | Soft cap for streaming ROI slabs during resample (~500 MB) |
| `multires.registration.reference_channel` | manifest default | Laser/channel slug used for Elastix (multichannel manifests) |
| `multires.registration.apply_transform_to` | all non-reference channels | Additional channels that receive the saved transform without re-running Elastix |

Example configs in the repository:

| Config | Microscope | Geometry mode |
|--------|------------|---------------|
| [`JulieBuron_multires.yaml`](../examples/config/multiresolution/JulieBuron_multires.yaml) | mesoSPIM | `metadata`, `lateral_flip: [1, -1]` |
| [`JulieBuron_multires_manifest.yaml`](../examples/config/multiresolution/JulieBuron_multires_manifest.yaml) | mesoSPIM | alias of `JulieBuron_multires.yaml` |
| [`marianna_multires.yaml`](../examples/config/multiresolution/marianna_multires.yaml) | mesoSPIM (0.8× ↔ 2.5×, dual channel) | `metadata`, `lateral_flip: [-1, -1]` |
| [`OP39M2_multires_manifest.yaml`](../examples/config/multiresolution/OP39M2_multires_manifest.yaml) | mesoSPIM (stitched + dual channel) | `metadata` |
| [`Multi_RES_SCANs_cortex_9x_manifest.yaml`](../examples/config/multiresolution/Multi_RES_SCANs_cortex_9x_manifest.yaml) | SmartSPIM | `hybrid` |
| [`Multi_RES_SCANs_single_fov_manifest.yaml`](../examples/config/multiresolution/Multi_RES_SCANs_single_fov_manifest.yaml) | SmartSPIM | `hybrid` |
| [`Multi_RES_SCANs_cerebellum_9x_manifest.yaml`](../examples/config/multiresolution/Multi_RES_SCANs_cerebellum_9x_manifest.yaml) | SmartSPIM | `metadata` (geometry QC example) |

---

## Step 3 — Validate the config

```bash
uv run lightsuite multires validate-config -c my_multires.yaml
```

This loads the YAML, resolves the pair manifest, and confirms both volume paths exist. Fix any `FileNotFoundError` before continuing.

---

## Step 3b — Inspect mesoSPIM geometry (`lateral_flip`)

Skip when metadata overlap already looks correct in `check-geometry` / `export-preview`.

For mesoSPIM pairs declared via `multires.channels`, stage `*_meta.txt` files record motor positions but **not** how TIFF voxel rows/columns map onto those motors. LightSuite bridges that gap with `multires.mesospim_geometry.lateral_flip` (default `[1, -1]`, matching the JulieBuron reference dataset).

```bash
uv sync --extra registration --extra gui
uv run lightsuite multires inspect-geometry -c my_multires.yaml
```

The Napari viewer shows the overview overlap crop beside the ROI **physically resampled** onto the same grid (not a naive array resize). Toggle **Flip X** / **Flip Y** on **both** volumes together, watch the physical NCC readout, then click **Apply to YAML** to patch `multires.mesospim_geometry` in place.

Headless scoring (no GUI):

```bash
uv run lightsuite multires inspect-geometry -c my_multires.yaml --headless
uv run lightsuite multires inspect-geometry -c my_multires.yaml --headless --write-config
```

### mesoSPIM `lateral_flip` axis convention

Each `lateral_flip: [fx, fy]` entry sets the sign of the ITK direction cosines for the X and Y image axes derived from mesoSPIM stage metadata (`x_pos`, `y_pos`, pixel size, `stage_xy_is_center`). In short: it answers “does increasing column index move in the +motor direction?” for each lateral axis.

| Dataset | Typical `lateral_flip` | Notes |
|---------|---------------------|-------|
| JulieBuron (reference) | `[1, -1]` | Default in `MesospimGeometryConfig` |
| Marianna CMU (0.8× ↔ 2.5×) | `[-1, -1]` | Both axes inverted vs default; metadata-only registration NCC ~0.92 |

**Why Marianna needs a Y flip as well as X:** the meta sidecar only stores absolute stage coordinates. It does not record sample handedness, which motor is wired as ITK dimension 0, or whether the camera readout order matches the historical mesoSPIM “increasing row = decreasing motor Y” convention. With the default `[1, -1]`, the predicted FOV overlap is close enough that `check-geometry` slice NCC looks plausible (~0.7), but the ROI appears **mirrored along Y** relative to the overview. Flipping X alone (`[-1, 1]`) fixes left–right but leaves cranio–caudal reversed; only `[-1, -1]` aligns both lateral axes so SimpleITK can resample the 2.5× stack into 0.8× space before Elastix.

**Rules of thumb:**

- **Overview and ROI must use identical `lateral_flip`.** Mismatched direction cosines break `stream_resample_to_reference` and yield near-zero registration NCC even when FOV boxes intersect.
- Prefer **`inspect-geometry` physical NCC** over `check-geometry --level slice-qc` alone when diagnosing mirrors: slice QC resizes arrays without mapping through physical space and can score two different flip settings similarly.
- After changing `lateral_flip`, delete stale `elastix_roi_to_overview/<experiment_name>/` outputs (or use a new `experiment_name`) before re-running `register`, so checkpoints and TIFFs are not mixed across conventions.

---

## Step 4 — Place landmarks (hybrid mode only)

Skip this step when `geometry_mode: metadata`.

```bash
uv run lightsuite multires match-points -c my_multires.yaml
```

The Napari GUI shows side-by-side **overview** and **ROI** crops (derived from metadata overlap). For each anatomical feature:

1. Click a point on the overview panel (yellow).
2. Click the corresponding point on the ROI panel (cyan).
3. Repeat until you have at least `min_pairs` matched pairs (3 by default).

Use **Save && Close** when done. Landmarks are written to the session path (default: `save_path/multires_landmarks_<pair_label>.json`).

To initialise a session file without opening Napari (e.g. on a headless node before copying to a workstation):

```bash
uv run lightsuite multires match-points -c my_multires.yaml --headless
```

---

## Step 5 — Check geometry

```bash
uv run lightsuite multires check-geometry -c my_multires.yaml
```

By default this runs a **full** check: loads overlap crops, writes QC plots, and saves a checkpoint. Use lighter levels for quick iteration:

```bash
# FOV schematic only (no voxel IO)
uv run lightsuite multires check-geometry -c my_multires.yaml --level metadata-only

# One mid-plane slice + normalized cross-correlation
uv run lightsuite multires check-geometry -c my_multires.yaml --level slice-qc

# Full overlap stacks + cropped overview preview TIFF
uv run lightsuite multires check-geometry -c my_multires.yaml --level full
```

### Geometry outputs

Written under `save_path/geometry/<pair_label>/`:

| File | Description |
|------|-------------|
| `geometry_report.json` | Overlap box, alignment metrics, optional landmark fit |
| `fov_overlap.png` | Schematic FOV diagram |
| `geometry_overlap_qc.png` | Slice or stack overlay QC |
| `<overview>_overview_crop.tif` | Overview cropped to the overlap box *(full level only)* |

A checkpoint is always updated at `save_path/multires_regopts.json`.

**Review before registering:** open `fov_overlap.png` and confirm the ROI box sits inside the overview. In hybrid mode, check the reported landmark RMS error (aim for low tens of µm or better depending on resolution).

---

## Step 6 — (Optional) Export alignment previews

Export small TIFF crops at several Z planes for offline inspection (Fiji, etc.):

```bash
uv run lightsuite multires export-preview -c my_multires.yaml
```

Optional flags: `--output-dir`, `--n-slices`, `--margin-um`, `--projection slices|max`.

Default output: `save_path/geometry/alignment_preview/<pair_label>/`.

---

## Step 7 — Register

```bash
uv run lightsuite multires register -c my_multires.yaml
```

This:

1. Streams overview and ROI overlap crops (using metadata or landmark geometry).
2. Optionally bins volumes for Elastix (`registration_bin`).
3. Runs the requested Elastix stages (default: translation then rigid).
4. Writes registered TIFFs and transform parameter files.

### Registration outputs

Under `save_path/elastix_roi_to_overview/<experiment_name>/`:

| File | Description |
|------|-------------|
| `<experiment>_<overview>_overview_crop.tif` | Fixed (overview) overlap crop |
| `<experiment>_<roi>_registered_to_<overview>.tif` | ROI warped into overview overlap space |
| `<experiment>_<roi>_registered_to_<overview>_in_full_overview.tif` | ROI embedded in full overview grid *(when `write_full_overview_canvas: true`)* |
| `TransformParameters.*.txt` | Elastix transform chain |
| `registration_overlay_qc.png` | Mid-plane overlay QC |

The checkpoint `save_path/multires_regopts.json` is updated with paths, overlap box, and optional NCC score.

> **Renamed in 2026-08:** the overview crop was previously `<...>_cropped_overlap.tif`, which read as "the overlap, cropped" rather than "the overview, cropped to the overlap". Re-run `check-geometry` / `register` to regenerate under the new name; old checkpoints still point at the old path.

---

## Step 8 — Inspect the registration in Napari

```bash
uv run lightsuite multires inspect-registration -c my_multires.yaml
```

Loads the overview and every registered ROI channel as additive layers so you can scroll through and confirm the warp. By default it uses the **overlap crop**, where the overview and ROI share one small grid. Add `--full-overview` to load the full-overview canvas instead — much larger, and only the reference channel is available there.

Volumes are memory-mapped where the TIFF layout allows, so opening a multi-GB canvas does not read it all into RAM.

| Flag | Effect |
|------|--------|
| *(default)* | Overview crop + all registered channels (`registered_roi_path`, `additional_channel_paths`) |
| `--full-overview` | Full overview volume + `registered_roi_full_overview_path` |
| `--headless` | Resolve and load layers without opening Napari (tests) |

Requires `uv sync --extra gui`.

---

## Step 9 — (Optional) Import segmentation from the ROI

Segmentation is usually run on the **high-resolution ROI**, but quantification happens on the overview (which is what gets registered to an atlas). This command carries ROI-native annotations across using the transforms from `register`:

```bash
uv run lightsuite multires import-annotations -c my_multires.yaml
```

Add an `import` block to the multires YAML using the same schema as the brain pipeline:

```yaml
import:
  write_csv: true
  annotations:
    - format: points_csv
      path: /data/segmentation/arivis_2p5x_points.csv
      label: arivis_cells
    - format: mask_tiff
      path: /data/segmentation/arivis_2p5x_mask.tif
      label: arivis_mask
```

Inputs are **1-based ROI voxel indices** (`x,y,z` CSV header) and masks on the ROI native grid. Points outside the ROI grid are dropped before warping; the count that lands inside the overview grid is reported.

### Import outputs

Written under `save_path/annotations_in_overview/`:

| File | Description |
|------|-------------|
| `<label>_in_overview.csv` | 1-based overview voxel indices (`x,y,z` + any numeric feature columns) |
| `<label>_overview_coords.npz` | `overviewptcoords` + original `roiptcoords` |
| `<label>_in_overview_crop.tif` | Warped mask on the overlap crop (uint8) |
| `<label>_in_overview.tif` | Warped mask on the full overview grid *(unless `--crop-only`)* |
| `import_annotations_summary.json` | Per-source counts and output paths |

Masks are warped with **nearest-neighbour** interpolation so labels are not blended. Use `--crop-only` to skip the full-overview canvas when the mask is sparse and the overview is large.

### Handing off to the brain pipeline

The CSV and TIFF are already in LightSuite Sample Space v1 on the overview grid, so they drop straight into the brain config with no conversion:

```yaml
# brain pipeline YAML — overview is the 'sample' here
import:
  annotations:
    - format: points_csv
      path: <save_path>/annotations_in_overview/arivis_cells_in_overview.csv
      label: arivis_cells
    - format: mask_tiff
      path: <save_path>/annotations_in_overview/arivis_mask_in_overview.tif
      label: arivis_mask
```

```bash
uv run lightsuite brain import-annotations -c brain.yaml
uv run lightsuite brain view-registration    -c brain.yaml
```

Requires `transformix` on `PATH`.

---

## End-to-end command cheat sheet

### mesoSPIM — metadata-only

```bash
uv sync --extra registration

# 1. Publish manifest (notebook or Python helper)
uv run lightsuite multires validate-config -c examples/config/multiresolution/JulieBuron_multires.yaml
uv run lightsuite multires check-geometry      -c examples/config/multiresolution/JulieBuron_multires.yaml
uv run lightsuite multires register            -c examples/config/multiresolution/JulieBuron_multires.yaml
```

### mesoSPIM — Marianna CMU (0.8× overview ↔ 2.5× ROI, dual channel)

```bash
uv sync --extra registration --extra gui

uv run lightsuite multires validate-config      -c examples/config/multiresolution/marianna_multires.yaml
uv run lightsuite multires inspect-geometry     -c examples/config/multiresolution/marianna_multires.yaml   # lateral_flip [-1, -1]
uv run lightsuite multires check-geometry       -c examples/config/multiresolution/marianna_multires.yaml
uv run lightsuite multires register             -c examples/config/multiresolution/marianna_multires.yaml
uv run lightsuite multires inspect-registration -c examples/config/multiresolution/marianna_multires.yaml
uv run lightsuite multires import-annotations   -c examples/config/multiresolution/marianna_multires.yaml
```

Atlas registration for the 0.8× overview: [`marianna_perens.yaml`](../examples/config/mesoSPIM/marianna_perens.yaml). Point its `import.annotations` at the files written under `annotations_in_overview/` to quantify the 2.5× segmentation by atlas region.

### SmartSPIM — hybrid (overview ↔ 9× cortex mosaic)

```bash
uv sync --extra registration --extra gui

uv run python scripts/rebuild_multi_res_manifests.py
uv run lightsuite multires match-points   -c examples/config/multiresolution/Multi_RES_SCANs_cortex_9x_manifest.yaml
uv run lightsuite multires check-geometry -c examples/config/multiresolution/Multi_RES_SCANs_cortex_9x_manifest.yaml
uv run lightsuite multires register       -c examples/config/multiresolution/Multi_RES_SCANs_cortex_9x_manifest.yaml
```

---

## Troubleshooting

### `Pair manifest does not exist` / volume path errors

Run `validate-config` and confirm `pair_manifest` points to the JSON you published in Step 1. Relative `volume_path` entries in the manifest are resolved from the manifest file’s directory.

### `Landmark session not found` (hybrid mode)

Run `multires match-points` first, or set `landmarks.session_path` to an existing JSON file.

### Napari fails to launch

Install GUI extras: `uv sync --extra gui`. On Linux you may need a working Qt platform plugin (`QT_QPA_PLATFORM`).

### `SimpleITK is required`

Install registration extras: `uv sync --extra registration`.

### Poor overlap / high landmark RMS

- Verify the correct overview and ROI acquisitions are paired in the manifest.
- For mesoSPIM, run `inspect-geometry` and try alternate `lateral_flip` settings before adding landmarks (see [mesoSPIM lateral_flip](#mesospim-lateral_flip-axis-convention)).
- For SmartSPIM mosaics, confirm tile origins in `metadata.txt` / `metadata.json` match the stitched stack.
- Add more landmark pairs and try `fit_mode: affine` if similarity is too rigid.
- Use `export-preview` or `check-geometry --level slice-qc` to inspect alignment before registering.

### Overview and ROI look mirrored / registration NCC stays low

- Set `multires.mesospim_geometry.overview.lateral_flip` and `.roi.lateral_flip` to the **same** `[fx, fy]` (each entry must be `+1` or `-1`).
- Use `multires inspect-geometry` and confirm physical NCC improves when toggling flips.
- Remove old `save_path/elastix_roi_to_overview/<experiment_name>/` folders after changing geometry so `multires_regopts.json` refers to the latest run.

### Registration is slow or runs out of memory

- `check-geometry` and `register` warn when the estimated overlap-crop peak RAM exceeds available machine memory.
- Increase `registration_bin` (e.g. `2`) for faster, lower-resolution Elastix.
- Use a negative `overlap_margin_um` to shrink the registration crop when the full overlap is too large.
- Lower `max_slab_bytes` to reduce streaming slab size during ROI resampling.
- Set `write_full_overview_canvas: false` when you only need the overlap-registered ROI.
- Ensure `sample.scratch` points to fast local SSD space.

---

## Related documentation

- [Installation](installation.md) — Python, `uv`, Elastix
- [Brain lightsheet usage](usage_lightsheet_brain.md) — atlas registration for whole-brain lightsheet data
