# Registration output spaces

LightSuite supports two **complementary** registration output directions:

| Direction | Transform | Typical use |
|-----------|-----------|-------------|
| **Atlas space** (`sample → atlas`) | Warp sample channels into the atlas grid | Cross-subject comparison on a fixed CCF grid |
| **Sample space** (`atlas → sample`) | Warp atlas labels/template onto the registration grid | Native-grid QC, label transfer without resampling huge volumes, brainreg-style inspection |

Both directions use the **same** elastix transform chain; only the resampling direction differs.

## Grids

| Pipeline | Atlas-space grid | Sample-space grid |
|----------|------------------|-------------------|
| Brain | Atlas native resolution (e.g. 10 µm Allen) | Permuted **20 µm registration grid** (`regvolsize` in `transform_params.json`) |
| Spinal cord | Fiederling native (10×10×20 µm export layout) | **Straightened 20 µm** grid (same as `qc/registration_*.png`) |

Sample-space outputs do **not** upsample to the full native stitched volume. Segmentation should still run on the native grid (`sample_reference.json`); import writes both atlas- and sample-space coordinates.

Preprocess writes `chan_*_sample_register_*um.tif` on the **unpermuted** downsampled grid. Registration, import (`regptcoords`, `*_in_sample_20um.tif`), and sample-space atlas exports use the **permuted** grid from `registration.orientation` / `permute_sample_to_atlas`. Tools that overlay channels with imported annotations apply that permutation when loading registration TIFFs. Sample-space atlas TIFFs written by `brain export --space sample` are stored on the permuted grid (`atlas_volumes_permuted: true` in `sample_space_manifest.json`); older exports without that flag are permuted at load time for Napari QC.

## Configuration

```yaml
export:
  spaces: [atlas]              # or [sample] or [atlas, sample]
  save_registered_volume: true # atlas-space channel TIFFs
  save_sample_space_volume: true # warped labels under sample_space/

analysis:
  stats_spaces: [atlas]          # or [sample] or [atlas, sample]
```

CLI override: `--space atlas|sample|both` on `brain export`, `spinal export`, and `spinal region-stats`.

## Output layout (brain)

```
volume_registered/
├── chan_01_registered_atlas.tif      # atlas space
├── chan01_region_stats.csv           # atlas-space intensities
├── region_stats.csv                  # atlas-space combined stats
├── region_stats_top10.csv            # wide top-N regions per channel
├── sample_space/
│   ├── annotation_in_sample_20um.tif
│   ├── template_in_sample_20um.tif
│   ├── division_labels_in_sample_20um.tif
│   ├── sample_space_manifest.json
│   ├── chan01_region_stats_sample.csv
│   ├── region_stats_sample.csv
│   └── region_stats_sample_top10.csv
├── {label}_atlas_coords.npz
├── {label}_sample_coords.npz         # key: regptcoords
└── {label}_in_sample_20um.tif        # masks only
```

## Output layout (spinal cord)

```
volume_registered/
├── chan01_channel1.tiff              # atlas native export
├── annotation_registered.tiff
├── sample_space/
│   ├── annotation_in_sample_20um.tif
│   ├── template_in_sample_20um.tif
│   ├── segments_in_sample_20um.tif
│   ├── hemisphere_in_sample_20um.tif
│   ├── chan_01_sample_straight_20um.tif
│   ├── sample_space_manifest.json
│   └── region_stats_sample.csv
├── {label}_atlas_coords.npz
└── {label}_sample_coords.npz
```

## Comparing statistics across spaces

**Do not** directly compare `region_stats.csv` and `region_stats_sample.csv`:

- **Voxel size** differs (atlas µm vs 20 µm registration).
- **Region volumes** in sample space reflect **warped** label footprints (deformation), not canonical atlas volumes.
- **Cell densities** use different denominators.

Use atlas space for cross-subject comparison; use sample space for QC and native-resolution workflows that bin into warped labels.

Sample-space hemisphere splits use ``hemisphere_in_sample_20um.tif`` (warped like the
annotation). The atlas→sample export mirrors the lateral axis, so left/right assignment
applies the complementary ``flip`` to stay consistent with atlas-space stats.

## Re-running registration

Sample-space export requires forward transform files written at `register`:

- Brain: `bspline_atlas_to_samp_20um.txt` + `tform_affine_atlas_to_samp20um_px` in `transform_params.json`
- Spinal: `transforms/bspline_atlas_to_samp_20um.txt`

Existing samples need a **re-run of `register`** (not full preprocess) before first sample-space export.
