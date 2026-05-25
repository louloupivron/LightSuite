# Widefield slice analysis

The **Python pipeline does not yet include the coronal slice module.** Use the MATLAB workflow described below until this module is ported.

For the Python brain pipeline, see [Brain lightsheet analysis](usage_lightsheet_brain.md).

---

## MATLAB workflow (current)

The slice module registers individual 2D coronal (or arbitrary) brain sections to the Allen atlas. It is optimized for widefield microscopy data rather than full 3D lightsheet volumes.

### Supported inputs (MATLAB)

- Series of 2D TIFF files (one file per slice)
- AxioScan `.czi` exports

### Getting started

1. Install LightSuite for MATLAB (see repository `README.md`)
2. Download Allen CCF atlas files
3. Run `demos/ls_analyze_slice_volume.m`

### Workflow highlights

- Per-slice registration with manual control-point refinement recommended
- B-spline registration per slice
- Cell detection and atlas mapping analogous to the brain module

### Python status

Slice registration and the Napari-based slice GUI are not implemented in Python yet. Planned as a separate pipeline under `lightsuite slice` (name TBD) after brain export and cell detection are complete.
