"""Tooltip text for the GUI config form (by workflow)."""

from __future__ import annotations

_COMMON: dict[str, str] = {
    "name": "Sample identifier used in logs, export filenames, and result folders.",
    "scratch": (
        "Fast local disk for intermediate volumes and checkpoints during preprocessing "
        "and registration. Should have enough free space for downsampled working copies."
    ),
    "save_path": (
        "Output directory for checkpoints (regopts.json, transform_params.json) and "
        "registered volumes. Created automatically when stages run."
    ),
    "workers": "Parallel worker count for preprocessing and other CPU-heavy stages.",
}

_BRAIN: dict[str, str] = {
    **_COMMON,
    "tiff_type": (
        "channelperfile: one TIFF file per channel (typical whole-brain stitch). "
        "planeperfile: one folder per channel with one TIFF per slice (Terastitcher multi-channel)."
    ),
    "source_path": (
        "Folder containing stitched TIFF stacks when using channelperfile layout. "
        "Ignored when channel folders are listed below."
    ),
    "channel_folders": (
        "One folder per imaging channel (planeperfile). All channels must share the same "
        "slice dimensions and ordering."
    ),
    "voxel_um": (
        "Native sample voxel size in microns [x, y, z]. Required for preprocess; used to "
        "downsample to registration resolution."
    ),
    "atlas_source": (
        "Local NIfTI files: point atlas directory at Allen or Gubra downloads. "
        "BrainGlobe: auto-download from the registry (needs uv sync --extra atlas)."
    ),
    "atlas": (
        "Allen CCF or Perens/Gubra LSFM atlas when using local NIfTI files. "
        "Resolution defaults: Allen 10 µm, Perens 20 µm."
    ),
    "brainglobe_atlas": (
        "BrainGlobe registry atlas; native resolution is set automatically from the "
        "selected pack (Allen 10 µm, Gubra 20 µm, Princeton 20 µm, or rat 39 µm)."
    ),
    "atlas_resolution": (
        "Native atlas voxel size in µm. Allen local files are 10 µm; Perens/Gubra are 20 µm."
    ),
    "atlas_dir": (
        "Folder with atlas NIfTIs. Allen: average_template_10.nii.gz + annotation_10.nii.gz. "
        "Perens: gubra_template_olf.nii.gz + gubra_ano_olf.nii.gz."
    ),
    "channel_primary": (
        "1-based channel index for registration (usually autofluorescence / structural channel)."
    ),
    "channel_secondary": (
        "Optional second channel for dual mutual-information registration. Set to none to disable."
    ),
    "registration_resolution": (
        "Isotropic resolution (µm) for registration working volumes. Sample is downsampled "
        "to this resolution; often 20 µm even when the atlas is 10 µm."
    ),
    "bspline_spatial_scale_mm": (
        "B-spline control-point grid spacing in millimeters. Smaller values allow finer "
        "deformations but need more landmarks and smoothness tuning."
    ),
    "control_point_weight": (
        "Landmark metric weight in Elastix (0–1). Higher values pull the warp toward "
        "matched control points from match-points."
    ),
    "augment_points": (
        "When enabled, thinned auto-landmarks are added to user control points before "
        "B-spline registration."
    ),
    "dual_channel_mi_weight_primary": (
        "Relative weight of the primary registration channel in dual mutual-information "
        "registration (requires a secondary channel)."
    ),
    "dual_channel_mi_weight_secondary": (
        "Relative weight of the secondary channel in dual mutual-information registration."
    ),
    "orientation": (
        "Axis permutation as three integers, e.g. 1, -3, 2. Leave empty to use "
        "brain_orientation.txt from check-orientation."
    ),
    "canvas_mode": (
        "Adjust the Elastix working grid before registration: None keeps MATLAB parity "
        "(atlas warped to sample shape); Pad/Crop/Union reconcile sample and atlas extents."
    ),
    "import_annotations": (
        "External annotations to register after the main volume: points_csv (Imaris/Fiji "
        "CSV) or mask_tiff. Each entry needs a path; label sets the output filename stem."
    ),
    "intensity_metrics": (
        "Per-region intensity statistics written to region_stats.csv during export: "
        "median, mean, std, variance, and/or region volume. Imported points are "
        "always counted into region_stats when present."
    ),
    "detection": (
        "Built-in cell detection is not implemented in Python yet — leave disabled and use "
        "import-annotations with external points.csv / mask.tif."
    ),
}

_SPINAL: dict[str, str] = {
    **_COMMON,
    "tiff_type": (
        "channelperfile: single stitched folder per channel. "
        "planeperfile: one folder per channel with one TIFF per slice (typical Terastitcher output)."
    ),
    "source_path": (
        "Stitched TIFF folder for a single channel (channelperfile). "
        "Use channel folders below for multi-channel plane-per-file layouts."
    ),
    "channel_folders": (
        "One folder per channel (planeperfile). Each folder contains one TIFF per axial slice."
    ),
    "voxel_um": (
        "Native sample voxel size in microns [x, y, z] before straightening and registration."
    ),
    "atlas_dir": (
        "Fiederling spinal cord atlas folder with Template.tif, Annotation.tif, Segments.csv, "
        "Atlas_Regions.csv, and Hemisphere_Annotation.tif (Mendeley download)."
    ),
    "channel_primary": (
        "1-based channel index used for registration and atlas alignment after straightening."
    ),
    "registration_resolution": (
        "Isotropic resolution (µm) for registration working volumes after preprocess "
        "(typically 20 µm)."
    ),
    "control_point_weight": (
        "Landmark metric weight in Elastix (0–1). Higher values pull the warp toward "
        "matched control points from match-points."
    ),
    "import_annotations": (
        "External annotations to register after cord registration: points_csv or mask_tiff. "
        "Each entry needs a path; label sets the output filename stem."
    ),
    "segmentation_suite": (
        "Where segmentation came from. Native = list annotation layers below. Vendor suites run "
        "Convert annotations after preprocess. Imaris splits Component Name into one layer each "
        "(label is the filename prefix). Custom needs convert_to_lightsuite(source, output, *, reference)."
    ),
    "converter_source": "Raw vendor export path (JSON / CSV / XLSX / …).",
    "converter_label": (
        "Single-layer name for SmartSPIM/FIJI/Arivis, or Imaris filename prefix "
        "(components become <label>_<component>_points.csv)."
    ),
    "converter_voxel_um": (
        "Imaris/FIJI calibration [x, y, z] µm — same units as Position columns in the export "
        "(often 1,1,1 for Imaris; not always sample.voxel_um)."
    ),
    "converter_custom": (
        "Python file for suite=custom. Must define convert_to_lightsuite; output is "
        "always validated against sample_reference.json."
    ),
    "intensity_metrics": (
        "Per-region intensity statistics in region_stats.csv: median, mean, std, variance, "
        "and/or volume per segment. Imported points are always counted when present."
    ),
    "parcellate_intensities": (
        "Compute intensity region stats from exported channel volumes during export/region-stats."
    ),
}

_MULTIRES: dict[str, str] = {
    **_COMMON,
    "pair_label": (
        "Short label for this overview/ROI pair, used in manifest filenames and log messages."
    ),
    "pair_manifest": (
        "JSON manifest describing the overview and ROI image pair, motor positions, and "
        "geometry (written by multires import or conversion tools)."
    ),
    "reference_channel": (
        "Channel name (e.g. 488, 555) whose overview↔ROI transform is estimated first and "
        "applied to co-registered channels. Can also be set in the GUI shell."
    ),
    "multires_channels": (
        "Per-channel overview (low mag) and ROI (high mag) paths. Each may be a TIFF or a "
        "stitched folder. Browse offers Select file… or Select folder…. Optional metadata "
        "sidecars when not auto-discovered beside the volume."
    ),
    "multires_channel_overview": (
        "Low-mag overview volume: a single TIFF or a stitched plane-per-file folder. "
        "Browse → Select file… or Select folder…."
    ),
    "multires_channel_roi": (
        "High-mag ROI volume: a single TIFF or a stitched folder. "
        "Browse → Select file… or Select folder…."
    ),
    "multires_channel_overview_meta": (
        "Required for mesoSPIM stitched overview folders (anchor tile *_meta.txt). "
        "SmartSPIM: metadata.txt/json when not auto-discovered beside the stack. "
        "Optional for a TIFF whose sidecar sits next to the file."
    ),
    "multires_channel_roi_meta": (
        "Required for a stitched ROI folder. Optional for a TIFF whose sidecar sits "
        "next to the file, or SmartSPIM metadata.txt/json beside the stack."
    ),
    "geometry_mode": (
        "metadata: use manifest motor geometry only. hybrid: refine overview↔ROI alignment "
        "with landmarks from match-points (set fit mode below)."
    ),
    "landmark_fit_mode": (
        "Transform model for landmark-based geometry in hybrid mode: similarity, affine, "
        "or rigid."
    ),
    "overlap_margin_um": (
        "Extra margin (µm) around the overview/ROI overlap region. Use negative values "
        "to shrink the working canvas and save memory."
    ),
    "write_full_overview_canvas": (
        "When enabled, registration writes a full-overview registered canvas in addition "
        "to overlap crops (needed for some import-annotations workflows)."
    ),
    "import_annotations": (
        "External annotations to warp with the overview↔ROI transform: points_csv or "
        "mask_tiff. Each entry needs a path; label sets the output filename stem."
    ),
}

_TOOLTIPS_BY_WORKFLOW: dict[str, dict[str, str]] = {
    "brain": _BRAIN,
    "spinal": _SPINAL,
    "multires": _MULTIRES,
}


def tooltips_for_workflow(workflow: str | None) -> dict[str, str]:
    """Return tooltip strings keyed by form field id."""
    if workflow is None:
        return dict(_COMMON)
    return dict(_TOOLTIPS_BY_WORKFLOW.get(workflow.strip().lower(), _COMMON))
