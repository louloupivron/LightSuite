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
        "Per-channel overview (low mag) and ROI (high mag) image paths. Channel names "
        "become YAML keys under multires.channels."
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
