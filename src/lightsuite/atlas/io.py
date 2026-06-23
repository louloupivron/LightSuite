"""Load atlas template / annotation volumes (NIfTI or TIFF)."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import tifffile


def load_atlas_volume(path: Path) -> np.ndarray:
    """Load a 3D atlas volume from NIfTI (``.nii`` / ``.nii.gz``) or TIFF (``.tif`` / ``.tiff``)."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        msg = f"Atlas volume not found: {resolved}"
        raise FileNotFoundError(msg)

    name_lower = resolved.name.lower()
    if name_lower.endswith((".tif", ".tiff")):
        data = np.asanyarray(tifffile.imread(resolved))
    elif name_lower.endswith((".nii", ".nii.gz")):
        data = np.asanyarray(nib.load(str(resolved)).dataobj)
    else:
        msg = (
            f"Unsupported atlas volume format: {resolved}. "
            "Expected .nii, .nii.gz, .tif, or .tiff."
        )
        raise ValueError(msg)

    if data.ndim != 3:
        msg = f"Expected 3D atlas volume in {resolved}, got shape {data.shape}"
        raise ValueError(msg)
    return data
