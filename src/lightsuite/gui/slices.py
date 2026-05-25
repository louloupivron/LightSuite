"""Extract 2D slices from 3D volumes (volumeIdtoImage.m)."""

from __future__ import annotations

import numpy as np


def volume_index_to_image(volume: np.ndarray, chooserow: np.ndarray) -> np.ndarray:
    """Extract a 2D slice from a 3D volume given chooselist row."""
    slice_index = int(chooserow[0])
    axis = int(chooserow[1])  # 1-based MATLAB dim
    slices = [slice(None)] * 3
    slices[axis - 1] = slice_index - 1
    image = volume[tuple(slices)]
    return np.squeeze(image)
