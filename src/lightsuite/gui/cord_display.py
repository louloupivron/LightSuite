"""Shared display helpers for spinal cord Napari GUIs."""

from __future__ import annotations

import numpy as np
from skimage.exposure import equalize_adapthist

CORD_ATLAS_PROVIDER = "cord"


def atlas_cut_axis_size(atlas_shape: tuple[int, int, int], chooserow: np.ndarray) -> int:
    cut_axis = int(chooserow[1]) - 1
    return int(atlas_shape[cut_axis])


def normalize_cord_display(image: np.ndarray) -> np.ndarray:
    data = image.astype(np.float32)
    if data.max() <= 0:
        return data
    hi = float(np.quantile(data, 0.999))
    data = np.clip(data / max(hi, 1e-6), 0, 1)
    if data.ndim == 2:
        try:
            data = equalize_adapthist(data, clip_limit=0.01)
        except ValueError:
            pass
    return data
