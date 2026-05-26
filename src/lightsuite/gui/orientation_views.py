"""Average-projection views for orientation checking."""

from __future__ import annotations

import numpy as np
from skimage.exposure import equalize_adapthist
from skimage.transform import resize

PROJECTION_AXIS_LABELS = [
    "Projection 0 — mean over axis 0 (Y)",
    "Projection 1 — mean over axis 1 (X)",
    "Projection 2 — mean over axis 2 (Z)",
]


def normalize_projection(image: np.ndarray) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    if data.size == 0:
        return data
    if data.max() <= 0:
        return np.zeros_like(data)
    hi = float(np.quantile(data, 0.999))
    out = np.clip(data / max(hi, 1e-6), 0, 1)
    try:
        out = equalize_adapthist(out, clip_limit=0.01)
    except ValueError:
        pass
    return out


def mean_projections(volume: np.ndarray) -> list[np.ndarray]:
    """Return mean projections along axes 0, 1, 2 of a (Y, X, Z) volume."""
    if volume.ndim != 3:
        msg = f"Expected 3D volume, got shape {volume.shape}"
        raise ValueError(msg)
    return [
        np.mean(volume, axis=0),
        np.mean(volume, axis=1),
        np.mean(volume, axis=2),
    ]


def projection_pair(
    atlas_projections: list[np.ndarray],
    sample_projections: list[np.ndarray],
    axis: int,
    *,
    match_height: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Return normalized atlas and sample projections for one axis at full sample resolution."""
    if axis not in {0, 1, 2}:
        msg = f"axis must be 0, 1, or 2, got {axis}"
        raise ValueError(msg)

    atlas = normalize_projection(atlas_projections[axis])
    sample = normalize_projection(sample_projections[axis])

    if match_height and atlas.shape != sample.shape:
        # Upscale atlas to sample grid for comparison; never downscale the sample.
        atlas = resize(
            atlas,
            sample.shape,
            order=1,
            preserve_range=True,
            anti_aliasing=True,
        ).astype(np.float32)

    return atlas, sample


def dual_panel_translate(atlas_shape: tuple[int, int], gap: int = 24) -> tuple[float, float]:
    """Place the sample panel to the right of the atlas panel."""
    return float(atlas_shape[1] + gap), 0.0
