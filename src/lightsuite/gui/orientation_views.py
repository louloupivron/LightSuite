"""Average-projection views for orientation checking."""

from __future__ import annotations

import numpy as np
from skimage.exposure import equalize_adapthist
from skimage.transform import resize


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


def _stack_projections(images: list[np.ndarray], gap_y: int) -> np.ndarray:
    if not images:
        return np.zeros((1, 1), dtype=np.float32)
    col_width = max(img.shape[1] for img in images)
    total_h = sum(img.shape[0] for img in images) + gap_y * (len(images) - 1)
    canvas = np.zeros((total_h, col_width), dtype=np.float32)
    y = 0
    for img in images:
        canvas[y : y + img.shape[0], : img.shape[1]] = img
        y += img.shape[0] + gap_y
    return canvas


def build_dual_panel_projections(
    atlas_projections: list[np.ndarray],
    sample_projections: list[np.ndarray],
    *,
    gap_x: int = 24,
    gap_y: int = 12,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Stack all three axis projections into left (atlas) and right (sample) panels."""
    pairs = [projection_pair(atlas_projections, sample_projections, axis) for axis in range(3)]
    atlas_panel = _stack_projections([pair[0] for pair in pairs], gap_y)
    sample_panel = _stack_projections([pair[1] for pair in pairs], gap_y)
    sample_translate_x = float(atlas_panel.shape[1] + gap_x)
    return atlas_panel, sample_panel, sample_translate_x
