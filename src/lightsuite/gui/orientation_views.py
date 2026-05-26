"""Average-projection views for orientation checking."""

from __future__ import annotations

import numpy as np
from skimage.exposure import equalize_adapthist


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


def projection_layout(
    atlas_projections: list[np.ndarray],
    sample_projections: list[np.ndarray],
) -> tuple[list[tuple[np.ndarray, tuple[float, float], tuple[float, float]]], float]:
    """Compute napari layer data, translate, and scale for side-by-side projections."""
    atlas_norm = [normalize_projection(p) for p in atlas_projections]
    sample_norm = [normalize_projection(p) for p in sample_projections]

    row_gap = float(max(p.shape[0] for p in atlas_norm + sample_norm))
    x_atlas = 0.0
    atlas_layers: list[tuple[np.ndarray, tuple[float, float], tuple[float, float]]] = []
    for idx, proj in enumerate(atlas_norm):
        translate = (x_atlas, 0.0)
        atlas_layers.append((proj, translate, (1.0, 1.0)))
        x_atlas += proj.shape[1]

    sample_layers: list[tuple[np.ndarray, tuple[float, float], tuple[float, float]]] = []
    x_sample = row_gap
    for idx, (proj, ref) in enumerate(zip(sample_norm, atlas_norm, strict=True)):
        scale_y = ref.shape[0] / max(proj.shape[0], 1)
        scale_x = ref.shape[1] / max(proj.shape[1], 1)
        translate = (x_sample, 0.0)
        sample_layers.append((proj, translate, (scale_y, scale_x)))
        x_sample += ref.shape[1]

    return atlas_layers + sample_layers, row_gap
