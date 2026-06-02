"""Control point utilities for registration."""

from __future__ import annotations

import numpy as np
from scipy.spatial.distance import pdist, squareform


def thin_point_list(points: np.ndarray, min_distance: float) -> np.ndarray:
    """Return boolean mask of points to keep (thinPointList.m)."""
    if points.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    dmat = squareform(pdist(points.astype(np.float32)))
    kept = np.zeros(points.shape[0], dtype=bool)
    kept[0] = True
    for idx in range(1, points.shape[0]):
        if np.all(dmat[idx, kept] >= min_distance):
            kept[idx] = True
    return kept


def subsample_point_pairs(
    atlas_pts: np.ndarray,
    sample_pts: np.ndarray,
    *,
    max_points: int,
    seed: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Uniformly subsample paired landmarks (same count/order preserved)."""
    atlas_pts = np.asarray(atlas_pts, dtype=float)
    sample_pts = np.asarray(sample_pts, dtype=float)
    if atlas_pts.shape[0] != sample_pts.shape[0]:
        msg = f"Point arrays must match length, got {atlas_pts.shape[0]} vs {sample_pts.shape[0]}"
        raise ValueError(msg)
    n = atlas_pts.shape[0]
    if n <= max_points:
        return atlas_pts, sample_pts
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(n, size=max_points, replace=False))
    return atlas_pts[idx], sample_pts[idx]
