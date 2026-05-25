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
