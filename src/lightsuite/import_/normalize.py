"""Native sample-space coordinate validation."""

from __future__ import annotations

import numpy as np


def filter_in_bounds_points(
    points_xyz: np.ndarray,
    *,
    target_size_yxz: tuple[int, int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Keep 1-based [x,y,z] points inside native sample bounds."""
    pts = np.asarray(points_xyz, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] < 3:
        msg = f"Expected Nx3+ coordinates, got {pts.shape}"
        raise ValueError(msg)

    ny, nx, nz = target_size_yxz
    in_bounds = (
        (pts[:, 0] >= 1)
        & (pts[:, 0] <= nx)
        & (pts[:, 1] >= 1)
        & (pts[:, 1] <= ny)
        & (pts[:, 2] >= 1)
        & (pts[:, 2] <= nz)
    )
    return pts, in_bounds
