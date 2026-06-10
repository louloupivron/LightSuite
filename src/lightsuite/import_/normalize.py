"""Map external coordinates to LightSuite native sample space."""

from __future__ import annotations

import numpy as np


def reorder_coordinates(
    points: np.ndarray,
    *,
    axis_order: str,
) -> np.ndarray:
    """Return Nx3 array in [x, y, z] column order."""
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] < 3:
        msg = f"Expected Nx3+ coordinates, got {pts.shape}"
        raise ValueError(msg)
    order = axis_order.strip().lower()
    if order == "xyz":
        return pts[:, :3].copy()
    if order == "zyx":
        return pts[:, [2, 1, 0]]
    msg = f"Unsupported axis_order {axis_order!r}"
    raise ValueError(msg)


def to_lightsuite_sample_indices(
    points_xyz: np.ndarray,
    *,
    index_base: int,
    target_size_yxz: tuple[int, int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Convert external [x,y,z] coords to 1-based LightSuite sample indices."""
    pts = np.asarray(points_xyz, dtype=np.float64)
    if index_base == 0:
        pts = pts + 1.0
    elif index_base != 1:
        msg = f"index_base must be 0 or 1, got {index_base}"
        raise ValueError(msg)

    ny, nx, nz = target_size_yxz
    out = pts.copy()
    in_bounds = (
        (out[:, 0] >= 1)
        & (out[:, 0] <= nx)
        & (out[:, 1] >= 1)
        & (out[:, 1] <= ny)
        & (out[:, 2] >= 1)
        & (out[:, 2] <= nz)
    )
    return out, in_bounds


def scale_coordinates_between_voxels(
    points_xyz: np.ndarray,
    source_voxel_um: list[float],
    target_voxel_um: list[float],
) -> np.ndarray:
    """Rescale index coordinates when source and target voxel sizes differ."""
    src = np.asarray(source_voxel_um, dtype=np.float64)
    tgt = np.asarray(target_voxel_um, dtype=np.float64)
    if src.shape != (3,) or tgt.shape != (3,):
        msg = "source_voxel_um and target_voxel_um must be length-3"
        raise ValueError(msg)
    if np.allclose(src, tgt):
        return np.asarray(points_xyz, dtype=np.float64)
    scale = src / tgt
    pts = np.asarray(points_xyz, dtype=np.float64)
    centered = pts - 1.0
    scaled = centered * scale + 1.0
    return scaled
