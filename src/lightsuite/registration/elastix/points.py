"""Landmark point files for elastix (-fp / -mp)."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def write_landmark_file(path: Path, points_xyz: np.ndarray) -> None:
    """Write transformix-style landmark file in physical mm coordinates."""
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    pts = np.asarray(points_xyz, dtype=float)
    if pts.size == 0:
        pts = pts.reshape(0, 3)
    if pts.ndim != 2 or pts.shape[1] != 3:
        msg = f"points_xyz must be Nx3, got {pts.shape}"
        raise ValueError(msg)

    lines = ["point", str(pts.shape[0])]
    for row in pts:
        lines.append(" ".join(f"{v:.10f}" for v in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def volume_indices_to_elastix_physical(
    points_yxz: np.ndarray,
    spacing_mm: float,
    *,
    zero_based: bool = True,
) -> np.ndarray:
    """Map volume indices (Y, X, Z) to elastix physical points (x, y, z) in mm.

    Matches ``writePointsToText.m`` (columns x, y, z) and ``write_mhd`` ITK layout.
    """
    pts = np.asarray(points_yxz, dtype=float)
    if pts.size == 0:
        return pts.reshape(0, 3)
    if not zero_based:
        pts = pts - 1.0
    spacing = float(spacing_mm)
    return np.column_stack([pts[:, 1], pts[:, 0], pts[:, 2]]) * spacing


def voxel_points_to_physical(points_1based: np.ndarray, spacing_mm: float) -> np.ndarray:
    """Backward-compatible alias; input is 0-based (Y, X, Z) unless legacy 1-based xyz."""
    return volume_indices_to_elastix_physical(points_1based, spacing_mm, zero_based=False)
