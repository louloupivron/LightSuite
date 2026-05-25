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


def voxel_points_to_physical(points_1based: np.ndarray, spacing_mm: float) -> np.ndarray:
    """Convert 1-based voxel indices (x,y,z) to physical mm (elastix convention)."""
    if points_1based.size == 0:
        return points_1based.reshape(0, 3)
    return (np.asarray(points_1based, dtype=float) - 1.0) * spacing_mm
