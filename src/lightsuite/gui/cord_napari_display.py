"""Napari display helpers for spinal cord registration / import QC."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lightsuite.preprocess.cord_checkpoint import CordTransformParamsCheckpoint


def load_cord_tofliprc(save_path: Path) -> bool:
    """Return whether preprocess applied a rostrocaudal atlas flip (``tofliprc``)."""
    json_path = Path(save_path).expanduser() / "transform_params.json"
    if not json_path.is_file():
        return False
    return bool(CordTransformParamsCheckpoint.load(json_path).tofliprc)


def flip_registration_z_volume_yxz(volume: np.ndarray) -> np.ndarray:
    """Flip a (Y, X, Z) registration-grid volume along the longitudinal axis."""
    return np.flip(np.asarray(volume), axis=2)


def flip_registration_z_points_xyz(points_xyz: np.ndarray, *, nz: int) -> np.ndarray:
    """Flip 1-based ``x, y, z`` point coordinates on a registration grid of depth ``nz``."""
    pts = np.asarray(points_xyz, dtype=np.float64)
    if pts.size == 0:
        return pts.reshape(0, pts.shape[1] if pts.ndim == 2 else 3)
    out = pts.copy()
    out[:, 2] = nz - out[:, 2] + 1
    return out


def align_sample_space_for_atlas_qc(
    *,
    template: np.ndarray,
    annotation: np.ndarray,
    channels: dict[int, np.ndarray],
    point_layers: dict[str, np.ndarray],
    hemisphere: np.ndarray | None = None,
    tofliprc: bool,
) -> tuple[np.ndarray, np.ndarray, dict[int, np.ndarray], dict[str, np.ndarray], np.ndarray | None]:
    """Match sample-space Napari display to atlas-space rostrocaudal orientation.

    Atlas-space imports apply ``tofliprc`` when mapping points onto the native
    Fiederling grid. Sample-space ``regptcoords`` live on the straightened grid
    without that flip, so scrolling Z appears inverted when comparing inspect
    views. Flipping sample volumes and points here is display-only.
    """
    if not tofliprc:
        return template, annotation, channels, point_layers, hemisphere

    nz = int(annotation.shape[2])
    template_out = flip_registration_z_volume_yxz(template)
    annotation_out = flip_registration_z_volume_yxz(annotation)
    channels_out = {ich: flip_registration_z_volume_yxz(vol) for ich, vol in channels.items()}
    points_out = {
        label: flip_registration_z_points_xyz(coords, nz=nz) for label, coords in point_layers.items()
    }
    hemisphere_out = (
        flip_registration_z_volume_yxz(hemisphere) if hemisphere is not None else None
    )
    return template_out, annotation_out, channels_out, points_out, hemisphere_out


__all__ = [
    "align_sample_space_for_atlas_qc",
    "flip_registration_z_points_xyz",
    "flip_registration_z_volume_yxz",
    "load_cord_tofliprc",
]
