"""Coordinate offsets for Option A (working grid vs canonical atlas/sample space)."""

from __future__ import annotations

import numpy as np

from lightsuite.registration.content_bbox import ContentBox


def offset_xyz_yxz(offset_yxz: tuple[int, int, int] | list[int] | None) -> np.ndarray:
    """Return a (Y, X, Z) shift as a length-3 float vector."""
    if offset_yxz is None:
        return np.zeros(3, dtype=float)
    return np.asarray(offset_yxz, dtype=float).reshape(3)


def atlas_indices_trimmed_to_native(
    points_yxz: np.ndarray,
    crop_start_yxz: tuple[int, int, int] | list[int] | None,
) -> np.ndarray:
    """Map trimmed-atlas (Y, X, Z) indices to native full-atlas indices."""
    shift = offset_xyz_yxz(crop_start_yxz)
    if not np.any(shift):
        return np.asarray(points_yxz, dtype=float)
    out = np.asarray(points_yxz, dtype=float).copy()
    out[:, :3] += shift
    return out


def atlas_indices_native_to_trimmed(
    points_yxz: np.ndarray,
    crop_start_yxz: tuple[int, int, int] | list[int] | None,
) -> np.ndarray:
    shift = offset_xyz_yxz(crop_start_yxz)
    if not np.any(shift):
        return np.asarray(points_yxz, dtype=float)
    out = np.asarray(points_yxz, dtype=float).copy()
    out[:, :3] -= shift
    return out


def registration_indices_uncropped_to_cropped(
    points_yxz: np.ndarray,
    crop_start_yxz: tuple[int, int, int] | list[int] | None,
) -> np.ndarray:
    """Subtract sample content-crop origin on the unpermuted registration grid."""
    shift = offset_xyz_yxz(crop_start_yxz)
    if not np.any(shift):
        return np.asarray(points_yxz, dtype=float)
    out = np.asarray(points_yxz, dtype=float).copy()
    out[:, :3] -= shift
    return out


def registration_indices_cropped_to_uncropped(
    points_yxz: np.ndarray,
    crop_start_yxz: tuple[int, int, int] | list[int] | None,
) -> np.ndarray:
    shift = offset_xyz_yxz(crop_start_yxz)
    if not np.any(shift):
        return np.asarray(points_yxz, dtype=float)
    out = np.asarray(points_yxz, dtype=float).copy()
    out[:, :3] += shift
    return out


def affine_with_source_offset(
    tform_4x4: np.ndarray,
    source_offset_yxz: tuple[int, int, int] | list[int] | None,
) -> np.ndarray:
    """Adjust affine for warping a cropped/trimmed moving volume.

    With ``transform_points``, ``p_out = p_in @ L + t``. Trimmed indices
    ``p' = p - offset`` correspond to native ``p``, so use ``t' = t + offset @ L``.
    """
    matrix = np.asarray(tform_4x4, dtype=float).copy()
    shift = offset_xyz_yxz(source_offset_yxz)
    if not np.any(shift):
        return matrix
    linear = matrix[:3, :3]
    matrix[:3, 3] += shift @ linear
    return matrix


def native_crop_offset_yxz(
    crop_start_yxz: tuple[int, int, int],
    *,
    voxel_um: list[float],
    registres_um: float,
) -> tuple[float, float, float]:
    """Map registration crop origin to native sample (x, y, z) 0-based offsets."""
    vy, vx, vz = float(voxel_um[1]), float(voxel_um[0]), float(voxel_um[2])
    y0, x0, z0 = (int(v) for v in crop_start_yxz)
    return (
        x0 * registres_um / vx,
        y0 * registres_um / vy,
        z0 * registres_um / vz,
    )


def content_box_to_lists(box: ContentBox) -> tuple[list[int], list[int]]:
    return list(box.start_yxz), list(box.size_yxz)
