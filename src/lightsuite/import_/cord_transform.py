"""Transform native sample-space points into Fiederling atlas coordinates."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lightsuite.gui.affine import transform_points
from lightsuite.import_.transform import interpolate_bspline_displacement
from lightsuite.preprocess.cord_checkpoint import CordTransformParamsCheckpoint
from lightsuite.registration.elastix.runner import (
    run_transformix_deformation_field,
    run_transformix_points,
)
from lightsuite.registration.straightening import transform_cord_points_slices


def native_xyz_to_cord_registration_yxz(
    points_xyz_1based: np.ndarray,
    *,
    sampleres_um: list[float],
    registrationres_um: float,
) -> np.ndarray:
    """Map native sample ``x,y,z`` (1-based) to registration-grid ``y,x,z`` (0-based)."""
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    vx, vy, vz = (float(v) for v in sampleres_um)
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    reg_y = (y - 1.0) * vy / float(registrationres_um)
    reg_x = (x - 1.0) * vx / float(registrationres_um)
    reg_z = (z - 1.0) * vz / float(registrationres_um)
    return np.column_stack([reg_y, reg_x, reg_z])


def permute_cord_registration_indices_yxz(
    points_yxz: np.ndarray,
    perm_1based: list[int],
) -> np.ndarray:
    """Apply the same axis reordering as ``np.transpose(finvol, perm)``."""
    perm_order = [int(p) - 1 for p in perm_1based]
    pts = np.asarray(points_yxz, dtype=np.float64)
    return pts[:, perm_order]


def crop_cord_registration_points(
    points_yxz: np.ndarray,
    *,
    yrange: list[int],
    xrange: list[int],
    zrange: list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Keep points inside preprocess crop ranges; return cropped ``y,x,z`` (0-based)."""
    y0, y1 = int(yrange[0]) - 1, int(yrange[1]) - 1
    x0, x1 = int(xrange[0]) - 1, int(xrange[1]) - 1
    z0, z1 = int(zrange[0]) - 1, int(zrange[1]) - 1
    pts = np.asarray(points_yxz, dtype=np.float64)
    keep = (
        (pts[:, 0] >= y0)
        & (pts[:, 0] <= y1)
        & (pts[:, 1] >= x0)
        & (pts[:, 1] <= x1)
        & (pts[:, 2] >= z0)
        & (pts[:, 2] <= z1)
    )
    cropped = pts[keep].copy()
    cropped[:, 0] -= y0
    cropped[:, 1] -= x0
    cropped[:, 2] -= z0
    return cropped, keep


def registration_yxz_to_cord_cloud_xyz(points_yxz_0based: np.ndarray) -> np.ndarray:
    """Convert registration ``y,x,z`` (0-based) to cord straightening ``x,y,z`` (1-based)."""
    pts = np.asarray(points_yxz_0based, dtype=np.float64)
    return np.column_stack([pts[:, 1] + 1.0, pts[:, 0] + 1.0, pts[:, 2] + 1.0])


def cord_cloud_xyz_to_registration_yxz(points_xyz_1based: np.ndarray) -> np.ndarray:
    """Inverse of :func:`registration_yxz_to_cord_cloud_xyz`."""
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    return np.column_stack([pts[:, 1] - 1.0, pts[:, 0] - 1.0, pts[:, 2] - 1.0])


def registration_atlas_yxz_to_native_cloud_xyz(
    points_yxz: np.ndarray,
    *,
    atlassize: tuple[int, int, int],
    template_native_shape: tuple[int, int, int],
) -> np.ndarray:
    """Upsample registration-grid atlas ``y,x,z`` (0-based) to native Fiederling ``x,y,z`` (1-based)."""
    pts = np.asarray(points_yxz, dtype=np.float64)
    ay, ax, az = (int(v) for v in atlassize)
    nz, ny, nx = (int(v) for v in template_native_shape)

    def _scale(value: np.ndarray, src: int, dst: int) -> np.ndarray:
        if src <= 1:
            return np.zeros_like(value)
        return value * (dst - 1) / (src - 1)

    nat_z = _scale(pts[:, 2], az, nz)
    nat_y = _scale(pts[:, 0], ay, ny)
    nat_x = _scale(pts[:, 1], ax, nx)
    return np.column_stack([nat_x + 1.0, nat_y + 1.0, nat_z + 1.0])


def filter_points_for_cord_registration(
    points_xyz_1based: np.ndarray,
    *,
    transform_params: CordTransformParamsCheckpoint,
) -> tuple[np.ndarray, np.ndarray]:
    """Downsample, permute, and crop-filter native points. Returns filtered coords and keep mask."""
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    if pts.size == 0:
        return pts.reshape(0, pts.shape[1] if pts.ndim == 2 else 3), np.zeros(0, dtype=bool)

    reg_um = float(transform_params.registrationres_um[0])
    reg_yxz = native_xyz_to_cord_registration_yxz(
        pts,
        sampleres_um=transform_params.sampleres_um,
        registrationres_um=reg_um,
    )
    reg_yxz = permute_cord_registration_indices_yxz(reg_yxz, transform_params.how_to_perm)
    _, keep = crop_cord_registration_points(
        reg_yxz,
        yrange=transform_params.samp_ikeepy,
        xrange=transform_params.samp_ikeepx,
        zrange=transform_params.samp_ikeeplong,
    )
    return pts[keep], keep


def sample_points_to_straightened_grid(
    points_xyz_1based: np.ndarray,
    *,
    transform_params: CordTransformParamsCheckpoint,
    tforms: list[np.ndarray],
) -> np.ndarray:
    """Map native sample points to the straightened registration grid (``y,x,z`` 0-based)."""
    reg_um = float(transform_params.registrationres_um[0])
    reg_yxz = native_xyz_to_cord_registration_yxz(
        points_xyz_1based,
        sampleres_um=transform_params.sampleres_um,
        registrationres_um=reg_um,
    )
    reg_yxz = permute_cord_registration_indices_yxz(reg_yxz, transform_params.how_to_perm)
    reg_yxz, _ = crop_cord_registration_points(
        reg_yxz,
        yrange=transform_params.samp_ikeepy,
        xrange=transform_params.samp_ikeepx,
        zrange=transform_params.samp_ikeeplong,
    )
    if reg_yxz.size == 0:
        return reg_yxz.reshape(0, 3)

    cloud = registration_yxz_to_cord_cloud_xyz(reg_yxz)
    straight_cloud = transform_cord_points_slices(cloud, tforms)
    return cord_cloud_xyz_to_registration_yxz(straight_cloud)


def transform_points_to_cord_atlas(
    points_xyz_1based: np.ndarray,
    *,
    transform_params: CordTransformParamsCheckpoint,
    elastix_affine_path: Path,
    transinit: np.ndarray,
    tforms: list[np.ndarray],
    spacing_mm: float,
    template_native_shape: tuple[int, int, int],
    temp_dir: Path,
) -> np.ndarray:
    """Transform native sample-space points to 1-based Fiederling atlas ``x,y,z`` indices.

    Chain: straighten → B-spline displacement → transformix(forward elastix) →
    ``inv(transinit)`` → optional rostrocaudal flip → native upsample.

    Volume export resamples with the *inverse* elastix file (``output(x)=input(T(x))``).
    Feature points must be pushed with the *forward* elastix transform so they land
    where the corresponding voxels appear. The composed 4×4 affine approximation
    introduces a Y–Z shear that tilts point clouds relative to the registered cord.
    """
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    if pts.size == 0:
        return pts.reshape(0, 3)

    straight_yxz = sample_points_to_straightened_grid(
        pts,
        transform_params=transform_params,
        tforms=tforms,
    )
    if straight_yxz.size == 0:
        return straight_yxz.reshape(0, 3)

    bspline_path = Path(transform_params.tform_bspline_samp20um_to_atlas_20um_px)
    displacement = run_transformix_deformation_field(
        transform_path=bspline_path,
        output_dir=temp_dir / "deformation_field",
    )
    regsize_mm = float(transform_params.registrationres_um[0]) * 1e-3
    displacement_vox = displacement / regsize_mm
    warped_yxz = straight_yxz + interpolate_bspline_displacement(straight_yxz, displacement_vox)

    after_affine = run_transformix_points(
        points_yxz=warped_yxz,
        transform_path=Path(elastix_affine_path),
        output_dir=temp_dir / "elastix_affine_points",
        spacing_mm=spacing_mm,
    )
    atlas_yxz = transform_points(after_affine, np.linalg.inv(np.asarray(transinit, dtype=float)))

    atlassize = tuple(int(v) for v in transform_params.atlassize)
    if transform_params.tofliprc:
        atlas_yxz = atlas_yxz.copy()
        atlas_yxz[:, 2] = (atlassize[2] - 1) - atlas_yxz[:, 2]

    atlas_xyz = registration_atlas_yxz_to_native_cloud_xyz(
        atlas_yxz,
        atlassize=atlassize,
        template_native_shape=template_native_shape,
    )

    if pts.shape[1] > 3:
        return np.column_stack([atlas_xyz, pts[:, 3:]])
    return atlas_xyz
