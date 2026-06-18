"""Transform sample-space annotations into atlas coordinates."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import interpn

from lightsuite.gui.affine import transform_points
from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.elastix.runner import run_transformix_deformation_field
from lightsuite.registration.points import volume_indices_to_cloud_xyz


def _unpermuted_registration_shape(
    permuted_shape: tuple[int, int, int],
    permute: list[int],
) -> tuple[int, int, int]:
    """Invert ``permute_brain_volume`` shape mapping (Y, X, Z)."""
    perm_order = [abs(int(v)) - 1 for v in permute]
    orig = [0, 0, 0]
    for new_axis, old_axis in enumerate(perm_order):
        orig[old_axis] = int(permuted_shape[new_axis])
    return tuple(orig)


def native_xyz_to_unpermuted_registration_yxz(
    points_xyz_1based: np.ndarray,
    *,
    ori_voxel_um: list[float],
    registres_um: float,
) -> np.ndarray:
    """Map native sample ``x,y,z`` (1-based) to unpermuted registration ``y,x,z`` (0-based)."""
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    vx, vy, vz = (float(v) for v in ori_voxel_um)
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    reg_y = (y - 1.0) * vy / float(registres_um)
    reg_x = (x - 1.0) * vx / float(registres_um)
    reg_z = (z - 1.0) * vz / float(registres_um)
    return np.column_stack([reg_y, reg_x, reg_z])


def permute_registration_indices_yxz(
    points_yxz: np.ndarray,
    shape_yxz: tuple[int, int, int],
    permute: list[int],
) -> np.ndarray:
    """Apply the same index remapping as :func:`permute_brain_volume` to ``y,x,z`` coords."""
    pts = np.asarray(points_yxz, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] != 3:
        msg = f"Expected Nx3 points, got {pts.shape}"
        raise ValueError(msg)

    perm_order = [abs(int(v)) - 1 for v in permute]
    comps = pts[:, [0, 1, 2]]
    remapped = comps[:, perm_order]
    out = remapped.copy()
    for dim, val in enumerate(permute):
        if int(val) < 0:
            axis_size = shape_yxz[perm_order[dim]]
            out[:, dim] = (axis_size - 1.0) - out[:, dim]
    return out


def sample_points_to_registration_voxels(
    points_xyz_1based: np.ndarray,
    transform_params: TransformParamsCheckpoint,
    *,
    registres_um: float,
) -> np.ndarray:
    """Map native sample ``x,y,z`` (1-based) to permuted registration ``y,x,z`` (0-based).

    Matches the registration / export volume grid: downsample native indices, then apply
    ``permute_brain_volume`` before B-spline warping.
    """
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] < 3:
        msg = f"Expected Nx3+ points, got {pts.shape}"
        raise ValueError(msg)

    permuted_shape = tuple(int(v) for v in transform_params.regvolsize)
    unperm_shape = _unpermuted_registration_shape(permuted_shape, transform_params.permute_sample_to_atlas)
    reg_yxz = native_xyz_to_unpermuted_registration_yxz(
        pts,
        ori_voxel_um=transform_params.ori_voxel_um,
        registres_um=registres_um,
    )
    pad = transform_params.warp_canvas_pad_before
    if pad is not None and any(int(p) for p in pad):
        reg_yxz = reg_yxz.copy()
        reg_yxz[:, 0] += float(pad[0])
        reg_yxz[:, 1] += float(pad[1])
        reg_yxz[:, 2] += float(pad[2])
    return permute_registration_indices_yxz(
        reg_yxz,
        unperm_shape,
        transform_params.permute_sample_to_atlas,
    )


def _displacement_field_registration_voxels(
    transform_params: TransformParamsCheckpoint,
    *,
    cache_dir: Path,
    registres_um: float,
) -> np.ndarray:
    """Load B-spline displacement field scaled to registration voxel units (Y,X,Z,3)."""
    del registres_um
    regsize_mm = float(transform_params.atlas_resolution_um) * 2.0 * 1e-3
    field_mm = run_transformix_deformation_field(
        transform_path=Path(transform_params.tform_bspline_samp20um_to_atlas_20um_px),
        output_dir=cache_dir,
    )
    return field_mm / regsize_mm


def interpolate_bspline_displacement(
    reg_points_yxz: np.ndarray,
    displacement_field_yxz: np.ndarray,
) -> np.ndarray:
    """Trilinear displacement sampling on the permuted registration grid (Y, X, Z)."""
    pts = np.asarray(reg_points_yxz, dtype=np.float64)
    field = np.asarray(displacement_field_yxz, dtype=np.float64)
    sy, sx, sz, _ = field.shape
    grid_y = np.arange(1, sy + 1)
    grid_x = np.arange(1, sx + 1)
    grid_z = np.arange(1, sz + 1)
    dx = interpn(
        (grid_y, grid_x, grid_z),
        field[:, :, :, 0],
        pts,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    dy = interpn(
        (grid_y, grid_x, grid_z),
        field[:, :, :, 1],
        pts,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    dz = interpn(
        (grid_y, grid_x, grid_z),
        field[:, :, :, 2],
        pts,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    disp = -np.column_stack([dx, dy, dz])
    disp[np.isnan(disp)] = 0.0
    return disp


def transform_points_to_atlas(
    points_xyz_1based: np.ndarray,
    transform_params: TransformParamsCheckpoint,
    *,
    registres_um: float,
    temp_dir: Path,
) -> np.ndarray:
    """Transform sample-space cell coordinates to 1-based atlas ``x,y,z`` voxel indices."""
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    reg_yxz = sample_points_to_registration_voxels(
        pts,
        transform_params,
        registres_um=registres_um,
    )
    displacement = _displacement_field_registration_voxels(
        transform_params,
        cache_dir=temp_dir / "deformation_field",
        registres_um=registres_um,
    )
    warped_yxz = reg_yxz + interpolate_bspline_displacement(reg_yxz, displacement)
    affine = np.asarray(transform_params.tform_affine_samp20um_to_atlas_10um_px, dtype=np.float64)
    atlas_yxz = transform_points(warped_yxz, affine)
    atlas_xyz = volume_indices_to_cloud_xyz(atlas_yxz)

    if pts.shape[1] > 3:
        return np.column_stack([atlas_xyz, pts[:, 3:]])
    return atlas_xyz
