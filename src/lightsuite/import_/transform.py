"""Transform sample-space annotations into atlas coordinates."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import interpn

from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.elastix.runner import run_transformix_deformation_field
from lightsuite.registration.warp import transform_points_affinetform


def sample_points_to_registration_voxels(
    points_xyz_1based: np.ndarray,
    transform_params: TransformParamsCheckpoint,
    *,
    registres_um: float,
) -> np.ndarray:
    """Map native sample [x,y,z] (1-based) to registration-resolution voxel indices.

    Port of ``coreTransform`` in transformPointsToAtlas.m (steps 1–4).
    """
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[1] < 3:
        msg = f"Expected Nx3+ points, got {pts.shape}"
        raise ValueError(msg)

    ori_px = np.asarray(transform_params.ori_voxel_um, dtype=np.float64)
    ori_size = np.asarray(transform_params.ori_size, dtype=np.float64)
    permute = transform_params.permute_sample_to_atlas
    regsize_mm = float(transform_params.atlas_resolution_um) * 2.0 * 1e-3

    coords = (pts[:, :3] - 1.0) * ori_px * 1e-3
    coords = coords[:, [1, 0, 2]]

    phys_size_yxz = (ori_size - 1.0) * ori_px[[1, 0, 2]] * 1e-3
    perm_order = [abs(int(p)) for p in permute]
    permuted = np.zeros_like(coords)
    for dim, pval in enumerate(permute):
        orig_dim = perm_order[dim] - 1
        permuted[:, dim] = coords[:, orig_dim]
        if int(pval) < 0:
            permuted[:, dim] = phys_size_yxz[orig_dim] - permuted[:, dim]

    reg_pts = permuted[:, [1, 0, 2]] / regsize_mm
    return reg_pts


def _displacement_field_registration_voxels(
    transform_params: TransformParamsCheckpoint,
    *,
    cache_dir: Path,
    registres_um: float,
) -> np.ndarray:
    """Load B-spline displacement field scaled to registration voxel units (Y,X,Z,3)."""
    regsize_mm = float(transform_params.atlas_resolution_um) * 2.0 * 1e-3
    field_mm = run_transformix_deformation_field(
        transform_path=Path(transform_params.tform_bspline_samp20um_to_atlas_20um_px),
        output_dir=cache_dir,
    )
    # MATLAB: permute(Dfield,[2 3 4 1]) / regsize_mm with Dfield x,y,z,components
    # Our field is Y,X,Z,components in mm; swap to x,y,z ordering then scale.
    field_xyz = field_mm[:, :, :, [1, 0, 2]] / regsize_mm
    return field_xyz


def interpolate_bspline_displacement(
    reg_points_xyz: np.ndarray,
    displacement_field_xyz: np.ndarray,
) -> np.ndarray:
    """Trilinear displacement sampling (transformPointsToAtlas.m interpn)."""
    pts = np.asarray(reg_points_xyz, dtype=np.float64)
    field = np.asarray(displacement_field_xyz, dtype=np.float64)
    sx, sy, sz, _ = field.shape
    grid_x = np.arange(1, sx + 1)
    grid_y = np.arange(1, sy + 1)
    grid_z = np.arange(1, sz + 1)
    dx = interpn(
        (grid_x, grid_y, grid_z),
        field[:, :, :, 0],
        pts,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    dy = interpn(
        (grid_x, grid_y, grid_z),
        field[:, :, :, 1],
        pts,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    dz = interpn(
        (grid_x, grid_y, grid_z),
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
    """Transform sample-space cell coordinates to atlas voxel indices."""
    pts = np.asarray(points_xyz_1based, dtype=np.float64)
    reg_pts = sample_points_to_registration_voxels(
        pts,
        transform_params,
        registres_um=registres_um,
    )
    pad = transform_params.warp_canvas_pad_before
    if pad is not None and any(int(p) for p in pad):
        reg_pts = reg_pts.copy()
        # B-spline displacement field is sampled on the padded registration grid.
        reg_pts[:, 0] += pad[1]
        reg_pts[:, 1] += pad[0]
        reg_pts[:, 2] += pad[2]
    displacement = _displacement_field_registration_voxels(
        transform_params,
        cache_dir=temp_dir / "deformation_field",
        registres_um=registres_um,
    )
    warped = reg_pts + interpolate_bspline_displacement(reg_pts, displacement)
    affine = np.asarray(transform_params.tform_affine_samp20um_to_atlas_10um_px, dtype=np.float64)
    atlas_pts = transform_points_affinetform(warped, affine)

    if pts.shape[1] > 3:
        return np.column_stack([atlas_pts, pts[:, 3:]])
    return atlas_pts
