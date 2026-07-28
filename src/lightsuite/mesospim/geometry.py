"""Physical geometry helpers for mesoSPIM overview / ROI volumes."""

from __future__ import annotations

import numpy as np
import SimpleITK as sitk

from lightsuite.mesospim.config_models import MesospimGeometryConfig
from lightsuite.multires.geometry import (
    crop_to_physical_box,
    embed_crop_in_full_overview,
    overlap_physical_bounds,
    physical_bounds,
    physical_center,
    physical_corners,
    prepare_registration_pair,
    resample_to_reference_grid,
    transform_physical_points,
    transformed_bounds_in_target_space,
    voxel_count_gb,
    zyx_points_to_physical,
)

__all__ = [
    "apply_image_geometry",
    "apply_voxel_geometry",
    "crop_to_physical_box",
    "embed_crop_in_full_overview",
    "geometry_report",
    "mesospim_geometry_fields",
    "overlap_physical_bounds",
    "physical_bounds",
    "physical_center",
    "physical_corners",
    "prepare_registration_pair",
    "resample_to_reference_grid",
    "stitched_mosaic_geometry_fields",
    "transform_physical_points",
    "transformed_bounds_in_target_space",
    "voxel_count_gb",
    "voxel_geometry_report",
    "zyx_points_to_physical",
]


def apply_image_geometry(
    image: sitk.Image,
    meta: dict[str, float | int | str],
    geometry: MesospimGeometryConfig,
) -> sitk.Image:
    """Set spacing, direction, and origin from mesoSPIM meta (in-place + returned)."""
    px = float(meta["Pixelsize in um"])
    z_start = float(meta["z_start"])
    z_end = float(meta["z_end"])
    z_step = float(meta["z_stepsize"])
    sx = float(meta["x_pos"])
    sy = float(meta["y_pos"])
    n0, n1, _ = image.GetSize()

    if geometry.itk_lateral_dim0_motor == "x":
        m0, m1 = sx, sy
    else:
        m0, m1 = sy, sx

    f0, f1 = geometry.lateral_flip
    fz = -1 if z_end < z_start else 1
    direction = np.diag([float(f0), float(f1), float(fz)])
    image.SetDirection(direction.ravel().tolist())
    image.SetSpacing((px, px, abs(z_step)))

    if geometry.stage_xy_is_center:
        ic0 = (n0 - 1) / 2.0
        ic1 = (n1 - 1) / 2.0
        ox = m0 - f0 * px * ic0
        oy = m1 - f1 * px * ic1
    else:
        ox, oy = m0, m1
    image.SetOrigin((ox, oy, z_start))
    return image


def mesospim_geometry_fields(
    shape_zyx: tuple[int, int, int],
    meta: dict[str, float | int | str],
    geometry: MesospimGeometryConfig,
) -> tuple[tuple[float, float, float], tuple[float, float, float], list[float]]:
    """Return spacing, origin, and direction for a ZYX shape without allocating voxels."""
    _z, y_size, x_size = shape_zyx
    px = float(meta["Pixelsize in um"])
    z_start = float(meta["z_start"])
    z_end = float(meta["z_end"])
    z_step = float(meta["z_stepsize"])
    sx = float(meta["x_pos"])
    sy = float(meta["y_pos"])

    if geometry.itk_lateral_dim0_motor == "x":
        m0, m1 = sx, sy
    else:
        m0, m1 = sy, sx

    f0, f1 = geometry.lateral_flip
    fz = -1 if z_end < z_start else 1
    direction = np.diag([float(f0), float(f1), float(fz)]).ravel().tolist()
    spacing = (px, px, abs(z_step))

    if geometry.stage_xy_is_center:
        ox = m0 - f0 * px * (x_size - 1) / 2.0
        oy = m1 - f1 * px * (y_size - 1) / 2.0
    else:
        ox, oy = m0, m1
    origin = (ox, oy, z_start)
    return spacing, origin, direction


def stitched_mosaic_geometry_fields(
    stitched_shape_zyx: tuple[int, int, int],
    anchor_meta: dict[str, float | int | str],
    geometry: MesospimGeometryConfig,
) -> tuple[tuple[float, float, float], tuple[float, float, float], list[float]]:
    """Geometry for a Y-stitched mesoSPIM mosaic using the northern anchor tile meta."""
    tile_shape = (
        int(anchor_meta["z_planes"]),
        int(anchor_meta["y_pixels"]),
        int(anchor_meta["x_pixels"]),
    )
    spacing, tile_origin, direction = mesospim_geometry_fields(tile_shape, anchor_meta, geometry)

    _nz, _ny, nx = stitched_shape_zyx
    px = spacing[0]
    f0 = direction[0]
    sx = float(anchor_meta["x_pos"])
    sy = float(anchor_meta["y_pos"])
    m0 = sx if geometry.itk_lateral_dim0_motor == "x" else sy
    ox = m0 - f0 * px * (nx - 1) / 2.0
    origin = (ox, float(tile_origin[1]), float(tile_origin[2]))
    return spacing, origin, direction


def apply_voxel_geometry(
    image: sitk.Image,
    voxel_um: tuple[float, float, float] | list[float],
) -> sitk.Image:
    """Set isotropic-ish spacing from YAML ``voxel_um`` [x, y, z] with origin at zero."""
    vx, vy, vz = (float(v) for v in voxel_um)
    image.SetDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    image.SetSpacing((vx, vy, vz))
    image.SetOrigin((0.0, 0.0, 0.0))
    return image


def voxel_geometry_report(
    label: str,
    shape_zyx: tuple[int, int, int],
    voxel_um: tuple[float, float, float] | list[float],
) -> dict[str, object]:
    """Build a diagnostic geometry report using YAML voxel sizes only."""
    from lightsuite.mesospim.io import empty_image_from_shape

    image = empty_image_from_shape(shape_zyx)
    apply_voxel_geometry(image, voxel_um)
    pmin, pmax = physical_bounds(image)
    center = physical_center(image)
    return {
        "label": label,
        "shape_zyx": shape_zyx,
        "spacing_um": tuple(float(s) for s in image.GetSpacing()),
        "origin_um": tuple(image.GetOrigin()),
        "phys_min": pmin,
        "phys_max": pmax,
        "phys_center": center,
    }


def geometry_report(
    label: str,
    meta: dict[str, float | int | str],
    shape_zyx: tuple[int, int, int],
    geometry: MesospimGeometryConfig,
) -> dict[str, object]:
    """Build a diagnostic geometry report without loading pixel data."""
    from lightsuite.mesospim.io import empty_image_from_shape

    image = empty_image_from_shape(shape_zyx)
    image = apply_image_geometry(image, meta, geometry)
    pmin, pmax = physical_bounds(image)
    center = physical_center(image)
    extent = tuple(float(x) for x in (pmax - pmin))
    z, y, x = shape_zyx
    size_xyz = (x, y, z)
    spacing = tuple(float(s) for s in image.GetSpacing())
    extent_nm1 = tuple((n - 1) * s for n, s in zip(size_xyz, spacing, strict=True))
    extent_n = tuple(n * s for n, s in zip(size_xyz, spacing, strict=True))
    meta_center = np.array([float(meta["x_pos"]), float(meta["y_pos"]), float(meta["z_start"])])
    return {
        "label": label,
        "shape_zyx": shape_zyx,
        "shape_xyz": size_xyz,
        "spacing_um": spacing,
        "origin_um": tuple(image.GetOrigin()),
        "phys_min": pmin,
        "phys_max": pmax,
        "phys_center": center,
        "extent_um": extent,
        "extent_n_minus_1_um": extent_nm1,
        "extent_n_um": extent_n,
        "meta_xy_z0": meta_center,
        "center_minus_meta": center - meta_center,
    }
