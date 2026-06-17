"""Physical geometry helpers for mesoSPIM overview / ROI volumes."""

from __future__ import annotations

import numpy as np
import SimpleITK as sitk

from lightsuite.mesospim.config_models import MesospimGeometryConfig


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


def zyx_points_to_physical(image: sitk.Image, points_zyx: np.ndarray) -> np.ndarray:
    """Convert remapped ZYX indices to physical coordinates (µm)."""
    if points_zyx.size == 0:
        return points_zyx.reshape(0, 3)
    if points_zyx.ndim != 2 or points_zyx.shape[1] != 3:
        msg = f"Expected Nx3 ZYX points, got shape {points_zyx.shape}"
        raise ValueError(msg)
    physical = []
    for z, y, x in points_zyx:
        index = (float(x), float(y), float(z))
        physical.append(image.TransformContinuousIndexToPhysicalPoint(index))
    return np.asarray(physical, dtype=float)


def physical_corners(pmin: np.ndarray, pmax: np.ndarray) -> np.ndarray:
    """Return 8 corners of an axis-aligned physical box as Nx3."""
    corners = []
    for x in (pmin[0], pmax[0]):
        for y in (pmin[1], pmax[1]):
            for z in (pmin[2], pmax[2]):
                corners.append((float(x), float(y), float(z)))
    return np.asarray(corners, dtype=float)


def transform_physical_points(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    if points.size == 0:
        return points.reshape(0, 3)
    hom = np.column_stack([points, np.ones(points.shape[0])])
    return (hom @ matrix.T)[:, :3]


def transformed_bounds_in_target_space(
    image: sitk.Image,
    source_to_target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map an image's physical bounding box into another volume's coordinate frame."""
    pmin, pmax = physical_bounds(image)
    transformed = transform_physical_points(physical_corners(pmin, pmax), source_to_target)
    return transformed.min(axis=0), transformed.max(axis=0)


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


def physical_bounds(image: sitk.Image) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned physical bounding box (µm), accounting for direction cosines."""
    size = image.GetSize()
    corners_phys = [
        image.TransformIndexToPhysicalPoint((i, j, k))
        for i in (0, size[0] - 1)
        for j in (0, size[1] - 1)
        for k in (0, size[2] - 1)
    ]
    corners_phys_arr = np.asarray(corners_phys, dtype=float)
    return corners_phys_arr.min(axis=0), corners_phys_arr.max(axis=0)


def physical_center(image: sitk.Image) -> np.ndarray:
    size = image.GetSize()
    idx = [(s - 1) / 2.0 for s in size]
    return np.array(image.TransformContinuousIndexToPhysicalPoint(idx), dtype=float)


def overlap_physical_bounds(
    image_a: sitk.Image,
    image_b: sitk.Image,
    *,
    margin_um: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Intersection of two volumes' physical bounding boxes."""
    min_a, max_a = physical_bounds(image_a)
    min_b, max_b = physical_bounds(image_b)
    overlap_min = np.maximum(min_a, min_b) - margin_um
    overlap_max = np.minimum(max_a, max_b) + margin_um
    if np.any(overlap_min >= overlap_max):
        msg = (
            f"No physical overlap between volumes.\n"
            f"  A: {min_a} .. {max_a}\n"
            f"  B: {min_b} .. {max_b}"
        )
        raise ValueError(msg)
    return overlap_min, overlap_max


def crop_to_physical_box(
    image: sitk.Image,
    phys_min: np.ndarray,
    phys_max: np.ndarray,
) -> tuple[sitk.Image, list[int]]:
    """Extract the index region that covers an axis-aligned physical box."""
    corners = []
    for x in (phys_min[0], phys_max[0]):
        for y in (phys_min[1], phys_max[1]):
            for z in (phys_min[2], phys_max[2]):
                corners.append((float(x), float(y), float(z)))

    indices = np.array([image.TransformPhysicalPointToContinuousIndex(p) for p in corners])
    idx_lo = np.floor(indices.min(axis=0)).astype(int)
    idx_hi = np.ceil(indices.max(axis=0)).astype(int)

    size_img = np.array(image.GetSize(), dtype=int)
    idx_lo = np.clip(idx_lo, 0, size_img - 1)
    idx_hi = np.clip(idx_hi, 0, size_img - 1)
    start = idx_lo.tolist()
    crop_size = (idx_hi - idx_lo + 1).tolist()
    cropped = sitk.RegionOfInterest(image, crop_size, start)
    return cropped, start


def embed_crop_in_full_overview(
    full_overview: sitk.Image,
    crop: sitk.Image,
    crop_start_index: list[int],
) -> sitk.Image:
    """Place a cropped overlap volume into a zero-filled full overview canvas."""
    canvas = sitk.Image(full_overview.GetSize(), sitk.sitkFloat32)
    canvas.CopyInformation(full_overview)
    return sitk.Paste(
        canvas,
        crop,
        crop.GetSize(),
        [0, 0, 0],
        crop_start_index,
    )


def resample_to_reference_grid(moving: sitk.Image, reference: sitk.Image) -> sitk.Image:
    """Resample *moving* onto *reference* spacing, origin, direction, and size."""
    return sitk.Resample(
        moving,
        reference,
        sitk.Transform(3, sitk.sitkIdentity),
        sitk.sitkLinear,
        0.0,
        moving.GetPixelID(),
    )


def prepare_registration_pair(
    fixed: sitk.Image,
    roi_full: sitk.Image,
    *,
    margin_um: float = 0.0,
) -> tuple[sitk.Image, sitk.Image, tuple[np.ndarray, np.ndarray], list[int]]:
    """Crop fixed to shared FOV and resample ROI onto that grid."""
    overlap_box = overlap_physical_bounds(fixed, roi_full, margin_um=margin_um)
    fixed_cropped, crop_start_index = crop_to_physical_box(fixed, *overlap_box)
    moving = resample_to_reference_grid(roi_full, fixed_cropped)
    return fixed_cropped, moving, overlap_box, crop_start_index


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


def voxel_count_gb(image: sitk.Image) -> float:
    count = np.prod(image.GetSize(), dtype=np.int64)
    return count * 4 / 1e9
