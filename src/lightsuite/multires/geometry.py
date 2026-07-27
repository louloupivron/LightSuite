"""Physical geometry helpers for multiresolution overview / ROI volumes."""

from __future__ import annotations

import numpy as np
import SimpleITK as sitk


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


def crop_index_range_from_continuous_indices(
    indices: np.ndarray,
    size_img: np.ndarray,
) -> tuple[list[int], list[int]]:
    """Return XYZ start index and crop size covering continuous index corners."""
    idx_lo = np.floor(indices.min(axis=0)).astype(int)
    idx_hi = np.ceil(indices.max(axis=0)).astype(int)
    idx_lo = np.clip(idx_lo, 0, size_img - 1)
    idx_hi = np.clip(idx_hi, 0, size_img - 1)
    start = idx_lo.tolist()
    crop_size = (idx_hi - idx_lo + 1).tolist()
    return start, crop_size


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
    size_img = np.array(image.GetSize(), dtype=int)
    start, crop_size = crop_index_range_from_continuous_indices(indices, size_img)
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


def voxel_count_gb(image: sitk.Image) -> float:
    count = np.prod(image.GetSize(), dtype=np.int64)
    return count * 4 / 1e9
