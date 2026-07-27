"""Manifest-only physical geometry (no voxel allocation)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import SimpleITK as sitk

from lightsuite.multires.geometry import crop_index_range_from_continuous_indices
from lightsuite.multires.volume import apply_manifest_geometry

if TYPE_CHECKING:
    from lightsuite.multires.models import ManifestVolumeSpec


def _direction_matrix(spec: ManifestVolumeSpec) -> np.ndarray:
    return np.array(spec.direction, dtype=float).reshape(3, 3)


def index_xyz_to_physical(spec: ManifestVolumeSpec, index_xyz: np.ndarray | tuple[float, ...]) -> np.ndarray:
    """Map continuous XYZ indices to physical coordinates (µm)."""
    origin = np.asarray(spec.origin_um, dtype=float)
    spacing = np.asarray(spec.spacing_um, dtype=float)
    direction = _direction_matrix(spec)
    index = np.asarray(index_xyz, dtype=float)
    return origin + direction @ (spacing * index)


def physical_to_continuous_index_xyz(
    spec: ManifestVolumeSpec,
    point_um: np.ndarray | tuple[float, ...],
) -> np.ndarray:
    """Map physical coordinates (µm) to continuous XYZ indices."""
    origin = np.asarray(spec.origin_um, dtype=float)
    spacing = np.asarray(spec.spacing_um, dtype=float)
    direction = _direction_matrix(spec)
    mat = direction @ np.diag(spacing)
    return np.linalg.solve(mat, np.asarray(point_um, dtype=float) - origin)


def shape_xyz(spec: ManifestVolumeSpec) -> tuple[int, int, int]:
    nz, ny, nx = (int(v) for v in spec.shape_zyx)
    return nx, ny, nz


def physical_bounds_from_spec(spec: ManifestVolumeSpec) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned physical bounding box from manifest fields only."""
    nx, ny, nz = shape_xyz(spec)
    corners = []
    for ix in (0, nx - 1):
        for iy in (0, ny - 1):
            for iz in (0, nz - 1):
                corners.append(index_xyz_to_physical(spec, (ix, iy, iz)))
    corners_arr = np.asarray(corners, dtype=float)
    return corners_arr.min(axis=0), corners_arr.max(axis=0)


def physical_center_from_spec(spec: ManifestVolumeSpec) -> np.ndarray:
    nx, ny, nz = shape_xyz(spec)
    return index_xyz_to_physical(spec, ((nx - 1) / 2.0, (ny - 1) / 2.0, (nz - 1) / 2.0))


def transformed_bounds_from_spec(
    spec: ManifestVolumeSpec,
    source_to_target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map a volume's physical bounding box into another coordinate frame."""
    from lightsuite.multires.geometry import physical_corners, transform_physical_points

    pmin, pmax = physical_bounds_from_spec(spec)
    transformed = transform_physical_points(physical_corners(pmin, pmax), source_to_target)
    return transformed.min(axis=0), transformed.max(axis=0)


def overlap_physical_bounds_from_specs(
    spec_a: ManifestVolumeSpec,
    spec_b: ManifestVolumeSpec,
    *,
    margin_um: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Intersection of two manifest volume bounding boxes."""
    min_a, max_a = physical_bounds_from_spec(spec_a)
    min_b, max_b = physical_bounds_from_spec(spec_b)
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


def crop_index_range_from_physical_box(
    spec: ManifestVolumeSpec,
    phys_min: np.ndarray,
    phys_max: np.ndarray,
) -> tuple[list[int], list[int]]:
    """Return XYZ start index and crop size covering a physical box."""
    nx, ny, nz = shape_xyz(spec)
    size_img = np.array([nx, ny, nz], dtype=int)
    corners = []
    for x in (float(phys_min[0]), float(phys_max[0])):
        for y in (float(phys_min[1]), float(phys_max[1])):
            for z in (float(phys_min[2]), float(phys_max[2])):
                corners.append((x, y, z))

    indices = np.array([physical_to_continuous_index_xyz(spec, point) for point in corners])
    start, crop_size = crop_index_range_from_continuous_indices(indices, size_img)
    return start, crop_size


def sitk_geometry_from_spec(spec: ManifestVolumeSpec) -> sitk.Image:
    """Minimal 1-voxel image carrying manifest spacing, origin, and direction."""
    image = sitk.Image([1, 1, 1], sitk.sitkFloat32)
    return apply_manifest_geometry(image, spec)


def overlap_box_from_landmark_specs(
    overview_spec: ManifestVolumeSpec,
    roi_spec: ManifestVolumeSpec,
    roi_to_overview: np.ndarray,
    *,
    margin_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Predict overlap crop box from landmark transform using manifest bounds only."""
    roi_min, roi_max = transformed_bounds_from_spec(roi_spec, roi_to_overview)
    overview_min, overview_max = physical_bounds_from_spec(overview_spec)
    overlap_min = np.maximum(roi_min, overview_min) - margin_um
    overlap_max = np.minimum(roi_max, overview_max) + margin_um
    if np.any(overlap_min >= overlap_max):
        msg = (
            "Landmark transform places ROI outside the overview bounds.\n"
            f"  ROI in overview space: {roi_min} .. {roi_max}\n"
            f"  Overview bounds: {overview_min} .. {overview_max}"
        )
        raise ValueError(msg)
    return overlap_min, overlap_max


def alignment_metrics_from_specs(
    overview_spec: ManifestVolumeSpec,
    roi_spec: ManifestVolumeSpec,
    *,
    overlap_min: np.ndarray | None = None,
    overlap_max: np.ndarray | None = None,
) -> dict[str, object]:
    """Quantify manifest-only alignment between overview and ROI volumes."""
    overview_center = physical_center_from_spec(overview_spec)
    roi_center = physical_center_from_spec(roi_spec)
    center_offset_um = roi_center - overview_center
    if overlap_min is None or overlap_max is None:
        overlap_min, overlap_max = overlap_physical_bounds_from_specs(overview_spec, roi_spec)
    overlap_extent = np.maximum(0.0, overlap_max - overlap_min)
    overlap_volume_um3 = float(np.prod(overlap_extent))
    roi_extent = physical_bounds_from_spec(roi_spec)[1] - physical_bounds_from_spec(roi_spec)[0]
    overview_extent = physical_bounds_from_spec(overview_spec)[1] - physical_bounds_from_spec(overview_spec)[0]
    roi_volume_um3 = float(np.prod(roi_extent))
    overview_volume_um3 = float(np.prod(overview_extent))
    return {
        "center_offset_um": center_offset_um.tolist(),
        "center_offset_norm_um": float(np.linalg.norm(center_offset_um)),
        "overlap_center_um": (0.5 * (overlap_min + overlap_max)).tolist(),
        "overlap_extent_um": overlap_extent.tolist(),
        "roi_overlap_fraction": overlap_volume_um3 / roi_volume_um3 if roi_volume_um3 > 0 else 0.0,
        "overview_overlap_fraction": overlap_volume_um3 / overview_volume_um3
        if overview_volume_um3 > 0
        else 0.0,
    }


def manifest_geometry_report_from_spec(label: str, spec: ManifestVolumeSpec) -> dict[str, object]:
    """Build a diagnostic geometry report without loading voxels."""
    pmin, pmax = physical_bounds_from_spec(spec)
    center = physical_center_from_spec(spec)
    z, y, x = (int(v) for v in spec.shape_zyx)
    spacing = tuple(float(s) for s in spec.spacing_um)
    return {
        "label": label,
        "shape_zyx": tuple(spec.shape_zyx),
        "shape_xyz": (x, y, z),
        "spacing_um": spacing,
        "origin_um": tuple(float(v) for v in spec.origin_um),
        "phys_min": pmin,
        "phys_max": pmax,
        "phys_center": center,
        "extent_um": tuple(float(v) for v in (pmax - pmin)),
    }
