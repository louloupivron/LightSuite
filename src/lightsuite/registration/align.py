"""Coarse similarity alignment (originalSimilarityTform.m / triageAndMatchClouds.m)."""

from __future__ import annotations

import numpy as np
import open3d as o3d
from scipy.spatial import cKDTree

from lightsuite.registration.warp import matlab_voxel_affine_from_icp


def _downsample_points(points: np.ndarray, target_count: int, seed: int = 1) -> np.ndarray:
    if points.shape[0] <= target_count:
        return points
    rng = np.random.default_rng(seed)
    idx = rng.choice(points.shape[0], size=target_count, replace=False)
    return points[idx]


def _voxel_grid_downsample(points: np.ndarray, max_points: int) -> np.ndarray:
    if points.shape[0] <= max_points:
        return points
    pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(points))
    voxel = max(1.0, (points.shape[0] / max_points) ** (1 / 3))
    down = pcd.voxel_down_sample(voxel)
    return np.asarray(down.points)


def _matlab_cloud_subset(points: np.ndarray, count_divisor: int) -> np.ndarray:
    """Match MATLAB triage/BCPD downsampling (use all points when target count < 6)."""
    target = int(np.ceil(points.shape[0] / count_divisor))
    if target >= 6:
        return _voxel_grid_downsample(points, target)
    return points


def _icp_cloud_subset(
    points: np.ndarray,
    count_divisor: int,
    *,
    cap: int,
) -> np.ndarray:
    """Choose an ICP subset; unlike MATLAB BCPD, Open3D needs thousands of points."""
    n = points.shape[0]
    target = int(np.ceil(n / count_divisor))
    if target < 6:
        count = min(n, cap)
    else:
        count = min(max(target, 500), cap)
    if n <= count:
        return points
    return _voxel_grid_downsample(points, count)


def similarity_scale(matrix: np.ndarray) -> float:
    """Uniform scale component of a 3D similarity transform."""
    return float(np.cbrt(np.linalg.det(np.asarray(matrix, dtype=float)[:3, :3])))


def _clamp_similarity_scale(
    transform: np.ndarray,
    *,
    lo: float = 0.75,
    hi: float = 1.35,
) -> np.ndarray:
    """Limit similarity scale to a plausible brain registration range."""
    out = np.asarray(transform, dtype=float).copy()
    scale = similarity_scale(out)
    if lo <= scale <= hi:
        return out
    clamped = float(np.clip(scale, lo, hi))
    linear = out[:3, :3] / scale * clamped
    out[:3, :3] = linear
    return out


def _similarity_init_from_centroids(
    sample_points: np.ndarray,
    atlas_points: np.ndarray,
) -> np.ndarray:
    """Rough similarity init mapping sample centroid/scale to atlas."""
    sample_centroid = sample_points.mean(axis=0)
    atlas_centroid = atlas_points.mean(axis=0)
    sample_scale = float(np.sqrt(np.mean(np.sum((sample_points - sample_centroid) ** 2, axis=1))))
    atlas_scale = float(np.sqrt(np.mean(np.sum((atlas_points - atlas_centroid) ** 2, axis=1))))
    scale = float(np.clip(atlas_scale / max(sample_scale, 1e-6), 0.75, 1.35))

    transform = np.eye(4, dtype=float)
    transform[:3, :3] *= scale
    transform[:3, 3] = atlas_centroid - scale * sample_centroid
    return transform


def downsample_point_cloud(points: np.ndarray, max_points: int) -> np.ndarray:
    """Reduce a point cloud to at most ``max_points`` via voxel downsampling."""
    return _voxel_grid_downsample(points, max_points)


def _run_icp(
    sample_points: np.ndarray,
    atlas_points: np.ndarray,
    *,
    init: np.ndarray,
    max_distance: float,
) -> np.ndarray:
    source = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(sample_points))
    target = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(atlas_points))
    result = o3d.pipelines.registration.registration_icp(
        source,
        target,
        max_correspondence_distance=max_distance,
        init=init,
        estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint(
            with_scaling=True
        ),
        criteria=o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=100),
    )
    return result.transformation.copy()


def estimate_similarity_transform(
    atlas_points: np.ndarray,
    sample_points: np.ndarray,
    *,
    max_distance: float = 25.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Estimate sample->atlas similarity transform.

    Returns ``(icp_transform, matlab_transform)`` where ``icp_transform`` operates
    on 0-based XYZ point indices and ``matlab_transform`` is stored in regopts.json.
    """
    sample_cap = min(sample_points.shape[0], 12_000)
    sample_use = _icp_cloud_subset(sample_points, 10_000, cap=sample_cap)
    atlas_cap = min(30_000, max(2_000, sample_use.shape[0] * 4))
    atlas_use = _icp_cloud_subset(atlas_points, 50_000, cap=atlas_cap)
    if atlas_use.shape[0] > atlas_cap:
        atlas_use = _voxel_grid_downsample(atlas_use, atlas_cap)

    init_transform = _similarity_init_from_centroids(sample_use, atlas_use)
    transform = init_transform
    for distance in (max(50.0, max_distance * 2.0), max_distance, max_distance * 0.6):
        transform = _run_icp(
            sample_use,
            atlas_use,
            init=transform,
            max_distance=distance,
        )

    sample_aligned = _transform_points(sample_use, transform)
    tree_a = cKDTree(atlas_use)
    tree_s = cKDTree(sample_aligned)
    dist_s, _ = tree_a.query(sample_aligned, k=1)
    dist_a, _ = tree_s.query(atlas_use, k=1)
    keep_s = dist_s <= max_distance
    keep_a = dist_a <= max_distance
    if keep_s.sum() >= 6 and keep_a.sum() >= 6:
        transform = _run_icp(
            sample_use[keep_s],
            atlas_use[keep_a],
            init=transform,
            max_distance=max_distance,
        )

    transform = _clamp_similarity_scale(transform)
    matlab_transform = matlab_voxel_affine_from_icp(transform)
    return transform, matlab_transform


def triage_and_match_clouds(
    sample_points: np.ndarray,
    atlas_points: np.ndarray,
    transform_icp: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Select candidate corresponding points after coarse alignment."""
    atlas_use = _matlab_cloud_subset(atlas_points, 50_000)
    sample_use = _matlab_cloud_subset(sample_points, 10_000)
    sample_fwd = _transform_points(sample_use, transform_icp)

    tree_a = cKDTree(atlas_use)
    tree_s = cKDTree(sample_fwd)
    dist_sf, _ = tree_a.query(sample_fwd, k=1)
    dist_fs, _ = tree_s.query(atlas_use, k=1)
    sample_fwd = sample_fwd[dist_sf <= 50]
    atlas_use = atlas_use[dist_fs <= 50]

    if sample_fwd.shape[0] < 6 or atlas_use.shape[0] < 6:
        return np.zeros((0, 3)), np.zeros((0, 3))

    tree_a = cKDTree(atlas_use)
    dists, nn_idx = tree_a.query(sample_fwd, k=1)
    keep = dists < 20
    if not np.any(keep):
        return np.zeros((0, 3)), np.zeros((0, 3))

    sample_kept = sample_fwd[keep]
    atlas_kept = atlas_use[nn_idx[keep]]
    sample_orig = _transform_points(sample_kept, np.linalg.inv(transform_icp))
    return sample_orig, atlas_kept


def _transform_points(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    if points.shape[0] == 0:
        return points
    hom = np.column_stack([points, np.ones(points.shape[0])])
    out = hom @ matrix.T
    return out[:, :3]
