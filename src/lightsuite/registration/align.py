"""Coarse similarity alignment (originalSimilarityTform.m / triageAndMatchClouds.m)."""

from __future__ import annotations

import numpy as np
import open3d as o3d
from scipy.spatial import cKDTree


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


def estimate_similarity_transform(
    atlas_points: np.ndarray,
    sample_points: np.ndarray,
    *,
    max_distance: float = 25.0,
) -> np.ndarray:
    """Estimate 4x4 transform mapping sample points toward atlas space."""
    n_down_ls = max(6, int(round(sample_points.shape[0] / 10_000)))
    n_down_tv = max(6, int(round(atlas_points.shape[0] / 50_000)))
    atlas_use = _voxel_grid_downsample(atlas_points, n_down_tv)
    sample_use = _voxel_grid_downsample(sample_points, n_down_ls)

    source = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(sample_use))
    target = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(atlas_use))

    result = o3d.pipelines.registration.registration_icp(
        source,
        target,
        max_correspondence_distance=max_distance,
        estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint(
            with_scaling=True
        ),
        criteria=o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=100),
    )
    transform = result.transformation.copy()

    # Outlier rejection pass (MATLAB originalSimilarityTform)
    sample_aligned = _transform_points(sample_use, transform)
    tree_a = cKDTree(atlas_use)
    tree_s = cKDTree(sample_aligned)
    dist_s, _ = tree_a.query(sample_aligned, k=1)
    dist_a, _ = tree_s.query(atlas_use, k=1)
    keep_s = dist_s <= max_distance
    keep_a = dist_a <= max_distance
    if keep_s.sum() >= 6 and keep_a.sum() >= 6:
        source2 = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(sample_use[keep_s]))
        target2 = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(atlas_use[keep_a]))
        result2 = o3d.pipelines.registration.registration_icp(
            source2,
            target2,
            max_correspondence_distance=max_distance,
            estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint(
                with_scaling=True
            ),
        )
        transform = result2.transformation

    # MATLAB stores the BCPD/ICP sample->atlas fit as original_trans (transinit).
    return transform


def triage_and_match_clouds(
    sample_points: np.ndarray,
    atlas_points: np.ndarray,
    transform: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Select candidate corresponding points after coarse alignment."""
    n_neigh_atlas = max(6, int(np.ceil(atlas_points.shape[0] / 50_000)))
    n_neigh_samp = max(6, int(np.ceil(sample_points.shape[0] / 10_000)))
    atlas_use = _voxel_grid_downsample(atlas_points, n_neigh_atlas)
    sample_use = _voxel_grid_downsample(sample_points, n_neigh_samp)
    sample_fwd = _transform_points(sample_use, transform)

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
    sample_orig = _transform_points(sample_kept, np.linalg.inv(transform))
    return sample_orig, atlas_kept


def _transform_points(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    if points.shape[0] == 0:
        return points
    hom = np.column_stack([points, np.ones(points.shape[0])])
    out = hom @ matrix.T
    return out[:, :3]
