"""MATLAB pointCloud downsampling ports (pcdownsample, pcdenoise)."""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def pcdownsample_random(
    points: np.ndarray,
    percentage: float,
    *,
    preserve_structure: bool = False,
    seed: int = 1,
) -> np.ndarray:
    """Port of ``pcdownsample(ptCloud, 'random', percentage, PreserveStructure=...)``."""
    points = np.asarray(points, dtype=np.float64)
    n = points.shape[0]
    if n == 0 or percentage >= 1.0:
        return points
    num_out = max(1, int(round(n * percentage)))
    if num_out >= n:
        return points
    rng = np.random.default_rng(seed)
    # MATLAB ``pcdownsample(..., 'random', p)`` draws ``round(Count*p)`` points.
    # ``PreserveStructure`` only affects organised clouds; ``extractSamplePoints.m``
    # builds an unorganised cloud, so the flag is a no-op there.
    del preserve_structure
    return points[rng.choice(n, size=num_out, replace=False)]


def pcdenoise(
    points: np.ndarray,
    *,
    num_neighbors: int = 4,
    std_ratio: float = 1.0,
) -> np.ndarray:
    """Remove outliers like MATLAB ``pcdenoise(ptcloud)`` / PCL statistical outlier removal.

  1. For each point, average distance to ``num_neighbors`` nearest neighbours
     (excluding the query point itself).
    2. Compute the global mean and sample standard deviation of those averages.
    3. Keep points with average distance ``<= mean + std_ratio * stddev``.

    This matches PCL's ``StatisticalOutlierRemoval`` and MATLAB's documented
    ``NumNeighbors`` / ``Threshold`` defaults (4 neighbours, 1σ). Open3D's
    implementation queries ``nb_neighbors`` points including the query itself,
    which is tighter and removes more points than MATLAB on dense clouds.
    """
    points = np.asarray(points, dtype=np.float64)
    n = points.shape[0]
    if n <= num_neighbors:
        return points

    tree = cKDTree(points)
    dists, _ = tree.query(points, k=num_neighbors + 1, workers=-1)
    mean_distances = dists[:, 1:].mean(axis=1)

    mu = float(mean_distances.mean())
    variance = float(
        ((mean_distances * mean_distances).sum() - mean_distances.sum() ** 2 / n)
        / max(n - 1, 1)
    )
    stddev = float(np.sqrt(max(variance, 0.0)))
    cutoff = mu + std_ratio * stddev
    return points[mean_distances <= cutoff]


def _kdtree_leaves(points: np.ndarray, max_num_points: int) -> list[np.ndarray]:
    """Split a cloud into ``2^ceil(log2(N/maxNumPoints))`` median-split kd-tree leaves.

    MATLAB's nonuniform grid methods take ``maxNumPoints`` per grid box rather
    than a grid step, so the leaf count is a power of two.
    """
    n = points.shape[0]
    levels = max(0, int(np.ceil(np.log2(n / max_num_points))))
    groups: list[np.ndarray] = [np.arange(n)]
    for _ in range(levels):
        next_groups: list[np.ndarray] = []
        for idx in groups:
            if idx.size <= 1:
                next_groups.append(idx)
                continue
            pts = points[idx]
            axis = int(np.argmax(pts.max(axis=0) - pts.min(axis=0)))
            order = np.argsort(pts[:, axis], kind="stable")
            half = idx.size // 2
            next_groups.append(idx[order[:half]])
            next_groups.append(idx[order[half:]])
        groups = next_groups
    return [idx for idx in groups if idx.size > 0]


def nonuniform_grid_sample(
    points: np.ndarray,
    max_num_points: int,
    *,
    seed: int = 1,
) -> np.ndarray:
    """Port of ``pcdownsample(..., 'nonuniformGridSample', maxNumPoints)``.

    Keeps one original point per kd-tree leaf (``rng(1)`` with ``PreserveStructure``).
    """
    points = np.asarray(points, dtype=np.float64)
    n = points.shape[0]
    if n == 0 or max_num_points < 6 or n <= max_num_points:
        return points

    rng = np.random.default_rng(seed)
    leaves = _kdtree_leaves(points, max_num_points)
    return np.vstack([points[idx[int(rng.integers(idx.size))]] for idx in leaves])


def nonuniform_grid(points: np.ndarray, max_num_points: int) -> np.ndarray:
    """Port of ``pcdownsample(..., 'nonuniformGrid', maxNumPoints)`` (leaf centroids)."""
    points = np.asarray(points, dtype=np.float64)
    n = points.shape[0]
    if n == 0 or max_num_points < 6 or n <= max_num_points:
        return points
    leaves = _kdtree_leaves(points, max_num_points)
    return np.vstack([points[idx].mean(axis=0) for idx in leaves])


def matlab_bcpd_grid_step(point_count: int, count_divisor: int) -> int:
    """``max(6, round(Count/count_divisor))`` from ``originalSimilarityTform.m``."""
    return max(6, int(round(point_count / count_divisor)))


def matlab_triage_grid_step(point_count: int, count_divisor: int) -> int:
    """``ceil(Count/count_divisor)`` from ``triageAndMatchClouds.m``."""
    return int(np.ceil(point_count / count_divisor))


def downsample_for_bcpd_similarity(points: np.ndarray, count_divisor: int) -> np.ndarray:
    """Downsample clouds before coarse BCPD (``nonuniformGridSample``)."""
    max_num_points = matlab_bcpd_grid_step(points.shape[0], count_divisor)
    if max_num_points < 6:
        return points
    return nonuniform_grid_sample(points, max_num_points)


def downsample_for_triage(points: np.ndarray, count_divisor: int) -> np.ndarray:
    """Downsample clouds before triage BCPD (``nonuniformGrid``)."""
    max_num_points = matlab_triage_grid_step(points.shape[0], count_divisor)
    if max_num_points < 6:
        return points
    return nonuniform_grid(points, max_num_points)
