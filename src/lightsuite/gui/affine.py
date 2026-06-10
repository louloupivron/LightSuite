"""Affine fitting for control-point alignment (fitAffineTrans3D.m)."""

from __future__ import annotations

import numpy as np


def fit_similarity_transform(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, float]:
    """Fit uniform-scale similarity mapping source -> target (fitSimilarityTrans3D.m)."""
    if source.shape != target.shape:
        msg = f"Point arrays must match shape, got {source.shape} vs {target.shape}"
        raise ValueError(msg)
    if source.shape[1] != 3:
        msg = "Points must be Nx3"
        raise ValueError(msg)
    n = source.shape[0]
    if n < 3:
        msg = "Need at least 3 point pairs for similarity fit"
        raise ValueError(msg)

    centroid_source = source.mean(axis=0)
    centroid_target = target.mean(axis=0)
    source_centered = source - centroid_source
    target_centered = target - centroid_target

    h = source_centered.T @ target_centered
    u, _, vt = np.linalg.svd(h)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T

    scale_num = float(np.sum((target_centered @ rotation) * source_centered))
    scale_den = float(np.sum(source_centered**2))
    scale = scale_num / scale_den if scale_den > 0 else 1.0
    translation = centroid_target - scale * rotation @ centroid_source

    matrix = np.eye(4)
    matrix[:3, :3] = scale * rotation
    matrix[:3, 3] = translation
    predicted = transform_points(source, matrix)
    mse = float(np.mean(np.sum((predicted - target) ** 2, axis=1)))
    return matrix, mse


def fit_affine_transform(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, float]:
    """Fit 4x4 affine mapping source -> target (least squares)."""
    if source.shape != target.shape:
        msg = f"Point arrays must match shape, got {source.shape} vs {target.shape}"
        raise ValueError(msg)
    if source.shape[1] != 3:
        msg = "Points must be Nx3"
        raise ValueError(msg)
    n = source.shape[0]
    if n < 4:
        msg = "Need at least 4 point pairs for affine fit"
        raise ValueError(msg)

    augmented = np.column_stack([source, np.ones(n)])
    coeffs, _, _, _ = np.linalg.lstsq(augmented, target, rcond=None)
    matrix = np.eye(4)
    matrix[:3, :] = coeffs.T
    predicted = augmented @ coeffs
    mse = float(np.mean(np.sum((predicted - target) ** 2, axis=1)))
    return matrix, mse


def transform_points(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    if points.size == 0:
        return points.reshape(0, 3)
    hom = np.column_stack([points, np.ones(points.shape[0])])
    return (hom @ matrix.T)[:, :3]


def transform_points_inverse(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    inv = np.linalg.inv(matrix)
    return transform_points(points, inv)


def affine_point_errors(
    source: np.ndarray,
    target: np.ndarray,
    matrix: np.ndarray,
) -> tuple[float, np.ndarray]:
    """Return (MSE, per-point Euclidean errors) for ``matrix`` mapping source -> target."""
    predicted = transform_points(source, matrix)
    diff = predicted - target
    per_point = np.linalg.norm(diff, axis=1)
    mse = float(np.mean(per_point**2))
    return mse, per_point


def summarize_point_errors(errors: np.ndarray) -> dict[str, float]:
    """Summary stats for a 1D error array (voxels)."""
    if errors.size == 0:
        return {"median": 0.0, "p95": 0.0, "max": 0.0, "mean": 0.0}
    return {
        "median": float(np.median(errors)),
        "p95": float(np.percentile(errors, 95)),
        "max": float(np.max(errors)),
        "mean": float(np.mean(errors)),
    }
