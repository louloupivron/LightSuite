"""Affine fitting for control-point alignment (fitAffineTrans3D.m)."""

from __future__ import annotations

import numpy as np


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
