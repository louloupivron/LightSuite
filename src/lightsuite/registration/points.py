"""Point cloud extraction for coarse registration."""

from __future__ import annotations

import numpy as np
from scipy import ndimage


def cloud_xyz_to_volume_indices(xyz: np.ndarray) -> np.ndarray:
    """Map point-cloud XYZ (x, y, z) to volume array indices (Y, X, Z)."""
    pts = np.asarray(xyz, dtype=float)
    if pts.size == 0:
        return pts.reshape(0, 3)
    return pts[:, [1, 0, 2]]


def volume_indices_to_cloud_xyz(indices: np.ndarray) -> np.ndarray:
    """Inverse of :func:`cloud_xyz_to_volume_indices`."""
    pts = np.asarray(indices, dtype=float)
    if pts.size == 0:
        return pts.reshape(0, 3)
    return pts[:, [1, 0, 2]]


def _random_subsample(points: np.ndarray, fraction: float, seed: int = 1) -> np.ndarray:
    if points.shape[0] == 0:
        return points
    if fraction >= 1.0:
        return points
    rng = np.random.default_rng(seed)
    keep = max(1, int(points.shape[0] * fraction))
    if keep >= points.shape[0]:
        return points
    idx = rng.choice(points.shape[0], size=keep, replace=False)
    return points[idx]


def extract_sample_points(
    volume: np.ndarray,
    threshold: float,
    *,
    subsample_fraction: float = 0.1,
) -> np.ndarray:
    """Extract gradient-based sample cloud (extractSamplePoints.m)."""
    rng = np.random.default_rng(1)
    grad = np.sqrt(
        ndimage.sobel(volume, axis=0) ** 2
        + ndimage.sobel(volume, axis=1) ** 2
        + ndimage.sobel(volume, axis=2) ** 2
    )
    overall_mode = float(np.median(volume[volume > 0])) if np.any(volume > 0) else 0.0

    sizevol = volume.shape
    batch = max(1, int(np.ceil(max(sizevol) / 3)))
    nb = [max(1, int(np.ceil(s / batch))) for s in sizevol]
    points_chunks: list[np.ndarray] = []

    for iby in range(nb[0]):
        ys = slice(iby * batch, min((iby + 1) * batch, sizevol[0]))
        for ibx in range(nb[1]):
            xs = slice(ibx * batch, min((ibx + 1) * batch, sizevol[1]))
            for ibz in range(nb[2]):
                zs = slice(ibz * batch, min((ibz + 1) * batch, sizevol[2]))
                vol_curr = ndimage.median_filter(volume[ys, xs, zs], size=1)
                grad_curr = grad[ys, xs, zs].copy()
                flat = vol_curr.ravel()
                sample_n = min(10_000, flat.size)
                idx = rng.choice(flat.size, size=sample_n, replace=False) if flat.size > sample_n else np.arange(flat.size)
                thres_init = max(float(np.quantile(flat[idx], 0.05)) * 2.0, overall_mode)
                grad_curr[vol_curr < thres_init] = 0
                mask = (grad_curr / np.maximum(vol_curr, 1e-6) > threshold) & (vol_curr > thres_init)
                if not np.any(mask):
                    continue
                rr, cc, dd = np.where(mask)
                y_idx = np.arange(ys.start, ys.stop)[rr]
                x_idx = np.arange(xs.start, xs.stop)[cc]
                z_idx = np.arange(zs.start, zs.stop)[dd]
                points_chunks.append(np.column_stack([x_idx, y_idx, z_idx]))

    if not points_chunks:
        pts = np.zeros((0, 3), dtype=np.float64)
    else:
        pts = np.vstack(points_chunks)
        n_trim = max(1, int(min(sizevol) / 100))
        keep = (
            (pts[:, 0] >= n_trim)
            & (pts[:, 0] <= sizevol[1] - n_trim)
            & (pts[:, 1] >= n_trim)
            & (pts[:, 1] <= sizevol[0] - n_trim)
            & (pts[:, 2] >= n_trim)
            & (pts[:, 2] <= sizevol[2] - n_trim)
        )
        pts = pts[keep]
        pts = _random_subsample(pts, subsample_fraction, seed=1)
    if pts.shape[0] == 0:
        flat_grad = grad.ravel()
        flat_vol = volume.ravel()
        n_keep = min(500, flat_grad.size)
        idx = np.argpartition(flat_grad, -n_keep)[-n_keep:]
        idx = idx[flat_vol[idx] > 0]
        if idx.size:
            rr, cc, dd = np.unravel_index(idx, volume.shape)
            pts = np.column_stack([cc, rr, dd])
    return pts


def extract_atlas_points_gradient(
    template: np.ndarray,
    annotation: np.ndarray,
    *,
    sigma: float = 20.0,
    threshold: float = 5.0,
) -> np.ndarray:
    """Extract atlas cloud masked by annotation (extractVolumePointsGradient.m)."""
    rng = np.random.default_rng(1)
    vol = template.astype(np.float32)
    if sigma > 0:
        div = ndimage.gaussian_filter(vol, sigma=sigma, mode="reflect")
    else:
        div = vol
    grad = np.sqrt(
        ndimage.sobel(vol, axis=0) ** 2
        + ndimage.sobel(vol, axis=1) ** 2
        + ndimage.sobel(vol, axis=2) ** 2
    )
    test = grad / np.maximum(div, 1e-6)
    vol_masked = vol.copy()
    vol_masked[annotation == 0] = 0
    flat = vol_masked.ravel()
    sample_n = min(10_000, flat.size)
    idx = rng.choice(flat.size, size=sample_n, replace=False) if flat.size > sample_n else np.arange(flat.size)
    thres_init = float(np.quantile(flat[idx], 0.05)) * 2.0
    test[vol_masked < thres_init] = 0
    mask = (test > threshold) & (vol_masked > thres_init)
    if not np.any(mask):
        return np.zeros((0, 3), dtype=np.float64)
    rr, cc, dd = np.where(mask)
    return np.column_stack([cc, rr, dd]).astype(np.float64)
