"""Point cloud extraction for coarse registration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage
from scipy.stats import mode as sp_mode

from lightsuite.registration.pc_downsample import pcdenoise, pcdownsample_random


@dataclass(frozen=True)
class SampleCloudExtractionStages:
    """Counts through ``extractSamplePoints.m`` post-processing."""

    mask_points: int
    trim_points: int
    downsample_points: int
    denoise_points: int


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


_MATLAB_IMGRADIENT3_SOBEL_Y = np.array(
    [
        [[1, 3, 1], [3, 6, 3], [1, 3, 1]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[-1, -3, -1], [-3, -6, -3], [-1, -3, -1]],
    ],
    dtype=np.float32,
)
_MATLAB_IMGRADIENT3_SOBEL_X = np.array(
    [
        [[1, 3, 1], [0, 0, 0], [-1, -3, -1]],
        [[3, 6, 3], [0, 0, 0], [-3, -6, -3]],
        [[1, 3, 1], [0, 0, 0], [-1, -3, -1]],
    ],
    dtype=np.float32,
)
_MATLAB_IMGRADIENT3_SOBEL_Z = np.array(
    [
        [[1, 0, -1], [3, 0, -3], [1, 0, -1]],
        [[3, 0, -3], [6, 0, -6], [3, 0, -3]],
        [[1, 0, -1], [3, 0, -3], [1, 0, -1]],
    ],
    dtype=np.float32,
)


def _imgradient3_magnitude(volume: np.ndarray) -> np.ndarray:
    """MATLAB ``imgradient3`` default Sobel magnitude.

    ``scipy.ndimage.sobel`` applies separable 1-D filters and does not match the
    full 3x3x3 kernels used by MATLAB (see ``imgradient3`` / ``imgradientxyz``).
    """
    vol = np.asarray(volume, dtype=np.float32, order="C")
    gy = ndimage.convolve(vol, _MATLAB_IMGRADIENT3_SOBEL_Y, mode="reflect")
    gx = ndimage.convolve(vol, _MATLAB_IMGRADIENT3_SOBEL_X, mode="reflect")
    gz = ndimage.convolve(vol, _MATLAB_IMGRADIENT3_SOBEL_Z, mode="reflect")
    return np.sqrt(gx * gx + gy * gy + gz * gz)


def _volume_mode_all(volume: np.ndarray) -> float:
    """MATLAB ``mode(voluse, 'all')`` on the Fortran-ordered flattening."""
    flat = np.asarray(volume, dtype=np.float64).ravel(order="F")
    if flat.size == 0:
        return 0.0
    return float(sp_mode(flat, keepdims=False).mode)


def _batch_intensity_threshold(
    vol_curr: np.ndarray,
    overall_mode: float,
    rng: np.random.Generator,
) -> float:
    """Per-batch ``thresinit`` from ``extractSamplePoints.m`` / ``extractVolumePointsGradient.m``."""
    flat = vol_curr.ravel(order="F")
    sample_n = min(10_000, flat.size)
    if flat.size > sample_n:
        idx = rng.choice(flat.size, size=sample_n, replace=False)
    else:
        idx = np.arange(flat.size)
    return max(float(np.quantile(flat[idx], 0.05)) * 2.0, overall_mode)


def _trim_border_points(
    points: np.ndarray,
    sizevol: tuple[int, ...],
) -> np.ndarray:
    n_trim = max(1, int(min(sizevol) / 100))
    keep = (
        (points[:, 0] >= n_trim)
        & (points[:, 0] <= sizevol[1] - n_trim)
        & (points[:, 1] >= n_trim)
        & (points[:, 1] <= sizevol[0] - n_trim)
        & (points[:, 2] >= n_trim)
        & (points[:, 2] <= sizevol[2] - n_trim)
    )
    return points[keep]


def _extract_sample_mask_points(
    volume: np.ndarray,
    threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Mask and border-trimmed sample points before downsample/denoise."""
    rng = np.random.default_rng(1)
    grad = _imgradient3_magnitude(volume)
    overall_mode = _volume_mode_all(volume)

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
                thres_init = _batch_intensity_threshold(vol_curr, overall_mode, rng)
                grad_curr[vol_curr < thres_init] = 0
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = grad_curr / vol_curr
                mask = (ratio > threshold) & (vol_curr > thres_init)
                if not np.any(mask):
                    continue
                rr, cc, dd = np.where(mask)
                y_idx = np.arange(ys.start, ys.stop)[rr]
                x_idx = np.arange(xs.start, xs.stop)[cc]
                z_idx = np.arange(zs.start, zs.stop)[dd]
                points_chunks.append(np.column_stack([x_idx, y_idx, z_idx]))

    if not points_chunks:
        mask_pts = np.zeros((0, 3), dtype=np.float64)
    else:
        mask_pts = np.vstack(points_chunks)
    return mask_pts, _trim_border_points(mask_pts, sizevol)


def extract_sample_points_stages(
    volume: np.ndarray,
    threshold: float,
    *,
    subsample_fraction: float = 0.1,
) -> tuple[np.ndarray, SampleCloudExtractionStages]:
    """Extract the sample cloud and return per-stage point counts."""
    grad = _imgradient3_magnitude(volume)
    mask_pts, trim_pts = _extract_sample_mask_points(volume, threshold)
    if trim_pts.shape[0] == 0:
        pts = _fallback_gradient_points(volume, grad)
        stages = SampleCloudExtractionStages(0, 0, 0, int(pts.shape[0]))
        return pts, stages
    down_pts = pcdownsample_random(
        trim_pts,
        subsample_fraction,
        preserve_structure=True,
        seed=1,
    )
    denoise_pts = pcdenoise(down_pts)
    if denoise_pts.shape[0] == 0:
        pts = _fallback_gradient_points(volume, grad)
        stages = SampleCloudExtractionStages(
            int(mask_pts.shape[0]),
            int(trim_pts.shape[0]),
            int(down_pts.shape[0]),
            int(pts.shape[0]),
        )
        return pts, stages
    stages = SampleCloudExtractionStages(
        mask_points=int(mask_pts.shape[0]),
        trim_points=int(trim_pts.shape[0]),
        downsample_points=int(down_pts.shape[0]),
        denoise_points=int(denoise_pts.shape[0]),
    )
    return denoise_pts, stages


def extract_sample_points(
    volume: np.ndarray,
    threshold: float,
    *,
    subsample_fraction: float = 0.1,
) -> np.ndarray:
    """Extract gradient-based sample cloud (extractSamplePoints.m)."""
    pts, _ = extract_sample_points_stages(
        volume,
        threshold,
        subsample_fraction=subsample_fraction,
    )
    return pts


def _fallback_gradient_points(volume: np.ndarray, grad: np.ndarray) -> np.ndarray:
    flat_grad = grad.ravel(order="F")
    flat_vol = volume.ravel(order="F")
    n_keep = min(500, flat_grad.size)
    idx = np.argpartition(flat_grad, -n_keep)[-n_keep:]
    idx = idx[flat_vol[idx] > 0]
    if idx.size:
        rr, cc, dd = np.unravel_index(idx, volume.shape, order="F")
        return np.column_stack([cc, rr, dd])
    return np.zeros((0, 3), dtype=np.float64)


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
    grad = _imgradient3_magnitude(vol)
    with np.errstate(divide="ignore", invalid="ignore"):
        test = grad / div
    vol_masked = vol.copy()
    vol_masked[annotation == 0] = 0
    flat = vol.ravel(order="F")
    sample_n = min(10_000, flat.size)
    if flat.size > sample_n:
        idx = rng.choice(flat.size, size=sample_n, replace=False)
    else:
        idx = np.arange(flat.size)
    thres_init = float(np.quantile(flat[idx], 0.05)) * 2.0
    test[vol_masked < thres_init] = 0
    mask = (test > threshold) & (vol_masked > thres_init)
    if not np.any(mask):
        return np.zeros((0, 3), dtype=np.float64)
    rr, cc, dd = np.where(mask)
    return np.column_stack([cc, rr, dd]).astype(np.float64)
