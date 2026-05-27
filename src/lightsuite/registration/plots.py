"""Registration preview figures (plotAnnotationComparison.m port)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import binary_dilation
from skimage.segmentation import find_boundaries

from lightsuite.gui.slices import volume_index_to_image
from lightsuite.registration.warp import warp_atlas_to_sample


def boundary_volume_from_annotation(annotation: np.ndarray) -> np.ndarray:
    """Fallback when annotation_boundary_*.nii.gz is unavailable."""
    labels = annotation.astype(np.int32)
    if not np.any(labels > 0):
        return np.zeros(labels.shape, dtype=np.uint8)
    edges = find_boundaries(labels, mode="inner")
    out = np.zeros(labels.shape, dtype=np.uint8)
    out[edges & (labels > 1)] = 255
    return out


def in_plane_pixel_size(
    dimplot: int,
    voxel_pxsize: tuple[float, float, float],
) -> tuple[float, float]:
    """Return (row spacing, col spacing) for a slice normal to ``dimplot``."""
    if dimplot not in {1, 2, 3}:
        msg = f"dimplot must be 1, 2, or 3, got {dimplot}"
        raise ValueError(msg)
    px = list(voxel_pxsize)
    px.pop(dimplot - 1)
    return float(px[0]), float(px[1])


def _prepare_boundary_slice(boundary_slice: np.ndarray, *, dilate: int = 1) -> np.ndarray:
    mask = np.asarray(boundary_slice) > 0
    if dilate > 0 and np.any(mask):
        mask = binary_dilation(mask, iterations=dilate)
    return mask


def plot_annotation_comparison(
    volume: np.ndarray,
    boundary: np.ndarray,
    dimplot: int,
    *,
    n_show: int = 8,
    voxel_pxsize: tuple[float, float, float] = (1.0, 1.0, 1.0),
    boundary_dilate: int = 1,
) -> plt.Figure:
    """Plot sample slices with warped atlas boundary mask overlaid (sample only)."""
    if dimplot not in {1, 2, 3}:
        msg = f"dimplot must be 1, 2, or 3, got {dimplot}"
        raise ValueError(msg)

    row_px, col_px = in_plane_pixel_size(dimplot, voxel_pxsize)
    ny = volume.shape[dimplot - 1]
    ishow = np.round(np.linspace(0.15 * ny, 0.85 * ny, n_show)).astype(int)
    ishow = np.clip(ishow, 1, ny)

    fig, axes = plt.subplots(2, n_show // 2, figsize=(17, 9))
    axes = np.atleast_1d(axes).ravel()

    for ii, islice in enumerate(ishow):
        ax = axes[ii]
        histim = volume_index_to_image(volume, np.array([islice, dimplot], dtype=int))
        boundim = volume_index_to_image(boundary, np.array([islice, dimplot], dtype=int))
        mask = _prepare_boundary_slice(boundim, dilate=boundary_dilate)

        height, width = histim.shape
        extent = (0.0, width * col_px, height * row_px, 0.0)
        ax.imshow(
            histim,
            cmap="gray",
            aspect="equal",
            origin="upper",
            extent=extent,
        )
        if np.any(mask):
            overlay = np.zeros((height, width, 4), dtype=float)
            overlay[mask] = (1.0, 0.8, 0.5, 0.9)
            ax.imshow(
                overlay,
                aspect="equal",
                origin="upper",
                extent=extent,
            )
        ax.set_title(str(islice))
        ax.axis("off")

    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.02, wspace=0.05, hspace=0.08)
    return fig


def save_initial_registration_previews(
    save_path: Path,
    sample: np.ndarray,
    boundary_atlas: np.ndarray,
    transform: np.ndarray,
    *,
    voxel_pxsize: tuple[float, float, float] = (1.0, 1.0, 1.0),
) -> int:
    """Write dim{1,2,3}_initial_registration.png using the coarse similarity transform."""
    save_path = Path(save_path)
    volmax = float(np.quantile(sample, 0.999)) or 1.0
    sample_u8 = np.clip(255.0 * sample / volmax, 0, 255).astype(np.uint8)
    boundary_warped = warp_atlas_to_sample(
        boundary_atlas.astype(np.float32),
        transform,
        sample.shape,
        order=0,
    )
    boundary_warped = (boundary_warped > 0).astype(np.uint8) * 255
    warped_count = int(np.count_nonzero(boundary_warped))

    for idim in range(1, 4):
        fig = plot_annotation_comparison(
            sample_u8,
            boundary_warped,
            idim,
            voxel_pxsize=voxel_pxsize,
        )
        out = save_path / f"dim{idim}_initial_registration.png"
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)

    return warped_count
