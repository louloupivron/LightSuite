"""Registration preview figures (plotAnnotationComparison.m port)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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


def mask_boundary_pixels(boundary_slice: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return row/column indices for a precomputed boundary mask (0/255)."""
    mask = np.asarray(boundary_slice)
    row, col = np.nonzero(mask > 0)
    return row, col


def annotation_boundary_pixels(annotation_slice: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return row/column indices from an annotation label slice (legacy fallback)."""
    atlasim = np.asarray(annotation_slice, dtype=np.float32)
    gy, gx = np.gradient(atlasim)
    boundaries = ((gx != 0) | (gy != 0)) & (atlasim > 1)
    row, col = np.nonzero(boundaries)
    return row, col


def plot_annotation_comparison(
    volume: np.ndarray,
    boundary: np.ndarray,
    dimplot: int,
    *,
    n_show: int = 8,
    pxsize: tuple[float, float] = (1.0, 1.0),
) -> plt.Figure:
    """Plot sample slices with warped atlas boundary mask overlaid (sample only)."""
    if dimplot not in {1, 2, 3}:
        msg = f"dimplot must be 1, 2, or 3, got {dimplot}"
        raise ValueError(msg)

    ny = volume.shape[dimplot - 1]
    ishow = np.round(np.linspace(0.15 * ny, 0.85 * ny, n_show)).astype(int)
    ishow = np.clip(ishow, 1, ny)

    fig, axes = plt.subplots(2, n_show // 2, figsize=(17, 9))
    axes = np.atleast_1d(axes).ravel()

    for ii, islice in enumerate(ishow):
        ax = axes[ii]
        histim = volume_index_to_image(volume, np.array([islice, dimplot], dtype=int))
        annotim = volume_index_to_image(boundary, np.array([islice, dimplot], dtype=int))
        if np.issubdtype(np.asarray(annotim).dtype, np.floating):
            row, col = annotation_boundary_pixels(annotim)
        else:
            row, col = mask_boundary_pixels(annotim)

        extent = (0, histim.shape[1], histim.shape[0], 0)
        ax.imshow(histim, cmap="gray", aspect="equal", origin="upper", extent=extent)
        if row.size:
            ax.plot(
                col * pxsize[1],
                row * pxsize[0],
                linestyle="none",
                marker=".",
                markersize=2.0,
                color=(1.0, 0.8, 0.5),
                alpha=0.95,
            )
        ax.set_title(str(islice))
        ax.axis("off")

    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.02, wspace=0.05, hspace=0.08)
    return fig


def save_initial_registration_previews(
    save_path: Path,
    sample: np.ndarray,
    annotation_atlas: np.ndarray,
    transform: np.ndarray,
    *,
    boundary_atlas: np.ndarray | None = None,
) -> int:
    """Write dim{1,2,3}_initial_registration.png using the coarse similarity transform.

    Overlays follow MATLAB ``plotAnnotationComparison``: gradient edges on the warped
    annotation label volume. An optional ``boundary_atlas`` mask is only used to report
    the warped voxel count.
    """
    save_path = Path(save_path)
    volmax = float(np.quantile(sample, 0.999)) or 1.0
    sample_u8 = np.clip(255.0 * sample / volmax, 0, 255).astype(np.uint8)
    annotation_warped = warp_atlas_to_sample(
        annotation_atlas.astype(np.float32),
        transform,
        sample.shape,
        order=0,
    )
    warped_count = int(np.count_nonzero(annotation_warped > 1))
    if boundary_atlas is not None:
        boundary_warped = warp_atlas_to_sample(
            boundary_atlas.astype(np.float32),
            transform,
            sample.shape,
            order=0,
        )
        warped_count = int(np.count_nonzero(boundary_warped > 0))

    for idim in range(1, 4):
        fig = plot_annotation_comparison(sample_u8, annotation_warped, idim)
        out = save_path / f"dim{idim}_initial_registration.png"
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)

    return warped_count
