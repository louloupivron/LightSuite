"""Registration preview figures (plotAnnotationComparison.m port)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import find_boundaries

from lightsuite.atlas.display import (
    canonical_view_name,
    canonical_view_slice,
    cut_axis_for_plot_dim,
)
from lightsuite.gui.slices import volume_index_to_image
from lightsuite.registration.warp import warp_atlas_to_sample

# Re-export for tests and callers that imported from plots historically.
reorient_registration_slice = canonical_view_slice
prepare_registration_slice = canonical_view_slice

__all__ = [
    "boundary_volume_from_annotation",
    "canonical_view_slice",
    "plot_annotation_comparison",
    "prepare_registration_slice",
    "reorient_registration_slice",
    "save_initial_registration_previews",
    "save_registration_stage_previews",
]


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
    """Return row/column indices of region edges (integer labels after transformix)."""
    labels = np.rint(annotation_slice).astype(np.int32)
    if not np.any(labels > 1):
        return np.array([], dtype=int), np.array([], dtype=int)
    edges = find_boundaries(labels, mode="inner")
    row, col = np.nonzero(edges & (labels > 1))
    return row, col


def _comparison_axis_limits(
    shape: tuple[int, int],
    pxsize: tuple[float, float],
    row: np.ndarray,
    col: np.ndarray,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Expand axis limits so warped atlas overlays remain visible when shifted."""
    height, width = shape
    row_px, col_px = pxsize
    x_min = 0.0
    x_max = float(width) * col_px
    y_min = 0.0
    y_max = float(height) * row_px
    if row.size:
        x_min = min(x_min, float(col.min()) * col_px)
        x_max = max(x_max, float(col.max()) * col_px)
        y_min = min(y_min, float(row.min()) * row_px)
        y_max = max(y_max, float(row.max()) * row_px)
    pad_x = 0.05 * max(x_max - x_min, 1e-6)
    pad_y = 0.05 * max(y_max - y_min, 1e-6)
    return (x_min - pad_x, x_max + pad_x), (y_max + pad_y, y_min - pad_y)


def plot_annotation_comparison(
    volume: np.ndarray,
    boundary: np.ndarray,
    dimplot: int,
    *,
    n_show: int = 8,
    pxsize: tuple[float, float] = (1.0, 1.0),
    atlas_provider: str = "allen",
) -> plt.Figure:
    """Plot sample slices with warped atlas boundary mask overlaid (sample only)."""
    if dimplot not in {1, 2, 3}:
        msg = f"dimplot must be 1, 2, or 3, got {dimplot}"
        raise ValueError(msg)

    cut_axis = cut_axis_for_plot_dim(atlas_provider, dimplot)
    ny = volume.shape[cut_axis - 1]
    ishow = np.round(np.linspace(0.15 * ny, 0.85 * ny, n_show)).astype(int)
    ishow = np.clip(ishow, 1, ny)

    fig, axes = plt.subplots(2, n_show // 2, figsize=(17, 9))
    axes = np.atleast_1d(axes).ravel()

    for ii, islice in enumerate(ishow):
        ax = axes[ii]
        chooserow = np.array([islice, cut_axis], dtype=int)
        histim = canonical_view_slice(
            volume_index_to_image(volume, chooserow),
            atlas_provider=atlas_provider,
            cut_axis=cut_axis,
        )
        annotim = canonical_view_slice(
            volume_index_to_image(boundary, chooserow),
            atlas_provider=atlas_provider,
            cut_axis=cut_axis,
        )
        if np.issubdtype(np.asarray(annotim).dtype, np.integer):
            row, col = annotation_boundary_pixels(annotim.astype(np.float32))
        elif np.issubdtype(np.asarray(annotim).dtype, np.floating):
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
                markersize=0.3,
                color=(1.0, 0.8, 0.5),
                alpha=0.95,
            )
        xlim, ylim = _comparison_axis_limits(histim.shape, pxsize, row, col)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        view = canonical_view_name(dimplot)
        ax.set_title(f"{islice} · {view}")
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
    atlas_provider: str = "allen",
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
        fig = plot_annotation_comparison(
            sample_u8,
            annotation_warped,
            idim,
            atlas_provider=atlas_provider,
        )
        out = save_path / f"dim{idim}_initial_registration.png"
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)

    return warped_count


def save_registration_stage_previews(
    save_path: Path,
    sample_name: str,
    sample_u8: np.ndarray,
    annotation_in_sample_space: np.ndarray,
    stage: str,
    *,
    atlas_provider: str = "allen",
) -> None:
    """Write ``{name}_dim{1,2,3}_{stage}.png`` (MATLAB ``plotAnnotationComparison`` style).

    ``stage`` is e.g. ``affine_registration`` or ``bspline_registration``. The annotation
    volume must already be in sample voxel space (``avaffine`` / ``avreg``).
    """
    save_path = Path(save_path)
    ann = np.asarray(annotation_in_sample_space, dtype=np.float32)
    for dimplot in range(1, 4):
        fig = plot_annotation_comparison(
            sample_u8, ann, dimplot, atlas_provider=atlas_provider
        )
        out = save_path / f"{sample_name}_dim{dimplot}_{stage}.png"
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)
