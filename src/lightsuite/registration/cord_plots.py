"""QC overlays for spinal cord registration (plotCordAnnotation.m port)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


def _validate_matching_shapes(volume: np.ndarray, annotation: np.ndarray) -> None:
    if volume.shape != annotation.shape:
        msg = (
            f"Volume and annotation shapes must match for QC plots; "
            f"got volume {volume.shape} vs annotation {annotation.shape}"
        )
        raise ValueError(msg)


def _tissue_threshold(volume: np.ndarray) -> int:
    """Background cutoff for straightened sample display volumes."""
    positive = volume[volume > 0]
    if positive.size == 0:
        return 2
    return max(2, int(np.percentile(positive, 5)))


def _annotation_boundaries(atlasim: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return row/col indices of label boundaries (plotCordAnnotation.m conv2 logic)."""
    labeled = atlasim.astype(np.float32)
    filled = ndimage.uniform_filter(labeled, size=3)
    boundary = np.round(filled) != labeled
    return np.where(boundary)


def _display_crop(
    histim: np.ndarray,
    rr: np.ndarray,
    cc: np.ndarray,
    *,
    margin: int = 12,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Crop axial/coronal/sagittal panels to the cord + annotation extent."""
    shape = histim.shape
    fg_rows: list[np.ndarray] = []
    fg_cols: list[np.ndarray] = []
    tissue = histim > 0
    if np.any(tissue):
        tr, tc = np.where(tissue)
        fg_rows.append(tr)
        fg_cols.append(tc)
    if rr.size:
        fg_rows.append(rr)
        fg_cols.append(cc)
    if not fg_rows:
        return histim, rr, cc

    rows = np.concatenate(fg_rows)
    cols = np.concatenate(fg_cols)
    r0 = max(0, int(rows.min()) - margin)
    r1 = min(shape[0], int(rows.max()) + margin + 1)
    c0 = max(0, int(cols.min()) - margin)
    c1 = min(shape[1], int(cols.max()) + margin + 1)
    cropped = histim[r0:r1, c0:c1]
    if rr.size:
        keep = (rr >= r0) & (rr < r1) & (cc >= c0) & (cc < c1)
        rr = rr[keep] - r0
        cc = cc[keep] - c0
    return cropped, rr, cc


def _select_axial_slices(
    volume: np.ndarray,
    annotation: np.ndarray,
    *,
    n_show: int,
    bg_threshold: int,
) -> np.ndarray:
    """Pick axial slices evenly spanning the full rostrocaudal extent."""
    n_length = volume.shape[2]
    if n_show <= 0:
        return np.array([], dtype=int)
    if n_show == 1:
        return np.array([n_length // 2], dtype=int)

    # Partition the full z range into bins and pick the best slice in each bin.
    edges = np.linspace(0, n_length, n_show + 1).astype(int)
    picked: list[int] = []
    for i in range(n_show):
        z0 = int(edges[i])
        z1 = int(edges[i + 1])
        if z1 <= z0:
            z1 = min(z0 + 1, n_length)
        candidates = np.arange(z0, min(z1, n_length))
        if candidates.size == 0:
            picked.append(min(z0, n_length - 1))
            continue

        best_z = int(candidates[len(candidates) // 2])
        best_score = -1
        for z in candidates:
            tissue = volume[:, :, z] > bg_threshold
            if not np.any(tissue):
                continue
            ann = annotation[:, :, z] > 0
            overlap = int(np.count_nonzero(tissue & ann))
            score = overlap if overlap > 0 else int(np.count_nonzero(tissue))
            if score > best_score:
                best_score = score
                best_z = int(z)
        picked.append(best_z)
    return np.array(picked, dtype=int)


def _cord_center_row(volume: np.ndarray, bg_threshold: int) -> int:
    """Median tissue row across sampled axial slices."""
    z_samples = np.linspace(0, volume.shape[2] - 1, min(50, volume.shape[2])).astype(int)
    rows: list[float] = []
    for z in z_samples:
        tissue = volume[:, :, z] > bg_threshold
        if np.any(tissue):
            rows.append(float(np.mean(np.where(tissue)[0])))
    if rows:
        return int(round(float(np.median(rows))))
    return volume.shape[0] // 2


def _show_slice(
    ax,
    histim: np.ndarray,
    rr: np.ndarray,
    cc: np.ndarray,
    ann_color: tuple[float, float, float],
    *,
    crop: bool = True,
) -> None:
    """Display one overlay panel with MATLAB-like axis equal + tight limits."""
    if crop:
        histim, rr, cc = _display_crop(histim, rr, cc)
    gray = histim.astype(np.float32)
    if gray.max() > 0:
        gray = gray / gray.max()
    rgb = np.stack([gray, gray, gray], axis=-1)
    if rr.size:
        # Paint boundary voxels directly so axial/coronal/sagittal overlays match.
        rgb[rr, cc, 0] = ann_color[0]
        rgb[rr, cc, 1] = ann_color[1]
        rgb[rr, cc, 2] = ann_color[2]
    ax.imshow(rgb, origin="upper", vmin=0.0, vmax=1.0)
    h, w = histim.shape
    ax.set_xlim(-0.5, w - 0.5)
    ax.set_ylim(h - 0.5, -0.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()


def save_cord_annotation_preview(
    volume: np.ndarray,
    annotation: np.ndarray,
    output_path: Path,
    *,
    title: str = "",
) -> Path:
    """Save multi-panel overlay PNG comparing sample and warped annotation."""
    _validate_matching_shapes(volume, annotation)

    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    atlas_size = annotation.shape
    n_show = 12
    n_long = 4
    bg_threshold = _tissue_threshold(volume)

    fig = plt.figure(figsize=(15, 8.5), dpi=120)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1.2], wspace=0.08, hspace=0.12)
    ax_coronal = gs[0, 0].subgridspec(1, n_long, wspace=0.05)
    ax_sagittal = gs[0, 1].subgridspec(1, n_long, wspace=0.05)
    ax_slices = gs[1, :].subgridspec(3, 4, wspace=0.02, hspace=0.02)

    ann_color = (1.0, 0.8, 0.5)
    show_slices = _select_axial_slices(volume, annotation, n_show=n_show, bg_threshold=bg_threshold)

    for ii, islice in enumerate(show_slices):
        row, col = divmod(ii, 4)
        ax = fig.add_subplot(ax_slices[row, col])
        histim = volume[:, :, islice]
        rr, cc = _annotation_boundaries(annotation[:, :, islice])
        _show_slice(ax, histim, rr, cc, ann_color)
        ax.text(
            0.95,
            0.95,
            str(islice),
            transform=ax.transAxes,
            color="white",
            ha="right",
            va="top",
            fontsize=7,
        )

    center_row = _cord_center_row(volume, bg_threshold)
    half_span = max(20, int(0.12 * atlas_size[0]))
    coronal_rows = np.clip(
        np.round(np.linspace(center_row - half_span, center_row + half_span, n_long)).astype(int),
        0,
        atlas_size[0] - 1,
    )
    for ii, islice in enumerate(coronal_rows):
        ax = fig.add_subplot(ax_coronal[ii])
        histim = np.squeeze(volume[islice, :, :]).T
        rr, cc = _annotation_boundaries(np.squeeze(annotation[islice, :, :]).T)
        _show_slice(ax, histim, rr, cc, ann_color)
        ax.set_title(str(islice), fontsize=8)

    n2 = atlas_size[1]
    sagittal_cols = np.round(np.linspace(0.25 * n2, 0.75 * n2, n_long)).astype(int)
    for ii, islice in enumerate(sagittal_cols):
        ax = fig.add_subplot(ax_sagittal[ii])
        histim = np.squeeze(volume[:, islice, :]).T
        rr, cc = _annotation_boundaries(np.squeeze(annotation[:, islice, :]).T)
        _show_slice(ax, histim, rr, cc, ann_color)
        ax.set_title(str(islice), fontsize=8)

    if title:
        fig.suptitle(title, fontsize=10)
    fig.savefig(output_path, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    return output_path
