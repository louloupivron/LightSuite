"""Per-segment anatomical heatmaps on spinal cord cross-sections.

Two backends:

* **slice** — paints registered ``annotation_registered.tiff`` slices with
  per-region metric values (fast, no extra dependencies).
* **bgh** — uses ``brainglobe-heatmap`` / brainrender meshes for publication-style
  region outlines (optional ``[viz]`` extra).
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.colors import Normalize
from skimage.segmentation import find_boundaries

from lightsuite.analysis.cord_hemisphere import HEMI_ACTIVE_VALUE
from lightsuite.analysis.viz.cord_io import CORD_METRICS, filter_cord_stats
from lightsuite.analysis.viz.cord_plots import _METRIC_LABELS, _percentile_limits

_BG_ATLAS_NAME = "allen_cord_20um"
_Z_VOXEL_UM = 20.0
_HEMISPHERE_SIDES = ("left", "right")


def _subplot_grid(n_panels: int, *, ncol_max: int = 4) -> tuple[int, int]:
    ncol = min(max(1, ncol_max), max(1, n_panels))
    nrow = (n_panels + ncol - 1) // ncol
    return nrow, ncol


def _metric_label(metric: str) -> str:
    return _METRIC_LABELS.get(metric, metric)


def load_segments_table(segments_csv: Path) -> pd.DataFrame:
    path = Path(segments_csv).expanduser().resolve()
    if not path.is_file():
        msg = f"Segments.csv not found: {path}"
        raise FileNotFoundError(msg)
    segments = pd.read_csv(path)
    if "Segment" not in segments.columns or "Start" not in segments.columns or "End" not in segments.columns:
        msg = f"Segments.csv must contain Segment, Start, End columns: {path}"
        raise ValueError(msg)
    segments["Segment"] = segments["Segment"].astype(str)
    segments["z_mid"] = ((segments["Start"].astype(int) + segments["End"].astype(int)) // 2).astype(int)
    return segments


def segment_z_index(segments_df: pd.DataFrame, segment: str) -> int:
    row = segments_df.loc[segments_df["Segment"] == str(segment)]
    if row.empty:
        msg = f"Unknown segment {segment!r}. Expected a label from Segments.csv."
        raise ValueError(msg)
    return int(row.iloc[0]["z_mid"])


def segment_z_um(segments_df: pd.DataFrame, segment: str, *, z_voxel_um: float = _Z_VOXEL_UM) -> float:
    return float(segment_z_index(segments_df, segment) * z_voxel_um)


@lru_cache(maxsize=1)
def _brainglobe_valid_acronyms() -> set[str]:
    from brainglobe_atlasapi.bg_atlas import BrainGlobeAtlas

    atlas = BrainGlobeAtlas(_BG_ATLAS_NAME)
    return {str(item["acronym"]) for item in atlas.structures_list if item.get("acronym")}


def segment_region_values(
    stats_df: pd.DataFrame,
    segment: str,
    *,
    channel: int | str,
    metric: str,
    rollup_level: str = "region",
    hemisphere: str | None = None,
) -> tuple[dict[int, float], dict[str, float]]:
    """Return per-region values for one segment as ``{id: value}`` and ``{acronym: value}``."""
    if metric not in CORD_METRICS:
        msg = f"Unknown metric {metric!r}. Expected one of {CORD_METRICS}."
        raise ValueError(msg)

    work = filter_cord_stats(
        stats_df,
        channel=channel,
        metric=metric,
        rollup_level=rollup_level,
        hemisphere=hemisphere,
    )
    work = work[work["segment"].astype(str) == str(segment)].copy()
    if work.empty and hemisphere is None:
        work = filter_cord_stats(
            stats_df,
            channel=channel,
            metric=metric,
            rollup_level=rollup_level,
        )
        work = work[work["segment"].astype(str) == str(segment)].copy()
        if not work.empty and "hemisphere" in work.columns and work["hemisphere"].astype(str).nunique() > 1:
            group_cols = ["parcellation_index", "acronym"]
            if metric == "cell_count":
                work = work.groupby(group_cols, sort=False)["value"].sum().reset_index()
            else:
                work = work.groupby(group_cols, sort=False)["value"].mean().reset_index()

    if work.empty:
        return {}, {}

    by_id = {
        int(row.parcellation_index): float(row.value)
        for row in work.itertuples()
        if pd.notna(row.value) and float(row.value) != 0.0
    }
    by_acronym = {
        str(row.acronym): float(row.value)
        for row in work.itertuples()
        if pd.notna(row.value) and float(row.value) != 0.0
    }
    return by_id, by_acronym


def _collect_bilateral_segment_values(
    stats_df: pd.DataFrame,
    segments: list[str],
    *,
    channel: int | str,
    metric: str,
    rollup_level: str,
) -> dict[str, dict[str, dict[int, float]]]:
    """Per-segment left/right region-id value maps."""
    out: dict[str, dict[str, dict[int, float]]] = {}
    for segment in segments:
        side_values: dict[str, dict[int, float]] = {}
        for side in _HEMISPHERE_SIDES:
            by_id, _by_acronym = segment_region_values(
                stats_df,
                segment,
                channel=channel,
                metric=metric,
                rollup_level=rollup_level,
                hemisphere=side,
            )
            side_values[side] = by_id
        if side_values["left"] or side_values["right"]:
            out[segment] = side_values
    return out


def _shared_color_limits_from_ids(
    value_maps: list[dict[int, float]],
    *,
    vmin: float | None,
    vmax: float | None,
) -> tuple[float, float]:
    if vmin is not None and vmax is not None:
        return float(vmin), float(vmax)
    all_vals = [value for value_map in value_maps for value in value_map.values() if np.isfinite(value)]
    if not all_vals:
        return 0.0, 1.0
    auto_vmin, auto_vmax = _percentile_limits(np.asarray(all_vals, dtype=float))
    return float(vmin if vmin is not None else auto_vmin), float(vmax if vmax is not None else auto_vmax)


def _hemisphere_side_masks(
    hemisphere_slice: np.ndarray,
    *,
    flip: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Return boolean left/right masks using the Fiederling convention (255=right by default)."""
    active = hemisphere_slice == HEMI_ACTIVE_VALUE
    inactive = hemisphere_slice == 0
    if flip:
        return active, inactive
    return inactive, active


def _require_hemisphere_volume(hemisphere_volume_path: Path | None) -> None:
    if hemisphere_volume_path is None or not Path(hemisphere_volume_path).is_file():
        msg = (
            "hemisphere_registered.tiff is required for left/right anatomy layouts. "
            "Run region-stats with split_hemispheres: true."
        )
        raise FileNotFoundError(msg)


def _draw_slice_panel(
    ax: plt.Axes,
    annotation_slice: np.ndarray,
    values_by_id: dict[int, float],
    *,
    cmap: str,
    vmin: float,
    vmax: float,
    title: str | None = None,
    draw_outlines: bool = True,
) -> None:
    rgba = _region_value_rgba(
        annotation_slice,
        values_by_id,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.imshow(rgba, origin="upper", aspect="equal")
    if draw_outlines and np.any(annotation_slice > 0):
        boundaries = find_boundaries(annotation_slice > 0, mode="outer")
        ys, xs = np.where(boundaries)
        ax.scatter(xs, ys, s=0.15, c="k", alpha=0.35, linewidths=0)
    if title is not None:
        ax.set_title(title, fontsize=10, fontweight="bold")
    ax.axis("off")


def _region_value_rgba_bilateral(
    annotation_slice: np.ndarray,
    hemisphere_slice: np.ndarray,
    values_left: dict[int, float],
    values_right: dict[int, float],
    *,
    cmap: str,
    vmin: float,
    vmax: float,
    flip: bool = False,
) -> np.ndarray:
    """Paint left hemivoxels from left stats and right hemivoxels from right stats."""
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    colormap = plt.get_cmap(cmap)
    left_mask, right_mask = _hemisphere_side_masks(hemisphere_slice, flip=flip)
    rgba = np.zeros((*annotation_slice.shape, 4), dtype=float)
    tissue_ids = np.unique(annotation_slice[annotation_slice > 0])
    for region_id in tissue_ids:
        rid = int(region_id)
        region_mask = annotation_slice == rid
        for side_mask, values in ((left_mask, values_left), (right_mask, values_right)):
            if rid not in values:
                continue
            mask = region_mask & side_mask
            if not np.any(mask):
                continue
            rgba[mask] = colormap(norm(float(values[rid])))
    return rgba


def _add_shared_colorbar(
    fig: plt.Figure,
    axes: np.ndarray,
    *,
    cmap: str,
    vmin: float,
    vmax: float,
    metric_label: str,
) -> None:
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    cbar.set_label(metric_label, fontsize=9)


def _add_row_segment_labels(
    axes: np.ndarray,
    segments: list[str],
    *,
    label_col: int = 0,
) -> None:
    """Place segment names to the left of each row (ylabels are hidden when axes are off)."""
    for row_idx, segment in enumerate(segments):
        axes[row_idx, label_col].text(
            -0.14,
            0.5,
            segment,
            transform=axes[row_idx, label_col].transAxes,
            ha="right",
            va="center",
            fontsize=11,
            fontweight="bold",
            clip_on=False,
        )


def plot_cord_segment_anatomy_slice_hemisphere_panel(
    stats_df: pd.DataFrame,
    *,
    segments: list[str],
    annotation_path: Path,
    segments_csv: Path,
    hemisphere_volume_path: Path,
    channel: int | str,
    metric: str = "median_intensity",
    rollup_level: str = "region",
    hemisphere_flip: bool = False,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    cmap: str = "inferno",
    vmin: float | None = None,
    vmax: float | None = None,
    draw_outlines: bool = True,
) -> tuple[plt.Figure, dict[str, dict[str, dict[int, float]]]]:
    """Left/right columns per segment using registered annotation slices."""
    if not segments:
        msg = "At least one segment is required."
        raise ValueError(msg)
    _require_hemisphere_volume(hemisphere_volume_path)

    segments_df = load_segments_table(segments_csv)
    segment_values = _collect_bilateral_segment_values(
        stats_df,
        segments,
        channel=channel,
        metric=metric,
        rollup_level=rollup_level,
    )
    if not segment_values:
        msg = f"No left/right {metric!r} data for segments {segments!r}."
        raise ValueError(msg)

    ordered_segments = [segment for segment in segments if segment in segment_values]
    all_maps = [
        segment_values[segment][side]
        for segment in ordered_segments
        for side in _HEMISPHERE_SIDES
    ]
    color_vmin, color_vmax = _shared_color_limits_from_ids(all_maps, vmin=vmin, vmax=vmax)

    annotation = tifffile.memmap(annotation_path)
    hemisphere_volume = tifffile.memmap(hemisphere_volume_path)

    nrow = len(ordered_segments)
    fig, axes = plt.subplots(nrow, 2, figsize=(6.8, 3.4 * nrow), squeeze=False)
    if nrow == 1:
        axes = axes.reshape(1, 2)

    plotted: dict[str, dict[str, dict[int, float]]] = {}
    for row_idx, segment in enumerate(ordered_segments):
        z_idx = segment_z_index(segments_df, segment)
        ann_slice = np.asarray(annotation[:, :, z_idx])
        hem_slice = np.asarray(hemisphere_volume[:, :, z_idx])
        for col_idx, side in enumerate(_HEMISPHERE_SIDES):
            ax = axes[row_idx, col_idx]
            masked_ann = _mask_annotation_for_side(
                ann_slice,
                hem_slice,
                side,
                flip=hemisphere_flip,
            )
            values_by_id = segment_values[segment][side]
            _draw_slice_panel(
                ax,
                masked_ann,
                values_by_id,
                cmap=cmap,
                vmin=color_vmin,
                vmax=color_vmax,
                draw_outlines=draw_outlines,
            )
        plotted[segment] = segment_values[segment]

    _add_row_segment_labels(axes, ordered_segments)

    if nrow > 0:
        axes[0, 0].set_title("Left", fontsize=10, fontweight="bold", color="#4e79a7")
        axes[0, 1].set_title("Right", fontsize=10, fontweight="bold", color="#e15759")

    metric_label = _metric_label(metric)
    fig.suptitle(
        title or f"Segment anatomy — {metric_label} (left vs right)",
        fontweight="bold",
        y=0.995,
    )
    _add_shared_colorbar(fig, axes, cmap=cmap, vmin=color_vmin, vmax=color_vmax, metric_label=metric_label)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)
    return fig, plotted


def plot_cord_segment_anatomy_slice_hemisphere_composite(
    stats_df: pd.DataFrame,
    *,
    segments: list[str],
    annotation_path: Path,
    segments_csv: Path,
    hemisphere_volume_path: Path,
    channel: int | str,
    metric: str = "median_intensity",
    rollup_level: str = "region",
    hemisphere_flip: bool = False,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    cmap: str = "inferno",
    vmin: float | None = None,
    vmax: float | None = None,
    ncol_max: int = 4,
    draw_outlines: bool = True,
) -> tuple[plt.Figure, dict[str, dict[str, dict[int, float]]]]:
    """One cross-section per segment with left/right stats on the same slice."""
    if not segments:
        msg = "At least one segment is required."
        raise ValueError(msg)
    _require_hemisphere_volume(hemisphere_volume_path)

    segments_df = load_segments_table(segments_csv)
    segment_values = _collect_bilateral_segment_values(
        stats_df,
        segments,
        channel=channel,
        metric=metric,
        rollup_level=rollup_level,
    )
    if not segment_values:
        msg = f"No left/right {metric!r} data for segments {segments!r}."
        raise ValueError(msg)

    ordered_segments = [segment for segment in segments if segment in segment_values]
    all_maps = [
        segment_values[segment][side]
        for segment in ordered_segments
        for side in _HEMISPHERE_SIDES
    ]
    color_vmin, color_vmax = _shared_color_limits_from_ids(all_maps, vmin=vmin, vmax=vmax)

    annotation = tifffile.memmap(annotation_path)
    hemisphere_volume = tifffile.memmap(hemisphere_volume_path)
    nrow, ncol = _subplot_grid(len(ordered_segments), ncol_max=ncol_max)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.6 * nrow), squeeze=False)

    plotted: dict[str, dict[str, dict[int, float]]] = {}
    for ax, segment in zip(axes.ravel(), ordered_segments, strict=False):
        z_idx = segment_z_index(segments_df, segment)
        ann_slice = np.asarray(annotation[:, :, z_idx])
        hem_slice = np.asarray(hemisphere_volume[:, :, z_idx])
        values_left = segment_values[segment]["left"]
        values_right = segment_values[segment]["right"]
        rgba = _region_value_rgba_bilateral(
            ann_slice,
            hem_slice,
            values_left,
            values_right,
            cmap=cmap,
            vmin=color_vmin,
            vmax=color_vmax,
            flip=hemisphere_flip,
        )
        ax.imshow(rgba, origin="upper", aspect="equal")
        if draw_outlines and np.any(ann_slice > 0):
            boundaries = find_boundaries(ann_slice > 0, mode="outer")
            ys, xs = np.where(boundaries)
            ax.scatter(xs, ys, s=0.15, c="k", alpha=0.35, linewidths=0)
        ax.set_title(segment, fontsize=10, fontweight="bold")
        ax.axis("off")
        plotted[segment] = segment_values[segment]

    for ax in axes.ravel()[len(ordered_segments) :]:
        ax.axis("off")

    metric_label = _metric_label(metric)
    fig.suptitle(
        title or f"Segment anatomy — {metric_label} (bilateral composite)",
        fontweight="bold",
        y=0.98,
    )
    _add_shared_colorbar(fig, axes, cmap=cmap, vmin=color_vmin, vmax=color_vmax, metric_label=metric_label)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)
    return fig, plotted


def _collect_segment_values(
    stats_df: pd.DataFrame,
    segments: list[str],
    *,
    channel: int | str,
    metric: str,
    rollup_level: str,
    hemisphere: str | None,
) -> dict[str, tuple[dict[int, float], dict[str, float]]]:
    out: dict[str, tuple[dict[int, float], dict[str, float]]] = {}
    for segment in segments:
        by_id, by_acronym = segment_region_values(
            stats_df,
            segment,
            channel=channel,
            metric=metric,
            rollup_level=rollup_level,
            hemisphere=hemisphere,
        )
        if by_id or by_acronym:
            out[segment] = (by_id, by_acronym)
    return out


def _mask_annotation_for_side(
    annotation_slice: np.ndarray,
    hemisphere_slice: np.ndarray,
    side: str,
    *,
    flip: bool = False,
) -> np.ndarray:
    left_mask, right_mask = _hemisphere_side_masks(hemisphere_slice, flip=flip)
    keep = left_mask if side == "left" else right_mask
    masked = annotation_slice.copy()
    masked[~keep] = 0
    return masked


def _shared_color_limits(
    segment_values: dict[str, tuple[dict[int, float], dict[str, float]]],
    *,
    vmin: float | None,
    vmax: float | None,
) -> tuple[float, float]:
    if vmin is not None and vmax is not None:
        return float(vmin), float(vmax)
    all_vals = [
        value
        for by_id, _by_acronym in segment_values.values()
        for value in by_id.values()
        if np.isfinite(value)
    ]
    if not all_vals:
        return 0.0, 1.0
    auto_vmin, auto_vmax = _percentile_limits(np.asarray(all_vals, dtype=float))
    return float(vmin if vmin is not None else auto_vmin), float(vmax if vmax is not None else auto_vmax)


def _apply_hemisphere_mask(
    annotation_slice: np.ndarray,
    hemisphere_slice: np.ndarray | None,
    *,
    hemisphere: str | None,
) -> np.ndarray:
    if hemisphere is None or hemisphere_slice is None:
        return annotation_slice
    side = str(hemisphere).strip().lower()
    if side not in {"left", "right"}:
        msg = f"hemisphere must be 'left' or 'right', got {hemisphere!r}"
        raise ValueError(msg)
    keep_value = 255 if side == "right" else 0
    masked = annotation_slice.copy()
    masked[hemisphere_slice != keep_value] = 0
    return masked


def _region_value_rgba(
    annotation_slice: np.ndarray,
    values_by_id: dict[int, float],
    *,
    cmap: str,
    vmin: float,
    vmax: float,
) -> np.ndarray:
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    colormap = plt.get_cmap(cmap)
    rgba = np.zeros((*annotation_slice.shape, 4), dtype=float)
    rgba[..., 3] = 0.0
    for region_id, value in values_by_id.items():
        if not np.isfinite(value):
            continue
        mask = annotation_slice == int(region_id)
        if not np.any(mask):
            continue
        color = colormap(norm(float(value)))
        rgba[mask] = color
    return rgba


def plot_cord_segment_anatomy_slice(
    stats_df: pd.DataFrame,
    *,
    segments: list[str],
    annotation_path: Path,
    segments_csv: Path,
    hemisphere_volume_path: Path | None = None,
    channel: int | str,
    metric: str = "median_intensity",
    rollup_level: str = "region",
    hemisphere: str | None = None,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    cmap: str = "inferno",
    vmin: float | None = None,
    vmax: float | None = None,
    ncol_max: int = 4,
    draw_outlines: bool = True,
) -> tuple[plt.Figure, dict[str, dict[int, float]]]:
    """Anatomical heatmap panels from registered annotation slices.

    Each panel shows one user-selected segment. Regions are filled with the
    requested metric value (``median_intensity``, ``cell_count``, etc.) on the
    mid-slice of the segment's Z range in ``Segments.csv``.
    """
    if not segments:
        msg = "At least one segment is required."
        raise ValueError(msg)

    segments_df = load_segments_table(segments_csv)
    segment_values = _collect_segment_values(
        stats_df,
        segments,
        channel=channel,
        metric=metric,
        rollup_level=rollup_level,
        hemisphere=hemisphere,
    )
    if not segment_values:
        msg = f"No {metric!r} data for segments {segments!r}."
        raise ValueError(msg)

    color_vmin, color_vmax = _shared_color_limits(segment_values, vmin=vmin, vmax=vmax)
    annotation = tifffile.memmap(annotation_path)
    hemisphere_volume = tifffile.memmap(hemisphere_volume_path) if hemisphere_volume_path else None

    ordered_segments = [segment for segment in segments if segment in segment_values]
    nrow, ncol = _subplot_grid(len(ordered_segments), ncol_max=ncol_max)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.6 * nrow), squeeze=False)

    plotted: dict[str, dict[int, float]] = {}
    for ax, segment in zip(axes.ravel(), ordered_segments, strict=False):
        z_idx = segment_z_index(segments_df, segment)
        ann_slice = np.asarray(annotation[:, :, z_idx])
        hem_slice = np.asarray(hemisphere_volume[:, :, z_idx]) if hemisphere_volume is not None else None
        ann_slice = _apply_hemisphere_mask(ann_slice, hem_slice, hemisphere=hemisphere)
        values_by_id, _values_by_acronym = segment_values[segment]
        rgba = _region_value_rgba(
            ann_slice,
            values_by_id,
            cmap=cmap,
            vmin=color_vmin,
            vmax=color_vmax,
        )
        ax.imshow(rgba, origin="upper", aspect="equal")
        if draw_outlines and np.any(ann_slice > 0):
            boundaries = find_boundaries(ann_slice > 0, mode="outer")
            ys, xs = np.where(boundaries)
            ax.scatter(xs, ys, s=0.15, c="k", alpha=0.35, linewidths=0)
        ax.set_title(segment, fontsize=10, fontweight="bold")
        ax.axis("off")
        plotted[segment] = values_by_id

    for ax in axes.ravel()[len(ordered_segments) :]:
        ax.axis("off")

    metric_label = _metric_label(metric)
    fig.suptitle(
        title or f"Segment anatomy — {metric_label} (annotation slices)",
        fontweight="bold",
        y=0.98,
    )
    norm = Normalize(vmin=color_vmin, vmax=color_vmax)
    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    cbar.set_label(metric_label, fontsize=9)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, plotted


def plot_cord_segment_anatomy_bgh(
    stats_df: pd.DataFrame,
    *,
    segments: list[str],
    segments_csv: Path,
    channel: int | str,
    metric: str = "median_intensity",
    rollup_level: str = "region",
    hemisphere: str | None = None,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    cmap: str = "inferno",
    vmin: float | None = None,
    vmax: float | None = None,
    ncol_max: int = 3,
    thickness_um: float = 400.0,
    atlas_name: str = _BG_ATLAS_NAME,
) -> tuple[plt.Figure, dict[str, dict[str, float]]]:
    """Anatomical heatmap panels via brainglobe-heatmap (brainrender meshes).

    Requires the optional ``viz`` extra: ``uv sync --extra viz``.
    """
    try:
        import brainglobe_heatmap as bgh
    except ImportError as exc:
        msg = (
            "brainglobe-heatmap is not installed. "
            "Install with: uv sync --extra viz"
        )
        raise ImportError(msg) from exc

    if not segments:
        msg = "At least one segment is required."
        raise ValueError(msg)

    segments_df = load_segments_table(segments_csv)
    segment_values = _collect_segment_values(
        stats_df,
        segments,
        channel=channel,
        metric=metric,
        rollup_level=rollup_level,
        hemisphere=hemisphere,
    )
    if not segment_values:
        msg = f"No {metric!r} data for segments {segments!r}."
        raise ValueError(msg)

    valid_acronyms = _brainglobe_valid_acronyms()
    color_vmin, color_vmax = _shared_color_limits(segment_values, vmin=vmin, vmax=vmax)
    ordered_segments = [segment for segment in segments if segment in segment_values]
    nrow, ncol = _subplot_grid(len(ordered_segments), ncol_max=ncol_max)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.0 * nrow), squeeze=False)

    plotted: dict[str, dict[str, float]] = {}
    bgh_hemisphere = hemisphere if hemisphere in {"left", "right"} else "both"

    for ax, segment in zip(axes.ravel(), ordered_segments, strict=False):
        _values_by_id, values_by_acronym = segment_values[segment]
        bgh_values = {
            acronym: value
            for acronym, value in values_by_acronym.items()
            if acronym in valid_acronyms
        }
        if not bgh_values:
            ax.text(0.5, 0.5, "No BG-compatible regions", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(segment, fontsize=10)
            ax.axis("off")
            continue

        position_um = segment_z_um(segments_df, segment)
        heatmap = bgh.Heatmap(
            bgh_values,
            position=position_um,
            orientation="frontal",
            thickness=thickness_um,
            atlas_name=atlas_name,
            format="2D",
            cmap=cmap,
            vmin=color_vmin,
            vmax=color_vmax,
            hemisphere=bgh_hemisphere,
            interactive=False,
            check_latest=False,
        )
        heatmap.plot_subplot(fig, ax, hide_axes=True, show_cbar=False)
        ax.set_title(segment, fontsize=10, fontweight="bold")
        plotted[segment] = bgh_values

    for ax in axes.ravel()[len(ordered_segments) :]:
        ax.axis("off")

    metric_label = _metric_label(metric)
    fig.suptitle(
        title or f"Segment anatomy — {metric_label} (brainrender)",
        fontweight="bold",
        y=0.98,
    )
    norm = Normalize(vmin=color_vmin, vmax=color_vmax)
    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap), norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    cbar.set_label(metric_label, fontsize=9)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, plotted


__all__ = [
    "load_segments_table",
    "plot_cord_segment_anatomy_bgh",
    "plot_cord_segment_anatomy_slice",
    "plot_cord_segment_anatomy_slice_hemisphere_composite",
    "plot_cord_segment_anatomy_slice_hemisphere_panel",
    "segment_region_values",
    "segment_z_index",
    "segment_z_um",
]
