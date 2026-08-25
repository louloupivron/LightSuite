"""Paper-style segment × lamina/WM heatmaps from spinal cord region stats."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import volume_registered_dir

MetricName = Literal[
    "median_intensity",
    "relative_median_intensity",
    "cell_count",
    "cell_density",
    "std",
]

PAPER_STRUCTURE_ACRONYMS: tuple[str, ...] = (
    "Lamina_I",
    "Lamina_II",
    "Lamina_III",
    "Lamina_IV",
    "Lamina_V",
    "Lamina_VI",
    "Lamina_VII",
    "Lamina_VIII",
    "Lamina_IX",
    "Lamina_X",
    "df",
    "lf",
    "vf",
)

_ROMAN_BY_ACRONYM: dict[str, str] = {
    "Lamina_I": "I",
    "Lamina_II": "II",
    "Lamina_III": "III",
    "Lamina_IV": "IV",
    "Lamina_V": "V",
    "Lamina_VI": "VI",
    "Lamina_VII": "VII",
    "Lamina_VIII": "VIII",
    "Lamina_IX": "IX",
    "Lamina_X": "X",
    "df": "df",
    "lf": "lf",
    "vf": "vf",
}

_METRIC_LABELS: dict[str, str] = {
    "median_intensity": "Median intensity",
    "relative_median_intensity": "Relative median intensity",
    "cell_count": "Cell count",
    "cell_density": "Cell density (cells/mm³)",
    "std": "Intensity std",
}

NO_DATA_COLOR = "#b8b8b8"

_SUM_HEMISPHERE_METRICS = frozenset({"cell_count", "volume_mm3"})
_WEIGHTED_HEMISPHERE_METRICS = frozenset(
    {"median_intensity", "relative_median_intensity", "std", "cell_density"}
)


def metric_label(metric: str) -> str:
    return _METRIC_LABELS.get(metric, metric.replace("_", " "))


def discover_cord_plot_options(df: pd.DataFrame) -> dict[str, object]:
    """Summarize metrics, channels, and hemispheres available in a stats table."""
    metrics = sorted(df["metric"].astype(str).unique())
    channels_by_metric: dict[str, list[str]] = {}
    for metric in metrics:
        subset = df[df["metric"].astype(str) == metric]
        channels_by_metric[metric] = sorted(subset["channel"].astype(str).unique())
    rollup_levels = sorted(df["rollup_level"].astype(str).unique())
    return {
        "metrics": metrics,
        "channels_by_metric": channels_by_metric,
        "hemispheres": resolve_hemisphere_options(df),
        "rollup_levels": rollup_levels,
    }


def resolve_hemisphere_options(df: pd.DataFrame) -> list[str]:
    """Return hemisphere choices for plotting, including synthetic ``whole`` when split."""
    present = sorted({str(item) for item in df["hemisphere"].astype(str).unique()})
    if "whole" in present:
        return ["whole", *[item for item in present if item != "whole"]]
    if "left" in present and "right" in present:
        return ["whole", "left", "right"]
    return present


def _metric_lookup(
    frame: pd.DataFrame,
    *,
    group_cols: list[str],
) -> dict[tuple[object, ...], float]:
    lookup: dict[tuple[object, ...], float] = {}
    for _, row in frame.iterrows():
        key = tuple(row[col] for col in group_cols)
        lookup[key] = float(row["value"])
    return lookup


def _combine_hemisphere_value(
    metric: str,
    values: dict[str, float],
    volumes: dict[str, float],
) -> float:
    present = {side: value for side, value in values.items() if np.isfinite(value)}
    if not present:
        return np.nan
    if metric in _SUM_HEMISPHERE_METRICS:
        return float(sum(present.values()))
    if metric in _WEIGHTED_HEMISPHERE_METRICS:
        weights = {
            side: volumes[side]
            for side in present
            if side in volumes and np.isfinite(volumes[side]) and volumes[side] > 0
        }
        if weights:
            total_weight = float(sum(weights.values()))
            return float(
                sum(present[side] * weights[side] for side in weights) / total_weight
            )
    return float(np.mean(list(present.values())))


def merge_split_hemisphere_stats(
    df: pd.DataFrame,
    *,
    metric: str,
    rollup_level: str,
    channel: int | str,
) -> pd.DataFrame:
    """Combine left/right tidy rows into synthetic ``whole`` rows for one metric."""
    channel_str = str(channel)
    group_cols = [
        "sample",
        "channel",
        "atlas",
        "parcellation_index",
        "acronym",
        "name",
        "structure",
        "division",
        "segment",
        "rollup_level",
    ]
    metric_df = df[
        (df["metric"] == metric)
        & (df["rollup_level"] == rollup_level)
        & (df["channel"].astype(str) == channel_str)
        & (df["hemisphere"].isin(["left", "right"]))
    ]
    if metric_df.empty:
        return pd.DataFrame(columns=df.columns)

    volume_df = df[
        (df["metric"] == "volume_mm3")
        & (df["rollup_level"] == rollup_level)
        & (df["channel"].astype(str) == channel_str)
        & (df["hemisphere"].isin(["left", "right"]))
    ]
    volume_lookup = _metric_lookup(
        volume_df,
        group_cols=[*group_cols, "hemisphere"],
    )

    records: list[dict[str, object]] = []
    for group_key, group in metric_df.groupby(group_cols, dropna=False):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        base = dict(zip(group_cols, group_key, strict=True))
        values = {
            str(row.hemisphere): float(row.value)
            for row in group.itertuples(index=False)
        }
        volumes = {
            side: volume_lookup.get((*group_key, side), np.nan)
            for side in ("left", "right")
        }
        combined = _combine_hemisphere_value(metric, values, volumes)
        records.append(
            {
                **base,
                "hemisphere": "whole",
                "metric": metric,
                "value": combined,
            }
        )
    return pd.DataFrame(records).reindex(columns=df.columns)


def _stats_rows_for_hemisphere(
    df: pd.DataFrame,
    *,
    metric: str,
    rollup_level: str,
    hemisphere: str,
    channel: int | str,
) -> pd.DataFrame:
    channel_str = str(channel)
    hemi = str(hemisphere)
    if hemi == "whole":
        whole_rows = df[
            (df["metric"] == metric)
            & (df["rollup_level"] == rollup_level)
            & (df["channel"].astype(str) == channel_str)
            & (df["hemisphere"].astype(str) == "whole")
        ]
        if not whole_rows.empty:
            return whole_rows
        return merge_split_hemisphere_stats(
            df,
            metric=metric,
            rollup_level=rollup_level,
            channel=channel,
        )
    return df[
        (df["metric"] == metric)
        & (df["rollup_level"] == rollup_level)
        & (df["channel"].astype(str) == channel_str)
        & (df["hemisphere"].astype(str) == hemi)
    ]


def load_segment_order(atlas_dir: Path) -> list[str]:
    """Return rostrocaudal segment labels from ``Segments.csv``."""
    segments_path = atlas_dir.expanduser() / "Segments.csv"
    if not segments_path.is_file():
        msg = f"Missing {segments_path}"
        raise FileNotFoundError(msg)
    segments = pd.read_csv(segments_path)["Segment"].astype(str).tolist()
    return segments


def default_paper_segments(atlas_dir: Path, *, through: str = "Co2") -> list[str]:
    """Segments from C1 through ``through`` (inclusive), in rostrocaudal order."""
    segments = load_segment_order(atlas_dir)
    if through not in segments:
        msg = f"Segment {through!r} not found in Segments.csv"
        raise ValueError(msg)
    end = segments.index(through) + 1
    return segments[:end]


def structure_acronym_to_label(acronym: str) -> str:
    """Map structure-roll-up acronyms to paper-style column labels."""
    return _ROMAN_BY_ACRONYM.get(str(acronym), str(acronym))


def parse_segment_range(text: str, *, atlas_dir: Path) -> list[str]:
    """Parse ``C1:Co2`` or comma-separated segment lists."""
    raw = str(text).strip()
    if not raw:
        return default_paper_segments(atlas_dir)
    if ":" in raw:
        start, end = raw.split(":", 1)
        order = load_segment_order(atlas_dir)
        if start not in order or end not in order:
            msg = f"Invalid segment range {raw!r}; expected labels from Segments.csv"
            raise ValueError(msg)
        i0, i1 = order.index(start), order.index(end)
        if i0 > i1:
            i0, i1 = i1, i0
        return order[i0 : i1 + 1]
    labels = [part.strip() for part in raw.split(",") if part.strip()]
    return labels


def resolve_cord_rollup_columns(
    rollup_level: str,
    work_df: pd.DataFrame,
    *,
    regions_df: pd.DataFrame | None = None,
) -> tuple[list[str], list[str]]:
    """Return (column_acronyms, display_labels) for a given rollup_level."""
    norm = str(rollup_level).strip().lower()
    if norm == "structure":
        return list(PAPER_STRUCTURE_ACRONYMS), [structure_acronym_to_label(col) for col in PAPER_STRUCTURE_ACRONYMS]
    if norm == "division":
        cols = ["GM", "WM"]
        return cols, cols
    if norm == "horn":
        cols = ["DH", "VH", "C"]
        return cols, cols
    # For 'region' or other levels:
    if regions_df is not None and "acronym" in regions_df.columns:
        present = set(work_df["acronym"].unique()) if "acronym" in work_df.columns else set()
        ordered = [str(a) for a in regions_df["acronym"] if str(a) in present]
        remaining = [str(a) for a in work_df["acronym"].unique() if str(a) not in ordered] if "acronym" in work_df.columns else []
        cols = ordered + remaining
    elif "acronym" in work_df.columns:
        cols = sorted(str(a) for a in work_df["acronym"].unique())
    else:
        cols = []
    return cols, [structure_acronym_to_label(c) for c in cols]


def ensure_cord_rollups_in_dataframe(
    df: pd.DataFrame,
    *,
    atlas_dir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Ensure division, structure, and horn rollups exist in the DataFrame."""
    rollups_present = set(df["rollup_level"].astype(str).unique()) if "rollup_level" in df.columns else set()
    needed = {"structure", "division", "horn"}
    regions_df: pd.DataFrame | None = None

    if atlas_dir is not None:
        candidate = atlas_dir.expanduser() / "Atlas_Regions.csv"
        if candidate.is_file():
            regions_df = pd.read_csv(candidate)

    if not needed.issubset(rollups_present) and regions_df is not None:
        from lightsuite.analysis.cord_rollup import apply_cord_rollups

        df = apply_cord_rollups(df, regions_df, ["division", "structure", "horn"])

    return df, regions_df


def resolve_region_stats_path(
    config: SpinalCordPipelineConfig,
    *,
    input_path: Path | None = None,
) -> Path:
    if input_path is not None:
        path = input_path.expanduser()
        if not path.is_file():
            msg = f"Region stats file not found: {path}"
            raise FileNotFoundError(msg)
        return path
    save_path = config.sample.save_path.expanduser()
    for candidate in (
        save_path / "stats" / "region_stats.csv",
        save_path / "volume_registered" / "region_stats.csv",
    ):
        if candidate.is_file():
            return candidate
    msg = (
        f"Missing region_stats.csv under {save_path / 'stats'} or "
        f"{save_path / 'volume_registered'}. Run export (and import-annotations for counts) first, "
        "or pass --input."
    )
    raise FileNotFoundError(msg)


def _resolve_channel(
    df: pd.DataFrame,
    *,
    metric: str,
    channel: int | str | None,
) -> int | str:
    subset = df[df["metric"] == metric]
    if channel is not None:
        matches = subset[subset["channel"].astype(str) == str(channel)]
        if matches.empty:
            available = sorted(subset["channel"].astype(str).unique())
            msg = f"Channel {channel!r} not found for metric {metric!r}; available: {available}"
            raise ValueError(msg)
        return matches["channel"].iloc[0]

    channels = sorted(subset["channel"].astype(str).unique())
    if not channels:
        msg = f"No rows for metric {metric!r}"
        raise ValueError(msg)
    if len(channels) > 1:
        msg = (
            f"Multiple channels available for {metric!r}: {channels}. "
            "Pass --channel (int for intensity, label for counts)."
        )
        raise ValueError(msg)
    value = subset["channel"].iloc[0]
    if isinstance(value, (int, np.integer)) or str(value).isdigit():
        return int(value)
    return str(value)


def cord_stats_to_matrix(
    df: pd.DataFrame,
    *,
    metric: MetricName,
    channel: int | str | None = None,
    rollup_level: str = "structure",
    hemisphere: str = "whole",
    segments: list[str] | None = None,
    region_acronyms: list[str] | None = None,
    regions_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Pivot tidy cord stats into a segment × region matrix."""
    if segments is None:
        msg = "segments must be provided"
        raise ValueError(msg)

    resolved_channel = _resolve_channel(df, metric=metric, channel=channel)
    work = _stats_rows_for_hemisphere(
        df,
        metric=metric,
        rollup_level=rollup_level,
        hemisphere=hemisphere,
        channel=resolved_channel,
    )
    if work.empty:
        msg = (
            f"No data for metric={metric!r}, rollup_level={rollup_level!r}, "
            f"hemisphere={hemisphere!r}, channel={resolved_channel!r}"
        )
        raise ValueError(msg)

    if region_acronyms is not None:
        regions = list(region_acronyms)
    else:
        regions, _ = resolve_cord_rollup_columns(rollup_level, work, regions_df=regions_df)

    matrix = work.pivot_table(
        index="segment",
        columns="acronym",
        values="value",
        aggfunc="first",
    )
    matrix = matrix.reindex(index=segments, columns=regions)
    matrix.index.name = "segment"
    return matrix


def _apply_normalization(
    values: np.ndarray,
    *,
    mode: str,
) -> np.ndarray:
    normalized = values.astype(np.float64, copy=True)
    if mode == "none":
        return normalized
    if mode == "row":
        for row in range(normalized.shape[0]):
            row_vals = normalized[row]
            finite = row_vals[np.isfinite(row_vals)]
            if finite.size == 0:
                continue
            denom = float(np.nanmax(finite))
            if denom > 0:
                normalized[row] = row_vals / denom
        return normalized
    if mode == "column":
        for col in range(normalized.shape[1]):
            col_vals = normalized[:, col]
            finite = col_vals[np.isfinite(col_vals)]
            if finite.size == 0:
                continue
            denom = float(np.nanmax(finite))
            if denom > 0:
                normalized[:, col] = col_vals / denom
        return normalized
    msg = f"normalize must be 'none', 'row', or 'column'; got {mode!r}"
    raise ValueError(msg)


COMPARE_HEMISPHERES: tuple[str, ...] = ("left", "right", "whole")
DIFFERENCE_CMAP = "RdBu_r"


def default_segment_range_bounds(atlas_dir: Path, *, through: str = "Co2") -> tuple[str, str]:
    """Return (start, end) labels for the default rostrocaudal plot range."""
    order = load_segment_order(atlas_dir)
    if not order:
        msg = f"Segments.csv in {atlas_dir} has no Segment values"
        raise ValueError(msg)
    end = through if through in order else order[-1]
    return order[0], end


def _heatmap_column_labels(matrix: pd.DataFrame, column_labels: list[str] | None) -> list[str]:
    n_cols = int(matrix.shape[1])
    if column_labels is not None and len(column_labels) == n_cols:
        return list(column_labels)
    return [structure_acronym_to_label(col) for col in matrix.columns]


def _heatmap_colormap(cmap: str):
    colormap = plt.get_cmap(cmap).copy()
    colormap.set_bad(color=NO_DATA_COLOR)
    return colormap


def _finite_value_limits(*arrays: np.ndarray) -> tuple[float, float]:
    chunks = [arr[np.isfinite(arr)] for arr in arrays if getattr(arr, "size", 0)]
    finite = np.concatenate(chunks) if chunks else np.array([], dtype=np.float64)
    if finite.size == 0:
        return 0.0, 1.0
    lo = float(np.min(finite))
    hi = float(np.max(finite))
    if lo == hi:
        hi = lo + 1.0
    return lo, hi


def difference_color_limits(matrix: pd.DataFrame) -> tuple[float, float]:
    """Symmetric color limits for a left−right difference heatmap."""
    values = matrix.to_numpy(dtype=np.float64)
    finite = values[np.isfinite(values)]
    peak = float(np.max(np.abs(finite))) if finite.size else 1.0
    if peak == 0.0:
        peak = 1.0
    return -peak, peak


def _style_heatmap_axis(
    ax,
    *,
    n_rows: int,
    n_cols: int,
    column_labels: list[str],
    row_labels: list[str],
    title: str = "",
    x_label: str | None = None,
    show_ylabel: bool = True,
    show_yticklabels: bool = True,
) -> None:
    rotation = 0
    ha = "center"
    fontsize = 9
    if n_cols > 20:
        rotation = 45
        ha = "right"
        fontsize = 7
    elif n_cols > 10:
        fontsize = 8
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(column_labels, rotation=rotation, ha=ha, fontsize=fontsize)
    ax.set_yticks(np.arange(n_rows))
    if show_yticklabels:
        ax.set_yticklabels(list(row_labels), fontsize=8)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel(x_label or "Region")
    if show_ylabel:
        ax.set_ylabel("Segment")
    if title:
        ax.set_title(title, fontsize=11, pad=10)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.4, alpha=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)


def _imshow_heatmap(
    ax,
    values: np.ndarray,
    *,
    cmap: str,
    vmin: float | None,
    vmax: float | None,
    log_scale: bool = False,
):
    masked = np.ma.masked_invalid(values)
    colormap = _heatmap_colormap(cmap)
    positive = masked.compressed()
    positive = positive[positive > 0] if positive.size else positive
    norm = None
    if log_scale:
        if positive.size == 0:
            msg = "log_scale requires positive values in the heatmap"
            raise ValueError(msg)
        lo = vmin if vmin is not None else float(positive.min())
        hi = vmax if vmax is not None else float(positive.max())
        if lo <= 0:
            lo = float(positive.min())
        norm = LogNorm(vmin=lo, vmax=hi)
    return ax.imshow(
        masked,
        aspect="auto",
        origin="upper",
        cmap=colormap,
        vmin=None if norm is not None else vmin,
        vmax=None if norm is not None else vmax,
        norm=norm,
        interpolation="nearest",
    )


def build_cord_heatmap_figure(
    matrix: pd.DataFrame,
    *,
    column_labels: list[str] | None = None,
    title: str = "",
    cmap: str = "hot",
    vmin: float | None = None,
    vmax: float | None = None,
    log_scale: bool = False,
    normalize: str = "none",
    figsize: tuple[float, float] | None = None,
    x_label: str | None = None,
) -> Figure:
    """Build a publication-style heatmap figure (segments on Y, laminae/WM on X)."""
    values = _apply_normalization(matrix.to_numpy(dtype=np.float64), mode=normalize)
    n_rows, n_cols = values.shape
    labels = _heatmap_column_labels(matrix, column_labels)
    if figsize is None:
        col_width = 0.45 if n_cols <= 20 else 0.25
        figsize = (max(6.0, col_width * n_cols + 2.0), max(8.0, 0.22 * n_rows + 2.0))

    fig = Figure(figsize=figsize)
    ax = fig.add_subplot(111)
    im = _imshow_heatmap(ax, values, cmap=cmap, vmin=vmin, vmax=vmax, log_scale=log_scale)
    _style_heatmap_axis(
        ax,
        n_rows=n_rows,
        n_cols=n_cols,
        column_labels=labels,
        row_labels=list(matrix.index),
        title=title,
        x_label=x_label,
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    return fig


def cord_hemisphere_matrices(
    df: pd.DataFrame,
    *,
    metric: MetricName,
    channel: int | str | None,
    rollup_level: str,
    segments: list[str],
    regions_df: pd.DataFrame | None = None,
    hemispheres: tuple[str, ...] = COMPARE_HEMISPHERES,
) -> dict[str, pd.DataFrame]:
    """Build aligned segment × region matrices for each requested hemisphere."""
    matrices: dict[str, pd.DataFrame] = {}
    for hemisphere in hemispheres:
        matrices[hemisphere] = cord_stats_to_matrix(
            df,
            metric=metric,
            channel=channel,
            rollup_level=rollup_level,
            hemisphere=hemisphere,
            segments=segments,
            regions_df=regions_df,
        )
    return matrices


def difference_heatmap_matrix(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    """Return left − right on a shared segment × region index."""
    aligned = right.reindex(index=left.index, columns=left.columns)
    return left.subtract(aligned)


def build_cord_heatmap_comparison_figure(
    matrices: dict[str, pd.DataFrame],
    *,
    column_labels: list[str] | None = None,
    title: str = "",
    cmap: str = "hot",
    normalize: str = "none",
    x_label: str | None = None,
    order: tuple[str, ...] = COMPARE_HEMISPHERES,
) -> Figure:
    """Build side-by-side heatmaps for left, right, and whole with a shared color scale."""
    panels = [(key, matrices[key]) for key in order if key in matrices]
    if not panels:
        msg = "No hemisphere matrices to compare"
        raise ValueError(msg)
    prepared = [
        (_apply_normalization(mat.to_numpy(dtype=np.float64), mode=normalize), mat, key)
        for key, mat in panels
    ]
    vmin, vmax = _finite_value_limits(*(values for values, _mat, _key in prepared))
    first_values, first_matrix, _ = prepared[0]
    n_rows, n_cols = first_values.shape
    labels = _heatmap_column_labels(first_matrix, column_labels)
    n_panels = len(prepared)
    fig = Figure(
        figsize=(max(10.0, 0.32 * n_cols * n_panels + 2.8), max(8.0, 0.22 * n_rows + 2.4)),
        layout="constrained",
    )
    axes = fig.subplots(1, n_panels, sharey=True)
    if n_panels == 1:
        axes = [axes]
    panel_titles = {"left": "Left", "right": "Right", "whole": "Whole"}
    image = None
    for index, (values, mat, key) in enumerate(prepared):
        ax = axes[index]
        image = _imshow_heatmap(ax, values, cmap=cmap, vmin=vmin, vmax=vmax)
        _style_heatmap_axis(
            ax,
            n_rows=n_rows,
            n_cols=n_cols,
            column_labels=labels,
            row_labels=list(mat.index),
            title=panel_titles.get(key, key),
            x_label=x_label,
            show_ylabel=index == 0,
            show_yticklabels=index == 0,
        )
    if image is not None:
        cbar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
        cbar.ax.tick_params(labelsize=8)
    if title:
        fig.suptitle(title, fontsize=12)
    return fig


def write_cord_heatmap_matrix_csv(matrix: pd.DataFrame, output_path: Path) -> Path:
    """Write a segment × region heatmap matrix to CSV."""
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out = matrix.copy()
    out.index.name = "segment"
    out.to_csv(output_path)
    return output_path


def write_cord_heatmap_compare_csv(
    matrices: dict[str, pd.DataFrame],
    output_path: Path,
    *,
    order: tuple[str, ...] = COMPARE_HEMISPHERES,
) -> Path:
    """Write a tidy CSV with one column per hemisphere (segment, acronym, left, right, whole)."""
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frames: list[pd.Series] = []
    for key in order:
        if key not in matrices:
            continue
        stacked = matrices[key].stack()
        stacked.name = key
        frames.append(stacked)
    if not frames:
        msg = "No hemisphere matrices to export"
        raise ValueError(msg)
    combined = pd.concat(frames, axis=1)
    combined.index.names = ["segment", "acronym"]
    combined.to_csv(output_path)
    return output_path


def plot_cord_heatmap(
    matrix: pd.DataFrame,
    output_path: Path,
    *,
    column_labels: list[str] | None = None,
    title: str = "",
    cmap: str = "hot",
    vmin: float | None = None,
    vmax: float | None = None,
    log_scale: bool = False,
    normalize: str = "none",
    dpi: int = 300,
    figsize: tuple[float, float] | None = None,
    x_label: str | None = None,
) -> Path:
    """Save a publication-style heatmap (segments on Y, laminae/WM on X)."""
    output_path = output_path.expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig = build_cord_heatmap_figure(
        matrix,
        column_labels=column_labels,
        title=title,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        log_scale=log_scale,
        normalize=normalize,
        figsize=figsize,
        x_label=x_label,
    )
    FigureCanvasAgg(fig)
    fig.savefig(output_path, bbox_inches="tight", pad_inches=0.08, dpi=dpi)
    plt.close(fig)
    return output_path


def default_heatmap_output_path(
    config: SpinalCordPipelineConfig,
    *,
    metric: str,
    channel: int | str,
    view: str = "single",
    suffix: str = ".png",
) -> Path:
    plots_dir = config.sample.save_path.expanduser() / "plots"
    safe_channel = re.sub(r"[^\w.-]+", "_", str(channel))
    view_part = "" if view in {"", "single"} else f"_{view}"
    return plots_dir / f"cord_heatmap_{metric}_{safe_channel}{view_part}{suffix}"


def run_cord_heatmap(
    config: SpinalCordPipelineConfig,
    *,
    metric: MetricName,
    channel: int | str | None = None,
    rollup_level: str = "structure",
    hemisphere: str = "whole",
    segments: list[str] | None = None,
    region_acronyms: list[str] | None = None,
    input_path: Path | None = None,
    output_path: Path | None = None,
    cmap: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    log_scale: bool = False,
    normalize: str = "none",
) -> Path:
    """Load region stats and write a paper-style heatmap PNG."""
    stats_path = resolve_region_stats_path(config, input_path=input_path)
    df, regions_df = ensure_cord_rollups_in_dataframe(
        pd.read_csv(stats_path),
        atlas_dir=config.atlas.atlas_dir,
    )

    segment_order = segments or default_paper_segments(config.atlas.atlas_dir)
    matrix = cord_stats_to_matrix(
        df,
        metric=metric,
        channel=channel,
        rollup_level=rollup_level,
        hemisphere=hemisphere,
        segments=segment_order,
        region_acronyms=region_acronyms,
        regions_df=regions_df,
    )

    _, column_labels = resolve_cord_rollup_columns(
        rollup_level,
        df[df["rollup_level"] == rollup_level] if "rollup_level" in df.columns else df,
        regions_df=regions_df,
    )
    x_labels_by_rollup = {
        "structure": "Lamina / white matter",
        "division": "Division",
        "horn": "Horn",
        "region": "Atlas region",
    }
    x_label = x_labels_by_rollup.get(rollup_level.lower(), "Region")

    resolved_channel = _resolve_channel(df, metric=metric, channel=channel)
    if output_path is None:
        output_path = default_heatmap_output_path(config, metric=metric, channel=resolved_channel)

    metric_title = metric_label(metric)
    title = f"{config.sample.name} — {metric_title} (ch {resolved_channel})"
    if hemisphere != "whole":
        title = f"{title}, {hemisphere}"

    if cmap is None:
        cmap = "hot" if "intensity" in metric or metric == "std" else "viridis"

    return plot_cord_heatmap(
        matrix,
        output_path,
        column_labels=column_labels,
        title=title,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        log_scale=log_scale,
        normalize=normalize,
        x_label=x_label,
    )
