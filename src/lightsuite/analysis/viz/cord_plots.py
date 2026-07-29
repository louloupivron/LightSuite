"""Matplotlib plots for spinal cord region statistics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lightsuite.analysis.viz.cord_io import (
    division_profile_table,
    filter_cord_stats,
    segment_grouped_totals_table,
    segment_totals_table,
    structure_heatmap_matrix,
    top_regions_table,
)

_DIVISION_COLORS = {
    "GM": "#2ca02c",
    "WM": "#9467bd",
}

_GROUPED_BAR_COLORS = [
    "#e15759",
    "#4e79a7",
    "#59a14f",
    "#f28e2b",
    "#76b7b2",
    "#edc949",
]

_METRIC_LABELS = {
    "median_intensity": "Median intensity",
    "relative_median_intensity": "Relative median intensity",
    "cell_count": "Cell count",
    "cell_density": "Cell density (per mm³)",
    "std": "Intensity std",
    "volume_mm3": "Volume (mm³)",
}


def plot_cord_structure_heatmap(
    df: pd.DataFrame,
    *,
    segment_order: list[str] | None = None,
    title: str | None = None,
    metric: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
    x_tick_stride: int = 2,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Heatmap: structure (rows) × rostrocaudal segment (columns)."""
    matrix, row_labels, col_labels = structure_heatmap_matrix(df, segment_order=segment_order)
    if matrix.empty:
        msg = "No structure-level data to plot."
        raise ValueError(msg)

    data = matrix.to_numpy(dtype=float)
    fig_h = max(6.0, len(row_labels) * 0.35)
    fig_w = max(10.0, len(col_labels) * 0.25)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(data, aspect="auto", cmap="hot", origin="upper")
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)

    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    tick_idx = list(range(0, len(col_labels), max(1, x_tick_stride)))
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([col_labels[i] for i in tick_idx], rotation=0, fontsize=8)
    ax.set_xlabel("Segment")
    ax.set_ylabel("Structure")

    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_title(title or f"Structure × segment ({metric_label})", fontweight="bold")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            matrix.to_csv(output_path.with_suffix(".csv"))

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, matrix


def plot_cord_division_profile(
    df: pd.DataFrame,
    segments_df: pd.DataFrame,
    *,
    title: str | None = None,
    metric: str | None = None,
    z_voxel_um: float = 20.0,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
    x_tick_stride: int = 2,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Line plot: division signal vs rostrocaudal position (mm)."""
    profile = division_profile_table(df, segments_df, z_voxel_um=z_voxel_um)
    if profile.empty:
        msg = "No division-level profile data to plot."
        raise ValueError(msg)

    fig, ax = plt.subplots(figsize=(12, 5))
    divisions = sorted(profile["division"].astype(str).unique())
    for division in divisions:
        sub = profile[profile["division"] == division]
        color = _DIVISION_COLORS.get(division, None)
        ax.plot(
            sub["center_mm"],
            sub["value"],
            marker="o",
            markersize=3,
            linewidth=1.5,
            label=division,
            color=color,
        )

    centers = sorted(profile["center_mm"].unique())
    tick_centers = centers[:: max(1, x_tick_stride)]
    tick_labels = []
    center_to_segment = profile.drop_duplicates(subset="center_mm").set_index("center_mm")["segment"]
    for center in tick_centers:
        tick_labels.append(str(center_to_segment.get(center, "")))

    ax.set_xticks(tick_centers)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_xlim(left=0)
    ax.set_xlabel("Rostrocaudal position (mm)")
    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_ylabel(metric_label)
    ax.set_title(title or f"Division profile ({metric_label})", fontweight="bold")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.25, linestyle="--")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            profile.to_csv(output_path.with_suffix(".csv"), index=False)

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, profile


def plot_cord_segment_bars(
    df: pd.DataFrame,
    *,
    segment_order: list[str] | None = None,
    title: str | None = None,
    metric: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Bar chart: total metric per rostrocaudal segment (summed over regions)."""
    totals = segment_totals_table(df)
    if totals.empty:
        msg = "No segment totals to plot."
        raise ValueError(msg)

    if segment_order:
        order = [s for s in segment_order if s in set(totals["segment"])]
        order.extend(s for s in totals["segment"] if s not in order)
        totals = totals.set_index("segment").reindex(order).reset_index()
    totals = totals.dropna(subset=["total"])

    fig, ax = plt.subplots(figsize=(max(10, len(totals) * 0.35), 5))
    x = np.arange(len(totals))
    ax.bar(x, totals["total"], color="#4C72B0", alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(totals["segment"], rotation=45, ha="right", fontsize=8)
    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_ylabel(metric_label)
    ax.set_xlabel("Segment")
    ax.set_title(title or f"Total {metric_label.lower()} per segment", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            totals.to_csv(output_path.with_suffix(".csv"), index=False)

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, totals


def _legend_label(channel: str) -> str:
    text = str(channel)
    for prefix in ("imaris_", "imaris "):
        if text.lower().startswith(prefix):
            text = text[len(prefix) :]
            break
    return text.replace("_", " ")


def plot_cord_segment_grouped_bars(
    df: pd.DataFrame,
    *,
    channels: list[str],
    segment_order: list[str] | None = None,
    title: str | None = None,
    metric: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
    min_total: float = 0.0,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Grouped bar chart: compare import labels per rostrocaudal segment."""
    pivot = segment_grouped_totals_table(
        df,
        channels=channels,
        segment_order=segment_order,
    )
    label_cols = [col for col in pivot.columns if col != "segment"]
    if min_total > 0:
        row_totals = pivot[label_cols].sum(axis=1)
        pivot = pivot.loc[row_totals >= min_total].reset_index(drop=True)
    if pivot.empty:
        msg = "No segment totals to plot after filtering."
        raise ValueError(msg)

    n_segments = len(pivot)
    n_labels = len(label_cols)
    bar_width = 0.8 / max(n_labels, 1)
    x = np.arange(n_segments)

    fig_w = max(10.0, n_segments * 0.35)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    for idx, channel in enumerate(label_cols):
        offset = (idx - (n_labels - 1) / 2.0) * bar_width
        color = _GROUPED_BAR_COLORS[idx % len(_GROUPED_BAR_COLORS)]
        ax.bar(
            x + offset,
            pivot[channel],
            width=bar_width,
            label=_legend_label(channel),
            color=color,
            alpha=0.9,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(pivot["segment"], rotation=45, ha="right", fontsize=8)
    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_ylabel(metric_label)
    ax.set_xlabel("Segment")
    ax.set_title(title or f"Colocalization {metric_label.lower()} per segment", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            pivot.to_csv(output_path.with_suffix(".csv"), index=False)

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, pivot


def _percentile_limits(values: np.ndarray, *, lo: float = 2.0, hi: float = 98.0) -> tuple[float, float]:
    flat = values[np.isfinite(values)]
    if flat.size == 0:
        return 0.0, 1.0
    vmin, vmax = np.percentile(flat, [lo, hi])
    if vmax <= vmin:
        vmax = vmin + 1.0
    return float(vmin), float(vmax)


def plot_cord_structure_panel(
    stats_df: pd.DataFrame,
    *,
    channels: list[str],
    metric: str = "cell_count",
    rollup_level: str = "structure",
    segment_order: list[str] | None = None,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
    x_tick_stride: int = 2,
) -> tuple[plt.Figure, dict[str, pd.DataFrame]]:
    """Side-by-side structure heatmaps for several import labels."""
    if not channels:
        msg = "At least one channel/label is required."
        raise ValueError(msg)

    matrices: dict[str, pd.DataFrame] = {}
    arrays: list[np.ndarray] = []
    for channel in channels:
        table = filter_cord_stats(
            stats_df,
            channel=channel,
            metric=metric,
            rollup_level=rollup_level,
        )
        matrix, _row_labels, _col_labels = structure_heatmap_matrix(
            table,
            segment_order=segment_order,
        )
        matrices[str(channel)] = matrix
        arrays.append(matrix.to_numpy(dtype=float))

    stacked = np.concatenate([arr.ravel() for arr in arrays if arr.size])
    vmin, vmax = _percentile_limits(stacked)

    n_cols = len(channels)
    sample_rows = max(len(next(iter(matrices.values())).index), 1)
    fig_h = max(6.0, sample_rows * 0.35)
    fig_w = max(10.0, n_cols * 4.5)
    fig, axes = plt.subplots(1, n_cols, figsize=(fig_w, fig_h), squeeze=False, constrained_layout=True)

    for idx, channel in enumerate(channels):
        ax = axes[0, idx]
        matrix = matrices[str(channel)]
        data = matrix.to_numpy(dtype=float)
        im = ax.imshow(data, aspect="auto", cmap="hot", origin="upper", vmin=vmin, vmax=vmax)
        ax.set_title(_legend_label(channel), fontsize=9, fontweight="bold")
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels(matrix.index.tolist(), fontsize=7)
        col_labels = matrix.columns.tolist()
        tick_idx = list(range(0, len(col_labels), max(1, x_tick_stride)))
        ax.set_xticks(tick_idx)
        ax.set_xticklabels([col_labels[i] for i in tick_idx], rotation=0, fontsize=7)
        if idx == 0:
            ax.set_ylabel("Structure")
        ax.set_xlabel("Segment")

    metric_label = _METRIC_LABELS.get(metric, metric)
    fig.suptitle(title or f"Structure × segment ({metric_label})", fontweight="bold")
    fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            for channel, matrix in matrices.items():
                safe = str(channel).replace("/", "_")
                matrix.to_csv(output_path.with_name(f"{output_path.stem}_{safe}.csv"))

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, matrices


def plot_cord_top_regions(
    df: pd.DataFrame,
    *,
    top_n: int = 15,
    segment: str | None = None,
    title: str | None = None,
    metric: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Horizontal bar chart of the top region × segment combinations."""
    top = top_regions_table(df, top_n=top_n, segment=segment)
    if top.empty:
        msg = "No regions with positive values to plot."
        raise ValueError(msg)

    top = top.sort_values("value", ascending=True)
    fig, ax = plt.subplots(figsize=(10, max(4.0, len(top) * 0.35)))
    ax.barh(range(len(top)), top["value"], color="#4C72B0", alpha=0.9)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["plot_label"], fontsize=8)
    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_xlabel(metric_label)
    ax.set_title(title or f"Top {len(top)} regions ({metric_label.lower()})", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            top.to_csv(output_path.with_suffix(".csv"), index=False)

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, top


def plot_cord_coloc_overlap(
    summary: pd.DataFrame,
    *,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Bar chart of pairwise colocalization fractions (source → target)."""
    if summary.empty:
        msg = "No colocalization overlap data to plot."
        raise ValueError(msg)

    work = summary.copy()
    pairwise = work[work.get("comparison", "pairwise").astype(str) != "triple"].copy()
    if pairwise.empty:
        pairwise = work.copy()
    pairwise["pair"] = pairwise["source"].astype(str) + " → " + pairwise["target"].astype(str)
    pairwise = pairwise.sort_values("frac_of_source", ascending=True)

    fig, ax = plt.subplots(figsize=(10, max(4.0, len(pairwise) * 0.5)))
    y = np.arange(len(pairwise))
    ax.barh(y, pairwise["frac_of_source"], color="#e15759", alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [
            f"{row.pair} ({int(row.n_overlap_source_to_target)}/{int(row.n_source)})"
            for row in pairwise.itertuples(index=False)
        ],
        fontsize=8,
    )
    xmax = float(pairwise["frac_of_source"].max(skipna=True)) if len(pairwise) else 1.0
    ax.set_xlim(0, min(1.05, max(1.0, xmax * 1.1)))
    ax.set_xlabel("Overlap fraction")
    ax.set_title(title or "Pairwise spot colocalization overlap", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            work.to_csv(output_path.with_suffix(".csv"), index=False)

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, work


__all__ = [
    "plot_cord_coloc_overlap",
    "plot_cord_division_profile",
    "plot_cord_segment_bars",
    "plot_cord_segment_grouped_bars",
    "plot_cord_structure_heatmap",
    "plot_cord_structure_panel",
    "plot_cord_top_regions",
]
