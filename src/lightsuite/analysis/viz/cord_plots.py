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
    segment_level_class,
    segment_totals_table,
    structure_heatmap_matrix,
    top_regions_table,
)

_DIVISION_COLORS = {
    "GM": "#2ca02c",
    "WM": "#9467bd",
}

_LEVEL_COLORS = {
    "C": "#4e79a7",
    "T": "#f28e2b",
    "L": "#59a14f",
    "S": "#e15759",
    "Co": "#b07aa1",
    "?": "#9c9c9c",
}

_LEVEL_LABELS = {
    "C": "Cervical",
    "T": "Thoracic",
    "L": "Lumbar",
    "S": "Sacral",
    "Co": "Coccygeal",
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
    x_tick_stride: int | None = None,
    crop_empty_segments: bool = True,
    cmap: str = "magma",
) -> tuple[plt.Figure, pd.DataFrame]:
    """Heatmap: structure (rows) × rostrocaudal segment (columns).

    Missing structure×segment cells are shown in light grey (not as zero).
    Empty leading/trailing segments are cropped by default so the figure
    focuses on the imaged cord span.
    """
    matrix, row_labels, col_labels = structure_heatmap_matrix(
        df,
        segment_order=segment_order,
        crop_empty_segments=crop_empty_segments,
    )
    if matrix.empty:
        msg = "No structure-level data to plot."
        raise ValueError(msg)

    data = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    fig_h = max(5.5, len(row_labels) * 0.38)
    fig_w = max(8.0, len(col_labels) * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("#d9d9d9")
    im = ax.imshow(data, aspect="auto", cmap=cmap_obj, origin="upper")
    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label(metric_label, fontsize=9)

    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    if x_tick_stride is None:
        stride = 1 if len(col_labels) <= 24 else 2
    else:
        stride = max(1, int(x_tick_stride))
    tick_idx = list(range(0, len(col_labels), stride))
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([col_labels[i] for i in tick_idx], rotation=0, fontsize=8)
    ax.set_xlabel("Segment")
    ax.set_ylabel("Structure")
    ax.set_title(title or f"Structure × segment ({metric_label})", fontweight="bold")
    ax.text(
        1.0,
        -0.12,
        "grey = no data",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color="#666666",
    )
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
    x_tick_stride: int | None = None,
    crop_empty_segments: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Line plot: division signal vs rostrocaudal position (mm).

    Crops leading/trailing empty segments by default and omits non-positive
    values so lines do not falsely dive to zero outside coverage. A secondary
    top axis shows position in mm while bottom ticks use segment names.
    """
    profile = division_profile_table(
        df,
        segments_df,
        z_voxel_um=z_voxel_um,
        crop_empty_segments=crop_empty_segments,
        drop_nonpositive=True,
    )
    if profile.empty or not np.isfinite(profile["value"].to_numpy(dtype=float)).any():
        msg = "No division-level profile data to plot."
        raise ValueError(msg)

    fig, ax = plt.subplots(figsize=(11, 5))
    divisions = sorted(profile["division"].astype(str).unique())
    plotted: dict[str, pd.DataFrame] = {}
    for division in divisions:
        sub = profile[profile["division"] == division].dropna(subset=["value"])
        if sub.empty:
            continue
        plotted[division] = sub
        color = _DIVISION_COLORS.get(division, None)
        ax.plot(
            sub["center_mm"],
            sub["value"],
            marker="o",
            markersize=4,
            linewidth=1.8,
            label=division,
            color=color,
        )

    # Soft fill between GM and WM when both are present (same x grid).
    if "GM" in plotted and "WM" in plotted:
        gm = plotted["GM"].set_index("center_mm")["value"]
        wm = plotted["WM"].set_index("center_mm")["value"]
        shared = gm.index.intersection(wm.index)
        if len(shared) >= 2:
            x = np.asarray(shared, dtype=float)
            y_gm = gm.loc[shared].to_numpy(dtype=float)
            y_wm = wm.loc[shared].to_numpy(dtype=float)
            ax.fill_between(x, y_wm, y_gm, color="#2ca02c", alpha=0.12, linewidth=0)

    tick_meta = (
        profile.drop_duplicates(subset="center_mm")
        .sort_values("center_mm")[["center_mm", "segment"]]
        .reset_index(drop=True)
    )
    centers = tick_meta["center_mm"].tolist()
    if x_tick_stride is None:
        stride = 1 if len(centers) <= 20 else 2
    else:
        stride = max(1, int(x_tick_stride))
    tick_centers = centers[::stride]
    tick_labels = [str(s) for s in tick_meta["segment"].iloc[::stride].tolist()]

    ax.set_xticks(tick_centers)
    ax.set_xticklabels(tick_labels, fontsize=8)
    if centers:
        pad = max(0.15, 0.04 * (centers[-1] - centers[0]))
        ax.set_xlim(centers[0] - pad, centers[-1] + pad)
    ax.set_xlabel("Segment")
    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_ylabel(metric_label)
    ax.set_title(title or f"Division profile ({metric_label})", fontweight="bold")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.25, linestyle="--")
    ax.spines["top"].set_visible(False)

    ax_mm = ax.secondary_xaxis("top")
    ax_mm.set_xticks(tick_centers)
    ax_mm.set_xticklabels([f"{c:.1f}" for c in tick_centers], fontsize=7, color="#555555")
    ax_mm.set_xlabel("Rostrocaudal position (mm)", fontsize=9, color="#555555")

    finite_vals = profile["value"].to_numpy(dtype=float)
    finite_vals = finite_vals[np.isfinite(finite_vals)]
    if finite_vals.size:
        ymin = float(finite_vals.min())
        ymax = float(finite_vals.max())
        pad_y = max(1.0, 0.08 * (ymax - ymin))
        ax.set_ylim(max(0.0, ymin - pad_y), ymax + pad_y)

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
    min_total: float = 0.0,
    annotate: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Bar chart: total metric per rostrocaudal segment (summed over regions).

    Bars are colored by cord level (C/T/L/S/Co) and optionally annotated with
    the numeric total so a dominant segment (e.g. L5) does not hide smaller counts.
    """
    totals = segment_totals_table(df, min_total=min_total)
    if totals.empty:
        msg = "No segment totals to plot."
        raise ValueError(msg)

    if segment_order:
        order = [s for s in segment_order if s in set(totals["segment"].astype(str))]
        order.extend(s for s in totals["segment"].astype(str) if s not in order)
        totals = totals.set_index("segment").reindex(order).reset_index()
    totals = totals.dropna(subset=["total"]).reset_index(drop=True)
    totals["level"] = totals["segment"].map(segment_level_class)
    totals["pct"] = 100.0 * totals["total"] / float(totals["total"].sum())

    colors = [_LEVEL_COLORS.get(level, _LEVEL_COLORS["?"]) for level in totals["level"]]
    fig, ax = plt.subplots(figsize=(max(9.0, len(totals) * 0.55), 5.2))
    x = np.arange(len(totals))
    bars = ax.bar(x, totals["total"], color=colors, alpha=0.92, edgecolor="white", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(totals["segment"], rotation=0, fontsize=9)
    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_ylabel(metric_label)
    ax.set_xlabel("Segment")
    grand = float(totals["total"].sum())
    default_title = f"Total {metric_label.lower()} per segment (n={grand:g})"
    ax.set_title(title or default_title, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ymax = float(totals["total"].max()) if len(totals) else 1.0
    ax.set_ylim(0, ymax * 1.14)

    if annotate:
        for bar, total, pct in zip(bars, totals["total"], totals["pct"], strict=True):
            height = float(total)
            label = f"{height:g}" if pct < 8 else f"{height:g}\n({pct:.0f}%)"
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height + ymax * 0.015,
                label,
                ha="center",
                va="bottom",
                fontsize=7.5,
                color="#333333",
            )

    present_levels = [level for level in ("C", "T", "L", "S", "Co") if level in set(totals["level"])]
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=_LEVEL_COLORS[level], label=_LEVEL_LABELS[level])
        for level in present_levels
    ]
    if handles:
        ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.92, title="Level")

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
    parts = [part for part in text.replace("_", " ").split() if part]
    pretty: list[str] = []
    for part in parts:
        low = part.lower()
        if low == "coloc":
            pretty.append("Coloc")
        else:
            pretty.append(part.capitalize())
    return " ".join(pretty)


def _color_for_import_label(channel: str, fallback_idx: int) -> str:
    """Pick a color from coloc color-name tokens when present."""
    low = str(channel).lower().replace("-", " ").replace("_", " ")
    has_pink = "pink" in low
    has_cyan = "cyan" in low or "blue" in low
    has_yellow = "yellow" in low or "green" in low
    if has_pink and has_yellow:
        return "#e15759"
    if has_cyan and has_pink:
        return "#4e79a7"
    if has_cyan and has_yellow:
        return "#edc949"
    if has_pink:
        return "#e15759"
    if has_cyan:
        return "#4e79a7"
    if has_yellow:
        return "#edc949"
    return _GROUPED_BAR_COLORS[fallback_idx % len(_GROUPED_BAR_COLORS)]


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
    show_composition: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Grouped bar chart: compare import labels per rostrocaudal segment.

    Uses semantic colors from label color-names when possible. Optionally adds a
    second panel with 100% stacked bars so relative composition stays readable
    when one segment (e.g. L5) dominates absolute counts.
    """
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
    colors = [_color_for_import_label(channel, idx) for idx, channel in enumerate(label_cols)]
    pretty = [_legend_label(channel) for channel in label_cols]
    totals_per_label = {channel: float(pivot[channel].sum()) for channel in label_cols}

    fig_w = max(10.0, n_segments * 0.7)
    if show_composition:
        fig, (ax, ax_comp) = plt.subplots(
            2,
            1,
            figsize=(fig_w, 7.2),
            sharex=True,
            gridspec_kw={"height_ratios": [2.2, 1.2], "hspace": 0.12},
            constrained_layout=True,
        )
    else:
        fig, ax = plt.subplots(figsize=(fig_w, 5.0))
        ax_comp = None

    for idx, channel in enumerate(label_cols):
        offset = (idx - (n_labels - 1) / 2.0) * bar_width
        label = f"{pretty[idx]} (n={totals_per_label[channel]:g})"
        ax.bar(
            x + offset,
            pivot[channel],
            width=bar_width,
            label=label,
            color=colors[idx],
            alpha=0.92,
            edgecolor="white",
            linewidth=0.4,
        )

    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_ylabel(metric_label)
    ax.set_title(title or f"Colocalization {metric_label.lower()} per segment", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ymax = float(pivot[label_cols].to_numpy(dtype=float).max()) if n_segments else 1.0
    ax.set_ylim(0, ymax * 1.12)

    # Annotate the tallest bar in each segment group.
    for seg_i, (_, row) in enumerate(pivot.iterrows()):
        vals = [float(row[channel]) for channel in label_cols]
        peak = max(vals) if vals else 0.0
        if peak <= 0:
            continue
        peak_idx = int(np.argmax(vals))
        offset = (peak_idx - (n_labels - 1) / 2.0) * bar_width
        ax.text(
            float(seg_i) + offset,
            peak + ymax * 0.015,
            f"{peak:g}",
            ha="center",
            va="bottom",
            fontsize=7,
            color="#333333",
        )

    if ax_comp is not None:
        bottoms = np.zeros(n_segments, dtype=float)
        row_sums = pivot[label_cols].sum(axis=1).to_numpy(dtype=float)
        row_sums = np.where(row_sums > 0, row_sums, 1.0)
        for idx, channel in enumerate(label_cols):
            shares = pivot[channel].to_numpy(dtype=float) / row_sums
            ax_comp.bar(
                x,
                shares,
                bottom=bottoms,
                width=0.72,
                color=colors[idx],
                alpha=0.92,
                edgecolor="white",
                linewidth=0.4,
                label=pretty[idx],
            )
            bottoms = bottoms + shares
        ax_comp.set_ylim(0, 1.0)
        ax_comp.set_ylabel("Share of segment")
        ax_comp.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax_comp.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=8)
        ax_comp.spines["top"].set_visible(False)
        ax_comp.spines["right"].set_visible(False)
        ax_comp.grid(axis="y", alpha=0.25, linestyle="--")
        ax_comp.set_xlabel("Segment")
        ax.tick_params(labelbottom=False)
        ax_set = ax_comp
    else:
        ax.set_xlabel("Segment")
        ax_set = ax

    ax_set.set_xticks(x)
    ax_set.set_xticklabels(pivot["segment"].tolist(), rotation=0, fontsize=9)
    if ax_comp is None:
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
            crop_empty_segments=True,
        )
        matrices[str(channel)] = matrix
        arrays.append(matrix.to_numpy(dtype=float))

    stacked = np.concatenate([arr[np.isfinite(arr)] for arr in arrays if arr.size])
    vmin, vmax = _percentile_limits(stacked)

    n_cols = len(channels)
    sample_rows = max(len(next(iter(matrices.values())).index), 1)
    fig_h = max(6.0, sample_rows * 0.35)
    fig_w = max(10.0, n_cols * 4.5)
    fig, axes = plt.subplots(1, n_cols, figsize=(fig_w, fig_h), squeeze=False, constrained_layout=True)
    cmap_obj = plt.get_cmap("magma").copy()
    cmap_obj.set_bad("#d9d9d9")

    for idx, channel in enumerate(channels):
        ax = axes[0, idx]
        matrix = matrices[str(channel)]
        data = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
        im = ax.imshow(data, aspect="auto", cmap=cmap_obj, origin="upper", vmin=vmin, vmax=vmax)
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
