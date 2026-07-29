"""Matplotlib plots for spinal cord region statistics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lightsuite.analysis.viz.cord_io import (
    DF_SUBREGION_ORDER,
    align_structure_heatmap_matrices,
    df_subregion_table,
    division_profile_table,
    filter_cord_stats,
    laminae_level_table,
    laminae_pct_gm_table,
    segment_grouped_totals_table,
    segment_level_class,
    segment_totals_table,
    structure_heatmap_matrix,
    structure_names_ordered,
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
    x_tick_stride: int | None = None,
    crop_empty_segments: bool = True,
) -> tuple[plt.Figure, dict[str, pd.DataFrame]]:
    """Side-by-side structure heatmaps for several import labels.

    Panels share the same structure rows and segment columns (cropped to the
    union data span) so coloc labels can be compared directly. Subplot titles
    use semantic colors and per-label totals.
    """
    if not channels:
        msg = "At least one channel/label is required."
        raise ValueError(msg)

    raw_matrices: dict[str, pd.DataFrame] = {}
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
            crop_empty_segments=crop_empty_segments,
        )
        raw_matrices[str(channel)] = matrix

    row_order = structure_names_ordered(
        stats_df,
        channels=channels,
        metric=metric,
        rollup_level=rollup_level,
    )
    matrices, row_labels, col_labels = align_structure_heatmap_matrices(
        raw_matrices,
        segment_order=segment_order,
        row_order=row_order or None,
        drop_empty_rows=True,
    )
    if not matrices or not row_labels or not col_labels:
        msg = "No structure-level panel data to plot."
        raise ValueError(msg)

    arrays = [m.to_numpy(dtype=float) for m in matrices.values() if m.size]
    stacked = np.concatenate([arr[np.isfinite(arr)] for arr in arrays if arr.size])
    vmin, vmax = _percentile_limits(stacked)

    n_cols = len(channels)
    fig_h = max(5.5, len(row_labels) * 0.38)
    fig_w = max(10.0, n_cols * 4.8)
    fig, axes = plt.subplots(1, n_cols, figsize=(fig_w, fig_h), squeeze=False, constrained_layout=True)
    cmap_obj = plt.get_cmap("magma").copy()
    cmap_obj.set_bad("#d9d9d9")

    if x_tick_stride is None:
        stride = 1 if len(col_labels) <= 20 else 2
    else:
        stride = max(1, int(x_tick_stride))
    tick_idx = list(range(0, len(col_labels), stride))

    im = None
    for idx, channel in enumerate(channels):
        ax = axes[0, idx]
        key = str(channel)
        matrix = matrices[key]
        data = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
        im = ax.imshow(data, aspect="auto", cmap=cmap_obj, origin="upper", vmin=vmin, vmax=vmax)
        label_total = float(matrix.fillna(0).to_numpy(dtype=float).sum())
        title_color = _color_for_import_label(channel, idx)
        ax.set_title(
            f"{_legend_label(channel)} (n={label_total:g})",
            fontsize=9,
            fontweight="bold",
            color=title_color,
        )
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=7)
        ax.set_xticks(tick_idx)
        ax.set_xticklabels([col_labels[i] for i in tick_idx], rotation=0, fontsize=7)
        if idx == 0:
            ax.set_ylabel("Structure")
        ax.set_xlabel("Segment")

    metric_label = _METRIC_LABELS.get(metric, metric)
    fig.suptitle(title or f"Structure × segment ({metric_label})", fontweight="bold")
    if im is not None:
        cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.02, pad=0.02)
        cbar.set_label(metric_label, fontsize=9)
    fig.text(
        0.99,
        0.01,
        "grey = no data",
        ha="right",
        va="bottom",
        fontsize=8,
        color="#666666",
    )

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
    annotate: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Horizontal bar chart of the top region × segment combinations.

    Labels lead with acronyms; when restricted to one segment the repeated
    ``@ segment`` suffix is omitted. Bars are colored by division (GM/WM) when
    available, with value and percent annotations.
    """
    top = top_regions_table(df, top_n=top_n, segment=segment)
    if top.empty:
        msg = "No regions with positive values to plot."
        raise ValueError(msg)

    top = top.sort_values("value", ascending=True).reset_index(drop=True)
    grand = float(top["value"].sum())
    top["pct"] = 100.0 * top["value"] / grand if grand > 0 else 0.0

    if "division" in top.columns:
        colors = [
            _DIVISION_COLORS.get(str(div), "#4C72B0") for div in top["division"].astype(str)
        ]
        use_div_legend = True
    else:
        colors = ["#4C72B0"] * len(top)
        use_div_legend = False

    fig, ax = plt.subplots(figsize=(11, max(4.2, len(top) * 0.42)))
    y = np.arange(len(top))
    bars = ax.barh(y, top["value"], color=colors, alpha=0.92, edgecolor="white", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(top["plot_label"], fontsize=8)
    metric_label = _METRIC_LABELS.get(metric or "", metric or "value")
    ax.set_xlabel(metric_label)
    if title is None:
        scope = f" @ {segment}" if segment else ""
        title = f"Top {len(top)} regions{scope} ({metric_label.lower()}, n={grand:g})"
    ax.set_title(title, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    xmax = float(top["value"].max()) if len(top) else 1.0
    ax.set_xlim(0, xmax * 1.18)

    if annotate:
        for bar, value, pct in zip(bars, top["value"], top["pct"], strict=True):
            width = float(value)
            ax.text(
                width + xmax * 0.015,
                bar.get_y() + bar.get_height() / 2.0,
                f"{width:g} ({pct:.0f}%)",
                va="center",
                ha="left",
                fontsize=8,
                color="#333333",
            )

    if use_div_legend:
        present = [d for d in ("GM", "WM") if d in set(top["division"].astype(str))]
        handles = [
            plt.Rectangle((0, 0), 1, 1, color=_DIVISION_COLORS[d], label=d) for d in present
        ]
        if handles:
            ax.legend(handles=handles, loc="lower right", fontsize=8, framealpha=0.92, title="Division")

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
    """Bar chart of pairwise colocalization fractions.

    Shows both directions side-by-side for each pair:
    - Left (solid): fraction of SOURCE spots that overlap a TARGET neighbour.
    - Right (hatched): fraction of TARGET spots that have a SOURCE neighbour.

    Labels use the semantic ``_legend_label`` cleaner (strips ``imaris_`` prefix).
    """
    if summary.empty:
        msg = "No colocalization overlap data to plot."
        raise ValueError(msg)

    work = summary.copy()
    pairwise = work[work["comparison"].astype(str) != "triple"].copy() if "comparison" in work.columns else work.copy()
    if pairwise.empty:
        pairwise = work.copy()

    pairwise = pairwise.sort_values("frac_of_source", ascending=True).reset_index(drop=True)

    def _clean(label: str) -> str:
        return _legend_label(label)

    pairwise["src_label"] = pairwise["source"].map(_clean)
    pairwise["tgt_label"] = pairwise["target"].map(_clean)
    pairwise["pair"] = pairwise["src_label"] + " → " + pairwise["tgt_label"]

    n_pairs = len(pairwise)
    bar_h = 0.35
    y = np.arange(n_pairs)

    _SRC_COLOR = "#e15759"
    _TGT_COLOR = "#4e79a7"

    fig, ax = plt.subplots(figsize=(10, max(4.2, n_pairs * 0.75)))

    bars_src = ax.barh(
        y + bar_h / 2,
        pairwise["frac_of_source"],
        height=bar_h,
        color=_SRC_COLOR,
        alpha=0.92,
        edgecolor="white",
        linewidth=0.5,
        label="frac of source",
    )
    bars_tgt = ax.barh(
        y - bar_h / 2,
        pairwise["frac_of_target"].fillna(0.0),
        height=bar_h,
        color=_TGT_COLOR,
        alpha=0.92,
        edgecolor="white",
        linewidth=0.5,
        label="frac of target",
    )

    for bar, row in zip(bars_src, pairwise.itertuples(index=False), strict=True):
        w = float(bar.get_width())
        ax.text(
            w + 0.01,
            bar.get_y() + bar.get_height() / 2.0,
            f"{int(row.n_overlap_source_to_target)}/{int(row.n_source)} ({w:.0%})",
            va="center",
            ha="left",
            fontsize=8,
            color=_SRC_COLOR,
        )

    for bar, row in zip(bars_tgt, pairwise.itertuples(index=False), strict=True):
        w = float(bar.get_width())
        if pd.notna(row.frac_of_target):
            ax.text(
                w + 0.01,
                bar.get_y() + bar.get_height() / 2.0,
                f"{int(row.n_overlap_target_to_source)}/{int(row.n_target)} ({w:.0%})",
                va="center",
                ha="left",
                fontsize=8,
                color=_TGT_COLOR,
            )

    ax.set_yticks(y)
    ax.set_yticklabels(pairwise["pair"], fontsize=9)
    xmax = max(
        float(pairwise["frac_of_source"].max(skipna=True)),
        float(pairwise["frac_of_target"].fillna(0).max()),
        1.0,
    )
    ax.set_xlim(0, min(1.45, xmax * 1.45))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    ax.set_xlabel("Overlap fraction")
    ax.set_title(title or "Pairwise spot colocalization overlap", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.92)
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


def plot_cord_laminae_pct_gm_bars(
    stats_df: pd.DataFrame,
    *,
    intensity_channel: int | str,
    cell_channel: int | str,
    segments: list[str] | None = None,
    intensity_metric: str = "median_intensity",
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Grouped bar chart: % GM share of intensity vs. cell counts across Rexed laminae I–X."""
    table = laminae_pct_gm_table(
        stats_df,
        intensity_channel=intensity_channel,
        cell_channel=cell_channel,
        segments=segments,
        intensity_metric=intensity_metric,
    )
    if table.empty:
        msg = "No laminae % GM data to plot."
        raise ValueError(msg)

    x = np.arange(len(table))
    width = 0.36
    fig, ax = plt.subplots(figsize=(max(10.0, len(table) * 0.85), 5.5))

    int_err = table["intensity_pct_gm_sem"].fillna(0.0).to_numpy(dtype=float)
    cell_err = table["cell_pct_gm_sem"].fillna(0.0).to_numpy(dtype=float)
    show_int_err = bool(np.any(int_err > 0))
    show_cell_err = bool(np.any(cell_err > 0))

    ax.bar(
        x - width / 2.0,
        table["intensity_pct_gm"],
        width=width,
        color="#7f7f7f",
        alpha=0.92,
        edgecolor="white",
        linewidth=0.6,
        label="Signal intensity",
        yerr=int_err if show_int_err else None,
        capsize=3,
        error_kw={"elinewidth": 1.0, "ecolor": "#555555"},
    )
    ax.bar(
        x + width / 2.0,
        table["cell_pct_gm"],
        width=width,
        color="#d62728",
        alpha=0.92,
        edgecolor="white",
        linewidth=0.6,
        label="Cell density",
        yerr=cell_err if show_cell_err else None,
        capsize=3,
        error_kw={"elinewidth": 1.0, "ecolor": "#8b1a1a"},
    )

    ax.set_xticks(x)
    ax.set_xticklabels(table["lamina"], fontsize=9)
    ax.set_xlabel("Rexed lamina")
    ax.set_ylabel("% of gray matter")
    ax.set_title(title or "% GM occupied by signal intensity vs. cell density", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ymax = float(
        max(
            (table["intensity_pct_gm"] + (int_err if show_int_err else 0)).max(),
            (table["cell_pct_gm"] + (cell_err if show_cell_err else 0)).max(),
            1.0,
        )
    )
    ax.set_ylim(0, ymax * 1.12)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.92)
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            table.to_csv(output_path.with_suffix(".csv"), index=False)

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, table


_LEVEL_BAR_STYLES = {
    "C": {"color": "#1b5e20", "hatch": "", "label": "Cervical"},
    "T": {"color": "#43a047", "hatch": "///", "label": "Thoracic"},
    "L": {"color": "#a5d6a7", "hatch": "...", "label": "Lumbar"},
    "S": {"color": "#c8e6c9", "hatch": "xx", "label": "Sacral"},
}


def plot_cord_laminae_level_bars(
    stats_df: pd.DataFrame,
    *,
    channel: int | str,
    metric: str = "median_intensity",
    segments: list[str] | None = None,
    levels: tuple[str, ...] = ("C", "T", "L"),
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Grouped bar chart: metric per Rexed lamina, averaged across cervical/thoracic/lumbar levels."""
    table = laminae_level_table(
        stats_df,
        channel=channel,
        metric=metric,
        segments=segments,
        levels=levels,
    )
    if table.empty:
        msg = "No laminae level data to plot."
        raise ValueError(msg)

    laminae = table["lamina"].astype(str).tolist()
    level_cols = [col for col in table.columns if col != "lamina"]
    x = np.arange(len(laminae))
    n_levels = len(level_cols)
    width = 0.8 / max(n_levels, 1)
    fig, ax = plt.subplots(figsize=(max(10.0, len(laminae) * 0.85), 5.5))

    for idx, level in enumerate(level_cols):
        style = _LEVEL_BAR_STYLES.get(level, {"color": _GROUPED_BAR_COLORS[idx % len(_GROUPED_BAR_COLORS)], "hatch": "", "label": level})
        offset = (idx - (n_levels - 1) / 2.0) * width
        ax.bar(
            x + offset,
            table[level],
            width=width,
            color=style["color"],
            hatch=style["hatch"],
            edgecolor="#333333",
            linewidth=0.5,
            alpha=0.95,
            label=style["label"],
        )

    metric_label = _METRIC_LABELS.get(metric, metric)
    ax.set_xticks(x)
    ax.set_xticklabels(laminae, fontsize=9)
    ax.set_xlabel("Rexed lamina")
    ax.set_ylabel(metric_label)
    ax.set_title(title or f"{metric_label} by lamina and cord level", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ymax = float(table[level_cols].to_numpy(dtype=float).max()) if level_cols else 1.0
    ax.set_ylim(0, ymax * 1.12)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.92, title="Level")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            table.to_csv(output_path.with_suffix(".csv"), index=False)

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, table


def plot_cord_df_subregion_heatmap(
    stats_df: pd.DataFrame,
    *,
    channel: int | str,
    metric: str = "median_intensity",
    segments: list[str] | None = None,
    segment_order: list[str] | None = None,
    include_parent_df: bool = True,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    save_csv: bool = True,
    x_tick_stride: int | None = None,
    crop_empty_segments: bool = True,
    cmap: str = "inferno",
) -> tuple[plt.Figure, pd.DataFrame]:
    """Heatmap: dorsal funiculus subregions (dcs, cu, gr, psdc, df) × rostrocaudal segment."""
    table = df_subregion_table(
        stats_df,
        channel=channel,
        metric=metric,
        segments=segments,
        include_parent_df=include_parent_df,
    )
    row_order = [acr for acr in DF_SUBREGION_ORDER if acr in set(table["acronym"].astype(str))]
    if not row_order:
        msg = "No dorsal funiculus subregion rows to plot."
        raise ValueError(msg)

    matrix, row_labels, col_labels = structure_heatmap_matrix(
        table,
        segment_order=segment_order,
        label_col="acronym",
        crop_empty_segments=crop_empty_segments,
    )
    ordered_rows = [row for row in row_order if row in matrix.index]
    matrix = matrix.reindex(ordered_rows)
    row_labels = ordered_rows
    if matrix.empty:
        msg = "No dorsal funiculus subregion data to plot."
        raise ValueError(msg)

    data = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    fig_h = max(4.5, len(row_labels) * 0.55)
    fig_w = max(8.0, len(col_labels) * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("#d9d9d9")
    im = ax.imshow(data, aspect="auto", cmap=cmap_obj, origin="upper")
    metric_label = _METRIC_LABELS.get(metric, metric)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label(metric_label, fontsize=9)

    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9)
    if x_tick_stride is None:
        stride = 1 if len(col_labels) <= 24 else 2
    else:
        stride = max(1, int(x_tick_stride))
    tick_idx = list(range(0, len(col_labels), stride))
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([col_labels[i] for i in tick_idx], rotation=0, fontsize=8)
    ax.set_xlabel("Segment")
    ax.set_ylabel("Dorsal funiculus subregion")
    ax.set_title(title or f"Dorsal funiculus subregions × segment ({metric_label})", fontweight="bold")
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


__all__ = [
    "plot_cord_coloc_overlap",
    "plot_cord_df_subregion_heatmap",
    "plot_cord_division_profile",
    "plot_cord_laminae_level_bars",
    "plot_cord_laminae_pct_gm_bars",
    "plot_cord_segment_bars",
    "plot_cord_segment_grouped_bars",
    "plot_cord_structure_heatmap",
    "plot_cord_structure_panel",
    "plot_cord_top_regions",
]
