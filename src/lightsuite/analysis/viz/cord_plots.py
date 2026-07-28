"""Matplotlib plots for spinal cord region statistics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lightsuite.analysis.viz.cord_io import (
    division_profile_table,
    segment_totals_table,
    structure_heatmap_matrix,
)

_DIVISION_COLORS = {
    "GM": "#2ca02c",
    "WM": "#9467bd",
}

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


__all__ = [
    "plot_cord_division_profile",
    "plot_cord_segment_bars",
    "plot_cord_structure_heatmap",
]
