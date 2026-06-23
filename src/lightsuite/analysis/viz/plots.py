"""Matplotlib plots for region statistics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from lightsuite.analysis.viz.io import aggregate_by_division, filter_divisions

# Default division palette (Allen-style divisions from notebook).
DIVISION_COLORS: dict[str, str] = {
    "Isocortex": "#4e79a7",
    "Thalamus": "#f28e2b",
    "Striatum": "#e15759",
    "Hypothalamus": "#76b7b2",
    "Hippocampal formation": "#59a14f",
    "Cortical subplate": "#edc949",
}


def plot_division_bars(
    df: pd.DataFrame,
    *,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    keep_divisions: list[str] | None = None,
    exclude_divisions: list[str] | None = None,
    save_csv: bool = True,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Horizontal grouped bars: left vs right mean per division."""
    filtered = filter_divisions(df, keep=keep_divisions, exclude=exclude_divisions)
    agg = aggregate_by_division(filtered)
    if agg.empty:
        msg = "No division data to plot after filtering."
        raise ValueError(msg)

    fig, ax = plt.subplots(figsize=(10, max(6, len(agg) * 0.35)))
    y_pos = range(len(agg))
    bar_h = 0.35

    ax.barh(
        [y - bar_h / 2 for y in y_pos],
        agg["left"],
        height=bar_h,
        color="#6A5ACD",
        alpha=0.85,
        label="Left",
    )
    ax.barh(
        [y + bar_h / 2 for y in y_pos],
        agg["right"],
        height=bar_h,
        color="#9370DB",
        alpha=0.85,
        label="Right",
    )

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(agg["division"], fontsize=9)
    ax.set_xlabel("Mean value (mean of regions in division)")
    ax.set_title(title or "Left vs right by division", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        if save_csv:
            agg.to_csv(output_path.with_suffix(".csv"), index=False)

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, agg


def plot_lr_scatter(
    df: pd.DataFrame,
    *,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    keep_divisions: list[str] | None = None,
    exclude_divisions: list[str] | None = None,
    axis_min: float | None = None,
    axis_max: float | None = None,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Scatter left vs right per region, coloured by division."""
    filtered = filter_divisions(df, keep=keep_divisions, exclude=exclude_divisions)
    work = filtered.dropna(subset=["left", "right"]).copy()
    if work.empty:
        msg = "No finite left/right values to plot."
        raise ValueError(msg)

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    divisions = sorted(work["division"].dropna().astype(str).unique())
    for div in divisions:
        sub = work[work["division"].astype(str) == div]
        color = DIVISION_COLORS.get(div, "#888888")
        ax.scatter(
            sub["left"],
            sub["right"],
            s=18,
            alpha=0.7,
            edgecolors="none",
            c=color,
            label=f"{div} (n={len(sub)})",
        )

    lo = float(min(work["left"].min(), work["right"].min()))
    hi = float(max(work["left"].max(), work["right"].max()))
    pad = 0.02 * (hi - lo) if hi > lo else 1.0
    lims = (
        (float(axis_min), float(axis_max))
        if axis_min is not None and axis_max is not None
        else (lo - pad, hi + pad)
    )
    ax.plot(lims, lims, "k--", lw=1.2, zorder=1, label="y = x")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Left")
    ax.set_ylabel("Right")
    ax.set_title(title or f"Left vs right by region (n={len(work)})", fontweight="bold")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig, work
