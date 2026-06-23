"""Simple cohort group comparison bar plots."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_group_division_bars(
    summary: pd.DataFrame,
    *,
    title: str | None = None,
    output_path: Path | None = None,
    dpi: int = 200,
    show: bool = False,
    value_col: str = "mean",
    error_col: str = "sem",
) -> plt.Figure:
    """Grouped horizontal bars per division, one bar pair per experimental group."""
    required = {"group", "division", value_col}
    missing = required - set(summary.columns)
    if missing:
        msg = f"Group summary missing columns: {sorted(missing)}"
        raise KeyError(msg)

    work = summary.copy()
    work["division"] = work["division"].astype(str)
    divisions = sorted(work["division"].unique())
    groups = sorted(work["group"].unique())
    if not divisions or not groups:
        msg = "Group summary has no division/group rows to plot."
        raise ValueError(msg)

    n_div = len(divisions)
    n_grp = len(groups)
    fig_h = max(5.0, 0.45 * n_div + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    y = range(n_div)
    total_h = 0.8
    bar_h = total_h / max(n_grp, 1)
    cmap = plt.get_cmap("tab10")

    for gi, group in enumerate(groups):
        sub = work[work["group"] == group].set_index("division")
        offsets = [i + (gi - (n_grp - 1) / 2) * bar_h for i in y]
        vals = [float(sub.loc[d, value_col]) if d in sub.index else 0.0 for d in divisions]
        errs = None
        if error_col in work.columns:
            errs = [float(sub.loc[d, error_col]) if d in sub.index else 0.0 for d in divisions]
        ax.barh(
            offsets,
            vals,
            height=bar_h * 0.9,
            xerr=errs,
            label=str(group),
            color=cmap(gi % 10),
            capsize=2,
            error_kw={"elinewidth": 0.8, "ecolor": "#333"},
        )

    ax.set_yticks(list(y))
    ax.set_yticklabels(divisions, fontsize=9)
    ax.set_xlabel(value_col)
    ax.set_title(title or "Group mean by division", fontweight="bold")
    ax.legend(loc="lower right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    elif output_path is not None:
        plt.close(fig)

    return fig
