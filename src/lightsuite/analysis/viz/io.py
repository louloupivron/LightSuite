"""Load region statistics tables for matplotlib plots."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import pandas as pd

from lightsuite.analysis.region_stats import METRICS, TIDY_COLUMNS

DivisionAggregate = Literal["mean", "sum"]

#: Legacy wide intensity column names.
_WIDE_INTENSITY = {
    "median_intensity": ("RightSideIntensity", "LeftSideIntensity"),
    "std": ("RightSideIntensityStd", "LeftSideIntensityStd"),
    "volume_mm3": ("RightSideVolume[mm3]", "LeftSideVolume[mm3]"),
    "cell_count": ("RightSideCellCount", "LeftSideCellCount"),
    "cell_density": ("RightSideCellDensity[per_mm3]", "LeftSideCellDensity[per_mm3]"),
}

#: Standard per-region frame used by plot functions.
PLOT_COLUMNS = [
    "parcellation_index",
    "name",
    "structure",
    "division",
    "left",
    "right",
]


def _normalize_channel(value: object) -> str:
    if isinstance(value, (int, float)) and value == int(value):
        return str(int(value))
    return str(value).strip()


def parse_plot_channel(value: str | int | None) -> int | str | None:
    """Parse CLI ``--channel``: numeric imaging channels or import label strings."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    if text.lstrip("-").isdigit():
        return int(text)
    return text


def load_region_plot_table(
    path: str | Path,
    *,
    channel: int | str | None = None,
    metric: str = "median_intensity",
) -> pd.DataFrame:
    """Load a tidy or legacy wide CSV into a standard per-region plotting table."""
    if metric not in METRICS:
        msg = f"Unknown metric {metric!r}. Expected one of {METRICS}."
        raise ValueError(msg)

    csv_path = Path(path).expanduser().resolve()
    if not csv_path.is_file():
        msg = f"Input not found: {csv_path}"
        raise FileNotFoundError(msg)

    df = pd.read_csv(csv_path)
    if _is_tidy(df):
        return _tidy_to_plot_table(df, channel=channel, metric=metric)
    return _wide_to_plot_table(df, metric=metric)


def resolve_region_stats_from_config(config_path: str | Path) -> Path:
    """Resolve ``volume_registered/region_stats.csv`` from a brain pipeline YAML."""
    from lightsuite.config.loader import load_config

    cfg = load_config(config_path)
    stats = cfg.sample.save_path.expanduser().resolve() / "volume_registered" / "region_stats.csv"
    if not stats.is_file():
        msg = f"region_stats not found: {stats}. Run 'lightsuite analysis region-stats' first."
        raise FileNotFoundError(msg)
    return stats


def filter_divisions(
    df: pd.DataFrame,
    *,
    keep: list[str] | None = None,
    exclude: list[str] | None = None,
) -> pd.DataFrame:
    """Filter rows by division name."""
    out = df.copy()
    if "division" not in out.columns:
        return out
    if keep:
        allowed = {d.strip() for d in keep}
        out = out[out["division"].astype(str).str.strip().isin(allowed)]
    if exclude:
        blocked = {d.strip() for d in exclude}
        out = out[~out["division"].astype(str).str.strip().isin(blocked)]
    return out


def default_division_aggregate(metric: str) -> DivisionAggregate:
    """Default division rollup: sum for counts, mean for continuous metrics."""
    return "sum" if metric == "cell_count" else "mean"


def select_top_regions(
    df: pd.DataFrame,
    *,
    top_n: int,
) -> pd.DataFrame:
    """Return the top *top_n* regions ranked by left + right total."""
    if top_n < 1:
        msg = f"top_n must be >= 1, got {top_n}."
        raise ValueError(msg)

    work = df.copy()
    work["total"] = work["left"].fillna(0) + work["right"].fillna(0)
    work = work[work["total"] > 0]
    if work.empty:
        return work.reindex(columns=[*PLOT_COLUMNS, "total"])

    ranked = work.sort_values("total", ascending=False).head(top_n)
    return ranked.sort_values("total", ascending=True).reset_index(drop=True)


def aggregate_by_division(
    df: pd.DataFrame,
    *,
    how: DivisionAggregate = "mean",
) -> pd.DataFrame:
    """Roll up fine regions to divisions (mean or sum of left/right per side)."""
    if df.empty:
        return pd.DataFrame(columns=["division", "left", "right", "count"])

    work = df.copy()
    work["division"] = work["division"].astype(str).str.strip()
    work = work[work["division"].notna() & (work["division"] != "")]
    if work.empty:
        return pd.DataFrame(columns=["division", "left", "right", "count"])

    reducer = "sum" if how == "sum" else "mean"
    agg = (
        work.groupby("division", sort=True)
        .agg(left=("left", reducer), right=("right", reducer), count=("division", "size"))
        .reset_index()
    )
    agg["total"] = agg["left"].fillna(0) + agg["right"].fillna(0)
    return agg.sort_values("total", ascending=True)


def _is_tidy(df: pd.DataFrame) -> bool:
    return set(TIDY_COLUMNS).issubset(set(df.columns))


def _tidy_to_plot_table(
    df: pd.DataFrame,
    *,
    channel: int | str | None,
    metric: str,
) -> pd.DataFrame:
    work = df[df["metric"] == metric].copy()
    if channel is not None:
        want = _normalize_channel(channel)
        work = work[work["channel"].map(_normalize_channel) == want]
    if work.empty:
        msg = f"No rows for channel={channel!r}, metric={metric!r}."
        raise ValueError(msg)

    meta_cols = ["parcellation_index", "name", "structure", "division"]
    meta = work.drop_duplicates(subset="parcellation_index")[meta_cols]
    pivot = (
        work.pivot_table(
            index="parcellation_index",
            columns="hemisphere",
            values="value",
            aggfunc="first",
        )
        .rename(columns={"right": "right", "left": "left"})
    )
    out = meta.set_index("parcellation_index").join(pivot, how="inner").reset_index()
    out = out[out["parcellation_index"].fillna(0).astype(int) != 0]
    return out.reindex(columns=PLOT_COLUMNS)


def _wide_to_plot_table(df: pd.DataFrame, *, metric: str) -> pd.DataFrame:
    right_col, left_col = _WIDE_INTENSITY[metric]
    for col in (right_col, left_col):
        if col not in df.columns:
            msg = f"Wide CSV missing column {col!r} for metric {metric!r}."
            raise KeyError(msg)

    out = df.copy()
    out["left"] = pd.to_numeric(out[left_col], errors="coerce")
    out["right"] = pd.to_numeric(out[right_col], errors="coerce")
    for col in ("name", "structure", "division"):
        if col not in out.columns:
            out[col] = pd.NA
    if "parcellation_index" not in out.columns:
        msg = "Wide CSV must include parcellation_index."
        raise KeyError(msg)
    out = out[out["parcellation_index"].fillna(0).astype(int) != 0]
    return out.reindex(columns=PLOT_COLUMNS)
