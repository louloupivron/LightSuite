"""Load region statistics tables for matplotlib plots."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from lightsuite.analysis.region_stats import METRICS, TIDY_COLUMNS

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


def aggregate_by_division(df: pd.DataFrame) -> pd.DataFrame:
    """Mean left/right per division (mean of fine-region values in each division)."""
    if df.empty:
        return pd.DataFrame(columns=["division", "left", "right", "count"])

    work = df.copy()
    work["division"] = work["division"].astype(str).str.strip()
    work = work[work["division"].notna() & (work["division"] != "")]
    if work.empty:
        return pd.DataFrame(columns=["division", "left", "right", "count"])

    agg = (
        work.groupby("division", sort=True)
        .agg(left=("left", "mean"), right=("right", "mean"), count=("division", "size"))
        .reset_index()
    )
    agg["total"] = agg["left"] + agg["right"]
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
