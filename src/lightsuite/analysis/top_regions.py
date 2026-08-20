"""Build a glanceable top-N regions table from tidy region_stats."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from lightsuite.analysis.intensity_metrics import DEFAULT_INTENSITY_METRICS
from lightsuite.analysis.region_stats import METRICS, write_region_stats_csv

_SUM_RANK_METRICS = frozenset({"cell_count", "cell_density", "volume_mm3"})
_LR_HEMISPHERES = frozenset({"left", "right"})
_WHOLE_HEMISPHERES = frozenset({"", "whole"})


def top_n_regions_csv_name(n: int, *, sample_space: bool = False) -> str:
    """Filename next to ``region_stats.csv`` / ``region_stats_sample.csv``."""
    if sample_space:
        return f"region_stats_sample_top{n}.csv"
    return f"region_stats_top{n}.csv"


def resolve_rank_metric(df: pd.DataFrame, rank_by: str | None) -> str | None:
    """Pick the ranking metric: preferred ``rank_by`` if present, else cell_count, else intensity."""
    if df.empty or "metric" not in df.columns:
        return None
    present = [str(m) for m in df["metric"].dropna().unique()]
    present_set = set(present)
    if rank_by is not None and str(rank_by).strip() in present_set:
        return str(rank_by).strip()
    if "cell_count" in present_set:
        return "cell_count"
    for metric in DEFAULT_INTENSITY_METRICS:
        if metric in present_set:
            return metric
    return present[0] if present else None


def top_n_regions_wide(
    tidy: pd.DataFrame,
    *,
    n: int,
    rank_by: str | None = None,
) -> pd.DataFrame:
    """Return a wide top-N table (one row per channel × region, plus segment for cord)."""
    if n <= 0 or tidy is None or tidy.empty:
        return pd.DataFrame()
    work = _filter_leaf_regions(tidy)
    if work.empty:
        return pd.DataFrame()
    keys = _region_keys(work)
    ranked = _rank_regions(work, keys=keys, rank_by=rank_by, n=n)
    if ranked.empty:
        return pd.DataFrame()

    subset = work.merge(ranked[keys + ["rank"]], on=keys, how="inner")
    return _to_wide(subset, keys=keys)


def maybe_write_top_n_regions_csv(
    directory: Path,
    tidy: pd.DataFrame,
    *,
    n: int,
    rank_by: str | None = None,
    sample_space: bool = False,
) -> Path | None:
    """Write ``region_stats_top{n}.csv`` when ``n > 0`` and ranking succeeds."""
    wide = top_n_regions_wide(tidy, n=n, rank_by=rank_by)
    if wide.empty:
        return None
    path = directory / top_n_regions_csv_name(n, sample_space=sample_space)
    return write_region_stats_csv(path, wide)


def _filter_leaf_regions(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "parcellation_index" in out.columns:
        index = pd.to_numeric(out["parcellation_index"], errors="coerce")
        out = out.loc[index.fillna(-1) != 0].copy()
    if "rollup_level" in out.columns:
        level = out["rollup_level"].astype("string").str.strip().str.lower()
        keep = level.isin(["region", ""]) | level.isna()
        out = out.loc[keep].copy()
    return out


def _region_keys(df: pd.DataFrame) -> list[str]:
    keys = ["channel", "parcellation_index"]
    if "segment" in df.columns:
        keys.append("segment")
    return keys


def _metric_column(metric: str, hemisphere: object) -> str:
    side = str(hemisphere).strip().lower()
    if side in _WHOLE_HEMISPHERES:
        return str(metric)
    return f"{metric}_{side}"


def _hemispheres_for_rank(sides: pd.Series) -> pd.Series:
    labels = sides.astype("string").str.strip().str.lower()
    if labels.isin(list(_LR_HEMISPHERES)).any():
        return labels.isin(list(_LR_HEMISPHERES))
    return pd.Series(True, index=sides.index)


def _rank_regions(
    df: pd.DataFrame,
    *,
    keys: list[str],
    rank_by: str | None,
    n: int,
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for _, group in df.groupby("channel", dropna=False, sort=False):
        metric = resolve_rank_metric(group, rank_by)
        if metric is None:
            continue
        ranked = _rank_one_channel(group, keys=keys, metric=metric, n=n)
        if not ranked.empty:
            parts.append(ranked)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def _rank_one_channel(
    df: pd.DataFrame,
    *,
    keys: list[str],
    metric: str,
    n: int,
) -> pd.DataFrame:
    rows = df.loc[df["metric"].astype(str) == metric].copy()
    if rows.empty:
        return pd.DataFrame()
    rows = rows.loc[_hemispheres_for_rank(rows["hemisphere"])].copy()
    rows["value"] = pd.to_numeric(rows["value"], errors="coerce")
    rows = rows.loc[rows["value"].notna()]
    if rows.empty:
        return pd.DataFrame()
    how = "sum" if metric in _SUM_RANK_METRICS else "mean"
    scores = (
        rows.groupby(keys, dropna=False)["value"]
        .agg(how)
        .rename("rank_value")
        .reset_index()
    )
    scores = scores.sort_values("rank_value", ascending=False, kind="mergesort")
    scores["rank"] = range(1, len(scores) + 1)
    return scores.loc[scores["rank"] <= n].copy()


def _to_wide(df: pd.DataFrame, *, keys: list[str]) -> pd.DataFrame:
    work = df.copy()
    work["_col"] = [
        _metric_column(str(metric), hemisphere)
        for metric, hemisphere in zip(work["metric"], work["hemisphere"])
    ]
    meta_cols = [col for col in ("acronym", "name", "division") if col in work.columns]
    meta = work.drop_duplicates(subset=keys)[keys + meta_cols + ["rank"]]
    wide = work.pivot_table(
        index=keys,
        columns="_col",
        values="value",
        aggfunc="first",
    )
    wide.columns.name = None
    out = meta.set_index(keys).join(wide).reset_index()

    id_cols = ["rank", "channel", "parcellation_index"]
    if "segment" in keys:
        id_cols.append("segment")
    id_cols.extend(col for col in ("acronym", "name", "division") if col in out.columns)

    metric_cols: list[str] = []
    present_cols = set(out.columns)
    hemispheres = sorted({str(h).strip().lower() for h in work["hemisphere"].dropna().unique()})
    for metric in METRICS:
        if metric in present_cols:
            metric_cols.append(metric)
        for side in hemispheres:
            col = _metric_column(metric, side)
            if col != metric and col in present_cols:
                metric_cols.append(col)
    ordered = [col for col in id_cols if col in out.columns]
    ordered.extend(col for col in metric_cols if col not in ordered)
    ordered.extend(col for col in out.columns if col not in ordered and col != "_channel_sort")
    out["_channel_sort"] = out["channel"].astype(str)
    return (
        out.sort_values(["_channel_sort", "rank"], kind="mergesort")[ordered]
        .reset_index(drop=True)
    )
