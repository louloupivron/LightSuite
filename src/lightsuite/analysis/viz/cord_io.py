"""Load and prepare spinal cord region_stats tables for plotting."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from lightsuite.analysis.cord_counts import CORD_TIDY_COLUMNS
from lightsuite.analysis.viz.io import _normalize_channel, parse_plot_channel

CORD_METRICS = (
    "median_intensity",
    "relative_median_intensity",
    "std",
    "volume_mm3",
    "cell_count",
    "cell_density",
)


def resolve_cord_region_stats_from_config(config_path: str | Path) -> Path:
    """Resolve ``volume_registered/region_stats.csv`` from a spinal cord YAML."""
    from lightsuite.config.loader import load_spinal_config

    cfg = load_spinal_config(config_path)
    stats = cfg.sample.save_path.expanduser().resolve() / "volume_registered" / "region_stats.csv"
    if not stats.is_file():
        msg = f"region_stats not found: {stats}. Run 'lightsuite spinal region-stats' first."
        raise FileNotFoundError(msg)
    return stats


def load_cord_stats_csv(path: str | Path) -> pd.DataFrame:
    csv_path = Path(path).expanduser().resolve()
    if not csv_path.is_file():
        msg = f"Input not found: {csv_path}"
        raise FileNotFoundError(msg)
    df = pd.read_csv(csv_path)
    required = {"segment", "metric", "value", "channel"}
    missing = required - set(df.columns)
    if missing:
        msg = f"Cord stats CSV missing columns: {sorted(missing)}"
        raise ValueError(msg)
    if "rollup_level" not in df.columns:
        df = df.copy()
        df["rollup_level"] = "region"
    return df


def filter_cord_stats(
    df: pd.DataFrame,
    *,
    channel: int | str | None = None,
    metric: str = "median_intensity",
    rollup_level: str = "structure",
    sample: str | None = None,
) -> pd.DataFrame:
    """Filter a cord tidy table to one channel, metric, rollup level, and optional sample."""
    if metric not in CORD_METRICS:
        msg = f"Unknown metric {metric!r}. Expected one of {CORD_METRICS}."
        raise ValueError(msg)

    work = df[df["metric"] == metric].copy()
    if channel is not None:
        want = _normalize_channel(channel)
        work = work[work["channel"].map(_normalize_channel) == want]
    if sample is not None:
        work = work[work["sample"].astype(str) == str(sample)]
    level = str(rollup_level).strip().lower()
    work = work[work["rollup_level"].astype(str).str.lower() == level]
    if work.empty:
        msg = (
            f"No rows for channel={channel!r}, metric={metric!r}, "
            f"rollup_level={rollup_level!r}."
        )
        raise ValueError(msg)
    return work.reindex(columns=[c for c in CORD_TIDY_COLUMNS if c in work.columns])


def load_segment_order(segments_csv: Path | None) -> list[str]:
    if segments_csv is None or not segments_csv.is_file():
        return []
    segments = pd.read_csv(segments_csv)
    if "Segment" not in segments.columns:
        return []
    return segments["Segment"].astype(str).tolist()


def segment_centers_mm(segments_df: pd.DataFrame, *, z_voxel_um: float = 20.0) -> pd.DataFrame:
    """Rostrocaudal segment center positions in mm (Fiederling Z spacing)."""
    out = segments_df.copy()
    out["Segment"] = out["Segment"].astype(str)
    centers = (out["Start"].astype(float) + out["End"].astype(float)) / 2.0
    out["center_mm"] = centers * float(z_voxel_um) * 1e-3
    return out


def structure_heatmap_matrix(
    df: pd.DataFrame,
    *,
    segment_order: list[str] | None = None,
    label_col: str = "name",
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Pivot structure-level stats to a (structures × segments) matrix."""
    work = df.copy()
    work[label_col] = work[label_col].astype(str).str.replace("_", " ", regex=False)
    labels = (
        work[[label_col, "parcellation_index"]]
        .drop_duplicates(subset=label_col)
        .sort_values("parcellation_index")
    )
    row_labels = labels[label_col].tolist()

    if segment_order:
        col_labels = [s for s in segment_order if s in set(work["segment"].astype(str))]
        col_labels.extend(sorted(set(work["segment"].astype(str)) - set(col_labels)))
    else:
        col_labels = sorted(work["segment"].astype(str).unique())

    matrix = (
        work.pivot_table(index=label_col, columns="segment", values="value", aggfunc="mean")
        .reindex(index=row_labels, columns=col_labels)
        .fillna(0.0)
    )
    return matrix, row_labels, col_labels


def division_profile_table(
    df: pd.DataFrame,
    segments_df: pd.DataFrame,
    *,
    z_voxel_um: float = 20.0,
) -> pd.DataFrame:
    """Long table with division, segment, center_mm, and value for line plots."""
    centers = segment_centers_mm(segments_df, z_voxel_um=z_voxel_um)
    work = df.merge(centers[["Segment", "center_mm"]], left_on="segment", right_on="Segment", how="inner")
    work["division"] = work["acronym"].astype(str)
    return work.sort_values(["division", "center_mm"])


def segment_totals_table(df: pd.DataFrame) -> pd.DataFrame:
    """Sum metric values across regions for each segment."""
    work = df.copy()
    totals = (
        work.groupby("segment", sort=False)["value"]
        .sum()
        .reset_index()
        .rename(columns={"value": "total"})
    )
    return totals.sort_values("total", ascending=False)


__all__ = [
    "CORD_METRICS",
    "division_profile_table",
    "filter_cord_stats",
    "load_cord_stats_csv",
    "load_segment_order",
    "resolve_cord_region_stats_from_config",
    "segment_centers_mm",
    "segment_totals_table",
    "structure_heatmap_matrix",
    "parse_plot_channel",
]
