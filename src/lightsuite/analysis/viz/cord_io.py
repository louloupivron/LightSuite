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

# Combined Rexed laminae (structure rollup targets 201–210).
LAMINAE_STRUCTURE_ACRONYMS: tuple[str, ...] = tuple(
    f"Lamina_{suffix}" for suffix in ("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X")
)
LAMINAE_DISPLAY_LABELS: dict[str, str] = {
    acronym: roman for acronym, roman in zip(LAMINAE_STRUCTURE_ACRONYMS, ("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"), strict=True)
}

# Dorsal funiculus subregions (finest atlas level).
DF_SUBREGION_ACRONYMS: tuple[str, ...] = ("dcs", "cu", "gr", "psdc")
DF_SUBREGION_ORDER: tuple[str, ...] = DF_SUBREGION_ACRONYMS + ("df",)


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
    hemisphere: str | None = None,
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
    if hemisphere is not None:
        work = work[work["hemisphere"].astype(str).str.lower() == str(hemisphere).lower()]
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
    fill_missing: float | None = None,
    crop_empty_segments: bool = False,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Pivot structure-level stats to a (structures × segments) matrix.

    Missing structure×segment combinations stay as NaN unless ``fill_missing``
    is set. When ``crop_empty_segments`` is True, leading/trailing columns with
    no finite non-zero values are dropped (span between first and last data).
    """
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

    matrix = work.pivot_table(
        index=label_col, columns="segment", values="value", aggfunc="mean"
    ).reindex(index=row_labels, columns=col_labels)
    if fill_missing is not None:
        matrix = matrix.fillna(float(fill_missing))

    if crop_empty_segments and not matrix.empty:
        nonempty = [
            col
            for col in matrix.columns
            if bool(matrix[col].notna().any()) and float(matrix[col].fillna(0.0).abs().sum()) > 0.0
        ]
        if nonempty:
            cols = list(matrix.columns)
            i0 = cols.index(nonempty[0])
            i1 = cols.index(nonempty[-1])
            matrix = matrix.iloc[:, i0 : i1 + 1]

    return matrix, list(matrix.index), list(matrix.columns)


def structure_names_ordered(
    stats_df: pd.DataFrame,
    *,
    channels: list[str],
    metric: str,
    rollup_level: str = "structure",
) -> list[str]:
    """Atlas structure names in parcellation order for panel alignment."""
    from lightsuite.analysis.viz.io import _normalize_channel

    work = stats_df[stats_df["metric"] == metric].copy()
    work = work[work["rollup_level"].astype(str).str.lower() == str(rollup_level).lower()]
    wanted = {_normalize_channel(channel) for channel in channels}
    work["channel_norm"] = work["channel"].map(_normalize_channel)
    work = work[work["channel_norm"].isin(wanted)]
    if work.empty:
        return []
    labels = work[["name", "parcellation_index"]].drop_duplicates(subset="name")
    labels["name"] = labels["name"].astype(str).str.replace("_", " ", regex=False)
    labels = labels.sort_values("parcellation_index")
    return labels["name"].tolist()


def align_structure_heatmap_matrices(
    matrices: dict[str, pd.DataFrame],
    *,
    segment_order: list[str] | None = None,
    row_order: list[str] | None = None,
    drop_empty_rows: bool = True,
) -> tuple[dict[str, pd.DataFrame], list[str], list[str]]:
    """Reindex panel matrices to shared rows/columns and crop to the data span."""
    if not matrices:
        return {}, [], []

    all_cols: set[str] = set()
    for matrix in matrices.values():
        all_cols.update(matrix.columns.astype(str))

    if segment_order:
        col_labels = [s for s in segment_order if s in all_cols]
        col_labels.extend(sorted(all_cols - set(col_labels)))
    else:
        col_labels = sorted(all_cols)

    nonempty_cols: list[str] = []
    for col in col_labels:
        for matrix in matrices.values():
            if col in matrix.columns and float(matrix[col].fillna(0).sum()) > 0:
                nonempty_cols.append(col)
                break
    if nonempty_cols:
        i0 = col_labels.index(nonempty_cols[0])
        i1 = col_labels.index(nonempty_cols[-1])
        col_labels = col_labels[i0 : i1 + 1]

    if row_order:
        all_row_set: set[str] = set()
        for matrix in matrices.values():
            all_row_set.update(str(r) for r in matrix.index)
        row_labels = [r for r in row_order if r in all_row_set]
        row_labels.extend(sorted(all_row_set - set(row_labels)))
    else:
        row_labels = sorted({str(r) for matrix in matrices.values() for r in matrix.index})

    aligned: dict[str, pd.DataFrame] = {}
    for key, matrix in matrices.items():
        aligned[key] = matrix.reindex(index=row_labels, columns=col_labels)

    if drop_empty_rows and row_labels:
        keep_rows = []
        for row in row_labels:
            if any(float(aligned[k].loc[row].fillna(0).sum()) > 0 for k in aligned):
                keep_rows.append(row)
        row_labels = keep_rows
        for key in aligned:
            aligned[key] = aligned[key].reindex(keep_rows)

    return aligned, row_labels, col_labels


def division_profile_table(
    df: pd.DataFrame,
    segments_df: pd.DataFrame,
    *,
    z_voxel_um: float = 20.0,
    crop_empty_segments: bool = False,
    drop_nonpositive: bool = False,
) -> pd.DataFrame:
    """Long table with division, segment, center_mm, and value for line plots.

    When ``crop_empty_segments`` is True, keep only the contiguous span between
    the first and last segment that has any positive value across divisions.
    When ``drop_nonpositive`` is True, replace non-positive values with NaN so
    lines break instead of diving to zero outside coverage.
    """
    centers = segment_centers_mm(segments_df, z_voxel_um=z_voxel_um)
    work = df.merge(centers[["Segment", "center_mm"]], left_on="segment", right_on="Segment", how="inner")
    work["division"] = work["acronym"].astype(str)
    work = work.sort_values(["division", "center_mm"]).reset_index(drop=True)

    if crop_empty_segments and not work.empty:
        positive = work[work["value"].astype(float) > 0.0]
        if not positive.empty:
            mm_min = float(positive["center_mm"].min())
            mm_max = float(positive["center_mm"].max())
            work = work[(work["center_mm"] >= mm_min) & (work["center_mm"] <= mm_max)].copy()

    if drop_nonpositive and not work.empty:
        work = work.copy()
        work.loc[work["value"].astype(float) <= 0.0, "value"] = float("nan")

    return work


def segment_totals_table(
    df: pd.DataFrame,
    *,
    min_total: float = 0.0,
) -> pd.DataFrame:
    """Sum metric values across regions for each segment (unsorted)."""
    work = df.copy()
    totals = (
        work.groupby("segment", sort=False)["value"]
        .sum()
        .reset_index()
        .rename(columns={"value": "total"})
    )
    if min_total > 0:
        totals = totals[totals["total"] >= float(min_total)].reset_index(drop=True)
    return totals


def segment_level_class(segment: str) -> str:
    """Map a Fiederling segment label to C/T/L/S/Co."""
    text = str(segment).strip()
    if text.startswith("Co"):
        return "Co"
    if text[:1] in {"C", "T", "L", "S"}:
        return text[:1]
    return "?"


def parse_plot_channels(value: str | None) -> list[str]:
    """Parse comma-separated CLI ``--channels`` into normalized label strings."""
    if value is None:
        return []
    return [part.strip() for part in str(value).split(",") if part.strip()]


def filter_cord_stats_multi(
    df: pd.DataFrame,
    *,
    channels: list[str],
    metric: str = "cell_count",
    rollup_level: str = "region",
    sample: str | None = None,
) -> pd.DataFrame:
    """Filter cord stats to several import labels (or channels) at one rollup level."""
    if not channels:
        msg = "At least one channel/label is required."
        raise ValueError(msg)
    if metric not in CORD_METRICS:
        msg = f"Unknown metric {metric!r}. Expected one of {CORD_METRICS}."
        raise ValueError(msg)

    wanted = {_normalize_channel(channel) for channel in channels}
    work = df[df["metric"] == metric].copy()
    work["channel_norm"] = work["channel"].map(_normalize_channel)
    work = work[work["channel_norm"].isin(wanted)]
    if sample is not None:
        work = work[work["sample"].astype(str) == str(sample)]
    level = str(rollup_level).strip().lower()
    work = work[work["rollup_level"].astype(str).str.lower() == level]
    if work.empty:
        msg = (
            f"No rows for channels={channels!r}, metric={metric!r}, "
            f"rollup_level={rollup_level!r}."
        )
        raise ValueError(msg)
    return work.reindex(columns=[c for c in CORD_TIDY_COLUMNS if c in work.columns])


def segment_grouped_totals_table(
    df: pd.DataFrame,
    *,
    channels: list[str],
    segment_order: list[str] | None = None,
) -> pd.DataFrame:
    """Pivot summed metric values to one row per segment, one column per label."""
    if not channels:
        msg = "At least one channel/label is required."
        raise ValueError(msg)

    work = df.copy()
    work["channel_norm"] = work["channel"].map(_normalize_channel)
    wanted = [_normalize_channel(channel) for channel in channels]
    work = work[work["channel_norm"].isin(set(wanted))]

    totals = (
        work.groupby(["segment", "channel_norm"], sort=False)["value"]
        .sum()
        .reset_index()
    )
    pivot = totals.pivot(index="segment", columns="channel_norm", values="value").fillna(0.0)
    col_order = [channel for channel in wanted if channel in pivot.columns]
    missing = [channel for channel in wanted if channel not in pivot.columns]
    if missing:
        msg = f"No data for channel(s): {missing}"
        raise ValueError(msg)
    pivot = pivot.reindex(columns=col_order)

    if segment_order:
        order = [segment for segment in segment_order if segment in pivot.index]
        order.extend(segment for segment in pivot.index if segment not in order)
        pivot = pivot.reindex(order)
    else:
        pivot = pivot.sort_index()

    return pivot.reset_index()




def _filter_segments(work: pd.DataFrame, segments: list[str] | None) -> pd.DataFrame:
    if not segments:
        return work
    wanted = {str(segment) for segment in segments}
    return work[work["segment"].astype(str).isin(wanted)].copy()


def _metric_pivot_for_laminae(
    stats_df: pd.DataFrame,
    *,
    channel: int | str,
    metric: str,
    segments: list[str] | None = None,
    sample: str | None = None,
) -> pd.DataFrame:
    work = filter_cord_stats(
        stats_df,
        channel=channel,
        metric=metric,
        rollup_level="structure",
        sample=sample,
    )
    work = _filter_segments(work, segments)
    work = work[work["acronym"].astype(str).isin(LAMINAE_STRUCTURE_ACRONYMS)].copy()
    return work.pivot_table(
        index=["sample", "segment"],
        columns="acronym",
        values="value",
        aggfunc="first",
    )


def laminae_pct_gm_table(
    stats_df: pd.DataFrame,
    *,
    intensity_channel: int | str,
    cell_channel: int | str,
    segments: list[str] | None = None,
    intensity_metric: str = "median_intensity",
) -> pd.DataFrame:
    """Compute % GM share per combined Rexed lamina for intensity and cell counts.

    Intensity share uses volume-weighted signal:
    ``median_intensity * volume_mm3`` summed per lamina, normalized to GM laminae.
    Cell share uses ``cell_count`` per lamina normalized to total GM laminae counts.
    When multiple samples are present, returns mean/std/sem across samples.
    """
    if intensity_metric not in CORD_METRICS:
        msg = f"Unknown intensity metric {intensity_metric!r}."
        raise ValueError(msg)

    samples = sorted(stats_df["sample"].astype(str).unique()) if "sample" in stats_df.columns else ["sample"]
    records: list[dict[str, object]] = []

    for sample in samples:
        try:
            intensity = _metric_pivot_for_laminae(
                stats_df,
                channel=intensity_channel,
                metric=intensity_metric,
                segments=segments,
                sample=sample if "sample" in stats_df.columns else None,
            )
            volume = _metric_pivot_for_laminae(
                stats_df,
                channel=intensity_channel,
                metric="volume_mm3",
                segments=segments,
                sample=sample if "sample" in stats_df.columns else None,
            )
            cells = _metric_pivot_for_laminae(
                stats_df,
                channel=cell_channel,
                metric="cell_count",
                segments=segments,
                sample=sample if "sample" in stats_df.columns else None,
            )
        except ValueError:
            continue

        if intensity.empty or cells.empty:
            continue

        intensity_signal = intensity.fillna(0.0) * volume.reindex_like(intensity).fillna(0.0)
        cell_totals = cells.fillna(0.0)
        lamina_totals: dict[str, tuple[float, float]] = {}
        for acronym in LAMINAE_STRUCTURE_ACRONYMS:
            int_val = float(intensity_signal[acronym].sum()) if acronym in intensity_signal.columns else 0.0
            cell_val = float(cell_totals[acronym].sum()) if acronym in cell_totals.columns else 0.0
            if int_val > 0.0 or cell_val > 0.0:
                lamina_totals[acronym] = (int_val, cell_val)

        total_intensity = sum(v[0] for v in lamina_totals.values())
        total_cells = sum(v[1] for v in lamina_totals.values())
        if total_intensity <= 0.0 and total_cells <= 0.0:
            continue

        for acronym, (int_val, cell_val) in lamina_totals.items():
            records.append(
                {
                    "sample": sample,
                    "acronym": acronym,
                    "lamina": LAMINAE_DISPLAY_LABELS[acronym],
                    "intensity_pct_gm": 100.0 * int_val / total_intensity if total_intensity > 0 else 0.0,
                    "cell_pct_gm": 100.0 * cell_val / total_cells if total_cells > 0 else 0.0,
                }
            )

    if not records:
        msg = "No laminae structure data for the requested channels/segments."
        raise ValueError(msg)

    work = pd.DataFrame.from_records(records)
    grouped = (
        work.groupby(["acronym", "lamina"], sort=False)[["intensity_pct_gm", "cell_pct_gm"]]
        .agg(["mean", "std", "sem"])
        .reset_index()
    )
    grouped.columns = [
        "acronym",
        "lamina",
        "intensity_pct_gm",
        "intensity_pct_gm_std",
        "intensity_pct_gm_sem",
        "cell_pct_gm",
        "cell_pct_gm_std",
        "cell_pct_gm_sem",
    ]
    order = {acr: idx for idx, acr in enumerate(LAMINAE_STRUCTURE_ACRONYMS)}
    grouped["order"] = grouped["acronym"].map(order)
    return grouped.sort_values("order").drop(columns="order").reset_index(drop=True)


def laminae_level_table(
    stats_df: pd.DataFrame,
    *,
    channel: int | str,
    metric: str = "median_intensity",
    segments: list[str] | None = None,
    levels: tuple[str, ...] = ("C", "T", "L"),
) -> pd.DataFrame:
    """Mean metric per combined Rexed lamina, averaged across segments in each cord level."""
    work = filter_cord_stats(
        stats_df,
        channel=channel,
        metric=metric,
        rollup_level="structure",
    )
    work = _filter_segments(work, segments)
    work = work[work["acronym"].astype(str).isin(LAMINAE_STRUCTURE_ACRONYMS)].copy()
    if work.empty:
        msg = f"No laminae data for channel={channel!r}, metric={metric!r}."
        raise ValueError(msg)

    work["level"] = work["segment"].map(segment_level_class)
    work = work[work["level"].isin(levels)].copy()
    if work.empty:
        msg = f"No laminae data for levels={levels!r}."
        raise ValueError(msg)

    grouped = (
        work.groupby(["acronym", "level"], sort=False)["value"]
        .mean()
        .reset_index()
    )
    pivot = grouped.pivot(index="acronym", columns="level", values="value")
    pivot = pivot.reindex(index=[acr for acr in LAMINAE_STRUCTURE_ACRONYMS if acr in pivot.index])
    pivot = pivot.reindex(columns=[level for level in levels if level in pivot.columns])
    pivot.index = [LAMINAE_DISPLAY_LABELS.get(str(acr), str(acr)) for acr in pivot.index]
    pivot.index.name = "lamina"
    return pivot.reset_index()


def df_subregion_table(
    stats_df: pd.DataFrame,
    *,
    channel: int | str,
    metric: str = "median_intensity",
    segments: list[str] | None = None,
    include_parent_df: bool = True,
    hemisphere: str | None = None,
) -> pd.DataFrame:
    """Long table for dorsal funiculus subregions at finest (region) rollup level."""
    work = filter_cord_stats(
        stats_df,
        channel=channel,
        metric=metric,
        rollup_level="region",
        hemisphere=hemisphere,
    )
    work = _filter_segments(work, segments)
    work = work[work["acronym"].astype(str).isin(DF_SUBREGION_ACRONYMS)].copy()

    frames = [work] if not work.empty else []
    if include_parent_df:
        try:
            parent = filter_cord_stats(
                stats_df,
                channel=channel,
                metric=metric,
                rollup_level="structure",
                hemisphere=hemisphere,
            )
            parent = _filter_segments(parent, segments)
            parent = parent[parent["acronym"].astype(str) == "df"].copy()
            if not parent.empty:
                frames.append(parent)
        except ValueError:
            pass

    if not frames:
        msg = f"No dorsal funiculus subregion data for channel={channel!r}, metric={metric!r}."
        raise ValueError(msg)
    return pd.concat(frames, ignore_index=True)


def top_regions_table(
    df: pd.DataFrame,
    *,
    top_n: int = 15,
    segment: str | None = None,
    min_value: float = 0.0,
) -> pd.DataFrame:
    """Select top region × segment rows by metric value."""
    work = df.copy()
    if segment is not None:
        work = work[work["segment"].astype(str) == str(segment)]
    work = work[work["value"] > float(min_value)]
    if work.empty:
        return work

    # Compact labels: acronym first; omit repeated "@ segment" when filtered.
    names = work["name"].astype(str)
    acronyms = work["acronym"].astype(str)
    segments = work["segment"].astype(str)
    if segment is not None:
        work["plot_label"] = acronyms + " — " + names
    else:
        work["plot_label"] = acronyms + " — " + names + " @ " + segments

    keep = ["plot_label", "name", "segment", "acronym", "value"]
    for col in ("structure", "division"):
        if col in work.columns:
            keep.append(col)
    top = work.nlargest(int(top_n), "value")
    return top[keep].reset_index(drop=True)


__all__ = [
    "CORD_METRICS",
    "DF_SUBREGION_ACRONYMS",
    "DF_SUBREGION_ORDER",
    "LAMINAE_DISPLAY_LABELS",
    "LAMINAE_STRUCTURE_ACRONYMS",
    "align_structure_heatmap_matrices",
    "df_subregion_table",
    "division_profile_table",
    "filter_cord_stats",
    "filter_cord_stats_multi",
    "laminae_level_table",
    "laminae_pct_gm_table",
    "load_cord_stats_csv",
    "load_segment_order",
    "parse_plot_channel",
    "parse_plot_channels",
    "resolve_cord_region_stats_from_config",
    "segment_centers_mm",
    "segment_grouped_totals_table",
    "segment_level_class",
    "segment_totals_table",
    "structure_heatmap_matrix",
    "structure_names_ordered",
    "top_regions_table",
]
