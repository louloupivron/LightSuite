"""Per-region, per-segment intensity statistics for registered spinal cord volumes."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from lightsuite.analysis.cord_counts import (
    CORD_ATLAS_ID,
    CORD_HEMISPHERE,
    CORD_NATIVE_VOXEL_UM_YXZ,
    CORD_TIDY_COLUMNS,
    assign_segment_names,
)
from lightsuite.analysis.cord_hemisphere import hemisphere_label_from_side
from lightsuite.analysis.cord_ontology import CordRegionTable

_META_COLUMNS = ["acronym", "name", "structure", "division"]


def _attach_cord_region_metadata(df: pd.DataFrame, region_table: CordRegionTable | None) -> pd.DataFrame:
    out = df.copy()
    if region_table is not None:
        meta = region_table.df[["parcellation_index", *_META_COLUMNS]].drop_duplicates(
            subset="parcellation_index"
        )
        out = out.drop(columns=[c for c in _META_COLUMNS if c in out.columns])
        out = out.merge(meta, on="parcellation_index", how="left")
    else:
        for col in _META_COLUMNS:
            if col not in out.columns:
                out[col] = pd.NA
    return out


def _voxel_mm3_yxz(voxel_um: tuple[float, float, float]) -> float:
    vy, vx, vz = voxel_um
    return (vy * 1e-3) * (vx * 1e-3) * (vz * 1e-3)


def _flatten_region_segment_values(
    registered_volume: np.ndarray,
    annotation: np.ndarray,
    segments: pd.DataFrame,
    *,
    hemisphere_side: np.ndarray | None = None,
) -> pd.DataFrame:
    """Return long table of (region_id, segment, value[, side]) for in-atlas voxels."""
    av = np.asarray(annotation)
    values = np.asarray(registered_volume)
    if av.shape != values.shape:
        msg = f"Annotation shape {av.shape} != intensity volume {values.shape}"
        raise ValueError(msg)
    if av.ndim != 3:
        msg = f"Expected 3D volumes, got shape {av.shape}"
        raise ValueError(msg)
    if hemisphere_side is not None and np.asarray(hemisphere_side).shape != av.shape:
        msg = f"Hemisphere side shape {hemisphere_side.shape} != annotation {av.shape}"
        raise ValueError(msg)

    _, _, nz = av.shape
    z_1based = np.arange(1, nz + 1, dtype=np.int64)
    seg_per_z = assign_segment_names(z_1based, segments)

    flat_av = av.ravel()
    flat_vals = values.ravel().astype(np.float64, copy=False)
    flat_z = np.arange(flat_av.size, dtype=np.int64) % nz
    flat_seg = seg_per_z[flat_z]

    valid = flat_seg != ""
    if not np.any(valid):
        return pd.DataFrame(columns=["parcellation_index", "segment", "value"])

    out = pd.DataFrame(
        {
            "parcellation_index": flat_av[valid].astype(np.int64),
            "segment": flat_seg[valid].astype(str),
            "value": flat_vals[valid],
        }
    )
    if hemisphere_side is not None:
        sides = np.asarray(hemisphere_side, dtype=np.int8).ravel()[valid]
        out["side"] = sides
    return out


def _background_median_by_segment(long: pd.DataFrame) -> pd.Series:
    bg = long[long["parcellation_index"] == 0]
    if bg.empty:
        return pd.Series(dtype=np.float64)
    return bg.groupby("segment", sort=False)["value"].median()


def _aggregate_intensity_groups(
    long: pd.DataFrame,
    *,
    voxel_mm3: float,
    rel_mode: str,
    background_source: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Aggregate voxel table to one row per region × segment × hemisphere."""
    if long.empty:
        return pd.DataFrame()

    grouped = (
        long.groupby(["parcellation_index", "segment", "hemisphere"], sort=True)
        .agg(
            median_intensity=("value", "median"),
            std=("value", lambda s: float(np.std(s.to_numpy(dtype=np.float64), ddof=0)) if len(s) > 1 else 0.0),
            n_voxels=("value", "count"),
        )
        .reset_index()
    )
    grouped["volume_mm3"] = grouped["n_voxels"].astype(np.float64) * voxel_mm3

    if rel_mode == "background":
        bg_medians = _background_median_by_segment(background_source if background_source is not None else long)
        if not bg_medians.empty:
            grouped["relative_median_intensity"] = grouped.apply(
                lambda row: _relative_to_background(
                    float(row["median_intensity"]),
                    float(bg_medians.get(row["segment"], np.nan)),
                ),
                axis=1,
            )
    return grouped


def parcellate_cord_intensities(
    registered_volume: np.ndarray,
    annotation: np.ndarray,
    segments: pd.DataFrame,
    *,
    sample: str,
    channel: int | str,
    region_table: CordRegionTable | None = None,
    voxel_um_yxz: tuple[float, float, float] = CORD_NATIVE_VOXEL_UM_YXZ,
    atlas_id: str = CORD_ATLAS_ID,
    relative_to: str = "none",
    drop_background_regions: bool = True,
    hemisphere_side: np.ndarray | None = None,
    split_hemispheres: bool = False,
    keep_whole: bool = False,
) -> pd.DataFrame:
    """Compute median intensity, std, and volume per Fiederling region and segment."""
    long = _flatten_region_segment_values(
        registered_volume,
        annotation,
        segments,
        hemisphere_side=hemisphere_side if split_hemispheres else None,
    )
    if long.empty:
        return pd.DataFrame(columns=CORD_TIDY_COLUMNS)

    rel_mode = str(relative_to).strip().lower()
    voxel_mm3 = _voxel_mm3_yxz(voxel_um_yxz)
    frames: list[pd.DataFrame] = []

    if split_hemispheres:
        if hemisphere_side is None:
            msg = "hemisphere_side is required when split_hemispheres=True."
            raise ValueError(msg)
        split_long = long.copy()
        split_long["hemisphere"] = split_long["side"].map(hemisphere_label_from_side)
        split_long = split_long[split_long["hemisphere"].notna()].drop(columns="side")
        if drop_background_regions:
            split_long = split_long[split_long["parcellation_index"] > 0].copy()
        if not split_long.empty:
            frames.append(
                _aggregate_intensity_groups(
                    split_long,
                    voxel_mm3=voxel_mm3,
                    rel_mode=rel_mode,
                    background_source=long,
                )
            )

    if not split_hemispheres or keep_whole:
        whole_long = long.drop(columns="side", errors="ignore").copy()
        whole_long["hemisphere"] = CORD_HEMISPHERE
        if drop_background_regions:
            whole_long = whole_long[whole_long["parcellation_index"] > 0].copy()
        if not whole_long.empty:
            frames.append(
                _aggregate_intensity_groups(
                    whole_long,
                    voxel_mm3=voxel_mm3,
                    rel_mode=rel_mode,
                    background_source=long,
                )
            )

    if not frames:
        return pd.DataFrame(columns=CORD_TIDY_COLUMNS)

    grouped = pd.concat(frames, ignore_index=True)
    metric_columns = ["median_intensity", "std", "volume_mm3"]
    if rel_mode == "background" and "relative_median_intensity" in grouped.columns:
        metric_columns.append("relative_median_intensity")

    records: list[dict] = []
    for row in grouped.itertuples(index=False):
        base = {
            "parcellation_index": int(row.parcellation_index),
            "segment": str(row.segment),
            "hemisphere": str(row.hemisphere),
        }
        for metric in metric_columns:
            value = float(getattr(row, metric))
            if not np.isfinite(value):
                continue
            records.append({**base, "metric": metric, "value": value})

    df = pd.DataFrame.from_records(records, columns=["parcellation_index", "segment", "hemisphere", "metric", "value"])
    df["sample"] = sample
    df["channel"] = channel
    df["atlas"] = atlas_id
    df["rollup_level"] = "region"
    df = _attach_cord_region_metadata(df, region_table)
    return df.reindex(columns=CORD_TIDY_COLUMNS)


def _relative_to_background(median: float, background: float) -> float:
    if not np.isfinite(background) or background == 0.0:
        return np.nan
    return (median - background) / background


__all__ = [
    "parcellate_cord_intensities",
]
