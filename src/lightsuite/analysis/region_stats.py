"""Canonical long-form (tidy) region-statistics schema and converters.

One tidy row = one (sample, channel, region, hemisphere, metric) measurement::

    sample, channel, atlas, parcellation_index, acronym, name, structure,
    division, hemisphere, metric, value

Keeping every metric as ``metric``/``value`` rows (rather than wide columns) means
new measurements — cell count, density, or anything from the future segmentation
import — are added without schema changes, and cross-atlas / cross-subject joins
stay trivial. :func:`tidy_to_wide` reproduces the legacy ``chanXX_intensities.csv``
layout (now enriched with region names) for backwards compatibility.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from lightsuite.analysis.hemisphere import SIDE_LABELS
from lightsuite.analysis.ontology import RegionTable

if TYPE_CHECKING:
    from lightsuite.export.parcellation import ParcellationResult

#: All metric names the tidy schema understands.
METRICS = (
    "median_intensity",
    "std",
    "volume_mm3",
    "cell_count",
    "cell_density",
)

#: Tidy column order.
TIDY_COLUMNS = [
    "sample",
    "channel",
    "atlas",
    "parcellation_index",
    "acronym",
    "name",
    "structure",
    "division",
    "hemisphere",
    "metric",
    "value",
]

_META_COLUMNS = ["acronym", "name", "structure", "division"]

#: ``ParcellationResult`` array field per intensity-derived metric.
_PARCELLATION_FIELDS = {
    "median_intensity": "median_over_areas",
    "std": "std_over_areas",
    "volume_mm3": "volume_over_areas",
}

#: (metric, hemisphere) → legacy wide column name.
_WIDE_NAME = {
    ("median_intensity", "right"): "RightSideIntensity",
    ("median_intensity", "left"): "LeftSideIntensity",
    ("std", "right"): "RightSideIntensityStd",
    ("std", "left"): "LeftSideIntensityStd",
    ("volume_mm3", "right"): "RightSideVolume[mm3]",
    ("volume_mm3", "left"): "LeftSideVolume[mm3]",
    ("cell_count", "right"): "RightSideCellCount",
    ("cell_count", "left"): "LeftSideCellCount",
    ("cell_density", "right"): "RightSideCellDensity[per_mm3]",
    ("cell_density", "left"): "LeftSideCellDensity[per_mm3]",
}


def attach_region_metadata(
    df: pd.DataFrame,
    region_table: RegionTable | None,
) -> pd.DataFrame:
    """Left-join region metadata onto a frame keyed by ``parcellation_index``."""
    out = df.copy()
    if region_table is not None:
        meta = (
            region_table.df[["parcellation_index", *_META_COLUMNS]]
            .drop_duplicates(subset="parcellation_index")
        )
        out = out.drop(columns=[c for c in _META_COLUMNS if c in out.columns])
        out = out.merge(meta, on="parcellation_index", how="left")
    else:
        for col in _META_COLUMNS:
            if col not in out.columns:
                out[col] = pd.NA
    return out


def parcellation_result_to_tidy(
    result: "ParcellationResult",
    region_table: RegionTable | None = None,
    *,
    sample: str,
    channel: int | str,
    atlas: str,
    drop_nan: bool = True,
) -> pd.DataFrame:
    """Expand a :class:`ParcellationResult` into tidy rows with region metadata."""
    area_ids = np.asarray(result.area_ids).astype("int64")
    frames: list[pd.DataFrame] = []
    for metric, field in _PARCELLATION_FIELDS.items():
        values = np.asarray(getattr(result, field), dtype=float)
        for side, hemisphere in enumerate(SIDE_LABELS):
            frames.append(
                pd.DataFrame(
                    {
                        "parcellation_index": area_ids,
                        "hemisphere": hemisphere,
                        "metric": metric,
                        "value": values[:, side],
                    }
                )
            )

    long = pd.concat(frames, ignore_index=True)
    long["sample"] = sample
    long["channel"] = channel
    long["atlas"] = atlas
    if drop_nan:
        long = long[np.isfinite(long["value"].to_numpy(dtype=float))].copy()
    long = attach_region_metadata(long, region_table)
    return long.reindex(columns=TIDY_COLUMNS)


def concat_tidy(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate tidy frames, preserving the canonical column order."""
    usable = [f for f in frames if f is not None and len(f) > 0]
    if not usable:
        return pd.DataFrame(columns=TIDY_COLUMNS)
    return pd.concat(usable, ignore_index=True).reindex(columns=TIDY_COLUMNS)


def tidy_to_wide(df: pd.DataFrame) -> pd.DataFrame:
    """Pivot tidy rows back to the legacy one-row-per-region wide layout."""
    work = df.copy()
    work["_col"] = [
        _WIDE_NAME.get((m, h)) for m, h in zip(work["metric"], work["hemisphere"])
    ]
    work = work[work["_col"].notna()]
    if work.empty:
        return pd.DataFrame(columns=["parcellation_index", *_META_COLUMNS])

    meta = (
        work.drop_duplicates(subset="parcellation_index")
        .set_index("parcellation_index")[_META_COLUMNS]
    )
    wide = work.pivot_table(
        index="parcellation_index",
        columns="_col",
        values="value",
        aggfunc="first",
    )
    wide = meta.join(wide).reset_index()
    wide.columns.name = None

    preferred = [
        "parcellation_index",
        "name",
        "structure",
        "division",
        "acronym",
        "RightSideIntensity",
        "LeftSideIntensity",
        "RightSideIntensityStd",
        "LeftSideIntensityStd",
        "RightSideVolume[mm3]",
        "LeftSideVolume[mm3]",
        "RightSideCellCount",
        "LeftSideCellCount",
        "RightSideCellDensity[per_mm3]",
        "LeftSideCellDensity[per_mm3]",
    ]
    ordered = [c for c in preferred if c in wide.columns]
    ordered += [c for c in wide.columns if c not in ordered]
    return wide[ordered]


def write_region_stats_csv(path: Path, df: pd.DataFrame) -> Path:
    """Write a tidy region-stats frame to CSV (creating parent dirs)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path
