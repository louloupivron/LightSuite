"""Volume-weighted rollups of spinal cord region stats to division/structure levels."""

from __future__ import annotations

import re
from typing import Literal

import numpy as np
import pandas as pd

from lightsuite.analysis.cord_counts import CORD_TIDY_COLUMNS

RollupLevel = Literal["division", "structure", "horn"]

_SUM_METRICS = frozenset({"cell_count", "volume_mm3"})
_WEIGHTED_METRICS = frozenset({"median_intensity", "relative_median_intensity", "std"})
_REGION_ROLLUP_LEVEL = "region"


def _parse_children_ids(value: object) -> list[int]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (int, np.integer)):
        return [int(value)]
    if isinstance(value, float):
        return [int(value)]
    text = str(value).strip()
    if not text:
        return []
    return [int(token) for token in re.findall(r"\d+", text)]


def get_descendants(parent_id: int, regions_df: pd.DataFrame) -> list[int]:
    """Return all descendant atlas ids for ``parent_id`` (recursive, MATLAB-compatible)."""
    parent_id = int(parent_id)
    by_parent = regions_df.groupby("parent_ID")["id"].apply(list).to_dict() if "parent_ID" in regions_df.columns else {}

    direct: list[int] = [int(v) for v in by_parent.get(parent_id, [])]
    if "children_IDs" in regions_df.columns:
        row = regions_df.loc[regions_df["id"] == parent_id, "children_IDs"]
        if not row.empty:
            direct.extend(_parse_children_ids(row.iloc[0]))
    direct = sorted(set(direct))

    descendants: list[int] = []
    for child_id in direct:
        descendants.append(child_id)
        descendants.extend(get_descendants(child_id, regions_df))
    return sorted(set(descendants))


def _resolve_target_row(regions_df: pd.DataFrame, target: str) -> pd.Series | None:
    if target.isdigit():
        rows = regions_df.loc[regions_df["id"] == int(target)]
        return rows.iloc[0] if not rows.empty else None

    rows = regions_df.loc[regions_df["acronym"].astype(str) == target]
    if rows.empty and target == "lf":
        rows = regions_df.loc[regions_df["acronym"].astype(str) == "lfc"]
    return rows.iloc[0] if not rows.empty else None


def resolve_rollup_targets(aggtype: RollupLevel, regions_df: pd.DataFrame) -> list[tuple[int, str, str]]:
    """Return rollup targets as ``(id, acronym, name)`` tuples."""
    agg = aggtype.lower()
    if agg == "division":
        target_keys = ["GM", "WM"]
    elif agg == "structure":
        target_keys = [str(rid) for rid in range(201, 211)] + ["df", "lf", "vf"]
    elif agg == "horn":
        # DH = Dorsal Horn, VH = Ventral Horn, C = Central (Lamina X + canal)
        target_keys = ["DH", "VH", "C"]
    else:
        msg = f"aggtype must be 'division', 'structure', or 'horn', got {aggtype!r}"
        raise ValueError(msg)

    targets: list[tuple[int, str, str]] = []
    for key in target_keys:
        row = _resolve_target_row(regions_df, key)
        if row is None:
            continue
        targets.append((int(row["id"]), str(row["acronym"]), str(row["name"])))
    return targets


def _aggregate_metric(
    metric: str,
    values: pd.Series,
    volumes: pd.Series,
) -> float:
    if values.empty:
        return np.nan

    if metric in _SUM_METRICS:
        return float(values.sum())

    if metric == "cell_density":
        if "cell_count" not in values.index.names and volumes.empty:
            return float(values.mean())
        # Recomputed upstream when both count and volume exist.
        return float(values.iloc[0]) if len(values) == 1 else float(np.nan)

    if metric in _WEIGHTED_METRICS:
        if volumes.empty or float(volumes.sum()) <= 0:
            return float(values.mean())
        weights = volumes / volumes.sum()
        aligned = values.reindex(weights.index)
        mask = aligned.notna() & weights.notna()
        if not mask.any():
            return np.nan
        return float((aligned[mask] * weights[mask]).sum())

    return float(values.mean())


def _region_metadata_row(regions_df: pd.DataFrame, region_id: int) -> dict[str, object]:
    row = regions_df.loc[regions_df["id"] == int(region_id)]
    if row.empty:
        return {
            "parcellation_index": int(region_id),
            "acronym": "",
            "name": "",
            "structure": "",
            "division": "",
        }
    item = row.iloc[0]
    parent_acronym = str(item.get("parent_acronym", "")) if "parent_acronym" in row.columns else ""
    return {
        "parcellation_index": int(region_id),
        "acronym": str(item.get("acronym", "")),
        "name": str(item.get("name", "")),
        "structure": parent_acronym,
        "division": parent_acronym if str(parent_acronym) in {"GM", "WM", "SC", "CNS"} else "",
    }


def rollup_cord_tidy(
    df: pd.DataFrame,
    regions_df: pd.DataFrame,
    aggtype: RollupLevel,
) -> pd.DataFrame:
    """Roll finest-level cord stats up to division or structure targets."""
    if df.empty:
        return pd.DataFrame(columns=CORD_TIDY_COLUMNS)

    work = df.copy()
    if "rollup_level" not in work.columns:
        work["rollup_level"] = _REGION_ROLLUP_LEVEL
    work = work[work["rollup_level"] == _REGION_ROLLUP_LEVEL].copy()
    if work.empty:
        return pd.DataFrame(columns=CORD_TIDY_COLUMNS)

    regions = regions_df.copy()
    regions["id"] = regions["id"].astype(int)
    targets = resolve_rollup_targets(aggtype, regions)
    if not targets:
        return pd.DataFrame(columns=CORD_TIDY_COLUMNS)

    group_cols = ["sample", "channel", "atlas", "hemisphere", "segment"]
    records: list[dict] = []

    for group_key, group_df in work.groupby(group_cols, sort=False):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        group_info = dict(zip(group_cols, group_key, strict=True))

        pivot = group_df.pivot_table(
            index="parcellation_index",
            columns="metric",
            values="value",
            aggfunc="first",
        )
        if pivot.empty:
            continue

        volume_series = pivot["volume_mm3"] if "volume_mm3" in pivot.columns else pd.Series(dtype=float)
        count_series = pivot["cell_count"] if "cell_count" in pivot.columns else pd.Series(dtype=float)

        for target_id, _acronym, _name in targets:
            member_ids = sorted(set(get_descendants(target_id, regions)) | {target_id})
            present = [rid for rid in member_ids if rid in pivot.index]
            if not present:
                continue

            subset = pivot.loc[present]
            vol_subset = volume_series.reindex(present).fillna(0.0)

            meta = _region_metadata_row(regions, target_id)
            base = {
                **group_info,
                **meta,
                "rollup_level": aggtype.lower(),
            }

            metrics = [m for m in subset.columns if m != "cell_density"]
            for metric in metrics:
                value = _aggregate_metric(metric, subset[metric], vol_subset)
                if not np.isfinite(value):
                    continue
                records.append({**base, "metric": metric, "value": value})

            if not count_series.empty and not vol_subset.empty:
                total_count = float(count_series.reindex(present).fillna(0.0).sum())
                total_vol = float(vol_subset.sum())
                if total_vol > 0 and total_count > 0:
                    records.append(
                        {
                            **base,
                            "metric": "cell_density",
                            "value": total_count / total_vol,
                        }
                    )

    if not records:
        return pd.DataFrame(columns=CORD_TIDY_COLUMNS)

    out = pd.DataFrame.from_records(records)
    return out.reindex(columns=CORD_TIDY_COLUMNS)


def apply_cord_rollups(
    df: pd.DataFrame,
    regions_df: pd.DataFrame,
    rollups: list[str],
) -> pd.DataFrame:
    """Append division/structure rollup rows to finest-level cord stats."""
    if not rollups:
        out = df.copy()
        if "rollup_level" not in out.columns:
            out["rollup_level"] = _REGION_ROLLUP_LEVEL
        return out.reindex(columns=CORD_TIDY_COLUMNS)

    finest = df.copy()
    if "rollup_level" not in finest.columns:
        finest["rollup_level"] = _REGION_ROLLUP_LEVEL
    frames = [finest.reindex(columns=CORD_TIDY_COLUMNS)]

    for level in rollups:
        normalized = str(level).strip().lower()
        if normalized not in ("division", "structure", "horn"):
            continue
        rolled = rollup_cord_tidy(finest, regions_df, normalized)  # type: ignore[arg-type]
        if len(rolled):
            frames.append(rolled)

    return pd.concat(frames, ignore_index=True).reindex(columns=CORD_TIDY_COLUMNS)


__all__ = [
    "apply_cord_rollups",
    "get_descendants",
    "resolve_rollup_targets",
    "rollup_cord_tidy",
]
