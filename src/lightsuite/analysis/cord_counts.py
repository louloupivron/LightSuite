"""Per-region and per-segment cell counts for Fiederling spinal cord atlas."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from lightsuite.analysis.cord_ontology import CordRegionTable
from lightsuite.analysis.counts import ATLAS_POINTS_KEY, atlas_points_to_voxel_indices, load_atlas_points

CORD_ATLAS_ID = "fiederling"
CORD_HEMISPHERE = "whole"

CORD_TIDY_COLUMNS = [
    "sample",
    "channel",
    "atlas",
    "parcellation_index",
    "acronym",
    "name",
    "structure",
    "division",
    "segment",
    "rollup_level",
    "hemisphere",
    "metric",
    "value",
]

# Native Fiederling export grid (Y, X, Z) voxel sizes in µm.
CORD_NATIVE_VOXEL_UM_YXZ = (10.0, 10.0, 20.0)

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


def assign_segment_names(z_1based: np.ndarray, segments: pd.DataFrame) -> np.ndarray:
    """Map 1-based length-axis indices to segment labels (``Segments.csv``)."""
    z = np.asarray(z_1based, dtype=np.int64)
    out = np.full(z.shape, "", dtype=object)
    for row in segments.itertuples(index=False):
        name = str(row.Segment)
        start = int(row.Start)
        end = int(row.End)
        out[(z >= start) & (z <= end)] = name
    return out


def _voxel_mm3_yxz(voxel_um: tuple[float, float, float]) -> float:
    vy, vx, vz = voxel_um
    return (vy * 1e-3) * (vx * 1e-3) * (vz * 1e-3)


def _region_segment_voxel_counts(
    annotation: np.ndarray,
    segments: pd.DataFrame,
) -> Counter[tuple[int, str]]:
    """Count annotation voxels per (region_id, segment) on the export (Y, X, Z) grid."""
    av = np.asarray(annotation)
    if av.ndim != 3:
        msg = f"Expected 3D annotation volume, got shape {av.shape}"
        raise ValueError(msg)

    _, _, nz = av.shape
    z_1based = np.arange(1, nz + 1, dtype=np.int64)
    seg_per_z = assign_segment_names(z_1based, segments)

    flat_av = av.ravel()
    flat_z = np.arange(flat_av.size, dtype=np.int64) % nz
    flat_seg = seg_per_z[flat_z]

    valid = (flat_av > 0) & (flat_seg != "")
    if not np.any(valid):
        return Counter()

    pairs = zip(flat_av[valid].astype(np.int64).tolist(), flat_seg[valid].tolist())
    return Counter(pairs)


def _segment_names_from_volume_indices(
    seg_indices: np.ndarray,
    segments_df: pd.DataFrame,
) -> np.ndarray:
    """Map segment row indices (1-based) from a label volume to segment names."""
    names = np.array(segments_df["Segment"].astype(str).tolist())
    out = np.full(seg_indices.shape, "", dtype=object)
    valid = seg_indices > 0
    if np.any(valid):
        idx = seg_indices[valid].astype(np.int64) - 1
        in_bounds = (idx >= 0) & (idx < len(names))
        out[valid] = names[np.clip(idx, 0, len(names) - 1)]
        out[valid & ~in_bounds] = ""
    return out


def _region_segment_voxel_counts_from_volume(
    annotation: np.ndarray,
    segments_volume: np.ndarray,
    segments_df: pd.DataFrame,
) -> Counter[tuple[int, str]]:
    av = np.asarray(annotation)
    seg_vol = np.asarray(segments_volume)
    valid = (av > 0) & (seg_vol > 0)
    if not np.any(valid):
        return Counter()
    seg_names = _segment_names_from_volume_indices(seg_vol[valid].astype(np.int64), segments_df)
    labels = av[valid].astype(np.int64)
    good = seg_names != ""
    pairs = zip(labels[good].tolist(), seg_names[good].tolist(), strict=True)
    return Counter(pairs)


def count_points_in_cord_regions(
    atlas_points_xyz: np.ndarray,
    annotation: np.ndarray,
    segments: pd.DataFrame,
    *,
    sample: str,
    channel: int | str,
    region_table: CordRegionTable | None = None,
    voxel_um_yxz: tuple[float, float, float] = CORD_NATIVE_VOXEL_UM_YXZ,
    atlas_id: str = CORD_ATLAS_ID,
    segments_volume: np.ndarray | None = None,
) -> pd.DataFrame:
    """Bin atlas-space or sample-space points into Fiederling regions and segments."""
    idx = atlas_points_to_voxel_indices(atlas_points_xyz, annotation.shape)
    if idx.size:
        labels = annotation[idx[:, 0], idx[:, 1], idx[:, 2]].astype(np.int64)
        if segments_volume is not None:
            seg_idx = segments_volume[idx[:, 0], idx[:, 1], idx[:, 2]].astype(np.int64)
            seg_names = _segment_names_from_volume_indices(seg_idx, segments)
        else:
            z_1based = idx[:, 2] + 1
            seg_names = assign_segment_names(z_1based, segments)
        valid = (labels > 0) & (seg_names != "")
        labels = labels[valid]
        seg_names = seg_names[valid]
    else:
        labels = np.zeros(0, dtype=np.int64)
        seg_names = np.zeros(0, dtype=object)

    point_counts: Counter[tuple[int, str]] = Counter()
    for rid, seg in zip(labels.tolist(), seg_names.tolist(), strict=True):
        point_counts[(int(rid), str(seg))] += 1

    if segments_volume is not None:
        voxel_counts = _region_segment_voxel_counts_from_volume(annotation, segments_volume, segments)
    else:
        voxel_counts = _region_segment_voxel_counts(annotation, segments)
    voxel_mm3 = _voxel_mm3_yxz(voxel_um_yxz)

    records: list[dict] = []
    keys = sorted(set(point_counts) | set(voxel_counts))
    for rid, seg in keys:
        count = float(point_counts.get((rid, seg), 0))
        if count <= 0:
            continue
        records.append(
            {
                "parcellation_index": int(rid),
                "segment": seg,
                "hemisphere": CORD_HEMISPHERE,
                "metric": "cell_count",
                "value": count,
            }
        )
        n_voxels = int(voxel_counts.get((rid, seg), 0))
        region_volume_mm3 = n_voxels * voxel_mm3
        if region_volume_mm3 > 0:
            records.append(
                {
                    "parcellation_index": int(rid),
                    "segment": seg,
                    "hemisphere": CORD_HEMISPHERE,
                    "metric": "cell_density",
                    "value": count / region_volume_mm3,
                }
            )

    df = pd.DataFrame.from_records(
        records,
        columns=["parcellation_index", "segment", "hemisphere", "metric", "value"],
    )
    df["sample"] = sample
    df["channel"] = channel
    df["atlas"] = atlas_id
    df["rollup_level"] = "region"

    meta_table = region_table
    df = _attach_cord_region_metadata(df, meta_table)
    return df.reindex(columns=CORD_TIDY_COLUMNS)


def write_cord_region_stats_csv(path: Path, df: pd.DataFrame) -> Path:
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


__all__ = [
    "ATLAS_POINTS_KEY",
    "CORD_ATLAS_ID",
    "CORD_TIDY_COLUMNS",
    "assign_segment_names",
    "count_points_in_cord_regions",
    "load_atlas_points",
    "write_cord_region_stats_csv",
]
