"""Per-region cell counts and densities from atlas-space point clouds.

``lightsuite brain import-annotations`` writes warped cell coordinates to
``<label>_atlas_coords.npz`` (key ``atlasptcoords``, 1-based atlas ``x, y, z``).
This module bins those points into atlas regions — using the *same* left/right
split as the intensity statistics (:mod:`lightsuite.analysis.hemisphere`) — and
emits tidy ``cell_count`` / ``cell_density`` rows that line up with the intensity
tidy table. Density is ``count / region_volume_mm3`` computed from the annotation
voxels on that side, so it never depends on a prior intensity export.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from lightsuite.analysis.hemisphere import SIDE_LABELS, hemisphere_side_volume
from lightsuite.analysis.ontology import RegionTable
from lightsuite.analysis.region_stats import TIDY_COLUMNS, attach_region_metadata
from lightsuite.registration.points import cloud_xyz_to_volume_indices

ATLAS_POINTS_KEY = "atlasptcoords"
SAMPLE_POINTS_KEY = "regptcoords"


def load_atlas_points(npz_path: Path, *, key: str = ATLAS_POINTS_KEY) -> np.ndarray:
    """Load an ``Nx3+`` atlas-space point array from an import ``.npz`` file."""
    with np.load(npz_path) as data:
        if key not in data:
            available = ", ".join(sorted(data.files))
            msg = f"{npz_path} has no array '{key}'. Available: {available}"
            raise KeyError(msg)
        return np.asarray(data[key], dtype=float)


def atlas_points_to_voxel_indices(
    atlas_points_xyz: np.ndarray,
    shape: tuple[int, int, int],
) -> np.ndarray:
    """Convert 1-based atlas ``x, y, z`` points to in-bounds 0-based ``(y, x, z)`` indices."""
    pts = np.asarray(atlas_points_xyz, dtype=float)
    if pts.ndim != 2 or pts.shape[1] < 3:
        msg = f"Expected Nx3+ atlas points, got shape {pts.shape}"
        raise ValueError(msg)
    yxz = cloud_xyz_to_volume_indices(pts[:, :3])
    idx = np.rint(yxz).astype(np.int64) - 1  # documented 1-based → 0-based
    bounds = np.asarray(shape)
    in_bounds = np.all((idx >= 0) & (idx < bounds), axis=1)
    return idx[in_bounds]


def count_points_in_regions(
    atlas_points_xyz: np.ndarray,
    annotation: np.ndarray,
    *,
    atlas_id: str,
    atlas_resolution_um: float,
    sample: str,
    channel: int | str,
    region_table: RegionTable | None = None,
    ml_axis: int = 2,
    brainglobe_name: str | None = None,
) -> pd.DataFrame:
    """Bin atlas-space points into regions and return tidy count/density rows."""
    av = np.asanyarray(annotation)
    side_volume = hemisphere_side_volume(
        av,
        atlas_id,
        ml_axis=ml_axis,
        brainglobe_name=brainglobe_name,
    )
    voxel_mm3 = (float(atlas_resolution_um) * 1e-3) ** 3

    idx = atlas_points_to_voxel_indices(atlas_points_xyz, av.shape)
    if idx.size:
        labels = av[idx[:, 0], idx[:, 1], idx[:, 2]].astype(np.int64)
        sides = side_volume[idx[:, 0], idx[:, 1], idx[:, 2]]
    else:
        labels = np.zeros(0, dtype=np.int64)
        sides = np.zeros(0, dtype=np.int64)

    records: list[dict] = []
    for side, hemisphere in enumerate(SIDE_LABELS):
        point_labels = labels[sides == side]
        point_labels = point_labels[point_labels > 0]
        if point_labels.size == 0:
            continue
        counts = Counter(point_labels.tolist())

        side_voxels = av[(side_volume == side) & (av > 0)].astype(np.int64)
        voxel_counts = np.bincount(side_voxels) if side_voxels.size else np.zeros(1, dtype=np.int64)

        for pidx, count in counts.items():
            records.append(
                {
                    "parcellation_index": int(pidx),
                    "hemisphere": hemisphere,
                    "metric": "cell_count",
                    "value": float(count),
                }
            )
            n_voxels = int(voxel_counts[pidx]) if pidx < voxel_counts.size else 0
            region_volume_mm3 = n_voxels * voxel_mm3
            if region_volume_mm3 > 0:
                records.append(
                    {
                        "parcellation_index": int(pidx),
                        "hemisphere": hemisphere,
                        "metric": "cell_density",
                        "value": float(count) / region_volume_mm3,
                    }
                )

    df = pd.DataFrame.from_records(records, columns=["parcellation_index", "hemisphere", "metric", "value"])
    df["sample"] = sample
    df["channel"] = channel
    df["atlas"] = atlas_id
    df = attach_region_metadata(df, region_table)
    return df.reindex(columns=TIDY_COLUMNS)
