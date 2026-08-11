"""Shared helpers for parcellation and per-sample region tables.

Used by ``brain export``, ``spinal region-stats``, and annotation import — not a
standalone CLI module.
"""

from __future__ import annotations

from lightsuite.analysis.counts import count_points_in_regions, load_atlas_points
from lightsuite.analysis.hemisphere import (
    SIDE_LABELS,
    hemisphere_side_masks,
    hemisphere_side_volume,
)
from lightsuite.analysis.ontology import RegionTable, load_region_table
from lightsuite.analysis.region_stats import (
    METRICS,
    TIDY_COLUMNS,
    concat_tidy,
    parcellation_result_to_tidy,
    tidy_to_wide,
    write_region_stats_csv,
)

__all__ = [
    "METRICS",
    "SIDE_LABELS",
    "TIDY_COLUMNS",
    "RegionTable",
    "concat_tidy",
    "count_points_in_regions",
    "hemisphere_side_masks",
    "hemisphere_side_volume",
    "load_atlas_points",
    "load_region_table",
    "parcellation_result_to_tidy",
    "tidy_to_wide",
    "write_region_stats_csv",
]
