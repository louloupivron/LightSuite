"""Post-registration analysis: region statistics, cell counts, and taxonomy.

This package consumes the atlas-space artifacts produced by ``lightsuite brain
export`` and ``lightsuite brain import-annotations`` and turns them into tidy,
cross-atlas region tables (median intensity, std, volume, cell count, density).

Foundation milestone modules:

- :mod:`lightsuite.analysis.ontology` — shared region metadata table keyed on the
  Allen ontology (works for both the Allen ABC and Perens atlases).
- :mod:`lightsuite.analysis.hemisphere` — single source of truth for the
  left/right voxel split used by both intensity stats and cell counts.
- :mod:`lightsuite.analysis.region_stats` — canonical long-form (tidy) schema and
  converters to/from the legacy wide ``chanXX_intensities.csv`` layout.
- :mod:`lightsuite.analysis.counts` — per-region cell counts and densities from
  imported atlas-space point clouds.
"""

from __future__ import annotations

from lightsuite.analysis.cohort_models import CohortConfig, CohortSampleEntry, GroupAnalysisConfig
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
    "CohortConfig",
    "CohortSampleEntry",
    "GroupAnalysisConfig",
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
