"""Matplotlib visualization for region and cohort statistics."""

from __future__ import annotations

from lightsuite.analysis.viz.cohort_plots import plot_group_division_bars
from lightsuite.analysis.viz.cord_plots import (
    plot_cord_coloc_overlap,
    plot_cord_df_subregion_heatmap,
    plot_cord_division_profile,
    plot_cord_laminae_level_bars,
    plot_cord_laminae_pct_gm_bars,
    plot_cord_segment_bars,
    plot_cord_segment_grouped_bars,
    plot_cord_structure_heatmap,
    plot_cord_structure_hemisphere_panel,
    plot_cord_structure_panel,
    plot_cord_top_regions,
)
from lightsuite.analysis.viz.io import (
    load_region_plot_table,
    resolve_region_stats_from_config,
)
from lightsuite.analysis.viz.plots import plot_division_bars, plot_lr_scatter, plot_top_region_bars

__all__ = [
    "load_region_plot_table",
    "plot_cord_coloc_overlap",
    "plot_cord_df_subregion_heatmap",
    "plot_cord_division_profile",
    "plot_cord_laminae_level_bars",
    "plot_cord_laminae_pct_gm_bars",
    "plot_cord_segment_bars",
    "plot_cord_segment_grouped_bars",
    "plot_cord_structure_heatmap",
    "plot_cord_structure_hemisphere_panel",
    "plot_cord_structure_panel",
    "plot_cord_top_regions",
    "plot_division_bars",
    "plot_group_division_bars",
    "plot_lr_scatter",
    "plot_top_region_bars",
    "resolve_region_stats_from_config",
]
