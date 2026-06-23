"""Matplotlib visualization for region and cohort statistics."""

from __future__ import annotations

from lightsuite.analysis.viz.cohort_plots import plot_group_division_bars
from lightsuite.analysis.viz.io import (
    load_region_plot_table,
    resolve_region_stats_from_config,
)
from lightsuite.analysis.viz.plots import plot_division_bars, plot_lr_scatter

__all__ = [
    "load_region_plot_table",
    "plot_division_bars",
    "plot_group_division_bars",
    "plot_lr_scatter",
    "resolve_region_stats_from_config",
]
