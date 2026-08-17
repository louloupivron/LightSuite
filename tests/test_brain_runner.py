"""Tests for brain region-stats assembly."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

from lightsuite.analysis.brain_runner import finalize_brain_region_stats
from lightsuite.analysis.region_stats import TIDY_COLUMNS, parcellation_result_to_tidy
from lightsuite.config.loader import load_config
from lightsuite.export.parcellation import ParcellationResult


def _brain_config(tmp_path: Path):
    cfg_path = tmp_path / "brain.yaml"
    cfg_path.write_text(
        f"""
sample:
  name: mouse
  source:
    format: tiff_stack
    path: {tmp_path / "src"}
    tiff_type: channelperfile
  scratch: {tmp_path / "scratch"}
  save_path: {tmp_path / "results"}
  voxel_um: [5, 5, 5]
atlas:
  provider: allen
  resolution_um: 10
  atlas_dir: {tmp_path / "atlas"}
registration:
  resolution_um: 20
  channel_primary: 1
analysis:
  intensity_metrics: [median_intensity, volume_mm3]
  count_points: true
export:
  spaces: [atlas]
""",
        encoding="utf-8",
    )
    (tmp_path / "src").mkdir()
    (tmp_path / "scratch").mkdir()
    (tmp_path / "results").mkdir()
    (tmp_path / "atlas").mkdir()
    return load_config(cfg_path)


def test_finalize_brain_region_stats_merges_point_counts(tmp_path: Path) -> None:
    config = _brain_config(tmp_path)
    register_path = tmp_path / "results" / "volume_registered"
    register_path.mkdir(parents=True)

    result = ParcellationResult(
        area_ids=np.array([7, 9], dtype=np.int64),
        median_over_areas=np.array([[10.0, 20.0], [30.0, 40.0]], dtype=np.float32),
        mean_over_areas=np.array([[10.0, 20.0], [30.0, 40.0]], dtype=np.float32),
        std_over_areas=np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        variance_over_areas=np.array([[1.0, 4.0], [9.0, 16.0]], dtype=np.float32),
        volume_over_areas=np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
    )
    tidy = parcellation_result_to_tidy(
        result,
        None,
        sample="mouse",
        channel=1,
        atlas="allen",
        intensity_metrics=config.analysis.intensity_metrics,
    )

    annotation = np.zeros((4, 4, 4), dtype=np.int32)
    annotation[:, :, :2] = 7
    annotation[:, :, 2:] = 9
    np.savez_compressed(
        register_path / "cells_atlas_coords.npz",
        atlasptcoords=np.array([[1.0, 1.0, 1.0], [4.0, 4.0, 4.0]], dtype=np.float32),
    )

    atlas = MagicMock()
    atlas.brain_atlas = "allen"
    atlas.brainglobe_name = None
    transform_params = MagicMock()
    transform_params.atlas_resolution_um = 10.0

    stats = finalize_brain_region_stats(
        config,
        register_path,
        [tidy],
        annotation=annotation,
        region_table=None,
        atlas=atlas,
        transform_params=transform_params,
    )
    assert stats.combined_path is not None
    combined = pd.read_csv(stats.combined_path)
    assert set(combined["metric"]) >= {"median_intensity", "volume_mm3", "cell_count"}
    assert list(combined.columns) == TIDY_COLUMNS
