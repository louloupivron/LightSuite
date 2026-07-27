"""Tests for mesoSPIM landmark match-points GUI helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_mesospim_config
from lightsuite.gui.match_points_shared import (
    layer_xy_from_zyx,
    pair_status,
    zyx_from_layer_xy,
)
from lightsuite.gui.mesospim_data import MesospimSliceSource, prepare_mesospim_match_points_session
from lightsuite.mesospim.config_models import MesospimTiffRemapConfig
from lightsuite.mesospim.landmark_session import default_landmark_session_path


def test_layer_xy_roundtrip() -> None:
    points = [[2.0, 10.0, 20.0], [5.0, 30.0, 40.0]]
    xy = layer_xy_from_zyx(points, z_index=2)
    assert xy.shape == (1, 2)
    assert np.allclose(xy[0], [10.0, 20.0])

    updated = zyx_from_layer_xy(np.array([[11.0, 21.0], [31.0, 41.0]]), 2, points)
    assert len(updated) == 3
    assert updated[0] == [5.0, 30.0, 40.0]
    assert updated[1] == [2.0, 11.0, 21.0]
    assert updated[2] == [2.0, 31.0, 41.0]


def test_pair_status_messages() -> None:
    assert "ROI point #2" in pair_status(2, 1)
    assert "overview point #2" in pair_status(1, 2)
    assert "matched" in pair_status(3, 3)


def test_prepare_mesospim_match_points_session_headless(tmp_path: Path) -> None:

    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    tifffile.imwrite(overview, np.zeros((8, 16, 16), dtype=np.uint16), imagej=True)
    tifffile.imwrite(roi, np.zeros((8, 16, 16), dtype=np.uint16), imagej=True)

    save_path = tmp_path / "out"
    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "t", "save_path": str(save_path)},
                "mesospim": {
                    "geometry_mode": "landmarks",
                    "overview": {"path": str(overview), "voxel_um": [1.0, 1.0, 1.0]},
                    "roi": {"path": str(roi), "voxel_um": [1.0, 1.0, 1.0]},
                },
            }
        ),
        encoding="utf-8",
    )
    cfg = load_mesospim_config(config_path)
    path = prepare_mesospim_match_points_session(cfg)
    assert path == default_landmark_session_path(save_path)
    assert path.is_file()


def test_mesospim_slice_source_lazy_read_and_cache(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    arr = np.arange(4 * 8 * 8, dtype=np.uint16).reshape(4, 8, 8)
    tifffile.imwrite(overview, arr, imagej=True)
    tifffile.imwrite(roi, arr, imagej=True)

    source = MesospimSliceSource.from_path(
        overview,
        overview_path=overview,
        roi_path=roi,
        remap=MesospimTiffRemapConfig(),
    )
    sl0 = source.read_display_slice(0)
    sl0_cached = source.read_display_slice(0)
    assert sl0.shape == (8, 8)
    assert sl0 is sl0_cached
    sl2 = source.read_display_slice(2)
    assert sl2.shape == (8, 8)
    assert not np.allclose(sl0, sl2)
