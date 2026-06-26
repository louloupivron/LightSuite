"""Tests for spinal cord volume I/O and streaming downsampling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from lightsuite.config.loader import load_spinal_config
from lightsuite.config.models import CordTiffLayout
from lightsuite.io.cord_reader import read_spinal_cord_sample


def _ensure_fixtures(root: Path) -> None:
    if not (root / "sample" / "planeperfile").is_dir():
        from tests.fixtures.spinal_cord.build_fixtures import build_fixtures

        build_fixtures(root)


def test_plane_per_file_streaming_downsample(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "results"
    save.mkdir()
    config_data = {
        "sample": {
            "name": "plane_fixture",
            "source": {
                "path": str(fixture_root / "sample" / "planeperfile"),
                "tiff_type": "planeperfile",
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [10.0, 10.0, 10.0],
        },
        "atlas": {"atlas_dir": str(fixture_root / "atlas")},
        "registration": {"resolution_um": 20},
    }
    config_path = tmp_path / "spinal_plane.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_spinal_config(config_path)

    sample = read_spinal_cord_sample(cfg)
    assert sample.layout == CordTiffLayout.PLANE_PER_FILE
    assert sample.native_orisize == (40, 30, 16)
    assert sample.volume.shape[:3] == (20, 15, 8)
    assert sample.n_channels == 1
    assert sample.volume.dtype == np.uint16


def test_auto_detect_plane_per_file(tmp_path: Path) -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    _ensure_fixtures(fixture_root)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "results"
    save.mkdir()
    config_data = {
        "sample": {
            "name": "auto_plane",
            "source": {
                "path": str(fixture_root / "sample" / "planeperfile"),
                "tiff_type": "auto",
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [20.0, 20.0, 20.0],
        },
        "atlas": {"atlas_dir": str(fixture_root / "atlas")},
    }
    config_path = tmp_path / "spinal_auto.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_spinal_config(config_path)

    sample = read_spinal_cord_sample(cfg)
    assert sample.layout == CordTiffLayout.PLANE_PER_FILE
    assert sample.volume.shape == (40, 30, 16, 1)
