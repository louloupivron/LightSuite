"""Tests for pydantic config models."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from lightsuite.config.loader import load_config, save_orientation_to_config
from lightsuite.config.models import BrainPipelineConfig, DetectionConfig


def test_detection_threshold_order() -> None:
    with pytest.raises(ValueError, match="threshold"):
        DetectionConfig(thresholds=[0.3, 0.5])


def test_source_channels_requires_planeperfile(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    ch2 = tmp_path / "ch2"
    ch2.mkdir()
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "results"
    save.mkdir()

    config_data = {
        "sample": {
            "name": "test_mouse",
            "source": {
                "format": "tiff_stack",
                "tiff_type": "channelperfile",
                "channels": [str(data_dir), str(ch2)],
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [5.0, 5.0, 5.0],
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")

    with pytest.raises(ValueError, match="planeperfile"):
        load_config(config_path)


def test_source_channels_sets_path_from_first_channel(tmp_path: Path) -> None:
    ch1 = tmp_path / "ch1"
    ch2 = tmp_path / "ch2"
    ch1.mkdir()
    ch2.mkdir()
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "results"
    save.mkdir()

    config_data = {
        "sample": {
            "name": "test_mouse",
            "source": {
                "format": "tiff_stack",
                "tiff_type": "planeperfile",
                "channels": [str(ch1), str(ch2)],
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [5.0, 5.0, 5.0],
        },
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")

    cfg = load_config(config_path)
    assert cfg.sample.source.path == ch1.resolve()
    assert cfg.sample.source.channel_roots == (ch1.resolve(), ch2.resolve())


def test_load_config_from_yaml(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "results"
    save.mkdir()

    config_data = {
        "sample": {
            "name": "test_mouse",
            "source": {
                "format": "tiff_stack",
                "path": str(data_dir),
                "tiff_type": "channelperfile",
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [5.0, 5.0, 5.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")

    cfg = load_config(config_path)
    assert isinstance(cfg, BrainPipelineConfig)
    assert cfg.sample.name == "test_mouse"
    assert cfg.registration.resolution_um == 20.0


def test_save_orientation_to_config_updates_existing_field(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "out"
    save.mkdir()
    config_text = f"""# pipeline config
sample:
  name: mouse
  source:
    path: {data_dir}
    tiff_type: channelperfile
  scratch: {scratch}
  save_path: {save}
  voxel_um: [5.0, 5.0, 5.0]
registration:
  channel_primary: 1
  orientation: [1, 2, 3]
"""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(config_text, encoding="utf-8")

    save_orientation_to_config(config_path, [2, 1, -3])

    updated_text = config_path.read_text(encoding="utf-8")
    assert "# pipeline config" in updated_text
    assert "orientation: [2, 1, -3]" in updated_text
    assert load_config(config_path).registration.orientation == [2, 1, -3]


def test_save_orientation_to_config_inserts_missing_field(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "out"
    save.mkdir()
    config_text = f"""sample:
  name: mouse
  source:
    path: {data_dir}
    tiff_type: channelperfile
  scratch: {scratch}
  save_path: {save}
  voxel_um: [5.0, 5.0, 5.0]
registration:
  channel_primary: 1
"""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(config_text, encoding="utf-8")

    save_orientation_to_config(config_path, [-1, 3, 2])

    assert "orientation: [-1, 3, 2]" in config_path.read_text(encoding="utf-8")
    assert load_config(config_path).registration.orientation == [-1, 3, 2]
