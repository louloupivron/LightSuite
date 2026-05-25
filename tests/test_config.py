"""Tests for pydantic config models."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from lightsuite.config.loader import load_config
from lightsuite.config.models import BrainPipelineConfig, DetectionConfig


def test_detection_threshold_order() -> None:
    with pytest.raises(ValueError, match="threshold"):
        DetectionConfig(thresholds=[0.3, 0.5])


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
