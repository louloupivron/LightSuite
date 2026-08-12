"""Tests for multires channel selection helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_multires_config, save_reference_channel_to_multires_config
from lightsuite.mesospim.meta import meta_path_for_tiff
from lightsuite.multires.channels import (
    default_apply_transform_to,
    multires_channel_names,
    resolved_apply_transform_to,
)


def _write_stack(path: Path) -> None:
    tifffile.imwrite(path, np.ones((5, 16, 16), dtype=np.uint16), imagej=True)


def _write_meta(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "[Pixelsize in um] 5.0",
                "[x_pos] 100.0",
                "[y_pos] 500.0",
                "[z_start] 0.0",
                "[z_end] 8.0",
                "[z_stepsize] 2.0",
                "[z_planes] 5",
                "[x_pixels] 16",
                "[y_pixels] 16",
            ]
        ),
        encoding="utf-8",
    )


def _write_multichannel_config(tmp_path: Path) -> Path:
    overview_488 = tmp_path / "overview_488.tif"
    roi_488 = tmp_path / "roi_488.tif"
    overview_555 = tmp_path / "overview_555.tif"
    roi_555 = tmp_path / "roi_555.tif"
    for path in (overview_488, roi_488, overview_555, roi_555):
        _write_stack(path)
    _write_meta(meta_path_for_tiff(overview_488))
    _write_meta(meta_path_for_tiff(roi_488))

    config_path = tmp_path / "multires.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "sample", "save_path": str(tmp_path / "out")},
                "multires": {
                    "pair_label": "test",
                    "channels": {
                        "488": {"overview": str(overview_488), "roi": str(roi_488)},
                        "555": {"overview": str(overview_555), "roi": str(roi_555)},
                    },
                    "registration": {
                        "reference_channel": "488",
                        "apply_transform_to": ["555"],
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    return config_path


def test_multires_channel_helpers(tmp_path: Path) -> None:
    cfg = load_multires_config(_write_multichannel_config(tmp_path))
    assert multires_channel_names(cfg) == ["488", "555"]
    assert default_apply_transform_to(cfg, "555") == ["488"]
    assert resolved_apply_transform_to(cfg) == ["555"]


def test_save_reference_channel_updates_yaml(tmp_path: Path) -> None:
    config_path = _write_multichannel_config(tmp_path)
    saved = save_reference_channel_to_multires_config(config_path, "555")
    text = saved.read_text(encoding="utf-8")
    assert 'reference_channel: "555"' in text
    assert 'apply_transform_to: ["488"]' in text
    cfg = load_multires_config(saved)
    assert cfg.multires.registration.reference_channel == "555"
    assert cfg.multires.registration.apply_transform_to == ["488"]
