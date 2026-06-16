"""Tests for mesoSPIM pipeline configuration."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_mesospim_config
from lightsuite.mesospim.config_models import MesospimPipelineConfig
from lightsuite.mesospim.meta import meta_path_for_tiff


def _write_stack(path: Path, shape_zyx: tuple[int, int, int]) -> None:
    arr = np.zeros(shape_zyx, dtype=np.uint16)
    tifffile.imwrite(path, arr, imagej=True)


def _write_meta(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "[Pixelsize in um] 5.0",
                "[x_pos] 0.0",
                "[y_pos] 0.0",
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


def test_load_mesospim_config(tmp_path: Path) -> None:
    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    _write_stack(overview, (5, 16, 16))
    _write_stack(roi, (5, 16, 16))
    _write_meta(meta_path_for_tiff(overview))
    _write_meta(meta_path_for_tiff(roi))

    config_path = tmp_path / "mesospim.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {
                    "name": "test_sample",
                    "save_path": str(tmp_path / "out"),
                },
                "mesospim": {
                    "overview": {"path": str(overview)},
                    "roi": {"path": str(roi)},
                },
            }
        ),
        encoding="utf-8",
    )

    cfg = load_mesospim_config(config_path)
    assert isinstance(cfg, MesospimPipelineConfig)
    assert cfg.sample.name == "test_sample"
    assert cfg.mesospim.registration.elastix_stages == ["translation", "rigid"]
