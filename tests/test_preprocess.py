"""Tests for brain lightsheet preprocessing."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_config
from lightsuite.preprocess.brain import preprocess_lightsheet_volume


def _write_channel_stack(path: Path, slices: list[np.ndarray]) -> None:
    tifffile.imwrite(path, np.stack(slices, axis=0), photometric="minisblack")


def test_preprocess_channel_per_file(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    save = tmp_path / "results"
    save.mkdir()

    slice_a = (np.arange(24, dtype=np.uint16).reshape(4, 6) * 10 + 100)
    slice_b = slice_a + 50
    _write_channel_stack(data_dir / "ch1.tif", [slice_a, slice_b])
    _write_channel_stack(data_dir / "ch2.tif", [slice_a + 1, slice_b + 1])

    config_data = {
        "sample": {
            "name": "test",
            "source": {
                "format": "tiff_stack",
                "path": str(data_dir),
                "tiff_type": "channelperfile",
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [10.0, 10.0, 10.0],
        },
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "detection": {"enabled": False},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")

    cfg = load_config(config_path)
    result = preprocess_lightsheet_volume(cfg)

    assert result.checkpoint.nchans == 2
    assert result.checkpoint.nz == 2
    reg_path = Path(result.checkpoint.regvolpath)
    assert reg_path.is_file()

    with tifffile.TiffFile(reg_path) as tif:
        pages = len(tif.pages)
        assert pages >= 1
        assert tif.pages[0].shape == (2, 3)  # halved 4x6 with ceil

    regopts = save / "regopts.json"
    assert regopts.is_file()
