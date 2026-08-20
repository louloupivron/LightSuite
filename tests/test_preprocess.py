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


def test_preprocess_planeperfile_multi_channel(tmp_path: Path) -> None:
    ch1 = tmp_path / "channel_488"
    ch2 = tmp_path / "channel_561"
    ch1.mkdir()
    ch2.mkdir()
    scratch = tmp_path / "scratch"
    save = tmp_path / "results"
    save.mkdir()

    for z in range(3):
        plane = (np.arange(12, dtype=np.uint16).reshape(3, 4) + z * 10)
        tifffile.imwrite(ch1 / f"plane_{z:03d}.tif", plane)
        tifffile.imwrite(ch2 / f"plane_{z:03d}.tif", plane + 50)

    config_data = {
        "sample": {
            "name": "test",
            "source": {
                "format": "tiff_stack",
                "tiff_type": "planeperfile",
                "channels": [str(ch1), str(ch2)],
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [10.0, 10.0, 10.0],
        },
        "registration": {
            "resolution_um": 20,
            "channel_primary": 1,
            "channel_secondary": 2,
        },
        "detection": {"enabled": False},
        "compute": {"workers": 2},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")

    cfg = load_config(config_path)
    assert cfg.sample.source.path == ch1.resolve()
    result = preprocess_lightsheet_volume(cfg)

    assert result.checkpoint.nchans == 2
    assert result.checkpoint.nz == 3
    assert result.checkpoint.regvolpath_secondary is not None
    assert (save / "chan_1_sample_register_20um.tif").is_file()
    assert (save / "chan_2_sample_register_20um.tif").is_file()


def test_preprocess_skips_cached_tiffs_when_yaml_metadata_changes(tmp_path: Path) -> None:
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
        "registration": {"resolution_um": 20, "channel_primary": 1, "channel_secondary": 2},
        "detection": {"enabled": False},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")

    cfg = load_config(config_path)
    preprocess_lightsheet_volume(cfg)

    ch1_mtime = (save / "chan_1_sample_register_20um.tif").stat().st_mtime
    ch2_mtime = (save / "chan_2_sample_register_20um.tif").stat().st_mtime

    config_data["registration"]["channel_primary"] = 2
    config_data["registration"]["channel_secondary"] = 1
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_config(config_path)
    result = preprocess_lightsheet_volume(cfg)

    assert Path(result.checkpoint.regvolpath).name == "chan_2_sample_register_20um.tif"
    assert (save / "chan_1_sample_register_20um.tif").stat().st_mtime == ch1_mtime
    assert (save / "chan_2_sample_register_20um.tif").stat().st_mtime == ch2_mtime


def test_preprocess_reuses_matching_volume_when_fingerprint_differs(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    save = tmp_path / "results"
    save.mkdir()

    slice_a = (np.arange(24, dtype=np.uint16).reshape(4, 6) * 10 + 100)
    _write_channel_stack(data_dir / "ch1.tif", [slice_a, slice_a + 1])

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
    preprocess_lightsheet_volume(cfg)
    cached = save / "chan_1_sample_register_20um.tif"
    mtime = cached.stat().st_mtime

    import json

    regopts = json.loads((save / "regopts.json").read_text(encoding="utf-8"))
    regopts["preprocess_fingerprint"]["tiff_type"] = "planeperfile"
    (save / "regopts.json").write_text(json.dumps(regopts), encoding="utf-8")

    result = preprocess_lightsheet_volume(load_config(config_path))
    assert cached.stat().st_mtime == mtime
    assert Path(result.checkpoint.regvolpath) == cached


def test_preprocess_force_redoes_downsample(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    save = tmp_path / "results"
    save.mkdir()

    slice_a = (np.arange(24, dtype=np.uint16).reshape(4, 6) * 10 + 100)
    _write_channel_stack(data_dir / "ch1.tif", [slice_a, slice_a + 1])

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
    preprocess_lightsheet_volume(cfg)
    first_mtime = (save / "chan_1_sample_register_20um.tif").stat().st_mtime

    preprocess_lightsheet_volume(cfg, force=True)
    second_mtime = (save / "chan_1_sample_register_20um.tif").stat().st_mtime
    assert second_mtime >= first_mtime
