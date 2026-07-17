"""Tests for content-box probing and YAML updates."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import nibabel as nib
import tifffile
import yaml

from lightsuite.config.loader import load_config, save_content_box_to_config
from lightsuite.config.models import ContentTrimMode
from lightsuite.registration.content_bbox import ContentBox
from lightsuite.registration.content_probe import (
    ContentProbeTarget,
    auto_detect_atlas_box,
    format_content_box_report,
    load_content_probe_data,
    preview_box_for_display,
)


def _write_minimal_config(tmp_path: Path, *, extra: dict | None = None) -> Path:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (40, 50, 30)
    tpl = np.zeros(shape, dtype=np.float32)
    ann = np.zeros(shape, dtype=np.uint16)
    tpl[8:32, 10:40, 5:25] = 1.0
    ann[8:32, 10:40, 5:25] = 7
    tpl_path = atlas_dir / "gubra_template_olf.nii.gz"
    ann_path = atlas_dir / "gubra_ano_olf.nii.gz"
    nib.save(nib.Nifti1Image(np.ascontiguousarray(tpl), np.eye(4)), str(tpl_path))
    nib.save(nib.Nifti1Image(np.ascontiguousarray(ann), np.eye(4)), str(ann_path))

    save_path = tmp_path / "out"
    save_path.mkdir()
    source_path = tmp_path / "source"
    source_path.mkdir()
    reg_path = save_path / "chan_1_sample_register_39um.tif"
    sample = np.zeros((20, 24, 16), dtype=np.uint16)
    sample[4:16, 5:18, 3:12] = 500
    with tifffile.TiffWriter(reg_path) as writer:
        for z in range(sample.shape[2]):
            writer.write(sample[:, :, z], compression="lzw", photometric="minisblack")

    config_data = {
        "sample": {
            "name": "test",
            "source": {"path": str(source_path), "tiff_type": "channelperfile"},
            "scratch": str(tmp_path / "scratch"),
            "save_path": str(save_path),
            "voxel_um": [8.0, 8.0, 5.0],
        },
        "atlas": {
            "provider": "perens",
            "resolution_um": 39,
            "atlas_dir": str(atlas_dir),
            "source": "files",
            "content_trim": "auto",
            "content_margin_vox": 0,
        },
        "registration": {
            "resolution_um": 39,
            "channel_primary": 1,
            "sample_content_margin_vox": 1,
            "sample_content_trim_z": False,
        },
    }
    if extra:
        for key, value in extra.items():
            if isinstance(value, dict) and key in config_data and isinstance(config_data[key], dict):
                config_data[key].update(value)
            else:
                config_data[key] = value

    regopts = {
        "sample_name": "test",
        "ny": 20,
        "nx": 24,
        "nz": 16,
        "nchans": 1,
        "voxel_um": [8.0, 8.0, 5.0],
        "registres_um": 39.0,
        "regvolpath": str(reg_path),
        "regvolpath_secondary": None,
        "regvolpaths": {"1": str(reg_path)},
        "tiff_type": "channelperfile",
        "channel_primary": 1,
        "channel_secondary": None,
    }
    (save_path / "regopts.json").write_text(
        __import__("json").dumps(regopts, indent=2),
        encoding="utf-8",
    )

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    return config_path


def test_auto_detect_atlas_box(tmp_path: Path) -> None:
    config_path = _write_minimal_config(tmp_path)
    cfg = load_config(config_path)
    box = auto_detect_atlas_box(cfg)
    assert box.y0 == 8
    assert box.x0 == 10
    assert box.z0 == 5


def test_load_content_probe_data_atlas(tmp_path: Path) -> None:
    config_path = _write_minimal_config(tmp_path)
    cfg = load_config(config_path)
    data = load_content_probe_data(cfg, target=ContentProbeTarget.ATLAS)
    assert data.full_shape == (40, 50, 30)
    assert data.box.y0 == 8
    report = format_content_box_report(data)
    assert "content_trim: manual" in report
    assert "content_box:" in report


def test_load_content_probe_data_sample(tmp_path: Path) -> None:
    config_path = _write_minimal_config(tmp_path)
    cfg = load_config(config_path)
    data = load_content_probe_data(cfg, target=ContentProbeTarget.SAMPLE)
    assert data.full_shape == (20, 24, 16)
    assert data.box.y0 == 2
    assert data.box.y1 == 17
    assert data.box.x0 == 3
    assert data.box.x1 == 19


def test_preview_box_for_display_maps_to_preview_shape() -> None:
    from lightsuite.registration.content_probe import ContentProbeData

    full = ContentBox(y0=0, y1=99, x0=0, x1=49, z0=0, z1=29)
    data = ContentProbeData(
        target=ContentProbeTarget.ATLAS,
        volume=np.zeros((10, 5, 3), dtype=np.float32),
        full_shape=(100, 50, 30),
        box=full,
        source_label="test",
        preview_downsampled=True,
    )
    preview_box = preview_box_for_display(data, full)
    assert preview_box.y1 == 9
    assert preview_box.x1 == 4
    assert preview_box.z1 == 2


def test_save_content_box_to_config_atlas(tmp_path: Path) -> None:
    config_path = _write_minimal_config(tmp_path)
    save_content_box_to_config(
        config_path,
        target="atlas",
        box=[1, 10, 2, 20, 3, 15],
    )
    text = config_path.read_text(encoding="utf-8")
    assert "content_trim: manual" in text
    assert "content_box: [1, 10, 2, 20, 3, 15]" in text
    cfg = load_config(config_path)
    assert cfg.atlas.content_trim == ContentTrimMode.MANUAL
    assert cfg.atlas.content_box == [1, 10, 2, 20, 3, 15]


def test_save_content_box_to_config_sample(tmp_path: Path) -> None:
    config_path = _write_minimal_config(tmp_path)
    save_content_box_to_config(
        config_path,
        target="sample",
        box=[4, 15, 5, 17, 3, 11],
    )
    text = config_path.read_text(encoding="utf-8")
    assert "sample_content_crop: manual" in text
    assert "sample_content_box: [4, 15, 5, 17, 3, 11]" in text
    cfg = load_config(config_path)
    assert cfg.registration.sample_content_box == [4, 15, 5, 17, 3, 11]
