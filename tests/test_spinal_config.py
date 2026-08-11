"""Tests for spinal cord config models."""

from __future__ import annotations

from pathlib import Path

import yaml

from lightsuite.config.loader import load_spinal_config


def test_load_spinal_config(tmp_path: Path) -> None:
    atlas = tmp_path / "atlas"
    atlas.mkdir()
    (atlas / "Template.tif").write_bytes(b"")
    (atlas / "Annotation.tif").write_bytes(b"")
    (atlas / "Segments.csv").write_text("Segment,Start,End\nC1,0,1\n", encoding="utf-8")
    (atlas / "Atlas_Regions.csv").write_text("id,name,children_IDs\n1,gm,1\n", encoding="utf-8")

    sample = tmp_path / "sample"
    sample.mkdir()
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    save = tmp_path / "results"
    save.mkdir()

    config_data = {
        "sample": {
            "name": "test_cord",
            "source": {"path": str(sample), "tiff_type": "channelperfile"},
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [20.0, 20.0, 20.0],
        },
        "atlas": {"atlas_dir": str(atlas)},
    }
    config_path = tmp_path / "spinal.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_spinal_config(config_path)
    assert cfg.sample.name == "test_cord"
    assert cfg.registration.resolution_um == 20.0


def test_load_spinal_config_multi_channel_folders(tmp_path: Path) -> None:
    atlas = tmp_path / "atlas"
    atlas.mkdir()
    (atlas / "Template.tif").write_bytes(b"")
    (atlas / "Annotation.tif").write_bytes(b"")
    (atlas / "Segments.csv").write_text("Segment,Start,End\nC1,0,1\n", encoding="utf-8")
    (atlas / "Atlas_Regions.csv").write_text("id,name,children_IDs\n1,gm,1\n", encoding="utf-8")

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
            "name": "multi_cord",
            "source": {
                "tiff_type": "planeperfile",
                "channels": [str(ch1), str(ch2)],
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [1.8, 1.8, 1.8],
        },
        "atlas": {"atlas_dir": str(atlas)},
        "registration": {"channel_primary": 2},
    }
    config_path = tmp_path / "spinal_multi.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_spinal_config(config_path)
    assert cfg.sample.source.path == ch1.resolve()
    assert cfg.sample.source.channel_roots == (ch1.resolve(), ch2.resolve())
    assert cfg.registration.channel_primary == 2


def test_spinal_config_channels_rejects_channelperfile(tmp_path: Path) -> None:
    import pytest
    from pydantic import ValidationError

    atlas = tmp_path / "atlas"
    atlas.mkdir()
    (atlas / "Template.tif").write_bytes(b"")
    (atlas / "Annotation.tif").write_bytes(b"")
    (atlas / "Segments.csv").write_text("Segment,Start,End\nC1,0,1\n", encoding="utf-8")
    (atlas / "Atlas_Regions.csv").write_text("id,name,children_IDs\n1,gm,1\n", encoding="utf-8")
    ch1 = tmp_path / "ch1"
    ch1.mkdir()

    config_data = {
        "sample": {
            "name": "bad",
            "source": {
                "tiff_type": "channelperfile",
                "channels": [str(ch1)],
            },
            "scratch": str(tmp_path / "scratch"),
            "save_path": str(tmp_path / "save"),
            "voxel_um": [20.0, 20.0, 20.0],
        },
        "atlas": {"atlas_dir": str(atlas)},
    }
    config_path = tmp_path / "bad.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    with pytest.raises(ValidationError, match="planeperfile"):
        load_spinal_config(config_path)


def _minimal_spinal_config(tmp_path: Path, *, analysis: dict | None = None) -> Path:
    atlas = tmp_path / "atlas"
    atlas.mkdir()
    (atlas / "Template.tif").write_bytes(b"")
    (atlas / "Annotation.tif").write_bytes(b"")
    (atlas / "Segments.csv").write_text("Segment,Start,End\nC1,0,1\n", encoding="utf-8")
    (atlas / "Atlas_Regions.csv").write_text("id,name,children_IDs\n1,gm,1\n", encoding="utf-8")

    sample = tmp_path / "sample"
    sample.mkdir()
    save = tmp_path / "results"
    save.mkdir()

    config_data: dict = {
        "sample": {
            "name": "test_cord",
            "source": {"path": str(sample), "tiff_type": "channelperfile"},
            "scratch": str(tmp_path / "scratch"),
            "save_path": str(save),
            "voxel_um": [20.0, 20.0, 20.0],
        },
        "atlas": {"atlas_dir": str(atlas)},
    }
    if analysis is not None:
        config_data["analysis"] = analysis
    config_path = tmp_path / "spinal.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    return config_path
