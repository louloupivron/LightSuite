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
