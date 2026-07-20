"""Tests for spinal cord import inspect discovery."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_spinal_config
from lightsuite.gui.inspect_cord_imports import (
    discover_cord_import_inspect_paths,
    load_cord_import_inspect_volumes,
)


def _write_minimal_config(tmp_path: Path) -> Path:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    (atlas_dir / "Template.tif").write_bytes(b"")
    (atlas_dir / "Annotation.tif").write_bytes(b"")
    (atlas_dir / "Segments.csv").write_text("Segment,Start,End\nC1,1,2\n", encoding="utf-8")
    (atlas_dir / "Atlas_Regions.csv").write_text(
        "id,name,acronym,parent_ID,parent_acronym\n1,a,1Sp,90,DH\n",
        encoding="utf-8",
    )

    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()
    tifffile.imwrite(sample_dir / "ch.tif", np.zeros((8, 8, 4), dtype=np.uint16))

    save = tmp_path / "registered"
    save.mkdir()
    (save / "transform_params.json").write_text("{}", encoding="utf-8")

    vr = save / "volume_registered"
    vr.mkdir()
    tifffile.imwrite(vr / "chan01_channel1.tiff", np.zeros((8, 8, 4), dtype=np.uint16))
    tifffile.imwrite(vr / "annotation_registered.tiff", np.ones((8, 8, 4), dtype=np.uint16))
    tifffile.imwrite(vr / "template_registered.tiff", np.zeros((8, 8, 4), dtype=np.float32))
    np.savez_compressed(
        vr / "cells_atlas_coords.npz",
        atlasptcoords=np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]),
        sampleptcoords=np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]),
    )

    config_data = {
        "sample": {
            "name": "fixture_cord",
            "source": {"path": str(sample_dir), "tiff_type": "channelperfile"},
            "scratch": str(tmp_path / "scratch"),
            "save_path": str(save),
            "voxel_um": [20.0, 20.0, 20.0],
        },
        "atlas": {"atlas_dir": str(atlas_dir)},
    }
    config_path = tmp_path / "spinal.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    return config_path


def test_discover_cord_import_inspect_paths(tmp_path: Path) -> None:
    config_path = _write_minimal_config(tmp_path)
    cfg = load_spinal_config(config_path)
    paths = discover_cord_import_inspect_paths(cfg)
    assert paths.volume_registered_dir.is_dir()
    assert "cells" in paths.point_npz_paths
    assert 1 in paths.registered_channels


def test_load_cord_import_inspect_volumes_headless(tmp_path: Path) -> None:
    config_path = _write_minimal_config(tmp_path)
    cfg = load_spinal_config(config_path)
    paths = discover_cord_import_inspect_paths(cfg)
    volumes = load_cord_import_inspect_volumes(cfg, paths=paths)
    assert volumes.annotation.shape == (8, 8, 4)
    assert "cells" in volumes.point_layers
    assert volumes.point_layers["cells"].shape == (2, 3)
