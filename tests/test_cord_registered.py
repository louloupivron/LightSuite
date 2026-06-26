"""Tests for registered spinal cord volume discovery helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile

from lightsuite.export.cord_registered import (
    REGISTERED_ANNOTATION_FILENAME,
    discover_registered_cord_paths,
    export_layout_from_native,
    load_registered_stack,
)


def test_export_layout_from_native_transposes_length_axis_last() -> None:
    native = np.zeros((40, 12, 16), dtype=np.uint8)
    export = export_layout_from_native(native)
    assert export.shape == (12, 16, 40)


def test_discover_registered_cord_paths(tmp_path: Path) -> None:
    vr = tmp_path / "volume_registered"
    vr.mkdir()
    tifffile.imwrite(vr / "chan01_channel1.tiff", np.zeros((10, 12, 20), dtype=np.uint16))
    tifffile.imwrite(vr / "chan02_channel2.tiff", np.zeros((10, 12, 20), dtype=np.uint16))

    from lightsuite.config.models import CordAtlasConfig, CordRegistrationConfig, CordSampleConfig, SpinalCordPipelineConfig

    cfg = SpinalCordPipelineConfig(
        sample=CordSampleConfig(
            name="t",
            source={"format": "tiff_stack", "path": str(tmp_path)},
            scratch=tmp_path / "scratch",
            save_path=tmp_path,
            voxel_um=[20.0, 20.0, 20.0],
        ),
        atlas=CordAtlasConfig(atlas_dir=tmp_path),
        registration=CordRegistrationConfig(),
    )
    paths = discover_registered_cord_paths(cfg)
    assert paths.registered_channels == {1: (vr / "chan01_channel1.tiff").resolve(), 2: (vr / "chan02_channel2.tiff").resolve()}
    assert paths.annotation_path.name == REGISTERED_ANNOTATION_FILENAME


def test_load_registered_stack(tmp_path: Path) -> None:
    path = tmp_path / "vol.tiff"
    tifffile.imwrite(path, np.ones((5, 6, 7), dtype=np.uint16))
    vol = load_registered_stack(path)
    assert vol.shape == (5, 6, 7)


def test_warp_output_to_uint16_clips_negative_values() -> None:
    from lightsuite.export.cord_registered import warp_output_to_uint16

    wrapped = warp_output_to_uint16(np.array([-10.0, 0.4, 100.6, 70000.0], dtype=np.float32))
    assert wrapped.tolist() == [0, 0, 101, 65535]


def test_load_registered_stack_rejects_non_3d(tmp_path: Path) -> None:
    path = tmp_path / "bad.tiff"
    tifffile.imwrite(path, np.ones((5, 6), dtype=np.uint16))
    with pytest.raises(ValueError, match="Expected 3D"):
        load_registered_stack(path)
