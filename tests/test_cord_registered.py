"""Tests for registered spinal cord volume discovery helpers."""

from __future__ import annotations

import json
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


def test_discover_cord_sample_space_paths(tmp_path: Path) -> None:
    from lightsuite.config.models import CordAtlasConfig, CordRegistrationConfig, CordSampleConfig, SpinalCordPipelineConfig
    from lightsuite.export.cord_sample_space import (
        ANNOTATION_IN_SAMPLE,
        MANIFEST_NAME,
        TEMPLATE_IN_SAMPLE,
        discover_cord_sample_space_paths,
        load_cord_sample_space_volumes,
    )
    from lightsuite.io.tiff_write import save_registration_volume

    save = tmp_path / "registered"
    out = save / "volume_registered" / "sample_space"
    out.mkdir(parents=True)
    shape = (12, 16, 20)
    save_registration_volume(np.zeros(shape, dtype=np.uint16), out / "chan_01_sample_straight_20um.tif")
    save_registration_volume(np.ones(shape, dtype=np.uint16), out / ANNOTATION_IN_SAMPLE)
    save_registration_volume(np.full(shape, 2, dtype=np.uint16), out / TEMPLATE_IN_SAMPLE)
    (out / MANIFEST_NAME).write_text(
        json.dumps(
            {
                "channel_paths": {"1": str(out / "chan_01_sample_straight_20um.tif")},
                "annotation_path": str(out / ANNOTATION_IN_SAMPLE),
                "template_path": str(out / TEMPLATE_IN_SAMPLE),
            }
        ),
        encoding="utf-8",
    )

    cfg = SpinalCordPipelineConfig(
        sample=CordSampleConfig(
            name="t",
            source={"format": "tiff_stack", "path": str(tmp_path)},
            scratch=tmp_path / "scratch",
            save_path=save,
            voxel_um=[20.0, 20.0, 20.0],
        ),
        atlas=CordAtlasConfig(atlas_dir=tmp_path),
        registration=CordRegistrationConfig(),
    )
    paths = discover_cord_sample_space_paths(cfg)
    assert paths.channel_paths == {1: (out / "chan_01_sample_straight_20um.tif").resolve()}
    assert paths.annotation_path.name == ANNOTATION_IN_SAMPLE

    volumes = load_cord_sample_space_volumes(cfg, paths=paths)
    assert volumes.channels[1].shape == shape
    assert volumes.annotation.shape == shape
    assert volumes.template.shape == shape
