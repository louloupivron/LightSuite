"""Tests for fast preprocess slice operations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
from skimage.transform import resize

from lightsuite.preprocess.slice_ops import (
    background_fill_fast,
    output_xy_shape,
    output_z_count,
    process_slice_job,
    resize_xy_fast,
    write_z_downsampled_volume,
    z_downsample_plane,
)
from lightsuite.preprocess.slice_ops import SliceLoadJob


def test_resize_xy_fast_output_shape() -> None:
    plane = np.arange(80 * 120, dtype=np.uint16).reshape(80, 120)
    scale = 0.263
    got = resize_xy_fast(plane, scale)
    assert got.shape == output_xy_shape(80, 120, scale)
    assert got.dtype == np.uint16


def test_background_fill_fast_zeros() -> None:
    plane = np.zeros((4, 4), dtype=np.uint16)
    plane[1, 1] = 42
    filled = background_fill_fast(plane)
    assert filled[0, 0] == 42


def test_z_downsample_streaming(tmp_path: Path) -> None:
    vol = np.arange(24, dtype=np.uint16).reshape(2, 3, 4)
    scale_z = 0.5
    expected = resize(
        vol,
        (2, 3, output_z_count(4, scale_z)),
        order=1,
        preserve_range=True,
        anti_aliasing=False,
    ).astype(np.uint16)

    out_path = tmp_path / "down.tif"
    write_z_downsampled_volume(vol, out_path, scale_z)
    with tifffile.TiffFile(out_path) as tif:
        pages = np.stack([p.asarray() for p in tif.pages], axis=0)
    got = np.moveaxis(pages, 0, -1)
    assert got.shape == expected.shape
    assert np.allclose(got, expected, atol=2)


def test_preprocess_plane_per_file(tmp_path: Path) -> None:
    from lightsuite.config.loader import load_config
    from lightsuite.preprocess.brain import preprocess_lightsheet_volume
    import yaml

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    save = tmp_path / "results"
    save.mkdir()

    for z in range(3):
        plane = (np.arange(12, dtype=np.uint16).reshape(3, 4) + z * 10)
        tifffile.imwrite(data_dir / f"plane_{z:03d}.tif", plane)

    config_data = {
        "sample": {
            "name": "test",
            "source": {
                "format": "tiff_stack",
                "path": str(data_dir),
                "tiff_type": "planeperfile",
            },
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [10.0, 10.0, 10.0],
        },
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "detection": {"enabled": False},
        "compute": {"workers": 2},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")

    cfg = load_config(config_path)
    result = preprocess_lightsheet_volume(cfg)

    assert result.checkpoint.nchans == 1
    assert result.checkpoint.nz == 3
    reg_path = Path(result.checkpoint.regvolpath)
    assert reg_path.is_file()
    assert (save / "regopts.json").is_file()
