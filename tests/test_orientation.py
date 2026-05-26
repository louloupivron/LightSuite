"""Tests for brain orientation utilities and GUI prep."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from lightsuite.gui.orientation_views import mean_projections, projection_layout
from lightsuite.registration.orientation import (
    DEFAULT_PERMVEC,
    permute_for_atlas,
    permvec_from_indices,
    save_orientation,
    validate_permvec,
)


def test_validate_permvec_rejects_duplicate_axes() -> None:
    with pytest.raises(ValueError, match="unique"):
        validate_permvec([1, 1, 3])


def test_permvec_from_indices_roundtrip() -> None:
    permvec = permvec_from_indices((0, 2, 4))
    assert permvec == [1, 2, 3]


def test_permute_for_atlas_flip() -> None:
    vol = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    out = permute_for_atlas(vol, [-1, 2, 3])
    assert out.shape == (2, 3, 4)
    assert np.allclose(out[0, :, :], vol[1, :, :])
    assert np.allclose(out[1, :, :], vol[0, :, :])


def test_mean_projections_shapes() -> None:
    vol = np.ones((5, 6, 7), dtype=np.float32)
    projs = mean_projections(vol)
    assert len(projs) == 3
    assert projs[0].shape == (6, 7)
    assert projs[1].shape == (5, 7)
    assert projs[2].shape == (5, 6)


def test_projection_layout_returns_six_layers() -> None:
    atlas = [np.ones((4, 8)), np.ones((6, 5)), np.ones((3, 9))]
    sample = [np.ones((2, 4)), np.ones((3, 5)), np.ones((3, 3))]
    layers, row_gap = projection_layout(atlas, sample)
    assert len(layers) == 6
    assert row_gap == 6


def test_save_orientation(tmp_path: Path) -> None:
    path = save_orientation(tmp_path, DEFAULT_PERMVEC)
    assert path.is_file()
    values = np.loadtxt(path, dtype=int).tolist()
    assert values == DEFAULT_PERMVEC


def test_prepare_orientation_session(tmp_path: Path) -> None:
    from lightsuite.config.loader import load_config
    import nibabel as nib
    import tifffile
    import yaml

    from lightsuite.gui.orientation_brain import prepare_orientation_session
    from lightsuite.preprocess.brain import preprocess_lightsheet_volume

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    save = tmp_path / "results"
    save.mkdir()
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (12, 12, 12)
    yy, xx, zz = np.mgrid[0 : shape[0], 0 : shape[1], 0 : shape[2]]
    tpl = (255 * np.exp(-((yy - 6) ** 2 + (xx - 6) ** 2 + (zz - 6) ** 2) / 20)).astype(np.float32)
    nib.save(nib.Nifti1Image(tpl, np.eye(4)), atlas_dir / "average_template_10.nii.gz")
    nib.save(nib.Nifti1Image(np.ones(shape, dtype=np.float32), np.eye(4)), atlas_dir / "annotation_10.nii.gz")

    yy2, xx2 = np.mgrid[0:16, 0:16]
    base = (100 + 50 * (((yy2 - 8) ** 2 + (xx2 - 8) ** 2) < 36)).astype(np.uint16)
    tifffile.imwrite(data_dir / "ch1.tif", np.stack([base + z for z in range(5)], axis=0))

    config_data = {
        "sample": {
            "name": "orient",
            "source": {"path": str(data_dir), "tiff_type": "channelperfile"},
            "scratch": str(tmp_path / "scratch"),
            "save_path": str(save),
            "voxel_um": [10.0, 10.0, 10.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": str(atlas_dir)},
        "registration": {"resolution_um": 20, "channel_primary": 1},
        "detection": {"enabled": False},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_config(config_path)
    preprocess_lightsheet_volume(cfg)
    path = prepare_orientation_session(cfg)
    assert path.is_file()
