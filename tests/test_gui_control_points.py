"""Tests for control-point GUI support modules."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_config
from lightsuite.gui.affine import fit_affine_transform
from lightsuite.gui.brain_data import estimate_atlas_plane_index, prepare_brain_match_points_session
from lightsuite.gui.match_points_brain import (
    _boundary_overlay,
    _chooselist_slice_label,
    _layer_xy_to_volume_point,
    _volume_points_to_layer_xy,
)
from lightsuite.registration.warp import warp_volume_affine
from lightsuite.gui.chooselist import generate_control_point_list
from lightsuite.gui.control_points import ControlPointSession
from lightsuite.gui.slices import volume_index_to_image
from lightsuite.preprocess.brain import preprocess_lightsheet_volume
from lightsuite.registration.init_brain import initialize_brain_registration


def test_estimate_atlas_plane_identity() -> None:
    vol = np.zeros((20, 20, 20), dtype=np.float32)
    vol[:, :, 11] = 100.0
    row = np.array([12, 3, 1, 1], dtype=int)
    plane = estimate_atlas_plane_index(vol, row, np.eye(4), vol.shape)
    assert plane == 12


def test_boundary_overlay_follows_affine_warp() -> None:
    ann = np.zeros((24, 24, 24), dtype=np.float32)
    ann[8:16, 8:16, 8:16] = 100.0
    row = np.array([12, 3, 1, 1], dtype=int)  # axial slice through the cube
    shape = ann.shape
    identity = np.eye(4)
    shifted = identity.copy()
    shifted[0, 3] = 6.0
    warped_id = warp_volume_affine(ann, identity, shape, order=0, point_coords="array")
    warped_shift = warp_volume_affine(ann, shifted, shape, order=0, point_coords="array")
    overlay_id = _boundary_overlay(warped_id, row)
    overlay_shift = _boundary_overlay(warped_shift, row)
    assert overlay_id.shape == overlay_shift.shape
    assert overlay_id.sum() > 0
    assert overlay_shift.sum() > 0
    assert not np.allclose(overlay_id, overlay_shift)


def test_control_point_layer_xy_roundtrip() -> None:
    chooserow = np.array([12, 2, 1, 1], dtype=int)  # cut along axis 2 (X)
    stored = [_layer_xy_to_volume_point((4.0, 7.0), chooserow, timestamp=1.0)]
    xy = _volume_points_to_layer_xy(stored, chooserow)
    assert xy.shape == (1, 2)
    assert np.allclose(xy[0], [4.0, 7.0])
    again = _layer_xy_to_volume_point((float(xy[0, 0]), float(xy[0, 1])), chooserow, timestamp=2.0)
    assert again[0] == stored[0][0]
    assert again[1] == stored[0][1]
    assert again[2] == stored[0][2]


def test_chooselist_slice_label() -> None:
    chooselist = np.array([[42, 3, 1, 1], [10, 1, 2, 2]], dtype=int)
    label = _chooselist_slice_label(chooselist, 1)
    assert "42" in label
    assert "Z" in label


def test_generate_control_point_list() -> None:
    cplist = generate_control_point_list((100, 120, 80))
    assert cplist.shape[1] == 4
    assert cplist.shape[0] > 0


def test_volume_index_to_image() -> None:
    vol = np.arange(24, dtype=np.float32).reshape(4, 3, 2)
    row = np.array([2, 3, 1, 1])
    sl = volume_index_to_image(vol, row)
    assert sl.shape == (4, 3)


def test_affine_fit_identity_with_noise() -> None:
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    matrix, mse = fit_affine_transform(pts, pts)
    assert mse < 1e-6
    assert matrix.shape == (4, 4)


def test_control_point_session_point_counts() -> None:
    session = ControlPointSession.empty(np.eye(4), n_slices=3)
    session.histology_control_points[0] = [[1, 2, 3, 0], [4, 5, 6, 0]]
    session.atlas_control_points[0] = [[1, 2, 3, 0], [7, 8, 9, 0]]
    session.histology_control_points[1] = [[1, 2, 3, 0]]
    session.atlas_control_points[1] = [[1, 2, 3, 0], [4, 5, 6, 0]]
    matched, total_s, total_a = session.point_counts()
    assert matched == 2
    assert total_s == 3
    assert total_a == 4


def test_control_point_session_roundtrip(tmp_path: Path) -> None:
    session = ControlPointSession.empty(np.eye(4), n_slices=3)
    session.histology_control_points[0] = [[1.0, 2.0, 3.0, float("nan")]]
    path = tmp_path / "cp.json"
    session.save(path)
    loaded = ControlPointSession.load(path)
    assert len(loaded.histology_control_points[0]) == 1


def _write_channel_stack(path: Path, slices: list[np.ndarray]) -> None:
    tifffile.imwrite(path, np.stack(slices, axis=0), photometric="minisblack")


def test_prepare_match_points_session(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    scratch = tmp_path / "scratch"
    save = tmp_path / "results"
    save.mkdir()
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (20, 20, 20)
    yy, xx, zz = np.mgrid[0 : shape[0], 0 : shape[1], 0 : shape[2]]
    template = (255 * np.exp(-((yy - 10) ** 2 + (xx - 10) ** 2 + (zz - 10) ** 2) / 30)).astype(np.float32)
    ann = np.zeros(shape, dtype=np.float32)
    ann[5:15, 5:15, 5:15] = 1
    nib.save(nib.Nifti1Image(template, np.eye(4)), atlas_dir / "average_template_10.nii.gz")
    nib.save(nib.Nifti1Image(ann, np.eye(4)), atlas_dir / "annotation_10.nii.gz")

    yy2, xx2 = np.mgrid[0:32, 0:32]
    ring = (((yy2 - 16) ** 2 + (xx2 - 16) ** 2) > 36) & (((yy2 - 16) ** 2 + (xx2 - 16) ** 2) < 100)
    base = (150 + 100 * ring).astype(np.uint16)
    _write_channel_stack(data_dir / "ch1.tif", [base + z * 5 for z in range(6)])

    config_data = {
        "sample": {
            "name": "test",
            "source": {"path": str(data_dir), "tiff_type": "channelperfile"},
            "scratch": str(scratch),
            "save_path": str(save),
            "voxel_um": [10.0, 10.0, 10.0],
        },
        "atlas": {"provider": "allen", "resolution_um": 10, "atlas_dir": str(atlas_dir)},
        "registration": {"resolution_um": 20, "channel_primary": 1, "orientation": [1, 2, 3]},
        "detection": {"enabled": False},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(config_data), encoding="utf-8")
    cfg = load_config(config_path)
    preprocess_lightsheet_volume(cfg)
    initialize_brain_registration(cfg)

    session_path = prepare_brain_match_points_session(cfg)
    assert session_path.is_file()
    session = ControlPointSession.load(session_path)
    assert session.chooselist is not None
    assert len(session.chooselist) > 0
