"""Tests for spinal cord match-points GUI helpers."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tifffile
import yaml

from lightsuite.config.loader import load_spinal_config
from lightsuite.gui.chooselist import generate_cord_control_point_list
from lightsuite.gui.control_points import ControlPointSession, through_axis_column_for_cut
from lightsuite.gui.cord_data import (
    default_cord_session_path,
    load_cord_match_points_data,
    prepare_cord_match_points_session,
    slice_pair,
)
from lightsuite.gui.match_points_cord import (
    _layer_xy_to_volume_point,
    _pair_status,
    _volume_points_to_layer_xy,
)
from lightsuite.gui.slices import volume_index_to_image


def test_generate_cord_control_point_list_shape_and_axis() -> None:
    chooselist = generate_cord_control_point_list(50, n_slices=100)
    assert chooselist.shape == (100, 4)
    assert np.all(chooselist[:, 1] == 3)
    assert chooselist[:, 0].min() >= 1
    assert chooselist[:, 0].max() <= 50


def test_generate_cord_control_point_list_is_deterministic() -> None:
    a = generate_cord_control_point_list(40)
    b = generate_cord_control_point_list(40)
    assert np.array_equal(a, b)


def test_cord_layer_xy_roundtrip() -> None:
    vol = np.zeros((24, 24, 30), dtype=np.float32)
    chooserow = np.array([12, 3, 1, 1], dtype=int)
    slice_shape = tuple(volume_index_to_image(vol, chooserow).shape)
    stored = [
        _layer_xy_to_volume_point(
            (4.0, 7.0),
            chooserow,
            slice_shape=slice_shape,
            timestamp=1.0,
        )
    ]
    xy = _volume_points_to_layer_xy(stored, chooserow, slice_shape=slice_shape)
    assert xy.shape == (1, 2)
    assert np.allclose(xy[0], [4.0, 7.0])


def test_pair_status_messages() -> None:
    assert "atlas point #2" in _pair_status(2, 1)
    assert "sample point #2" in _pair_status(1, 2)
    assert "matched" in _pair_status(3, 3)


def _write_minimal_cord_checkpoint(tmp_path: Path) -> Path:
    save_path = tmp_path / "out"
    cache = save_path / "cache"
    cache.mkdir(parents=True)
    shape = (20, 24, 30)
    straightvol = np.random.default_rng(0).random(shape).astype(np.float32)
    tv = np.random.default_rng(1).random(shape).astype(np.float32)
    av = np.zeros(shape, dtype=np.uint16)
    av[8:14, 8:16, :] = 42

    straightvol_path = cache / "straightvol.tif"
    tv_path = cache / "tv.tif"
    av_path = cache / "av.tif"
    tifffile.imwrite(straightvol_path, straightvol)
    tifffile.imwrite(tv_path, tv)
    tifffile.imwrite(av_path, av)

    regopts = {
        "sample_name": "t",
        "data_folder": str(tmp_path),
        "lsfolder": str(save_path),
        "orisize": list(shape),
        "nchans": 1,
        "sampleres_um": [20.0, 20.0, 20.0],
        "registrationres_um": [20.0, 20.0, 20.0],
        "reg_channel": 1,
        "sample_perm": [1, 2, 3],
        "tofliprc": False,
        "ikeeprange": [1, shape[2]],
        "xrange": [1, shape[1]],
        "yrange": [1, shape[0]],
        "regvol_path": str(straightvol_path),
        "tv_path": str(tv_path),
        "av_path": str(av_path),
        "smpts_path": str(tmp_path / "sm.npy"),
        "tvpts_path": str(tmp_path / "tv.npy"),
        "atlas_res_um": [20.0, 10.0, 10.0],
        "segments_path": str(tmp_path / "seg.csv"),
        "regions_path": str(tmp_path / "reg.csv"),
        "tiff_type": "channelperfile",
        "straightvol_path": str(straightvol_path),
        "affine_atlas_to_samp": np.eye(4).tolist(),
    }
    np.save(tmp_path / "sm.npy", np.zeros((0, 3), dtype=np.float32))
    np.save(tmp_path / "tv.npy", np.zeros((0, 3), dtype=np.float32))
    save_path.mkdir(parents=True, exist_ok=True)
    regopts_path = save_path / "regopts.json"
    regopts_path.write_text(json.dumps(regopts, indent=2), encoding="utf-8")

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {
                    "name": "t",
                    "save_path": str(save_path),
                    "scratch": str(tmp_path / "scratch"),
                    "voxel_um": [20.0, 20.0, 20.0],
                    "source": {"format": "tiff_stack", "path": str(data_dir)},
                },
                "atlas": {"atlas_dir": str(atlas_dir)},
                "registration": {"resolution_um": 20, "channel_primary": 1},
            }
        ),
        encoding="utf-8",
    )
    return config_path


def test_load_cord_match_points_data_matching_shapes(tmp_path: Path) -> None:
    config_path = _write_minimal_cord_checkpoint(tmp_path)
    cfg = load_spinal_config(config_path)
    data = load_cord_match_points_data(cfg)
    assert data.sample_volume.shape == data.atlas_template.shape == data.atlas_annotation.shape
    assert data.chooselist.shape[0] == len(data.session.histology_control_points) == 100
    sample, atlas = slice_pair(data, 1)
    assert sample.shape == atlas.shape


def test_prepare_cord_match_points_session_headless(tmp_path: Path) -> None:
    config_path = _write_minimal_cord_checkpoint(tmp_path)
    cfg = load_spinal_config(config_path)
    path = prepare_cord_match_points_session(cfg)
    assert path == default_cord_session_path(cfg.sample.save_path)
    session = ControlPointSession.load(path)
    assert len(session.histology_control_points) == 100


def test_through_axis_column_for_cut() -> None:
    assert through_axis_column_for_cut(3) == 2
    assert through_axis_column_for_cut(1) == 1
    assert through_axis_column_for_cut(2) == 0


def test_cord_manual_alignment_constrains_through_axis() -> None:
    """Z offsets from atlas-plane scrolling must not shear the overlay affine."""
    session = ControlPointSession.empty(np.eye(4), n_slices=16)
    sample_z = [861.0, 65.0, 54.0, 585.0, 362.0, 712.0, 596.0]
    atlas_z = [911.0, 226.0, 219.0, 634.0, 451.0, 762.0, 645.0]
    for i in range(16):
        z_sample = sample_z[i % len(sample_z)]
        z_atlas = atlas_z[i % len(atlas_z)]
        session.histology_control_points[i] = [[60.0 + i, 80.0 + i, z_sample, 1.0]]
        session.atlas_control_points[i] = [[61.0 + i, 81.0 + i, z_atlas, 1.0]]

    session.update_manual_alignment(min_pairs=16, constrain_cut_axis=3)
    constrained = np.asarray(session.atlas2histology_tform, dtype=float)
    assert abs(constrained[2, 3]) < 1.0


def test_paired_points_xyz_dim_order() -> None:
    session = ControlPointSession.empty(np.eye(4), n_slices=2)
    session.histology_control_points[0] = [[10.0, 20.0, 5.0, 1.0]]
    session.atlas_control_points[0] = [[11.0, 21.0, 6.0, 1.0]]
    atlas_pts, sample_pts = session.paired_points_xyz()
    assert atlas_pts.shape == (1, 3)
    assert sample_pts.shape == (1, 3)
    assert np.allclose(atlas_pts[0], [21.0, 11.0, 6.0])
    assert np.allclose(sample_pts[0], [20.0, 10.0, 5.0])


def test_paired_points_volume_yxz_zero_based() -> None:
    session = ControlPointSession.empty(np.eye(4), n_slices=2)
    session.histology_control_points[0] = [[10.0, 20.0, 5.0, 1.0]]
    session.atlas_control_points[0] = [[11.0, 21.0, 6.0, 1.0]]
    atlas_pts, sample_pts = session.paired_points_volume_yxz()
    assert np.allclose(sample_pts[0], [10.0, 20.0, 4.0])
    assert np.allclose(atlas_pts[0], [11.0, 21.0, 5.0])
