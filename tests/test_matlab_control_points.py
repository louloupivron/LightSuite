"""Tests for MATLAB atlas2histology_tform.mat import."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat, savemat

from lightsuite.gui.control_points import (
    ControlPointSession,
    default_session_path,
    load_registration_control_point_session,
)
from lightsuite.gui.match_points_brain import _volume_points_to_layer_xy
from lightsuite.gui.slices import prepare_display_slice, volume_index_to_image
from lightsuite.import_.matlab_control_points import (
    export_matlab_control_points,
    find_matlab_control_point_session,
    import_matlab_control_points,
    load_control_point_session_from_mat,
)


def _write_brain_mat(
    path: Path,
    *,
    with_ori_trans: bool = False,
) -> None:
    histology = np.empty((2, 1), dtype=object)
    atlas = np.empty((2, 1), dtype=object)
    histology[0, 0] = np.zeros((0, 3))
    atlas[0, 0] = np.zeros((0, 3))
    histology[1, 0] = np.array([[12.0, 34.0, 56.0, 0.0], [13.0, 35.0, 57.0, 0.0]])
    atlas[1, 0] = np.array([[100.0, 200.0, 50.0], [110.0, 210.0, 60.0]])
    payload: dict = {
        "atlas2histology_tform": np.eye(4),
        "histology_control_points": histology,
        "atlas_control_points": atlas,
    }
    if with_ori_trans:
        ori = np.eye(4)
        ori[0, 3] = 7.0
        payload["ori_trans"] = ori
    savemat(path, payload)


def test_load_control_point_session_from_mat(tmp_path: Path) -> None:
    mat_path = tmp_path / "atlas2histology_tform.mat"
    _write_brain_mat(mat_path)
    regopts_ori = np.eye(4)
    regopts_ori[1, 3] = 3.0

    session = load_control_point_session_from_mat(mat_path, original_trans=regopts_ori)
    matched, _, _ = session.point_counts()
    assert matched == 2
    assert np.allclose(session.ori_trans, regopts_ori)

    atlas_pts, sample_pts = session.paired_points_xyz()
    assert atlas_pts.shape == (2, 3)
    # Stored GUI coords; paired_points_xyz applies the registration [2 1 3] swap.
    assert np.allclose(atlas_pts[0], [200.0, 100.0, 50.0])
    assert np.allclose(sample_pts[0], [34.0, 12.0, 56.0])
    assert np.allclose(session.atlas_control_points[1][0][:3], [100.0, 200.0, 50.0])


def test_load_mat_uses_file_ori_trans_when_present(tmp_path: Path) -> None:
    mat_path = tmp_path / "atlas2histology_tform.mat"
    _write_brain_mat(mat_path, with_ori_trans=True)
    regopts_ori = np.eye(4)
    regopts_ori[1, 3] = 99.0

    session = load_control_point_session_from_mat(mat_path, original_trans=regopts_ori)
    assert session.ori_trans[0][3] == 7.0


def test_find_matlab_control_point_session_glob(tmp_path: Path) -> None:
    _write_brain_mat(tmp_path / "custom_brain_tform.mat")
    found = find_matlab_control_point_session(tmp_path)
    assert found is not None
    assert found.name == "custom_brain_tform.mat"


def test_load_registration_prefers_json_over_mat(tmp_path: Path) -> None:
    _write_brain_mat(tmp_path / "atlas2histology_tform.mat")
    json_session = ControlPointSession.empty(np.eye(4), n_slices=1)
    json_session.histology_control_points = [[[1.0, 2.0, 3.0]]]
    json_session.atlas_control_points = [[[4.0, 5.0, 6.0]]]
    json_session.save(default_session_path(tmp_path))

    loaded = load_registration_control_point_session(tmp_path, original_trans=np.eye(4))
    assert len(loaded.histology_control_points[0]) == 1


def test_load_registration_falls_back_to_mat(tmp_path: Path) -> None:
    _write_brain_mat(tmp_path / "atlas2histology_tform.mat")
    loaded = load_registration_control_point_session(tmp_path, original_trans=np.eye(4))
    matched, _, _ = loaded.point_counts()
    assert matched == 2


def test_import_matlab_control_points_writes_json(tmp_path: Path) -> None:
    mat_path = tmp_path / "atlas2histology_tform.mat"
    _write_brain_mat(mat_path)
    out = import_matlab_control_points(mat_path)
    assert out.is_file()
    assert out.suffix == ".json"
    roundtrip = ControlPointSession.load(out)
    matched, _, _ = roundtrip.point_counts()
    assert matched == 2


def _write_slice_mat(path: Path, cut_axis: int, slice_index: int) -> None:
    """MATLAB file with one annotated cell placed on ``cut_axis`` at ``slice_index``."""
    point = [45.3, 55.7, 62.4]
    point[cut_axis - 1] = float(slice_index)

    histology = np.empty((2, 1), dtype=object)
    atlas = np.empty((2, 1), dtype=object)
    histology[0, 0] = np.zeros((0, 4))
    atlas[0, 0] = np.zeros((0, 4))
    histology[1, 0] = np.array([[*point, 0.0]])
    atlas[1, 0] = np.array([[*point, 0.0]])
    savemat(
        path,
        {
            "atlas2histology_tform": np.eye(4),
            "histology_control_points": histology,
            "atlas_control_points": atlas,
        },
    )


def test_matlab_import_recovers_chooselist_from_points(tmp_path: Path) -> None:
    """MATLAB writes cpt(idim) = chooselist(slice, 1), so the row is recoverable."""
    mat_path = tmp_path / "atlas2histology_tform.mat"
    _write_slice_mat(mat_path, cut_axis=2, slice_index=136)

    session = load_control_point_session_from_mat(mat_path)
    assert session.chooselist is not None
    assert session.chooselist[1][:2] == [136, 2]
    # Blanking flags are not stored by MATLAB and must stay "unknown".
    assert session.chooselist[1][2:] == [0, 0]


def test_matlab_import_places_points_on_slice(tmp_path: Path) -> None:
    """Imported MATLAB GUI coords must map inside the displayed slice."""
    volume_shape = (120, 140, 100)
    mat_path = tmp_path / "atlas2histology_tform.mat"
    _write_slice_mat(mat_path, cut_axis=2, slice_index=70)

    session = load_control_point_session_from_mat(mat_path)
    assert session.point_coord_source == "matlab"
    assert session.coord_schema_version == 2
    chooserow = np.asarray(session.chooselist[1], dtype=int)
    vol = np.zeros(volume_shape, dtype=np.float32)
    slice_shape = tuple(volume_index_to_image(vol, chooserow).shape)

    xy = _volume_points_to_layer_xy(
        session.histology_control_points[1],
        chooserow,
        slice_shape=slice_shape,
        atlas_provider="perens",
    )
    assert xy.shape == (1, 2)
    displayed = prepare_display_slice(
        volume_index_to_image(vol, chooserow), int(chooserow[1]), "perens"
    )
    disp_row, disp_col = float(xy[0, 0]), float(xy[0, 1])
    assert 0 <= disp_row < displayed.shape[0]
    assert 0 <= disp_col < displayed.shape[1]


def test_missing_fields_raise(tmp_path: Path) -> None:
    mat_path = tmp_path / "bad.mat"
    savemat(mat_path, {"atlas2histology_tform": np.eye(4)})
    with pytest.raises(ValueError, match="histology_control_points"):
        load_control_point_session_from_mat(mat_path)


def test_export_matlab_roundtrip(tmp_path: Path) -> None:
    mat_path = tmp_path / "atlas2histology_tform.mat"
    _write_brain_mat(mat_path, with_ori_trans=True)

    json_path = import_matlab_control_points(mat_path)
    exported = export_matlab_control_points(json_path)
    assert exported.is_file()

    original = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    roundtrip = loadmat(exported, squeeze_me=True, struct_as_record=False)

    assert np.allclose(original["atlas2histology_tform"], roundtrip["atlas2histology_tform"])
    assert np.allclose(original["ori_trans"], roundtrip["ori_trans"])

    for key in ("histology_control_points", "atlas_control_points"):
        orig_cells = np.asarray(original[key], dtype=object).ravel()
        rt_cells = np.asarray(roundtrip[key], dtype=object).ravel()
        assert len(orig_cells) == len(rt_cells)
        for orig, rt in zip(orig_cells, rt_cells, strict=True):
            assert np.allclose(orig, rt, equal_nan=True)
