"""Unit tests for elastix helper modules."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lightsuite.registration.elastix.mhd import read_mhd_spacing, write_mhd
from lightsuite.registration.elastix.runner import (
    _discover_transformix_result,
    _patch_transformix_params,
)
from lightsuite.registration.elastix.params import build_bspline_params, write_parameter_file
from lightsuite.registration.elastix.points import voxel_points_to_physical, write_landmark_file
from lightsuite.registration.points_utils import thin_point_list
from lightsuite.registration.warp import warp_volume_affine


def test_write_and_read_mhd_roundtrip(tmp_path: Path) -> None:
    vol = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    base = tmp_path / "testvol"
    mhd = write_mhd(vol, base, 0.02)
    spacing = read_mhd_spacing(mhd)
    assert np.allclose(spacing, [0.02, 0.02, 0.02])
    assert mhd.with_suffix(".raw").is_file()


def test_build_bspline_params_single_channel() -> None:
    params = build_bspline_params(
        dual_channel=False,
        control_point_weight=0.2,
        n_histogram_bins=48,
        bspline_spatial_scale_mm=0.64,
        fixed_shape=(20, 20, 20),
        spacing_mm=0.02,
    )
    assert params["Registration"] == "MultiMetricMultiResolutionRegistration"
    assert len(params["Metric"]) == 2
    assert params["Metric1Weight"] == 0.2


def test_write_parameter_file(tmp_path: Path) -> None:
    params = build_bspline_params(
        dual_channel=False,
        control_point_weight=0.1,
        n_histogram_bins=32,
        bspline_spatial_scale_mm=0.64,
        fixed_shape=(10, 10, 10),
        spacing_mm=0.02,
        use_multistep=False,
    )
    path = tmp_path / "params.txt"
    write_parameter_file(path, params)
    text = path.read_text(encoding="utf-8")
    assert "MultiMetricMultiResolutionRegistration" in text
    assert "Metric1Weight" in text


def test_voxel_points_to_physical() -> None:
    pts = np.array([[1.0, 1.0, 1.0], [2.0, 3.0, 4.0]])
    phys = voxel_points_to_physical(pts, 0.02)
    assert np.allclose(phys[0], 0.0)
    assert np.allclose(phys[1], [0.02, 0.04, 0.06])


def test_write_landmark_file(tmp_path: Path) -> None:
    pts = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    path = tmp_path / "fixed.txt"
    write_landmark_file(path, pts)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "point"
    assert lines[1] == "2"


def test_thin_point_list() -> None:
    pts = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [5.0, 5.0, 5.0]])
    kept = thin_point_list(pts, min_distance=1.0)
    assert kept.tolist() == [True, False, True]


def test_patch_transformix_params_forces_mhd_output() -> None:
    params = "(FinalBSplineInterpolationOrder 3)\n"
    patched = _patch_transformix_params(params, nearest=True)
    assert '(WriteResultImage "true")' in patched
    assert '(ResultImageFormat "mhd")' in patched
    assert '(FinalBSplineInterpolationOrder "0")' in patched


def test_patch_transformix_params_replaces_duplicate_keys() -> None:
    params = '(ResultImageFormat "nii.gz")\n(ResultImagePixelType "short")\n'
    patched = _patch_transformix_params(params, nearest=False)
    assert patched.lower().count("resultimageformat") == 1
    assert '(ResultImageFormat "mhd")' in patched
    assert '(ResultImagePixelType "float")' in patched


def test_discover_transformix_result_nii_and_mhd(tmp_path: Path) -> None:
    (tmp_path / "result.nii.gz").write_bytes(b"")
    assert _discover_transformix_result(tmp_path).name == "result.nii.gz"
    (tmp_path / "result.nii.gz").unlink()
    mhd = tmp_path / "result.0.mhd"
    mhd.write_text("x", encoding="utf-8")
    assert _discover_transformix_result(tmp_path) == mhd


def test_warp_volume_affine_identity() -> None:
    vol = np.ones((4, 5, 6), dtype=np.float32)
    out = warp_volume_affine(vol, np.eye(4), vol.shape, order=0)
    assert out.shape == vol.shape
    assert np.allclose(out, 1.0)
