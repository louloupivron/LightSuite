"""Unit tests for elastix helper modules."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from lightsuite.registration.elastix.mhd import (
    read_mhd_spacing,
    read_mhd_volume,
    scale_volume_for_elastix_mi,
    write_mhd,
)
from lightsuite.registration.elastix.runner import (
    _discover_transformix_result,
    _patch_transformix_params,
    read_elastix_landmark_metric_mm,
    volume_shape_from_transform_params,
)
from lightsuite.registration.elastix.params import build_bspline_params, write_parameter_file
from lightsuite.registration.elastix.points import (
    volume_indices_to_elastix_physical,
    write_landmark_file,
)
from lightsuite.registration.points_utils import subsample_point_pairs, thin_point_list
from lightsuite.registration.warp import warp_volume_affine


def test_read_mhd_volume_met_short(tmp_path: Path) -> None:
    ny, nx, nz = 4, 5, 6
    data = np.arange(ny * nx * nz, dtype=np.int16).reshape(ny, nx, nz)
    flat = np.transpose(data, (1, 0, 2)).ravel(order="C")
    raw_path = tmp_path / "result.raw"
    flat.tofile(raw_path)
    mhd_path = tmp_path / "result.mhd"
    mhd_path.write_text(
        "\n".join(
            [
                "ObjectType = Image",
                "NDims = 3",
                "BinaryData = True",
                f"DimSize = {nx} {ny} {nz}",
                "ElementType = MET_SHORT",
                "ElementDataFile = result.raw",
            ]
        ),
        encoding="utf-8",
    )
    loaded = read_mhd_volume(mhd_path)
    assert loaded.shape == (ny, nx, nz)
    assert np.allclose(loaded, data.astype(np.float32))


def test_write_and_read_mhd_roundtrip(tmp_path: Path) -> None:
    vol = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    base = tmp_path / "testvol"
    mhd = write_mhd(vol, base, 0.02)
    spacing = read_mhd_spacing(mhd)
    assert np.allclose(spacing, [0.02, 0.02, 0.02])
    assert mhd.with_suffix(".raw").is_file()


def test_write_and_read_mhd_uint16_roundtrip(tmp_path: Path) -> None:
    vol = np.arange(24, dtype=np.uint16).reshape(2, 3, 4) * 100
    base = tmp_path / "testvol_u16"
    mhd = write_mhd(vol, base, 0.02)
    text = mhd.read_text(encoding="utf-8")
    assert "MET_USHORT" in text
    loaded = read_mhd_volume(mhd)
    assert loaded.shape == vol.shape
    assert np.allclose(loaded, vol.astype(np.float32))


def test_scale_volume_for_elastix_mi() -> None:
    sample = np.zeros((10, 10, 10), dtype=np.float32)
    sample[2:8, 2:8, 2:8] = 15000.0
    atlas = np.zeros((10, 10, 10), dtype=np.float32)
    atlas[2:8, 2:8, 2:8] = 500.0
    sample_u16 = scale_volume_for_elastix_mi(sample)
    atlas_u16 = scale_volume_for_elastix_mi(atlas)
    assert sample_u16.dtype == np.uint16
    assert atlas_u16.dtype == np.uint16
    assert int(sample_u16.max()) > 60000
    assert int(atlas_u16.max()) > 60000


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
    assert params["SP_a"] == 500


def test_build_bspline_params_auto_landmarks_only() -> None:
    params = build_bspline_params(
        dual_channel=False,
        control_point_weight=0.2,
        n_histogram_bins=48,
        bspline_spatial_scale_mm=0.64,
        fixed_shape=(20, 20, 20),
        spacing_mm=0.02,
        auto_landmarks_only=True,
    )
    assert params["Metric"] == [
        "AdvancedMattesMutualInformation",
        "CorrespondingPointsEuclideanDistanceMetric",
    ]
    assert params["Optimizer"] == "AdaptiveStochasticGradientDescent"
    assert params["Metric0Weight"] == 0.05
    assert params["Metric1Weight"] == 1.0
    assert params["NumberOfResolutions"] == 1
    assert params["MaximumNumberOfIterations"] == [2000]
    assert params["FinalGridSpacingInPhysicalUnits"] == [1.28, 1.28, 1.28]
    assert params["ImagePyramidSchedule"] == [1, 1, 1]


def test_build_bspline_params_auto_landmarks_identity_only() -> None:
    params = build_bspline_params(
        dual_channel=False,
        control_point_weight=0.2,
        n_histogram_bins=48,
        bspline_spatial_scale_mm=0.64,
        fixed_shape=(20, 20, 20),
        spacing_mm=0.02,
        auto_landmarks_only=True,
        identity_only=True,
    )
    assert params["MaximumNumberOfIterations"] == [0]
    assert params["Metric0Weight"] == 1.0
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
    assert "UseDirectionCosines" in text


def test_volume_indices_to_elastix_physical() -> None:
    # 0-based volume indices (Y, X, Z)
    pts = np.array([[10.0, 20.0, 30.0], [11.0, 22.0, 35.0]])
    phys = volume_indices_to_elastix_physical(pts, 0.02)
    assert np.allclose(phys[0], [0.40, 0.20, 0.60])
    assert np.allclose(phys[1], [0.44, 0.22, 0.70])


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


def test_subsample_point_pairs() -> None:
    atlas = np.arange(300, dtype=float).reshape(100, 3)
    sample = atlas + 1.0
    sub_a, sub_s = subsample_point_pairs(atlas, sample, max_points=20)
    assert sub_a.shape == (20, 3)
    assert np.allclose(sub_s, sub_a + 1.0)


def test_read_elastix_landmark_metric_mm(tmp_path: Path) -> None:
    info = tmp_path / "IterationInfo.0.R0.txt"
    info.write_text(
        "\n".join(
            [
                "1:ItNr\t2:Metric\t2:Metric0\t2:Metric1",
                "0\t0.200\t-0.100\t0.180",
                "1\t0.150\t-0.110\t0.122",
            ]
        ),
        encoding="utf-8",
    )
    assert read_elastix_landmark_metric_mm(tmp_path) == pytest.approx(0.122)


def test_read_elastix_landmark_metric_mm_single_metric(tmp_path: Path) -> None:
    info = tmp_path / "IterationInfo.0.R0.txt"
    info.write_text(
        "\n".join(
            [
                "1:ItNr\t2:Metric\t2:Metric0",
                "0\t0.167402\t0.167402",
                "1\t0.143975\t0.143975",
            ]
        ),
        encoding="utf-8",
    )
    assert read_elastix_landmark_metric_mm(tmp_path) == pytest.approx(0.143975)
    patched = _patch_transformix_params("(UseDirectionCosines true)\n", nearest=True)
    assert '(UseDirectionCosines "false")' in patched
    assert '(DefaultPixelValue "0")' in patched


def test_volume_shape_from_transform_params() -> None:
    text = "(Size 341 671 671)\n"
    assert volume_shape_from_transform_params(text) == (671, 341, 671)


def test_patch_transformix_params_forces_mhd_output() -> None:
    params = "(FinalBSplineInterpolationOrder 3)\n(ResampleInterpolator \"FinalBSplineInterpolator\")\n"
    patched = _patch_transformix_params(params, nearest=True)
    assert '(WriteResultImage "true")' in patched
    assert '(ResultImageFormat "mhd")' in patched
    assert '(FinalBSplineInterpolationOrder "0")' in patched
    assert '(ResultImagePixelType "double")' in patched
    assert "FinalNearestNeighborInterpolator" not in patched


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
