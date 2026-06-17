"""Tests for mesoSPIM landmark-based geometry."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk
import yaml

from lightsuite.config.loader import load_mesospim_config
from lightsuite.mesospim.config_models import MesospimGeometryMode
from lightsuite.mesospim.landmark_geometry import (
    fit_landmark_transform,
    overlap_box_from_landmark_transform,
    prepare_registration_pair_from_landmarks,
)
from lightsuite.mesospim.landmark_session import (
    MesospimLandmarkSession,
    default_landmark_session_path,
)
from lightsuite.mesospim.prepare import prepare_mesospim_registration_pair


def _image(shape_zyx: tuple[int, int, int], spacing: float = 1.0) -> sitk.Image:
    arr = np.zeros(shape_zyx, dtype=np.float32)
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((spacing, spacing, spacing))
    image.SetOrigin((0.0, 0.0, 0.0))
    return image


def _translation_session(offset: int = 5) -> MesospimLandmarkSession:
    return MesospimLandmarkSession(
        overview_points_zyx=[
            [offset, offset, offset],
            [offset + 9, offset + 9, offset + 9],
            [offset, offset + 9, offset],
        ],
        roi_points_zyx=[
            [0, 0, 0],
            [9, 9, 9],
            [0, 9, 0],
        ],
        fit_mode="rigid",
    )


def test_fit_landmark_translation() -> None:
    overview = _image((20, 20, 20))
    roi = _image((10, 10, 10))
    session = _translation_session(offset=5)

    fit = fit_landmark_transform(
        overview=overview,
        roi=roi,
        session=session,
        fit_mode="rigid",
        min_pairs=3,
    )
    assert fit.rms_error_um < 1e-6
    assert fit.roi_to_overview_tform[0, 3] == pytest.approx(5.0)
    assert fit.roi_to_overview_tform[1, 3] == pytest.approx(5.0)
    assert fit.roi_to_overview_tform[2, 3] == pytest.approx(5.0)


def test_prepare_registration_pair_from_landmarks() -> None:
    overview = _image((20, 20, 20))
    roi = _image((10, 10, 10))
    session = _translation_session(offset=5)

    fixed_cropped, moving, overlap_box, crop_start, fit = prepare_registration_pair_from_landmarks(
        overview,
        roi,
        session=session,
        fit_mode="rigid",
        min_pairs=3,
        margin_um=0.0,
    )
    assert fit.rms_error_um < 1e-6
    assert fixed_cropped.GetSize() == moving.GetSize()
    assert crop_start[0] == 5
    assert crop_start[1] == 5
    assert crop_start[2] in (4, 5)
    overlap_min, overlap_max = overlap_box
    assert overlap_min[0] == 5.0
    assert overlap_max[0] == 14.0


def test_overlap_box_rejects_out_of_bounds_placement() -> None:
    overview = _image((10, 10, 10))
    roi = _image((5, 5, 5))
    session = MesospimLandmarkSession(
        overview_points_zyx=[[9, 9, 9], [9, 9, 4], [9, 4, 9]],
        roi_points_zyx=[[0, 0, 0], [0, 0, 5], [0, 5, 0]],
        fit_mode="rigid",
    )
    fit = fit_landmark_transform(
        overview=overview,
        roi=roi,
        session=session,
        fit_mode="rigid",
        min_pairs=3,
    )
    try:
        overlap_box_from_landmark_transform(
            overview,
            roi,
            fit.roi_to_overview_tform,
            margin_um=0.0,
        )
    except ValueError as exc:
        assert "outside the overview bounds" in str(exc)
    else:
        raise AssertionError("expected ValueError for out-of-bounds landmark placement")


def test_landmarks_config_requires_voxel_um(tmp_path: Path) -> None:
    import tifffile

    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    tifffile.imwrite(overview, np.zeros((5, 8, 8), dtype=np.uint16), imagej=True)
    tifffile.imwrite(roi, np.zeros((5, 8, 8), dtype=np.uint16), imagej=True)

    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "t", "save_path": str(tmp_path / "out")},
                "mesospim": {
                    "geometry_mode": "landmarks",
                    "overview": {"path": str(overview)},
                    "roi": {"path": str(roi)},
                },
            }
        ),
        encoding="utf-8",
    )
    try:
        load_mesospim_config(config_path)
    except ValueError as exc:
        assert "voxel_um" in str(exc)
    else:
        raise AssertionError("expected validation error for missing voxel_um")


def test_landmark_roi_report_matches_overlap_box() -> None:
    """ROI FOV in QA reports must use the ROI volume, not the overview."""
    from lightsuite.mesospim.runner import _landmark_roi_report

    overview = _image((20, 20, 20))
    roi = _image((10, 10, 10))
    session = _translation_session(offset=5)
    fit = fit_landmark_transform(
        overview=overview,
        roi=roi,
        session=session,
        fit_mode="rigid",
        min_pairs=3,
    )

    rep_roi = _landmark_roi_report(roi, fit)
    overlap_min, overlap_max = overlap_box_from_landmark_transform(
        overview,
        roi,
        fit.roi_to_overview_tform,
        margin_um=0.0,
    )

    assert np.allclose(rep_roi["phys_min"], overlap_min)
    assert np.allclose(rep_roi["phys_max"], overlap_max)


def test_prepare_landmarks_pipeline(tmp_path: Path) -> None:
    import tifffile

    overview = tmp_path / "overview.tif"
    roi = tmp_path / "roi.tif"
    tifffile.imwrite(overview, np.zeros((20, 20, 20), dtype=np.uint16), imagej=True)
    tifffile.imwrite(roi, np.zeros((10, 10, 10), dtype=np.uint16), imagej=True)

    save_path = tmp_path / "out"
    session = _translation_session(offset=5)
    session.save(default_landmark_session_path(save_path))

    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample": {"name": "t", "save_path": str(save_path)},
                "mesospim": {
                    "geometry_mode": "landmarks",
                    "overview": {"path": str(overview), "voxel_um": [1.0, 1.0, 1.0]},
                    "roi": {"path": str(roi), "voxel_um": [1.0, 1.0, 1.0]},
                    "landmarks": {"fit_mode": "rigid", "min_pairs": 3},
                },
            }
        ),
        encoding="utf-8",
    )
    cfg = load_mesospim_config(config_path)
    assert cfg.mesospim.geometry_mode == MesospimGeometryMode.LANDMARKS

    prepared = prepare_mesospim_registration_pair(cfg)
    assert prepared.crop_start_index[0] == 5
    assert prepared.crop_start_index[1] == 5
    assert prepared.crop_start_index[2] in (4, 5)
    assert prepared.landmark_fit is not None
    assert prepared.landmark_fit.rms_error_um < 1e-6
