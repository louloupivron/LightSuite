"""Tests for mesoSPIM geometry helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import SimpleITK as sitk

from lightsuite.mesospim.config_models import MesospimGeometryConfig
from lightsuite.mesospim.geometry import (
    apply_image_geometry,
    overlap_physical_bounds,
    physical_bounds,
    prepare_registration_pair,
)


def _meta(*, x_pos: float, y_pos: float, z_start: float, z_end: float) -> dict[str, float]:
    return {
        "Pixelsize in um": 5.0,
        "x_pos": x_pos,
        "y_pos": y_pos,
        "z_start": z_start,
        "z_end": z_end,
        "z_stepsize": 2.0,
    }


def _volume(shape_zyx: tuple[int, int, int]) -> sitk.Image:
    arr = np.zeros(shape_zyx, dtype=np.float32)
    return sitk.GetImageFromArray(arr)


def test_overlap_and_prepare_pair() -> None:
    geometry = MesospimGeometryConfig(
        stage_xy_is_center=True,
        itk_lateral_dim0_motor="x",
        lateral_flip=[1, 1],
    )
    meta_a = _meta(x_pos=0.0, y_pos=0.0, z_start=0.0, z_end=10.0)
    meta_b = _meta(x_pos=20.0, y_pos=0.0, z_start=0.0, z_end=10.0)

    fixed = _volume((5, 20, 20))
    moving = _volume((5, 20, 20))
    apply_image_geometry(fixed, meta_a, geometry)
    apply_image_geometry(moving, meta_b, geometry)

    overlap_min, overlap_max = overlap_physical_bounds(fixed, moving, margin_um=0.0)
    assert overlap_min[0] < overlap_max[0]
    assert overlap_min[1] <= overlap_max[1]

    fixed_crop, moving_resampled, box, crop_start = prepare_registration_pair(
        fixed, moving, margin_um=0.0
    )
    assert fixed_crop.GetSize() == moving_resampled.GetSize()
    assert box[0].shape == (3,)
    assert box[1].shape == (3,)
    assert len(crop_start) == 3

    pmin, pmax = physical_bounds(fixed_crop)
    assert np.all(pmin <= pmax)


def test_embed_crop_in_full_overview() -> None:
    from lightsuite.mesospim.geometry import embed_crop_in_full_overview

    full_arr = np.zeros((10, 20, 20), dtype=np.float32)
    crop_arr = np.ones((4, 6, 8), dtype=np.float32)
    full = sitk.GetImageFromArray(full_arr)
    crop = sitk.GetImageFromArray(crop_arr)
    crop_start = [5, 7, 3]

    embedded = embed_crop_in_full_overview(full, crop, crop_start)
    out = sitk.GetArrayFromImage(embedded)
    assert out.shape == full_arr.shape
    assert np.count_nonzero(out) == crop_arr.size
    assert out[3:7, 7:13, 5:13].sum() == crop_arr.sum()


def test_save_geometry_overlap_qc_plot_metadata_only(tmp_path: Path) -> None:
    from lightsuite.mesospim.plots import save_geometry_overlap_qc_plot

    geometry = MesospimGeometryConfig(
        stage_xy_is_center=True,
        itk_lateral_dim0_motor="x",
        lateral_flip=[1, 1],
    )
    meta_a = _meta(x_pos=0.0, y_pos=0.0, z_start=0.0, z_end=10.0)
    meta_b = _meta(x_pos=20.0, y_pos=0.0, z_start=0.0, z_end=10.0)

    overview = _volume((5, 20, 20))
    roi = _volume((5, 20, 20))
    overview_arr = sitk.GetArrayFromImage(overview)
    overview_arr[2, 10, 10] = 100.0
    overview = sitk.GetImageFromArray(overview_arr)
    roi_arr = sitk.GetArrayFromImage(roi)
    roi_arr[2, 10, 10] = 100.0
    roi = sitk.GetImageFromArray(roi_arr)
    apply_image_geometry(overview, meta_a, geometry)
    apply_image_geometry(roi, meta_b, geometry)

    overlap_min, overlap_max = overlap_physical_bounds(overview, roi, margin_um=0.0)
    output_path = tmp_path / "geometry_overlap_qc.png"
    save_geometry_overlap_qc_plot(
        overview=overview,
        roi=roi,
        metadata_overlap_min=overlap_min,
        metadata_overlap_max=overlap_max,
        output_path=output_path,
        geometry_mode="metadata",
    )
    assert output_path.is_file()
    assert output_path.stat().st_size > 1000


def test_save_geometry_overlap_qc_plot_hybrid_comparison(tmp_path: Path) -> None:
    from lightsuite.mesospim.landmark_geometry import fit_landmark_transform, overlap_box_from_landmark_transform
    from lightsuite.mesospim.landmark_session import MesospimLandmarkSession
    from lightsuite.mesospim.plots import save_geometry_overlap_qc_plot

    overview = _volume((20, 20, 20))
    roi = _volume((10, 10, 10))
    session = MesospimLandmarkSession(
        overview_points_zyx=[[5, 5, 5], [5, 14, 14], [5, 5, 14]],
        roi_points_zyx=[[0, 0, 0], [0, 9, 9], [0, 0, 9]],
        fit_mode="rigid",
    )
    fit = fit_landmark_transform(
        overview=overview,
        roi=roi,
        session=session,
        fit_mode="rigid",
        min_pairs=3,
    )
    metadata_overlap = overlap_physical_bounds(overview, roi, margin_um=0.0)
    hybrid_overlap = overlap_box_from_landmark_transform(
        overview,
        roi,
        fit.roi_to_overview_tform,
        margin_um=0.0,
    )
    output_path = tmp_path / "geometry_overlap_qc.png"
    save_geometry_overlap_qc_plot(
        overview=overview,
        roi=roi,
        metadata_overlap_min=metadata_overlap[0],
        metadata_overlap_max=metadata_overlap[1],
        hybrid_overlap_min=hybrid_overlap[0],
        hybrid_overlap_max=hybrid_overlap[1],
        roi_to_overview=fit.roi_to_overview_tform,
        output_path=output_path,
        geometry_mode="hybrid",
    )
    assert output_path.is_file()
    assert output_path.stat().st_size > 1000
