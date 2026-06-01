"""Tests for registration preview plots."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np

from lightsuite.gui.slices import volume_index_to_image
from lightsuite.registration.plots import (
    annotation_boundary_pixels,
    mask_boundary_pixels,
    plot_annotation_comparison,
    save_initial_registration_previews,
    save_registration_stage_previews,
)
from lightsuite.gui.affine import fit_affine_transform, transform_points
from lightsuite.registration.points import cloud_xyz_to_volume_indices
from lightsuite.registration.volume import resize_atlas_volume
from lightsuite.registration.warp import (
    affinetform_rows_to_internal,
    imwarp_volume,
    matlab_voxel_affine_from_icp,
    matlab_voxel_affine_to_icp,
    transform_points_affinetform,
    warp_atlas_to_sample,
    warp_volume_affine,
)


def test_mask_boundary_pixels_detects_mask() -> None:
    mask = np.zeros((10, 10), dtype=np.uint8)
    mask[3, 3:7] = 255
    row, col = mask_boundary_pixels(mask)
    assert row.size == 4
    assert col.size == 4


def test_annotation_boundary_pixels_detects_edges() -> None:
    annot = np.zeros((10, 10), dtype=np.float32)
    annot[2:8, 2:8] = 5
    row, col = annotation_boundary_pixels(annot)
    assert row.size > 0
    assert col.size == row.size


def test_annotation_boundary_pixels_ignores_float_interpolation_speckle() -> None:
    annot = np.zeros((10, 10), dtype=np.float32)
    annot[2:8, 2:8] = 5.0
    annot[5, 5] = 5.3  # speckle inside region
    row, col = annotation_boundary_pixels(annot)
    assert row.size < 80


def test_plot_annotation_comparison_layout() -> None:
    volume = np.random.randint(0, 255, size=(24, 32, 20), dtype=np.uint8)
    boundary = np.zeros((24, 32, 20), dtype=np.uint8)
    boundary[4:20, 6:26, 4:16] = 255
    fig = plot_annotation_comparison(volume, boundary, dimplot=3)
    assert len(fig.axes) == 8
    fig.clf()


def test_save_registration_stage_previews(tmp_path) -> None:
    sample_u8 = np.random.randint(0, 255, size=(24, 32, 20), dtype=np.uint8)
    annotation = np.zeros((24, 32, 20), dtype=np.float32)
    annotation[4:20, 6:26, 4:16] = 5

    save_registration_stage_previews(
        tmp_path, "mouse1", sample_u8, annotation, "affine_registration"
    )

    for idim in range(1, 4):
        png = tmp_path / f"mouse1_dim{idim}_affine_registration.png"
        assert png.is_file()
        assert png.stat().st_size > 10_000


def test_save_initial_registration_previews(tmp_path) -> None:
    sample = np.random.rand(20, 24, 18).astype(np.float32)
    annotation = np.zeros((20, 24, 18), dtype=np.float32)
    annotation[4:16, 5:19, 3:15] = 5
    transform = np.eye(4)

    warped_count = save_initial_registration_previews(tmp_path, sample, annotation, transform)

    assert warped_count > 0
    for idim in range(1, 4):
        png = tmp_path / f"dim{idim}_initial_registration.png"
        assert png.is_file()
        assert png.stat().st_size > 10_000


def test_affinetform_rows_to_internal_matches_transform_points() -> None:
    affinetform = np.eye(4, dtype=float)
    affinetform[:3, :3] = np.array(
        [
            [0.95, -0.12, 0.03],
            [0.11, 0.97, 0.04],
            [-0.02, -0.03, 1.01],
        ]
    )
    affinetform[:3, 3] = [6.0, -4.0, 2.0]
    points = np.array([[0.0, 0.0, 0.0], [12.0, 5.0, 3.0], [40.0, 20.0, 10.0]])
    internal = affinetform_rows_to_internal(affinetform)
    from lightsuite.registration.align import _transform_points

    assert np.allclose(
        _transform_points(points, internal),
        transform_points_affinetform(points, affinetform),
        atol=1e-6,
    )


def test_warp_atlas_to_sample_with_rotation() -> None:
    boundary = np.zeros((20, 24, 18), dtype=np.uint8)
    boundary[4:16, 6:18, 4:14] = 255
    sample_shape = (40, 48, 32)

    affinetform = np.eye(4, dtype=float)
    affinetform[:3, :3] = np.array(
        [
            [0.98, -0.17, 0.0],
            [0.17, 0.98, 0.0],
            [0.0, 0.0, 1.0],
        ]
    ) * 0.9
    affinetform[:3, 3] = [4.0, -2.0, 3.0]
    stored = matlab_voxel_affine_from_icp(affinetform_rows_to_internal(affinetform))

    warped = warp_atlas_to_sample(boundary, stored, sample_shape, order=0)
    assert np.count_nonzero(warped) > 100


def test_matlab_voxel_affine_roundtrip_with_rotation() -> None:
    internal = np.eye(4, dtype=float)
    internal[:3, :3] = np.array(
        [
            [0.95, -0.12, 0.03],
            [0.11, 0.97, 0.04],
            [-0.02, -0.03, 1.01],
        ]
    ).T
    internal[:3, 3] = [6.0, -4.0, 2.0]
    stored = matlab_voxel_affine_from_icp(internal)
    recovered = matlab_voxel_affine_to_icp(stored)
    assert np.allclose(recovered[:3, :3], internal[:3, :3], atol=1e-6)
    assert np.allclose(recovered[:3, 3], internal[:3, 3], atol=1e-6)


def test_warp_atlas_to_sample_matches_matlab_imwarp() -> None:
    boundary = np.zeros((20, 24, 18), dtype=np.uint8)
    boundary[4:16, 6:18, 4:14] = 255
    sample_shape = (40, 48, 32)
    sample_to_atlas = matlab_voxel_affine_from_icp(np.eye(4) * 0.5)

    av_correct = warp_atlas_to_sample(boundary, sample_to_atlas, sample_shape, order=0)
    av_wrong = imwarp_volume(
        boundary,
        sample_to_atlas,
        sample_shape,
        order=0,
        point_coords="xyz",
    )

    assert np.count_nonzero(av_correct) > 500
    assert np.count_nonzero(av_wrong) < np.count_nonzero(av_correct) / 10


def test_affine_fit_and_warp_share_volume_axis_order() -> None:
    """Auto landmarks from clouds are XYZ; imwarp uses (Y, X, Z) array indices."""
    downfac = 0.5
    ann_full = np.zeros((40, 48, 36), dtype=np.float32)
    ann_full[8:32, 12:36, 6:30] = 5
    sample_shape = (80, 96, 64)

    atlas_xyz = np.array(
        [
            [12.0, 24.0, 18.0],
            [14.0, 26.0, 20.0],
            [16.0, 28.0, 22.0],
            [18.0, 30.0, 24.0],
        ]
    )
    sample_xyz = atlas_xyz + np.array([5.0, -3.0, 2.0])
    atlas_vol = cloud_xyz_to_volume_indices(atlas_xyz / downfac)
    sample_vol = cloud_xyz_to_volume_indices(sample_xyz)
    tform_wrong, _ = fit_affine_transform(atlas_xyz / downfac, sample_xyz)
    tform_right, _ = fit_affine_transform(atlas_vol, sample_vol)

    pred_wrong = transform_points(atlas_vol, tform_wrong)
    pred_right = transform_points(atlas_vol, tform_right)
    err_wrong = np.linalg.norm(pred_wrong - sample_vol, axis=1)
    err_right = np.linalg.norm(pred_right - sample_vol, axis=1)
    assert np.median(err_right) < np.median(err_wrong) / 5


def test_affine_volume_warp_must_use_full_atlas_not_downsampled() -> None:
    """Regression: tform_aff is fit in full-atlas coords (see multiobjRegistration.m)."""
    downfac = 0.5
    sample_shape = (80, 96, 64)
    ann_full = np.zeros((40, 48, 36), dtype=np.float32)
    ann_full[8:32, 12:36, 6:30] = 5
    ann_reg = resize_atlas_volume(ann_full, downfac, nearest=True)

    atlas_pts = np.array(
        [
            [20.0, 24.0, 18.0],
            [30.0, 30.0, 20.0],
            [15.0, 20.0, 12.0],
            [25.0, 28.0, 22.0],
        ],
        dtype=float,
    )
    sample_pts = np.array(
        [
            [40.0, 50.0, 30.0],
            [55.0, 60.0, 35.0],
            [35.0, 45.0, 25.0],
            [48.0, 58.0, 32.0],
        ],
        dtype=float,
    )
    af_atlas = atlas_pts / downfac
    tform_aff, _ = fit_affine_transform(af_atlas, sample_pts)

    warped_full = warp_volume_affine(ann_full, tform_aff, sample_shape, order=0)
    warped_reg = warp_volume_affine(ann_reg, tform_aff, sample_shape, order=0)

    mid = sample_shape[1] // 2
    sl_full = volume_index_to_image(warped_full, np.array([mid, 2], dtype=int))
    sl_reg = volume_index_to_image(warped_reg, np.array([mid, 2], dtype=int))
    n_full = mask_boundary_pixels(sl_full)[0].size
    n_reg = mask_boundary_pixels(sl_reg)[0].size
    assert n_full > 30
    assert n_reg < n_full / 5

    pred = transform_points(af_atlas, tform_aff)
    assert np.median(np.linalg.norm(pred - sample_pts, axis=1)) < 2.0


def test_preview_warp_produces_visible_boundaries() -> None:
    sample = np.full((40, 48, 32), 120, dtype=np.uint8)
    boundary = np.zeros((20, 24, 18), dtype=np.uint8)
    boundary[4:16, 6:18, 4:14] = 255
    sample_to_atlas = matlab_voxel_affine_from_icp(np.eye(4) * 0.5)

    boundary_warped = warp_atlas_to_sample(boundary, sample_to_atlas, sample.shape, order=0)
    mid = sample.shape[1] // 2
    sl = volume_index_to_image(boundary_warped, np.array([mid, 2], dtype=int))
    assert mask_boundary_pixels(sl)[0].size > 20
