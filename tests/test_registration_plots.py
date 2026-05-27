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
)
from lightsuite.registration.warp import imwarp_volume, matlab_voxel_affine_from_icp, warp_atlas_to_sample


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


def test_plot_annotation_comparison_layout() -> None:
    volume = np.random.randint(0, 255, size=(24, 32, 20), dtype=np.uint8)
    boundary = np.zeros((24, 32, 20), dtype=np.uint8)
    boundary[4:20, 6:26, 4:16] = 255
    fig = plot_annotation_comparison(volume, boundary, dimplot=3)
    assert len(fig.axes) == 8
    fig.clf()


def test_save_initial_registration_previews(tmp_path) -> None:
    sample = np.random.rand(20, 24, 18).astype(np.float32)
    boundary = np.zeros((20, 24, 18), dtype=np.uint8)
    boundary[4:16, 5:19, 3:15] = 255
    transform = np.eye(4)

    warped_count = save_initial_registration_previews(tmp_path, sample, boundary, transform)

    assert warped_count > 0
    for idim in range(1, 4):
        png = tmp_path / f"dim{idim}_initial_registration.png"
        assert png.is_file()
        assert png.stat().st_size > 10_000


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


def test_preview_warp_produces_visible_boundaries() -> None:
    sample = np.full((40, 48, 32), 120, dtype=np.uint8)
    boundary = np.zeros((20, 24, 18), dtype=np.uint8)
    boundary[4:16, 6:18, 4:14] = 255
    sample_to_atlas = matlab_voxel_affine_from_icp(np.eye(4) * 0.5)

    boundary_warped = warp_atlas_to_sample(boundary, sample_to_atlas, sample.shape, order=0)
    mid = sample.shape[1] // 2
    sl = volume_index_to_image(boundary_warped, np.array([mid, 2], dtype=int))
    assert mask_boundary_pixels(sl)[0].size > 20
