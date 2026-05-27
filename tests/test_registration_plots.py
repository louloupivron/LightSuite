"""Tests for registration preview plots."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np

from lightsuite.gui.slices import volume_index_to_image
from lightsuite.registration.plots import (
    annotation_boundary_pixels,
    plot_annotation_comparison,
    save_initial_registration_previews,
)
from lightsuite.registration.warp import warp_atlas_to_sample, warp_volume_affine


def test_annotation_boundary_pixels_detects_edges() -> None:
    annot = np.zeros((10, 10), dtype=np.float32)
    annot[2:8, 2:8] = 5
    row, col = annotation_boundary_pixels(annot)
    assert row.size > 0
    assert col.size == row.size


def test_plot_annotation_comparison_layout() -> None:
    volume = np.random.randint(0, 255, size=(24, 32, 20), dtype=np.uint8)
    annotation = np.zeros((24, 32, 20), dtype=np.float32)
    annotation[4:20, 6:26, 4:16] = 10
    fig = plot_annotation_comparison(volume, annotation, dimplot=3)
    assert len(fig.axes) == 8
    fig.clf()


def test_save_initial_registration_previews(tmp_path) -> None:
    sample = np.random.rand(20, 24, 18).astype(np.float32)
    annotation = np.zeros((20, 24, 18), dtype=np.float32)
    annotation[4:16, 5:19, 3:15] = 55
    transform = np.eye(4)

    save_initial_registration_previews(tmp_path, sample, annotation, transform)

    for idim in range(1, 4):
        png = tmp_path / f"dim{idim}_initial_registration.png"
        assert png.is_file()
        assert png.stat().st_size > 10_000


def test_warp_atlas_to_sample_uses_inverse_transform() -> None:
    atlas = np.zeros((20, 24, 18), dtype=np.float32)
    atlas[4:16, 6:18, 4:14] = 88
    sample_shape = (40, 48, 32)
    sample_to_atlas = np.eye(4) * 0.5
    sample_to_atlas[3, 3] = 1.0

    av_correct = warp_atlas_to_sample(atlas, sample_to_atlas, sample_shape, order=0)
    av_wrong = warp_volume_affine(atlas, sample_to_atlas, sample_shape, order=0)

    assert np.count_nonzero(av_correct > 1) > 500
    assert np.count_nonzero(av_wrong > 1) < np.count_nonzero(av_correct > 1) / 10
