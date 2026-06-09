"""Tests for correspondence-informed affine correction in register."""

from __future__ import annotations

import numpy as np

from lightsuite.gui.affine import transform_points, transform_points_inverse
from lightsuite.gui.slice_correspondence import SliceAnchor, SliceCorrespondence
from lightsuite.registration.correspondence_affine import (
    build_correspondence_anchor_pairs,
    compose_correspondence_affine_correction,
)


def test_build_correspondence_anchor_pairs_counts_confirmed_only() -> None:
    volume = np.zeros((40, 50, 30), dtype=np.float32)
    volume[10:30, 15:35, 8:22] = 200.0
    correspondence = SliceCorrespondence(
        original_trans=np.eye(4).tolist(),
        axes={
            1: [
                SliceAnchor(10, 11, True),
                SliceAnchor(30, 31, False),
            ],
            2: [
                SliceAnchor(12, 22, True),
                SliceAnchor(32, 42, True),
            ],
        },
    )
    atlas_pts, sample_pts = build_correspondence_anchor_pairs(
        correspondence,
        sample_warped=volume,
        original_trans_vol=np.eye(4),
        downfac=1.0,
    )
    assert atlas_pts.shape == (3, 3)
    assert sample_pts.shape == (3, 3)


def test_compose_correspondence_affine_reduces_anchor_residuals() -> None:
    rng = np.random.default_rng(0)
    corr_atlas = rng.uniform(10, 30, size=(8, 3))
    bias = np.array([5.0, -2.0, 1.5])
    tform_aff = np.eye(4)
    corr_sample = corr_atlas + bias

    before = np.linalg.norm(
        transform_points(corr_atlas, tform_aff) - corr_sample,
        axis=1,
    )
    tform_corrected, stats = compose_correspondence_affine_correction(
        tform_aff,
        corr_atlas,
        corr_sample,
    )
    after = np.linalg.norm(
        transform_points(corr_atlas, tform_corrected) - corr_sample,
        axis=1,
    )
    assert stats.applied
    assert float(np.median(after)) < float(np.median(before))
    assert float(np.max(after)) < 1e-6


def test_compose_skips_when_fewer_than_four_anchors() -> None:
    tform = np.eye(4)
    corr_atlas = np.array([[1.0, 2.0, 3.0]])
    corr_sample = np.array([[1.0, 2.0, 3.0]])
    tform_out, stats = compose_correspondence_affine_correction(tform, corr_atlas, corr_sample)
    assert np.allclose(tform_out, tform)
    assert not stats.applied
    assert stats.skip_reason == "fewer_than_4_anchor_pairs"


def test_compose_skips_when_already_aligned() -> None:
    atlas_pts = np.array(
        [
            [10.0, 20.0, 8.0],
            [12.0, 22.0, 10.0],
            [14.0, 24.0, 12.0],
            [16.0, 26.0, 14.0],
        ]
    )
    tform_aff = np.eye(4)
    tform_out, stats = compose_correspondence_affine_correction(
        tform_aff,
        atlas_pts,
        atlas_pts.copy(),
    )
    assert not stats.applied
    assert stats.skip_reason == "already_aligned"
    assert np.allclose(tform_out, tform_aff)
