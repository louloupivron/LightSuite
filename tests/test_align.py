"""Tests for coarse similarity alignment."""

from __future__ import annotations

from unittest.mock import patch

import numpy as np

from lightsuite.registration.align import (
    _icp_cloud_subset,
    _transform_points,
    estimate_similarity_transform,
    similarity_scale,
    triage_and_match_clouds,
)
from lightsuite.registration.bcpd import atlas_to_sample_affinetform
from lightsuite.registration.pc_downsample import downsample_for_bcpd_similarity
from lightsuite.registration.warp import matlab_voxel_affine_from_icp


def test_icp_cloud_subset_uses_full_small_sample_cloud() -> None:
    points = np.arange(120, dtype=float).reshape(40, 3)
    subset = _icp_cloud_subset(points, 10_000, cap=12_000)
    assert subset.shape[0] == 40


def test_downsample_for_bcpd_similarity_keeps_tiny_clouds() -> None:
    points = np.arange(15, dtype=float).reshape(5, 3)
    subset = downsample_for_bcpd_similarity(points, 10_000)
    assert subset.shape == points.shape


def test_atlas_to_sample_affinetform_fixes_translation_vs_raw_inv() -> None:
    """Regression for MATLAB premultiply inversion (not raw np.linalg.inv)."""
    hybrid = np.array(
        [
            [0.9242, 0.1880, 0.0625, -50.0859],
            [-0.1829, 0.9243, -0.0757, -40.2302],
            [-0.0761, 0.0620, 0.9401, -132.4498],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    sample_to_atlas = atlas_to_sample_affinetform(hybrid)
    wrong_affine = np.linalg.inv(hybrid)

    assert sample_to_atlas[1, 3] > 0.0
    assert abs(sample_to_atlas[1, 3] - 4.5) < abs(wrong_affine[1, 3] - 4.5)


@patch("lightsuite.registration.align.resolve_bcpd_executable", return_value=None)
def test_triage_uses_icp_transform_frame(_mock_bcpd: object) -> None:
    rng = np.random.default_rng(0)
    sample = rng.random((200, 3)) * 40.0
    atlas = sample + np.array([12.0, 8.0, 5.0])
    transform_icp = np.eye(4)
    transform_icp[:3, :3] *= 1.08
    transform_icp[:3, 3] = [12.0, 8.0, 5.0]
    transform_matlab = matlab_voxel_affine_from_icp(transform_icp)

    icp_pairs = triage_and_match_clouds(sample, atlas, transform_icp)
    assert icp_pairs[0].shape[0] > 0

    icp_err = np.linalg.norm(_transform_points(sample, transform_icp) - atlas, axis=1).mean()
    matlab_err = np.linalg.norm(_transform_points(sample, transform_matlab) - atlas, axis=1).mean()
    assert icp_err < matlab_err


def test_triage_finds_pairs_with_huge_atlas_cloud() -> None:
    rng = np.random.default_rng(2)
    sample = rng.random((4_000, 3)) * 35.0
    transform = np.eye(4)
    transform[:3, :3] *= 1.05
    transform[:3, 3] = [120.0, 80.0, 200.0]
    aligned = _transform_points(sample, transform)
    atlas_far = rng.random((600_000, 3)) * np.array([660.0, 400.0, 570.0])
    atlas_near = aligned + rng.normal(0.0, 2.5, aligned.shape)
    atlas = np.vstack([atlas_far, atlas_near])
    pairs = triage_and_match_clouds(sample, atlas, transform)
    assert pairs[0].shape[0] > 100


@patch("lightsuite.registration.align.resolve_bcpd_executable", return_value=None)
def test_estimate_similarity_transform_returns_both_frames(_mock_bcpd: object) -> None:
    rng = np.random.default_rng(1)
    sample = rng.random((500, 3)) * 30.0
    atlas = sample * 1.05 + np.array([4.0, 2.0, 1.0])
    transform_icp, transform_matlab, backend = estimate_similarity_transform(atlas, sample)
    assert backend == "icp"
    assert transform_icp.shape == (4, 4)
    assert transform_matlab.shape == (4, 4)
    assert not np.allclose(transform_icp, transform_matlab)
    assert 0.75 <= similarity_scale(transform_icp) <= 1.35

    pairs = triage_and_match_clouds(sample, atlas, transform_icp)
    assert pairs[0].shape[0] > 0
