"""Tests for coarse similarity alignment."""

from __future__ import annotations

import numpy as np

from lightsuite.registration.align import (
    _icp_cloud_subset,
    _matlab_cloud_subset,
    _transform_points,
    estimate_similarity_transform,
    similarity_scale,
    triage_and_match_clouds,
)
from lightsuite.registration.warp import matlab_voxel_affine_from_icp


def test_icp_cloud_subset_uses_full_small_sample_cloud() -> None:
    points = np.arange(120, dtype=float).reshape(40, 3)
    subset = _icp_cloud_subset(points, 10_000, cap=12_000)
    assert subset.shape[0] == 40


def test_matlab_cloud_subset_keeps_small_clouds() -> None:
    points = np.arange(30, dtype=float).reshape(10, 3)
    subset = _matlab_cloud_subset(points, 10_000)
    assert subset.shape == points.shape


def test_triage_uses_icp_transform_frame() -> None:
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


def test_estimate_similarity_transform_returns_both_frames() -> None:
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
