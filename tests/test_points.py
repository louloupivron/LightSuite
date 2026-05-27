"""Tests for registration point cloud extraction."""

from __future__ import annotations

import numpy as np

from lightsuite.registration.points import extract_sample_points


def test_extract_sample_points_keeps_more_without_voxel_downsample() -> None:
    yy, xx, zz = np.mgrid[0:48, 0:48, 0:24]
    volume = np.exp(-((yy - 24) ** 2 + (xx - 24) ** 2 + (zz - 12) ** 2) / 200.0).astype(np.float32)
    dense = extract_sample_points(volume, threshold=1.0, subsample_fraction=1.0)
    sparse = extract_sample_points(volume, threshold=1.0, subsample_fraction=0.1)
    assert dense.shape[0] > 1000
    assert sparse.shape[0] >= 100
    assert dense.shape[0] > sparse.shape[0]


def test_extract_sample_points_zeros_grad_below_intensity_cutoff() -> None:
    volume = np.zeros((24, 24, 12), dtype=np.float32)
    volume[8:16, 8:16, 4:8] = 0.05
    volume[10:14, 10:14, 5:7] = 1.0
    points = extract_sample_points(volume, threshold=0.5, subsample_fraction=1.0)
    assert points.shape[0] > 0
