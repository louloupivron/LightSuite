"""Tests for MATLAB point-cloud downsampling ports."""

from __future__ import annotations

import numpy as np

from lightsuite.registration.pc_downsample import (
    downsample_for_bcpd_similarity,
    downsample_for_triage,
    matlab_bcpd_grid_step,
    matlab_triage_grid_step,
    nonuniform_grid,
    nonuniform_grid_sample,
)


def test_matlab_bcpd_grid_step_uses_round() -> None:
    assert matlab_bcpd_grid_step(362_005, 10_000) == 36
    assert matlab_bcpd_grid_step(5_000, 10_000) == 6
    assert matlab_bcpd_grid_step(523_993, 10_000) == 52


def test_matlab_triage_grid_step_uses_ceil() -> None:
    assert matlab_triage_grid_step(362_005, 10_000) == 37
    assert matlab_triage_grid_step(5_000, 10_000) == 1


def test_nonuniform_grid_sample_matches_matlab_leaf_count() -> None:
    rng = np.random.default_rng(0)
    points = rng.random((523_993, 3)) * 200.0
    out = nonuniform_grid_sample(points, 52)
    assert out.shape[0] == 16_384


def test_nonuniform_grid_returns_bin_centroids() -> None:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [10.0, 10.0, 10.0],
        ],
        dtype=float,
    )
    out = nonuniform_grid(points, 2.0)
    assert out.shape[0] == 2
    assert np.allclose(out[0], [0.5, 0.0, 0.0])


def test_downsample_for_bcpd_similarity_reduces_large_cloud() -> None:
    rng = np.random.default_rng(1)
    points = rng.random((20_000, 3)) * 200.0
    out = downsample_for_bcpd_similarity(points, 10_000)
    assert out.shape[0] < points.shape[0]
    assert out.shape[0] > 100


def test_downsample_for_triage_keeps_small_cloud() -> None:
    points = np.arange(30, dtype=float).reshape(10, 3)
    out = downsample_for_triage(points, 10_000)
    assert out.shape[0] > 0
