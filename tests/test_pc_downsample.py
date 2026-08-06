"""Tests for MATLAB point-cloud downsampling ports."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from lightsuite.registration.pc_downsample import (
    downsample_for_bcpd_similarity,
    downsample_for_triage,
    matlab_bcpd_grid_step,
    matlab_triage_grid_step,
    nonuniform_grid,
    nonuniform_grid_sample,
    pcdenoise,
    pcdownsample_random,
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


def test_nonuniform_grid_returns_leaf_centroids() -> None:
    rng = np.random.default_rng(0)
    points = rng.random((10_691_773 // 100, 3)) * 200.0
    out = nonuniform_grid(points, 214 // 100 + 6)
    assert out.shape[0] == 2 ** int(np.ceil(np.log2(points.shape[0] / 8)))
    assert np.allclose(out.mean(axis=0), points.mean(axis=0), atol=1.0)


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


def test_pcdownsample_random_targets_rounded_fraction() -> None:
    rng = np.random.default_rng(1)
    points = rng.random((50_000, 3)) * 200.0
    out = pcdownsample_random(points, 0.1, preserve_structure=True, seed=1)
    assert out.shape[0] == 5_000


def test_pcdenoise_keeps_dense_cluster() -> None:
    rng = np.random.default_rng(0)
    core = rng.normal(size=(5_000, 3)) * 0.05
    noise = rng.uniform(-5, 5, size=(200, 3))
    points = np.vstack([core, noise])
    out = pcdenoise(points, num_neighbors=4, std_ratio=1.0)
    assert out.shape[0] > 4_500
    assert out.shape[0] < points.shape[0]


def test_pcdenoise_matches_matlab_marianna_export() -> None:
    export = Path(
        "/run/user/1002/gvfs/smb-share:server=h1data1.wysscenter.ch,share=computingdata"
        "/Alice/0003_CBT_WYSS_LIGHTSHEET/DATA/MESOSPIM/Marianna CMU/voluse_export"
    )
    pts_path = export / "downsample_pts.txt"
    outlier_path = export / "outlier_idx.txt"
    if not pts_path.is_file() or not outlier_path.is_file():
        return

    down = np.loadtxt(pts_path, delimiter=",", dtype=np.float64)
    raw = outlier_path.read_text().strip().replace("\n", ",")
    outlier_idx = np.array([int(x) for x in raw.split(",") if x], dtype=np.int64) - 1
    mat_outliers = np.zeros(down.shape[0], dtype=bool)
    mat_outliers[outlier_idx] = True

    kept = pcdenoise(down)
    assert kept.shape[0] == 523_993

    tree = cKDTree(down)
    dists, _ = tree.query(down, k=5, workers=-1)
    mean_sq = np.square(dists[:, 1:]).mean(axis=1)
    cutoff = float(mean_sq.mean() + mean_sq.std(ddof=1))
    py_outliers = mean_sq > cutoff
    assert np.array_equal(py_outliers, mat_outliers)
