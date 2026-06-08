"""Tests for init-registration diagnostics."""

from __future__ import annotations

import numpy as np

from lightsuite.registration.align import coarse_alignment_metrics
from lightsuite.registration.init_diagnostics import classify_init_registration_status


def test_coarse_alignment_metrics_bidirectional() -> None:
    rng = np.random.default_rng(0)
    sample = rng.random((200, 3)) * 40.0
    transform = np.eye(4)
    transform[:3, :3] *= 1.05
    transform[:3, 3] = [5.0, 3.0, 2.0]
    aligned = sample @ transform[:3, :3].T + transform[:3, 3]
    atlas = aligned + rng.normal(0.0, 2.0, aligned.shape)

    metrics = coarse_alignment_metrics(sample, atlas, transform)
    assert metrics["median_sample_to_atlas_vox"] < 5.0
    assert metrics["median_atlas_to_sample_vox"] < 5.0
    assert metrics["median_vox"] == max(
        metrics["median_sample_to_atlas_vox"],
        metrics["median_atlas_to_sample_vox"],
    )
    assert metrics["inlier_fraction"] > 0.8


def test_classify_init_registration_status_good() -> None:
    status, message, warnings = classify_init_registration_status(
        median_error_vox=8.0,
        auto_pairs=120,
        inlier_fraction=0.9,
        similarity_scale=1.02,
        warped_boundary_voxels=10_000,
        sample_cloud_points=5_000,
        atlas_cloud_points=50_000,
        alignment_backend="bcpd",
    )
    assert status == "good"
    assert "ready" in message.lower()
    assert warnings == []


def test_classify_init_registration_status_poor() -> None:
    status, _, _ = classify_init_registration_status(
        median_error_vox=30.0,
        auto_pairs=80,
        inlier_fraction=0.2,
        similarity_scale=1.0,
        warped_boundary_voxels=0,
        sample_cloud_points=5_000,
        atlas_cloud_points=50_000,
        alignment_backend="icp",
    )
    assert status == "poor"


def test_classify_init_registration_status_failed() -> None:
    status, message, _ = classify_init_registration_status(
        median_error_vox=5.0,
        auto_pairs=2,
        inlier_fraction=0.9,
        similarity_scale=1.0,
        warped_boundary_voxels=1_000,
        sample_cloud_points=5_000,
        atlas_cloud_points=50_000,
        alignment_backend="bcpd",
    )
    assert status == "failed"
    assert "too few" in message.lower()
