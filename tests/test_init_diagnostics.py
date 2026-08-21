"""Tests for init-registration diagnostics."""

from __future__ import annotations

from io import StringIO

import numpy as np
from rich.console import Console

from lightsuite.registration.align import coarse_alignment_metrics
from lightsuite.registration.init_diagnostics import (
    InitRegistrationDiagnostics,
    classify_init_registration_status,
)


def test_init_print_summary_shows_status_panel() -> None:
    diag = InitRegistrationDiagnostics(
        sample_shape=[10, 10, 10],
        atlas_shape=[8, 8, 8],
        orientation=[1, 2, 3],
        registration_resolution_um=20.0,
        cloud_threshold=5.0,
        sample_cloud_subsample=0.1,
        sample_cloud_points=2_307_391,
        atlas_cloud_points=10_523_389,
        alignment_backend="bcpd",
        similarity_scale=1.001,
        median_error_sample_to_atlas_vox=5.2,
        median_error_atlas_to_sample_vox=1.9,
        median_error_vox=5.2,
        p95_error_sample_to_atlas_vox=43.6,
        p95_error_atlas_to_sample_vox=5.5,
        inlier_fraction=0.89,
        inlier_threshold_vox=25.0,
        auto_pairs=13_598,
        warped_boundary_voxels=7_147_517,
        alignment_elapsed_s=32.4,
        triage_elapsed_s=39.1,
        preview_elapsed_s=9.1,
        status="good",
        status_message="Coarse alignment ready for match-points or register.",
        sample_mask_points=25_828_872,
        sample_trim_points=25_531_779,
        sample_downsample_points=2_553_178,
        sample_denoise_points=2_307_391,
    )
    buf = StringIO()
    diag.print_summary(console=Console(file=buf, force_terminal=False, width=120, color_system=None))
    text = buf.getvalue()
    assert "GOOD" in text
    assert "Coarse alignment ready" in text
    assert "Point clouds" in text
    assert "BCPD" in text

    diag.warnings = ["Sparse sample cloud — try lowering registration.cloud_threshold."]
    buf = StringIO()
    diag.print_summary(console=Console(file=buf, force_terminal=False, width=120, color_system=None))
    text = buf.getvalue()
    assert "Warning:" in text
    assert "Sparse sample cloud" in text


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
