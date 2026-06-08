"""Tests for register-step diagnostics."""

from __future__ import annotations

from lightsuite.registration.register_diagnostics import (
    classify_registration_status,
    landmark_mm_to_vox,
)


def test_landmark_mm_to_vox_at_20um() -> None:
    assert landmark_mm_to_vox(1.0, 20.0) == 50.0


def test_classify_registration_status_good() -> None:
    status, message, warnings = classify_registration_status(
        affine_median_error_vox=6.0,
        affine_p95_error_vox=12.0,
        n_manual=12,
        n_landmark_pairs=80,
        bspline_landmark_metric_vox=5.0,
        annotation_label_voxels=1_000_000,
        use_multistep=True,
    )
    assert status == "good"
    assert "export" in message.lower()
    assert warnings == []


def test_classify_registration_status_poor_affine() -> None:
    status, _, _ = classify_registration_status(
        affine_median_error_vox=30.0,
        affine_p95_error_vox=40.0,
        n_manual=0,
        n_landmark_pairs=50,
        bspline_landmark_metric_vox=8.0,
        annotation_label_voxels=500_000,
        use_multistep=True,
    )
    assert status == "poor"


def test_classify_registration_status_moderate_auto_only() -> None:
    status, _, warnings = classify_registration_status(
        affine_median_error_vox=10.0,
        affine_p95_error_vox=18.0,
        n_manual=0,
        n_landmark_pairs=120,
        bspline_landmark_metric_vox=8.0,
        annotation_label_voxels=500_000,
        use_multistep=True,
    )
    assert status == "moderate"
    assert any("auto-only" in w.lower() for w in warnings)
