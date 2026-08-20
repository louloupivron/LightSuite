"""Tests for register-step diagnostics."""

from __future__ import annotations

from io import StringIO

from rich.console import Console

from lightsuite.registration.register_diagnostics import (
    RegistrationDiagnostics,
    classify_registration_status,
    landmark_mm_to_vox,
)


def test_register_print_summary_is_status_only() -> None:
    diag = RegistrationDiagnostics(
        sample_shape=[10, 10, 10],
        atlas_shape=[8, 8, 8],
        orientation=[1, 2, 3],
        registration_resolution_um=20.0,
        n_manual_pairs=0,
        n_auto_pairs=32,
        n_landmark_pairs=32,
        control_point_weight=0.1,
        use_multistep=True,
        use_dual_channel_mi=False,
        bspline_spatial_scale_mm=0.64,
        status="poor",
        status_message="Registration quality is low — review previews and add manual landmarks.",
        warnings=["Auto-only landmarks — consider match-points for difficult samples."],
    )
    buf = StringIO()
    diag.print_summary(console=Console(file=buf, force_terminal=False, width=80, color_system=None))
    text = buf.getvalue()
    assert "POOR" in text
    assert "Registration quality is low" in text
    assert "Warning:" in text
    assert "Control points" not in text
    assert "Affine fit" not in text


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
