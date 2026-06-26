"""Tests for spinal cord B-spline elastix parameter generation."""

from __future__ import annotations

from lightsuite.registration.elastix.cord_bspline import build_cord_bspline_params


def test_build_cord_bspline_params_without_control_points_uses_mi_only() -> None:
    params = build_cord_bspline_params(
        control_point_weight=0.2,
        fixed_shape=(126, 161, 1052),
        use_control_points=False,
    )
    assert params["Registration"] == "MultiResolutionRegistration"
    assert params["Metric"] == "AdvancedMattesMutualInformation"
    assert "Metric1Weight" not in params


def test_build_cord_bspline_params_with_control_points_uses_dual_metric() -> None:
    params = build_cord_bspline_params(
        control_point_weight=0.2,
        fixed_shape=(126, 161, 1052),
        use_control_points=True,
    )
    assert params["Registration"] == "MultiMetricMultiResolutionRegistration"
    assert params["Metric"] == [
        "AdvancedMattesMutualInformation",
        "CorrespondingPointsEuclideanDistanceMetric",
    ]
    assert params["Metric1Weight"] == 0.2
