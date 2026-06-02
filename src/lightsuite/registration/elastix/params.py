"""Elastix parameter file generation (performMultObjBsplineRegistration.m)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def build_bspline_params(
    *,
    dual_channel: bool,
    control_point_weight: float,
    n_histogram_bins: int,
    bspline_spatial_scale_mm: float,
    fixed_shape: tuple[int, int, int],
    spacing_mm: float,
    use_multistep: bool = True,
    dual_weight_autofluor: float = 1.0,
    dual_weight_signal: float = 0.5,
) -> dict[str, Any]:
    """Build elastix parameter dict for multi-metric B-spline registration.

    Mirrors ``performMultObjBsplineRegistration.m`` exactly: the no-manual-landmark
    (auto-only) case uses the same multi-resolution schedule and image-dominant metric
    weights as the manual case. The caller simply halves the control-point weight for
    auto-only runs (MATLAB "revert to automated mode").
    """
    cpwt = control_point_weight
    ny, nx, nz = fixed_shape
    fixed_size_mm = np.array([ny, nx, nz], dtype=float) * spacing_mm
    max_sample_region = float(min(fixed_size_mm) / 3.0)
    base_sizes = [4.5, 4.0, 3.0, 2.0] if use_multistep else [2.0]
    base_sizes = [min(v, max_sample_region) for v in base_sizes]

    params: dict[str, Any] = {
        "Registration": "MultiMetricMultiResolutionRegistration",
        "Transform": "RecursiveBSplineTransform",
        "Optimizer": "AdaptiveStochasticGradientDescent",
        "AutomaticParameterEstimation": "true",
        "AutomaticScalesEstimation": "false",
        "BSplineInterpolationOrder": 3,
        "FinalBSplineInterpolationOrder": 3,
        "FixedImageDimension": 3,
        "MovingImageDimension": 3,
        "UseRandomSampleRegion": "true",
        "NewSamplesEveryIteration": "true",
        "NumberOfResolutions": 4,
        "NumberOfHistogramBins": n_histogram_bins,
        "SP_A": 20,
        "SP_a": 500,
        "SP_alpha": 0.6,
        "MaximumNumberOfIterations": [500, 1000, 1500, 2000],
        "NumberOfSpatialSamples": 5000,
        "ImagePyramidSchedule": [8] * 3 + [4] * 3 + [2] * 3 + [1] * 3,
        "FinalGridSpacingInPhysicalUnits": [bspline_spatial_scale_mm] * 3,
        "FixedLimitRangeRatio": 0.01,
        "MovingLimitRangeRatio": 0.01,
        "FixedImageBSplineInterpolationOrder": 1,
        # Elastix 5.1 defaults to true; our MHD headers use identity axes only (see elastix.log warning).
        "UseDirectionCosines": "false",
    }
    if dual_channel:
        params["Metric"] = [
            "AdvancedMattesMutualInformation",
            "AdvancedMattesMutualInformation",
            "CorrespondingPointsEuclideanDistanceMetric",
        ]
        params["FixedImagePyramid"] = ["FixedRecursiveImagePyramid"] * 3
        params["MovingImagePyramid"] = ["MovingRecursiveImagePyramid"] * 3
        params["ImageSampler"] = ["RandomCoordinate"] * 3
        params["Interpolator"] = ["BSplineInterpolator"] * 3
        params["Metric0Weight"] = dual_weight_autofluor
        params["Metric1Weight"] = dual_weight_signal
        params["Metric2Weight"] = cpwt
    else:
        params["Metric"] = [
            "AdvancedMattesMutualInformation",
            "CorrespondingPointsEuclideanDistanceMetric",
        ]
        params["FixedImagePyramid"] = "FixedRecursiveImagePyramid"
        params["MovingImagePyramid"] = "MovingRecursiveImagePyramid"
        params["ImageSampler"] = "RandomCoordinate"
        params["Interpolator"] = "BSplineInterpolator"
        params["Metric0Weight"] = 1.0
        params["Metric1Weight"] = cpwt

    if use_multistep:
        sample_regions: list[float] = []
        for size in base_sizes:
            sample_regions.extend([size] * 3)
        params["SampleRegionSize"] = sample_regions
    else:
        params["SampleRegionSize"] = [base_sizes[0]] * 3

    return params


def write_parameter_file(path: Path, params: dict[str, Any]) -> None:
    """Write elastix parameter file from a parameter dictionary."""
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for key, value in params.items():
        lines.append(_format_param(key, value))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _format_param(key: str, value: Any) -> str:
    if isinstance(value, str):
        return f'({key} "{value}")'
    if isinstance(value, bool):
        return f"({key} {'true' if value else 'false'})"
    if isinstance(value, (list, tuple)):
        if value and isinstance(value[0], str):
            quoted = " ".join(f'"{v}"' for v in value)
            return f"({key} {quoted})"
        nums = " ".join(_num(v) for v in value)
        return f"({key} {nums})"
    return f"({key} {_num(value)})"


def _num(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(int(value))
