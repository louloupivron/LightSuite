"""Elastix B-spline registration for spinal cord (performCordBsplineRegistration.m port)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from lightsuite.registration.elastix.params import write_parameter_file


def build_cord_bspline_params(
    *,
    control_point_weight: float,
    fixed_shape: tuple[int, int, int],
    use_control_points: bool = True,
) -> dict[str, Any]:
    """Cord-specific elastix parameters from performCordBsplineRegistration.m."""
    base: dict[str, Any] = {
        "Transform": "RecursiveBSplineTransform",
        "Optimizer": "AdaptiveStochasticGradientDescent",
        "ImageSampler": "RandomCoordinate",
        "AutomaticParameterEstimation": "true",
        "AutomaticScalesEstimation": "false",
        "BSplineInterpolationOrder": 3,
        "FinalBSplineInterpolationOrder": 3,
        "FixedImageDimension": 3,
        "MovingImageDimension": 3,
        "FixedImagePyramid": "FixedRecursiveImagePyramid",
        "MovingImagePyramid": "MovingRecursiveImagePyramid",
        "UseRandomSampleRegion": "true",
        "NewSamplesEveryIteration": "true",
        "NumberOfResolutions": 4,
        "NumberOfHistogramBins": 16,
        "SP_A": 20,
        "MaximumNumberOfIterations": [500, 1000, 1500, 2000],
        "NumberOfSpatialSamples": 5000,
        "ImagePyramidSchedule": [8] * 3 + [4] * 3 + [2] * 3 + [1] * 3,
        "FinalGridSpacingInPhysicalUnits": [0.96] * 3,
        "SampleRegionSize": [2, 2, 10, 2, 2, 8, 2, 2, 6, 2, 2, 4],
        "UseDirectionCosines": "false",
        "WriteResultImage": "true",
        "ResultImageFormat": "mhd",
        "FixedInternalImagePixelType": "float",
        "MovingInternalImagePixelType": "float",
        "RequiredRatioOfValidSamples": 0.1,
        "MaximumNumberOfSamplingAttempts": 10,
        "DefaultPixelValue": 0,
    }
    if use_control_points and control_point_weight > 0:
        return {
            **base,
            "Registration": "MultiMetricMultiResolutionRegistration",
            "Metric": ["AdvancedMattesMutualInformation", "CorrespondingPointsEuclideanDistanceMetric"],
            "Metric1Weight": control_point_weight,
            "Metric0Weight": 1.0,
        }
    return {
        **base,
        "Registration": "MultiResolutionRegistration",
        "Metric": "AdvancedMattesMutualInformation",
    }


def write_cord_bspline_parameter_file(path: Path, *, cpwt: float, fixed_shape: tuple[int, int, int]) -> Path:
    params = build_cord_bspline_params(control_point_weight=cpwt, fixed_shape=fixed_shape)
    write_parameter_file(path, params)
    return path
