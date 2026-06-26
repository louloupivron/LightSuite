"""Elastix affine registration helpers (performElastixAffineRegistration.m port)."""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lightsuite.registration.elastix.mhd import scale_volume_for_elastix_mi, write_mhd
from lightsuite.registration.elastix.params import write_parameter_file
from lightsuite.registration.elastix.runner import clear_elastix_workspace


@dataclass
class AffineRegistrationResult:
    output_dir: Path
    transform_path: Path
    copied_transform_path: Path


def build_affine_params() -> dict:
    """Parameter dict matching performElastixAffineRegistration.m."""
    return {
        "Registration": "MultiResolutionRegistration",
        "Metric": "AdvancedMattesMutualInformation",
        "Transform": "AffineTransform",
        "Optimizer": "AdaptiveStochasticGradientDescent",
        "ImageSampler": "RandomCoordinate",
        "AutomaticParameterEstimation": "true",
        "AutomaticScalesEstimation": "true",
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
        "MaximumNumberOfIterations": [500, 1000, 1500, 2000],
        "NumberOfSpatialSamples": 5000,
        "ImagePyramidSchedule": [8] * 3 + [4] * 3 + [2] * 3 + [1] * 3,
        "UseDirectionCosines": "false",
        "WriteResultImage": "true",
        "ResultImageFormat": "mhd",
        # Spinal straightvol extends beyond z-scaled atlas at rostral/caudal ends; default
        # RequiredRatioOfValidSamples (0.25) aborts at full resolution (~24% OOB samples).
        "RequiredRatioOfValidSamples": 0.1,
        "MaximumNumberOfSamplingAttempts": 10,
        "DefaultPixelValue": 0,
    }


def run_affine_registration(
    *,
    fixed_volume: np.ndarray,
    moving_volume: np.ndarray,
    save_path: Path,
    spacing_mm: float,
    output_name: str = "affine_atlas_to_samp_20um.txt",
    work_dir: Path | None = None,
    output_path: Path | None = None,
) -> AffineRegistrationResult:
    """Run elastix affine registration (moving → fixed)."""
    if shutil.which("elastix") is None:
        msg = "elastix not found on PATH"
        raise RuntimeError(msg)

    save_path = save_path.expanduser()
    output_dir = (work_dir or save_path / "elastix_temp_affine").expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    clear_elastix_workspace(output_dir)

    param_path = output_dir / "affine_parameters.txt"
    write_parameter_file(param_path, build_affine_params())

    sp = [spacing_mm, spacing_mm, spacing_mm]
    fixed_mhd = output_dir / "fixed"
    moving_mhd = output_dir / "moving"
    write_mhd(scale_volume_for_elastix_mi(fixed_volume), fixed_mhd, sp)
    write_mhd(scale_volume_for_elastix_mi(moving_volume), moving_mhd, sp)

    cmd = [
        "elastix",
        "-f",
        str(fixed_mhd.with_suffix(".mhd")),
        "-m",
        str(moving_mhd.with_suffix(".mhd")),
        "-out",
        str(output_dir),
        "-p",
        str(param_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        msg = f"elastix affine failed (exit {result.returncode}):\n{result.stdout}\n{result.stderr}"
        raise RuntimeError(msg)

    transforms = sorted(output_dir.glob("TransformParameters.*.txt"))
    if not transforms:
        msg = f"No TransformParameters.*.txt in {output_dir}"
        raise RuntimeError(msg)

    copied = (output_path or save_path / output_name).expanduser()
    copied.parent.mkdir(parents=True, exist_ok=True)
    copied.write_text(transforms[0].read_text(encoding="utf-8"), encoding="utf-8")
    return AffineRegistrationResult(
        output_dir=output_dir,
        transform_path=transforms[0],
        copied_transform_path=copied,
    )


def _get_param(text: str, name: str, count: int) -> np.ndarray:
    pat = rf"\(\s*{re.escape(name)}\s+([^)]+)\)"
    match = re.search(pat, text)
    if match is None:
        msg = f"Parameter {name} not found in elastix transform file."
        raise ValueError(msg)
    vals = np.fromstring(match.group(1), sep=" ")
    if vals.size < count:
        msg = f"Expected {count} values for {name}, got {vals.size}"
        raise ValueError(msg)
    return vals[:count]


def parse_elastix_affine(path: Path) -> np.ndarray:
    """Parse elastix affine to MATLAB affinetform3d 4x4 (moving → fixed, 1-based voxels)."""
    raw = path.expanduser().read_text(encoding="utf-8")
    params = _get_param(raw, "TransformParameters", 12)
    center = _get_param(raw, "CenterOfRotationPoint", 3)
    r_flat = params[:9]
    r_mat = r_flat.reshape(3, 3).T
    t_vec = params[9:12]
    t_eff = t_vec + center - (r_mat @ center)

    m_elastix = np.eye(4, dtype=float)
    m_elastix[:3, :3] = r_mat
    m_elastix[:3, 3] = t_eff
    m_matlab_phys = np.linalg.inv(m_elastix)

    shift_to_zero = np.eye(4)
    shift_to_zero[:3, 3] = -1.0
    shift_to_one = np.eye(4)
    shift_to_one[:3, 3] = 1.0
    # affinetform3d / imwarp_volume use translation in the 4th column (column-vector layout).
    return shift_to_one @ m_matlab_phys @ shift_to_zero


def compose_affinetform(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Compose affinetform3d matrices (MATLAB post-multiply: points @ A)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return a @ b
