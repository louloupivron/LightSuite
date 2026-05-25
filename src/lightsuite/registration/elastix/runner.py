"""Run elastix / transformix via subprocess."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from rich.console import Console

from lightsuite.registration.elastix.mhd import write_mhd
from lightsuite.registration.elastix.params import build_bspline_params, write_parameter_file
from lightsuite.registration.elastix.points import write_landmark_file

console = Console()


@dataclass
class BsplineRegistrationResult:
    output_dir: Path
    transform_path: Path
    copied_transform_path: Path | None = None


def clear_elastix_workspace(directory: Path) -> None:
    """Remove stale elastix outputs (clearElastixWorkspaceForNewRun.m)."""
    directory = directory.expanduser()
    if not directory.is_dir():
        return
    for pattern in ("TransformParameters.*.txt", "elastix.log"):
        for path in directory.glob(pattern):
            if path.is_file():
                path.unlink()


def run_bspline_registration(
    *,
    fixed_volume,
    moving_volume,
    fixed_secondary,
    moving_points_mm,
    fixed_points_mm,
    output_dir: Path,
    save_path: Path,
    spacing_mm: float,
    control_point_weight: float,
    n_histogram_bins: int,
    bspline_spatial_scale_mm: float,
    use_multistep: bool,
    dual_weight_autofluor: float,
    dual_weight_signal: float,
) -> BsplineRegistrationResult:
    """Run single- or dual-channel elastix B-spline registration."""
    if shutil.which("elastix") is None:
        msg = "elastix not found on PATH"
        raise RuntimeError(msg)

    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    clear_elastix_workspace(output_dir)

    fixpath = output_dir / "fixed.txt"
    movpath = output_dir / "moving.txt"
    write_landmark_file(fixpath, fixed_points_mm)
    write_landmark_file(movpath, moving_points_mm)

    dual = fixed_secondary is not None
    params = build_bspline_params(
        dual_channel=dual,
        control_point_weight=control_point_weight,
        n_histogram_bins=n_histogram_bins,
        bspline_spatial_scale_mm=bspline_spatial_scale_mm,
        fixed_shape=fixed_volume.shape,
        spacing_mm=spacing_mm,
        use_multistep=use_multistep,
        dual_weight_autofluor=dual_weight_autofluor,
        dual_weight_signal=dual_weight_signal,
    )
    param_path = output_dir / "bspline_parameters.txt"
    write_parameter_file(param_path, params)

    sp = [spacing_mm, spacing_mm, spacing_mm]
    if dual:
        stem = output_dir.name
        base_f0 = output_dir / f"{stem}_dual_f0"
        base_m0 = output_dir / f"{stem}_dual_m0"
        base_f1 = output_dir / f"{stem}_dual_f1"
        base_m1 = output_dir / f"{stem}_dual_m1"
        base_f2 = output_dir / f"{stem}_dual_f2"
        base_m2 = output_dir / f"{stem}_dual_m2"
        write_mhd(fixed_volume, base_f0, sp)
        write_mhd(moving_volume, base_m0, sp)
        write_mhd(fixed_secondary, base_f1, sp)
        write_mhd(moving_volume, base_m1, sp)
        write_mhd(fixed_volume, base_f2, sp)
        write_mhd(moving_volume, base_m2, sp)
        cmd = [
            "elastix",
            "-f0",
            str(base_f0.with_suffix(".mhd")),
            "-m0",
            str(base_m0.with_suffix(".mhd")),
            "-f1",
            str(base_f1.with_suffix(".mhd")),
            "-m1",
            str(base_m1.with_suffix(".mhd")),
            "-f2",
            str(base_f2.with_suffix(".mhd")),
            "-m2",
            str(base_m2.with_suffix(".mhd")),
            "-out",
            str(output_dir),
            "-fp",
            str(fixpath),
            "-mp",
            str(movpath),
            "-p",
            str(param_path),
        ]
    else:
        fixed_mhd = output_dir / "fixed.mhd"
        moving_mhd = output_dir / "moving.mhd"
        write_mhd(fixed_volume, fixed_mhd.with_suffix(""), sp)
        write_mhd(moving_volume, moving_mhd.with_suffix(""), sp)
        cmd = [
            "elastix",
            "-f",
            str(fixed_mhd),
            "-m",
            str(moving_mhd),
            "-out",
            str(output_dir),
            "-fp",
            str(fixpath),
            "-mp",
            str(movpath),
            "-p",
            str(param_path),
        ]

    (output_dir / "CMD.txt").write_text(" ".join(cmd), encoding="utf-8")
    console.print("Running elastix B-spline registration...")
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        msg = f"elastix failed (exit {result.returncode}):\n{result.stdout}\n{result.stderr}"
        raise RuntimeError(msg)

    transforms = sorted(output_dir.glob("TransformParameters.*.txt"))
    if not transforms:
        msg = f"No TransformParameters.*.txt found in {output_dir}"
        raise RuntimeError(msg)

    forward_copy = save_path / "bspline_atlas_to_samp_20um.txt"
    forward_copy.write_text(transforms[0].read_text(encoding="utf-8"), encoding="utf-8")
    return BsplineRegistrationResult(
        output_dir=output_dir,
        transform_path=transforms[0],
        copied_transform_path=forward_copy,
    )


def run_transformix(
    *,
    moving_volume,
    transform_path: Path,
    output_dir: Path,
    spacing_mm: float,
    nearest: bool = False,
):
    """Apply a transformix transform and load the result volume."""
    if shutil.which("transformix") is None:
        msg = "transformix not found on PATH"
        raise RuntimeError(msg)

    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("result.*", "TransformParameters.*.txt"):
        for path in output_dir.glob(pattern):
            if path.is_file():
                path.unlink()

    params = transform_path.read_text(encoding="utf-8")
    if nearest:
        params = params.replace(
            "FinalBSplineInterpolationOrder = 3",
            "FinalBSplineInterpolationOrder = 0",
        )
        params = params.replace(
            "(FinalBSplineInterpolationOrder 3)",
            "(FinalBSplineInterpolationOrder 0)",
        )
        if "ResultImagePixelType" not in params:
            params += '\n(ResultImagePixelType "float")\n'
    temp_param = output_dir / "transformix_params.txt"
    temp_param.write_text(params, encoding="utf-8")

    moving_mhd = output_dir / "moving.mhd"
    write_mhd(moving_volume, moving_mhd.with_suffix(""), [spacing_mm] * 3)
    cmd = [
        "transformix",
        "-in",
        str(moving_mhd),
        "-out",
        str(output_dir),
        "-tp",
        str(temp_param),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        msg = f"transformix failed (exit {proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
        raise RuntimeError(msg)

    result_mhd = next(output_dir.glob("result.*.mhd"), None)
    if result_mhd is None:
        msg = f"No transformix result MHD in {output_dir}"
        raise RuntimeError(msg)
    return _read_mhd_volume(result_mhd)


def _read_mhd_volume(mhd_path: Path):
    import numpy as np

    text = mhd_path.read_text(encoding="utf-8")
    dim_match = None
    for line in text.splitlines():
        if line.lower().startswith("dimsize"):
            dim_match = [int(v) for v in line.split("=", 1)[1].split()]
            break
    if dim_match is None:
        msg = f"DimSize missing in {mhd_path}"
        raise ValueError(msg)
    nx, ny, nz = dim_match
    raw_name = None
    for line in text.splitlines():
        if line.lower().startswith("elementdatafile"):
            raw_name = line.split("=", 1)[1].strip()
            break
    if raw_name is None:
        msg = f"ElementDataFile missing in {mhd_path}"
        raise ValueError(msg)
    raw_path = mhd_path.parent / raw_name
    flat = np.fromfile(raw_path, dtype=np.float32)
    vol = flat.reshape((nx, ny, nz), order="C")
    return np.transpose(vol, (1, 0, 2))
