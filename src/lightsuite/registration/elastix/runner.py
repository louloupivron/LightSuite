"""Run elastix / transformix via subprocess."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from rich.console import Console

from lightsuite.registration.elastix.mhd import read_mhd_volume, write_mhd
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


def _patch_transformix_params(params: str, *, nearest: bool) -> str:
    """Match transformAnnotationVolume.m: label resampling + always write a result image."""
    if nearest:
        for old, new in (
            ("FinalBSplineInterpolationOrder = 3", "FinalBSplineInterpolationOrder = 0"),
            ("(FinalBSplineInterpolationOrder 3)", "(FinalBSplineInterpolationOrder 0)"),
        ):
            params = params.replace(old, new)
    overrides = [
        '(WriteResultImage "true")',
        '(ResultImageFormat "mhd")',
        '(ResultImagePixelType "float")',
    ]
    if nearest:
        overrides.append('(FinalBSplineInterpolationOrder 0)')
    return params.rstrip() + "\n" + "\n".join(overrides) + "\n"


def _discover_transformix_result(output_dir: Path) -> Path | None:
    """Find transformix output volume (Elastix 4.x writes MHD; 5.x may use NIfTI)."""
    patterns = (
        "result*.mhd",
        "Result*.mhd",
        "result.nii.gz",
        "result*.nii.gz",
        "result.nii",
        "result*.nii",
    )
    for pattern in patterns:
        hits = sorted(output_dir.glob(pattern))
        if hits:
            return hits[0]
    for pattern in patterns:
        hits = sorted(output_dir.rglob(pattern))
        if hits:
            return hits[0]
    return None


def _read_transformix_result(result_path: Path):
    """Load transformix output into Y, X, Z numpy order."""
    import numpy as np

    name = result_path.name.lower()
    if name.endswith(".nii.gz") or name.endswith(".nii"):
        import nibabel as nib

        data = np.asanyarray(nib.load(str(result_path)).dataobj).astype(np.float32)
        if data.ndim != 3:
            msg = f"Expected 3D transformix result, got {data.shape} from {result_path}"
            raise ValueError(msg)
        return np.transpose(data, (1, 0, 2))
    return read_mhd_volume(result_path)


def _transformix_failure_message(output_dir: Path, proc: subprocess.CompletedProcess) -> str:
    listing = ", ".join(sorted(p.name for p in output_dir.iterdir())) or "(empty)"
    stderr_tail = (proc.stderr or "").strip()[-4000:]
    stdout_tail = (proc.stdout or "").strip()[-2000:]
    return (
        f"No transformix result image in {output_dir} (expected result.mhd or result.nii.gz). "
        f"Directory listing: {listing}. "
        f"transformix stderr (tail): {stderr_tail or '(empty)'} "
        f"stdout (tail): {stdout_tail or '(empty)'}"
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
    for pattern in ("result*", "Result*", "TransformParameters.*.txt", "transformix_params.txt"):
        for path in output_dir.glob(pattern):
            if path.is_file():
                path.unlink()

    params = _patch_transformix_params(transform_path.read_text(encoding="utf-8"), nearest=nearest)
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
    (output_dir / "transformix_cmd.txt").write_text(" ".join(cmd), encoding="utf-8")
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (output_dir / "transformix.stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
    (output_dir / "transformix.stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
    if proc.returncode != 0:
        msg = f"transformix failed (exit {proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
        raise RuntimeError(msg)

    result_path = _discover_transformix_result(output_dir)
    if result_path is None:
        raise RuntimeError(_transformix_failure_message(output_dir, proc))
    return _read_transformix_result(result_path)
