"""Run elastix / transformix via subprocess."""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from lightsuite.registration.elastix.mhd import (
    _mhd_element_dtype,
    _mhd_header_field,
    read_mhd_volume,
    scale_volume_for_elastix_mi,
    write_mhd,
)
from lightsuite.registration.elastix.params import build_bspline_params, write_parameter_file
from lightsuite.registration.elastix.points import (
    volume_indices_to_elastix_physical,
    write_landmark_file,
)


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
    bending_energy_weight: float = 0.0,
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
        bending_energy_weight=bending_energy_weight,
    )
    param_path = output_dir / "bspline_parameters.txt"
    write_parameter_file(param_path, params)

    fixed_u16 = scale_volume_for_elastix_mi(fixed_volume)
    moving_u16 = scale_volume_for_elastix_mi(moving_volume)
    fixed_secondary_u16 = (
        scale_volume_for_elastix_mi(fixed_secondary) if dual else None
    )

    sp = [spacing_mm, spacing_mm, spacing_mm]
    if dual:
        # One image pair per metric: AF+atlas, signal+atlas, then AF+atlas duplicates for
        # the landmark and (optional) bending-energy slots, which need valid inputs but
        # do not read the intensities.
        stem = output_dir.name
        n_metrics = len(params["Metric"])
        fixed_for_metric = [fixed_u16, fixed_secondary_u16] + [fixed_u16] * (n_metrics - 2)
        cmd = ["elastix"]
        for index in range(n_metrics):
            base_f = output_dir / f"{stem}_dual_f{index}"
            base_m = output_dir / f"{stem}_dual_m{index}"
            write_mhd(fixed_for_metric[index], base_f, sp)
            write_mhd(moving_u16, base_m, sp)
            cmd += [
                f"-f{index}",
                str(base_f.with_suffix(".mhd")),
                f"-m{index}",
                str(base_m.with_suffix(".mhd")),
            ]
        cmd += [
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
        write_mhd(fixed_u16, fixed_mhd.with_suffix(""), sp)
        write_mhd(moving_u16, moving_mhd.with_suffix(""), sp)
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


def read_elastix_landmark_metric_mm(output_dir: Path) -> float | None:
    """Return the final landmark metric (mm) from elastix IterationInfo, if present."""
    output_dir = output_dir.expanduser()
    candidates = sorted(output_dir.glob("IterationInfo.*.R*.txt"))
    if not candidates:
        return None
    text = candidates[-1].read_text(encoding="utf-8", errors="replace")
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        return None
    header = lines[0]
    metric_cols = [
        idx
        for idx, name in enumerate(header.split())
        if name.endswith(":Metric1") or name == "2:Metric1"
    ]
    if not metric_cols:
        metric_cols = [
            idx
            for idx, name in enumerate(header.split())
            if name.endswith(":Metric0") or name == "2:Metric0"
        ]
    if not metric_cols:
        return None
    col = metric_cols[0]
    data_lines = [line for line in lines[1:] if not line.startswith("1:")]
    if not data_lines:
        return None
    fields = data_lines[-1].split()
    if len(fields) <= col:
        return None
    try:
        return float(fields[col])
    except ValueError:
        return None


def _upsert_elastix_param(params: str, key: str, value: str) -> str:
    """Set one elastix parameter, removing any prior lines for the same key."""
    import re

    pat_new = re.compile(rf"^\s*\(\s*{re.escape(key)}\s+", re.IGNORECASE)
    pat_old = re.compile(rf"^\s*{re.escape(key)}\s*=", re.IGNORECASE)
    lines = [
        line
        for line in params.splitlines()
        if line.strip() and not pat_new.match(line) and not pat_old.match(line)
    ]
    lines.append(f'({key} "{value}")')
    return "\n".join(lines) + "\n"


def _patch_transformix_params(params: str, *, nearest: bool) -> str:
    """Match transformAnnotationVolume.m: label resampling + always write a result image."""
    params = _upsert_elastix_param(params, "WriteResultImage", "true")
    params = _upsert_elastix_param(params, "ResultImageFormat", "mhd")
    params = _upsert_elastix_param(params, "UseDirectionCosines", "false")
    params = _upsert_elastix_param(params, "DefaultPixelValue", "0")
    if nearest:
        params = _upsert_elastix_param(params, "ResultImagePixelType", "double")
        params = _upsert_elastix_param(params, "FinalBSplineInterpolationOrder", "0")
    else:
        params = _upsert_elastix_param(params, "ResultImagePixelType", "float")
    return params


def volume_shape_from_transform_params(params_text: str) -> tuple[int, int, int] | None:
    """Parse elastix (Size x y z) into numpy volume shape (Y, X, Z)."""
    import re

    match = re.search(r"(?im)^\(Size\s+(\d+)\s+(\d+)\s+(\d+)", params_text)
    if not match:
        return None
    nx, ny, nz = (int(v) for v in match.groups())
    return ny, nx, nz


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


def run_transformix_deformation_field(
    *,
    transform_path: Path,
    output_dir: Path,
) -> np.ndarray:
    """Load B-spline deformation field via ``transformix -def all`` (transformPointsToAtlas.m).

    Returns an array with shape (Y, X, Z, 3) in physical mm displacement components
    matching elastix / MetaImage vector layout.
    """
    if shutil.which("transformix") is None:
        msg = "transformix not found on PATH"
        raise RuntimeError(msg)

    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("deformationField*", "transformix_params.txt"):
        for path in output_dir.glob(pattern):
            if path.is_file():
                path.unlink()

    raw_params = transform_path.read_text(encoding="utf-8")
    params = _upsert_elastix_param(raw_params, "WriteDeformationField", "true")
    params = _upsert_elastix_param(params, "WriteResultImage", "false")
    temp_param = output_dir / "transformix_def_params.txt"
    temp_param.write_text(params, encoding="utf-8")

    cmd = [
        "transformix",
        "-def",
        "all",
        "-out",
        str(output_dir),
        "-tp",
        str(temp_param),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        msg = f"transformix -def all failed (exit {proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
        raise RuntimeError(msg)

    field_path = output_dir / "deformationField.mhd"
    if not field_path.is_file():
        hits = sorted(output_dir.glob("deformationField*.mhd"))
        if not hits:
            listing = ", ".join(sorted(p.name for p in output_dir.iterdir())) or "(empty)"
            msg = f"No deformationField.mhd in {output_dir} (listing: {listing})"
            raise FileNotFoundError(msg)
        field_path = hits[0]

    return _read_mhd_vector_field(field_path)


def _read_mhd_vector_field(mhd_path: Path) -> np.ndarray:
    """Read elastix vector deformation field into (Y, X, Z, 3) float32."""
    text = mhd_path.expanduser().read_text(encoding="utf-8")
    dim_field = _mhd_header_field(text, "DimSize")
    if dim_field is None:
        msg = f"DimSize missing in {mhd_path}"
        raise ValueError(msg)
    parts = [int(v) for v in dim_field.split()]
    if len(parts) == 3:
        nx, ny, nz = parts
        channels_field = _mhd_header_field(text, "ElementNumberOfChannels")
        if channels_field is None:
            msg = (
                f"Expected ElementNumberOfChannels for 3D vector DimSize in {mhd_path}, "
                f"got DimSize = {parts}"
            )
            raise ValueError(msg)
        nc = int(channels_field)
    elif len(parts) == 4:
        nx, ny, nz, nc = parts
    else:
        msg = f"Expected 3D or 4D vector DimSize in {mhd_path}, got {parts}"
        raise ValueError(msg)
    if nc != 3:
        msg = f"Expected 3-vector deformation field, got {nc} components"
        raise ValueError(msg)

    raw_name = _mhd_header_field(text, "ElementDataFile")
    if raw_name is None:
        msg = f"ElementDataFile missing in {mhd_path}"
        raise ValueError(msg)
    element_type = _mhd_header_field(text, "ElementType") or "MET_FLOAT"
    dtype = _mhd_element_dtype(element_type)
    raw_path = mhd_path.parent / raw_name
    flat = np.fromfile(raw_path, dtype=dtype)
    expected = nx * ny * nz * nc
    if flat.size != expected:
        msg = f"RAW size mismatch for vector field {raw_path}: {flat.size} vs {expected}"
        raise ValueError(msg)
    # MetaIO vector image: components vary fastest, then X, Y, Z (ITK x,y,z DimSize).
    vol = flat.reshape((nz, ny, nx, nc), order="C")
    return np.transpose(vol, (1, 2, 0, 3)).astype(np.float32, copy=False)


def run_transformix_points(
    *,
    points_yxz: np.ndarray,
    transform_path: Path,
    output_dir: Path,
    spacing_mm: float,
) -> np.ndarray:
    """Transform 0-based (Y, X, Z) points with ``transformix -def`` (physical mm I/O).

    Elastix transforms map fixed→moving for volume resampling. Applying the same
    parameter file to points with ``-def`` maps feature coordinates in the moving
    image to the fixed domain — the inverse relationship of volume resampling.
    """
    if shutil.which("transformix") is None:
        msg = "transformix not found on PATH"
        raise RuntimeError(msg)

    pts = np.asarray(points_yxz, dtype=float)
    if pts.size == 0:
        return pts.reshape(0, 3)
    if pts.ndim != 2 or pts.shape[1] != 3:
        msg = f"points_yxz must be Nx3, got {pts.shape}"
        raise ValueError(msg)

    output_dir = Path(output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("inputPoints.txt", "outputpoints.txt", "transformix.log"):
        path = output_dir / pattern
        if path.is_file():
            path.unlink()

    input_path = output_dir / "inputPoints.txt"
    phys = volume_indices_to_elastix_physical(pts, spacing_mm, zero_based=True)
    write_landmark_file(input_path, phys)
    cmd = [
        "transformix",
        "-def",
        str(input_path),
        "-out",
        str(output_dir),
        "-tp",
        str(Path(transform_path).expanduser()),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (output_dir / "transformix.stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
    (output_dir / "transformix.stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
    if proc.returncode != 0:
        msg = f"transformix -def points failed (exit {proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
        raise RuntimeError(msg)

    out_path = output_dir / "outputpoints.txt"
    if not out_path.is_file():
        msg = f"transformix did not write {out_path}"
        raise FileNotFoundError(msg)

    sp = float(spacing_mm)
    mapped: list[list[float]] = []
    for line in out_path.read_text(encoding="utf-8").splitlines():
        match = re.search(r"OutputPoint\s*=\s*\[([^\]]+)\]", line)
        if match is None:
            continue
        x_mm, y_mm, z_mm = (float(v) for v in match.group(1).split())
        # ITK physical (x,y,z) mm → volume indices (Y,X,Z)
        mapped.append([y_mm / sp, x_mm / sp, z_mm / sp])
    if len(mapped) != pts.shape[0]:
        msg = f"Expected {pts.shape[0]} transformed points, got {len(mapped)} in {out_path}"
        raise RuntimeError(msg)
    return np.asarray(mapped, dtype=float)


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

    raw_params = transform_path.read_text(encoding="utf-8")
    expected_shape = volume_shape_from_transform_params(raw_params)
    if expected_shape is not None and tuple(moving_volume.shape) != expected_shape:
        msg = (
            f"Moving volume shape {moving_volume.shape} does not match elastix transform "
            f"Size {expected_shape[1]} {expected_shape[0]} {expected_shape[2]} (ITK x,y,z). "
            "Re-run register on the same sample; do not mix TransformParameters from another run."
        )
        raise ValueError(msg)

    params = _patch_transformix_params(raw_params, nearest=nearest)
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
    volume = _read_transformix_result(result_path)
    if nearest:
        volume = np.rint(volume).astype(np.float32)
    return volume
