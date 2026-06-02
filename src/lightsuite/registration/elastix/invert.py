"""Invert elastix transforms (invertElastixTransformCP.m port)."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

from lightsuite.registration.elastix.mhd import read_mhd_spacing


def invert_elastix_transform(transform_dir: Path, output_dir: Path | None = None) -> Path:
    """Invert forward elastix B-spline transform; return inverted TransformParameters path."""
    transform_dir = transform_dir.expanduser()
    if output_dir is None:
        output_dir = transform_dir.parent / "elastix_inverse_temp"
    output_dir = output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = transform_dir / "elastix.log"
    if not log_path.is_file():
        msg = f"Missing elastix.log in {transform_dir}"
        raise FileNotFoundError(msg)

    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    fixed_mhd = _parse_fixed_image_path(log_text, transform_dir)
    param_paths = _parse_parameter_file_paths(log_text, transform_dir)
    if not param_paths:
        msg = f"Could not parse parameter file paths from {log_path}"
        raise RuntimeError(msg)

    coef_files = sorted(transform_dir.glob("TransformParameters.*.txt"), reverse=True)
    if not coef_files:
        msg = f"No TransformParameters.*.txt in {transform_dir}"
        raise FileNotFoundError(msg)

    _ = read_mhd_spacing(fixed_mhd)

    # Invert with the elastix-documented DisplacementMagnitudePenalty metric rather than
    # re-running the forward Mutual-Information registration of the fixed image against
    # itself. MI(self) is flat/noisy near the optimum, so at the finest resolution the
    # B-spline coefficients drift and fold (grid-period striping in the exported volume).
    # DisplacementMagnitudePenalty directly minimises ||T_inv(T_fwd(x)) - x||: no image
    # content, no landmarks, and a clean monotone convergence to the true inverse.
    inverse_param = _write_inverse_parameter_file(param_paths[0], output_dir)

    cmd = [
        "elastix",
        "-f",
        str(fixed_mhd),
        "-m",
        str(fixed_mhd),
        "-out",
        str(output_dir),
        "-t0",
        str(coef_files[0]),
        "-p",
        str(inverse_param),
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        msg = f"elastix inversion failed (exit {proc.returncode}):\n{proc.stdout}\n{proc.stderr}"
        raise RuntimeError(msg)

    inverted = sorted(output_dir.glob("TransformParameters.*.txt"))
    if not inverted:
        msg = f"No inverted TransformParameters in {output_dir}"
        raise RuntimeError(msg)

    text = inverted[0].read_text(encoding="utf-8")
    if "InitialTransformParametersFileName" not in text:
        text += '\n(InitialTransformParametersFileName "NoInitialTransform")\n'
        inverted[0].write_text(text, encoding="utf-8")
    return inverted[0]


def _param_line_key(line: str) -> str | None:
    match = re.match(r"\(\s*([A-Za-z0-9_]+)\b", line.strip())
    return match.group(1) if match else None


def _first_quoted_value(line: str) -> str | None:
    match = re.search(r'"([^"]+)"', line)
    return match.group(1) if match else None


def _write_inverse_parameter_file(forward_param: Path, output_dir: Path) -> Path:
    """Derive a DisplacementMagnitudePenalty inverse parameter file from the forward one.

    Keeps the transform/grid/pyramid/optimizer schedule but replaces the image+landmark
    metrics with a single displacement penalty (elastix manual, "Inverting a transform").
    """
    text = forward_param.read_text(encoding="utf-8")
    single_value_keys = {
        "FixedImagePyramid",
        "MovingImagePyramid",
        "ImageSampler",
        "Interpolator",
    }
    drop_keys = {"Metric0Weight", "Metric1Weight", "Metric2Weight"}
    out_lines: list[str] = []
    for line in text.splitlines():
        key = _param_line_key(line)
        if key in drop_keys:
            continue
        if key == "Registration":
            out_lines.append('(Registration "MultiResolutionRegistration")')
            continue
        if key == "Metric":
            out_lines.append('(Metric "DisplacementMagnitudePenalty")')
            continue
        if key in single_value_keys:
            value = _first_quoted_value(line)
            if value is not None:
                out_lines.append(f'({key} "{value}")')
                continue
        out_lines.append(line)
    inverse_path = output_dir / "inverse_parameters.txt"
    inverse_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return inverse_path


def write_inverted_transform_copy(source: Path, destination: Path) -> Path:
    destination = destination.expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    return destination


def _parse_fixed_image_path(log_text: str, transform_dir: Path) -> Path:
    patterns = [
        r"-f0\s+(.+?)\s+-m0\b",
        r"-f\s+(.+?)(?=\s+(?<![fm])-m\b)",
    ]
    for pattern in patterns:
        matches = re.findall(pattern, log_text, flags=re.DOTALL)
        if matches:
            candidate = _strip_quotes(matches[-1].strip())
            path = Path(candidate)
            if path.is_file():
                return path
    guesses = list(transform_dir.glob("*_dual_f0.mhd")) + list(transform_dir.glob("fixed.mhd"))
    if guesses:
        return guesses[0]
    msg = f"Could not determine fixed image path from elastix log in {transform_dir}"
    raise RuntimeError(msg)


def _parse_parameter_file_paths(log_text: str, transform_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for pattern in (
        r"end of ParameterFile:\s*([^\r\n]+)",
        r'(?<![fm])-p\s+"([^"]+)"',
        r"(?<![fm])-p\s+(\S+)",
    ):
        for match in re.findall(pattern, log_text):
            cleaned = _strip_quotes(str(match).strip().split("=")[0])
            if cleaned:
                paths.append(Path(cleaned))
    unique: list[Path] = []
    seen: set[str] = set()
    for path in paths:
        key = str(path)
        if key not in seen and path.is_file():
            unique.append(path)
            seen.add(key)
    if not unique:
        for path in transform_dir.glob("*parameters*.txt"):
            if path.is_file():
                unique.append(path)
    return unique


def _strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    return value
