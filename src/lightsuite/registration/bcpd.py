"""Bayesian Coherent Point Drift registration (pcregisterBCPD.m port)."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Literal

import numpy as np

TransformType = Literal[
    "similarity",
    "rigid",
    "affine",
    "affine_nonrigid",
    "similarity_nonrigid",
    "nonrigid",
]

_TRANSFORM_FLAGS: dict[TransformType, str] = {
    "similarity_nonrigid": "-Tsrn",
    "affine_nonrigid": "-Tan",
    "affine": "-Ta",
    "similarity": "-Tsr",
    "rigid": "-Tr",
    "nonrigid": "-Tn",
}

_SAVE_FLAGS: dict[TransformType, str] = {
    "rigid": "cyT",
    "similarity": "cyT",
    "affine": "cyT",
    "nonrigid": "cyTv",
    "similarity_nonrigid": "cyTv",
    "affine_nonrigid": "cyTv",
}


def find_bcpd_executable(path: str | Path | None = None) -> Path | None:
    """Locate the BCPD binary on PATH or at an explicit location."""
    if path is not None:
        candidate = Path(path).expanduser()
        if candidate.is_file():
            return candidate.resolve()
    for name in ("bcpd", "bcpd.exe"):
        found = shutil.which(name)
        if found:
            return Path(found).resolve()
    local = Path("./bcpd")
    if local.is_file():
        return local.resolve()
    return None


def _write_points(path: Path, points: np.ndarray) -> None:
    np.savetxt(path, np.asarray(points, dtype=float), fmt="%.10f")


def _read_matrix(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(path, dtype=float))


def _clean_rotation(rotation: np.ndarray) -> np.ndarray:
    """Orthogonalize rotation and remove reflections (pcregisterBCPD.m)."""
    matrix = np.asarray(rotation, dtype=float)
    u, _, vt = np.linalg.svd(matrix)
    rotation_clean = u @ vt
    if np.linalg.det(rotation_clean) < 0:
        vt = vt.copy()
        vt[-1, :] *= -1.0
        rotation_clean = u @ vt
    return rotation_clean


def _build_global_matrix(
    transform_type: TransformType,
    *,
    scale: np.ndarray | None = None,
    rotation: np.ndarray | None = None,
    translation: np.ndarray | None = None,
    affine: np.ndarray | None = None,
) -> np.ndarray:
    matrix = np.eye(4, dtype=float)
    if transform_type in {"similarity", "similarity_nonrigid"}:
        rotation_clean = _clean_rotation(rotation)
        uniform_scale = float(np.asarray(scale, dtype=float).ravel()[0])
        matrix[:3, :3] = uniform_scale * rotation_clean
        matrix[:3, 3] = np.asarray(translation, dtype=float).ravel()[:3]
    elif transform_type == "rigid":
        rotation_clean = _clean_rotation(rotation)
        matrix[:3, :3] = rotation_clean
        matrix[:3, 3] = np.asarray(translation, dtype=float).ravel()[:3]
    elif transform_type in {"affine", "affine_nonrigid"}:
        affine_matrix = np.asarray(affine, dtype=float)
        if affine_matrix.shape == (4, 4):
            matrix[:3, :3] = affine_matrix[:3, :3]
            matrix[:3, 3] = affine_matrix[:3, 3]
        else:
            matrix[:3, :3] = affine_matrix
            matrix[:3, 3] = np.asarray(translation, dtype=float).ravel()[:3]
    return matrix


def register_bcpd(
    moving: np.ndarray,
    fixed: np.ndarray,
    transform_type: TransformType = "similarity",
    *,
    bcpd_path: Path | str,
    outlier_ratio: float = 0.01,
    gamma: float = 1.0,
    lambda_: float = 50.0,
    beta: float = 2.0,
    convergence_tolerance: float = 1e-8,
    max_iterations: int = 500,
    normalize_common: bool = False,
    verbose: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Register ``moving`` onto ``fixed`` using the external BCPD binary.

    Returns ``(registered_moving, moving_to_fixed)`` where both point sets are
    Nx3 row vectors in the same XYZ frame as the inputs.
    """
    moving = np.asarray(moving, dtype=float)
    fixed = np.asarray(fixed, dtype=float)
    if moving.shape[0] < 3 or fixed.shape[0] < 3:
        msg = "BCPD requires at least three points in each cloud."
        raise ValueError(msg)

    executable = Path(bcpd_path)
    if not executable.is_file():
        msg = f"BCPD executable not found: {executable}"
        raise FileNotFoundError(msg)

    with tempfile.TemporaryDirectory(prefix="lightsuite_bcpd_") as tmpdir:
        tmp = Path(tmpdir)
        fixed_file = tmp / "fixed.txt"
        moving_file = tmp / "moving.txt"
        output_prefix = tmp / "output_"
        _write_points(fixed_file, fixed)
        _write_points(moving_file, moving)

        command = [
            str(executable),
            "-x",
            str(fixed_file),
            "-y",
            str(moving_file),
            "-o",
            str(output_prefix),
            _TRANSFORM_FLAGS[transform_type],
            "-r1",
            "-p",
            "-Db,5000,0.05",
            "-s",
            _SAVE_FLAGS[transform_type],
            f"-w{outlier_ratio:g}",
            f"-b{beta:g}",
            f"-l{lambda_:g}",
            f"-g{gamma:g}",
            f"-c{convergence_tolerance:g}",
            f"-n{max_iterations:d}",
        ]
        if not normalize_common:
            command.append("-ux")
        if not verbose:
            command.append("-q")

        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            msg = (
                "BCPD execution failed "
                f"(exit {result.returncode}).\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
            raise RuntimeError(msg)

        registered_path = Path(f"{output_prefix}y.txt")
        if not registered_path.is_file():
            msg = f"BCPD did not produce output file: {registered_path}"
            raise RuntimeError(msg)
        registered = np.atleast_2d(np.loadtxt(registered_path, dtype=float))

        if transform_type in {"similarity", "similarity_nonrigid"}:
            scale = _read_matrix(Path(f"{output_prefix}s.txt"))
            rotation = _read_matrix(Path(f"{output_prefix}R.txt"))
            translation = _read_matrix(Path(f"{output_prefix}t.txt"))
            moving_to_fixed = _build_global_matrix(
                transform_type,
                scale=scale,
                rotation=rotation,
                translation=translation,
            )
        elif transform_type == "rigid":
            rotation = _read_matrix(Path(f"{output_prefix}R.txt"))
            translation = _read_matrix(Path(f"{output_prefix}t.txt"))
            moving_to_fixed = _build_global_matrix(
                transform_type,
                rotation=rotation,
                translation=translation,
            )
        elif transform_type == "affine":
            affine = _read_matrix(Path(f"{output_prefix}A.txt"))
            translation = _read_matrix(Path(f"{output_prefix}t.txt"))
            moving_to_fixed = _build_global_matrix(
                transform_type,
                affine=affine,
                translation=translation,
            )
        else:
            affine = _read_matrix(Path(f"{output_prefix}A.txt"))
            translation = _read_matrix(Path(f"{output_prefix}t.txt"))
            moving_to_fixed = _build_global_matrix(
                "affine",
                affine=affine,
                translation=translation,
            )

    return registered, moving_to_fixed
