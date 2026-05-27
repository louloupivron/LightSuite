"""Bayesian Coherent Point Drift registration (pcregisterBCPD.m port)."""

from __future__ import annotations

import platform
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


def _native_bcpd_candidates(configured: Path) -> list[Path]:
    """Suggest native BCPD binaries when a Windows ``.exe`` was configured."""
    candidates: list[Path] = []
    if configured.name.lower() == "bcpd.exe":
        candidates.append(configured.with_name("bcpd"))
        if configured.parent.name.lower() == "win":
            candidates.append(configured.parent.parent / "bcpd")
    return candidates


def find_bcpd_executable(path: str | Path | None = None) -> Path | None:
    """Locate the BCPD binary on PATH or at an explicit location."""
    is_windows = platform.system().lower().startswith("win")

    if path is not None:
        configured = Path(path).expanduser()
        if configured.is_file():
            if not is_windows and configured.suffix.lower() == ".exe":
                for candidate in _native_bcpd_candidates(configured):
                    if candidate.is_file():
                        return candidate.resolve()
                found = shutil.which("bcpd")
                if found:
                    return Path(found).resolve()
            return configured.resolve()

    for name in ("bcpd.exe", "bcpd") if is_windows else ("bcpd", "bcpd.exe"):
        found = shutil.which(name)
        if found:
            resolved = Path(found).resolve()
            if not is_windows and resolved.suffix.lower() == ".exe":
                continue
            return resolved

    local = Path("./bcpd")
    if local.is_file():
        return local.resolve()
    return None


def to_matlab_voxel_points(points: np.ndarray) -> np.ndarray:
    """Convert 0-based array indices to MATLAB 1-based voxel coordinates."""
    return np.asarray(points, dtype=float) + 1.0


def from_matlab_voxel_points(points: np.ndarray) -> np.ndarray:
    """Convert MATLAB 1-based voxel coordinates to 0-based array indices."""
    return np.asarray(points, dtype=float) - 1.0


def _exec_format_hint(executable: Path) -> str:
    if platform.system().lower().startswith("win"):
        return ""
    if executable.suffix.lower() != ".exe":
        return ""
    native = _native_bcpd_candidates(executable)
    hints = [
        "On Linux/macOS, win/bcpd.exe is a Windows binary and cannot be executed.",
        "Build a native binary in the BCPD source tree: make OPT=-DUSE_OPENMP ENV=LINUX",
        "Then set registration.bcpd_path to that bcpd file, or remove bcpd_path and rely on PATH.",
    ]
    if native:
        hints.insert(1, f"Expected native binary alongside the source tree, e.g. {native[0]}.")
    found = shutil.which("bcpd")
    if found:
        hints.append(f"A native bcpd was found on PATH: {found}")
    return "\n".join(hints)


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


def _transform_rows(points: np.ndarray, linear: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """Apply ``points @ linear + translation`` (MATLAB affinetform3d / simtform3d rows)."""
    return points @ linear + translation.reshape(1, 3)


def _similarity_forward(
    points: np.ndarray,
    scale: float,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> np.ndarray:
    """Match MATLAB ``simtform3d`` forward: ``points * R * s + t``."""
    rotation_clean = _clean_rotation(rotation)
    return _transform_rows(points, rotation_clean * scale, translation)


def _similarity_inverse_matrix(
    scale: float,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> np.ndarray:
    """Invert a BCPD/MATLAB similarity (moving->fixed) to fixed->moving."""
    rotation_clean = _clean_rotation(rotation)
    inv_scale = 1.0 / scale
    inv_rotation = rotation_clean.T
    inv_translation = -inv_scale * (translation @ inv_rotation)
    matrix = np.eye(4, dtype=float)
    matrix[:3, :3] = inv_rotation * inv_scale
    matrix[:3, 3] = inv_translation
    return matrix


def _internal_from_affinetform(matrix: np.ndarray) -> np.ndarray:
    """Convert affinetform3d ``points @ A + t`` storage to ``points @ M.T + t`` form."""
    aff = np.asarray(matrix, dtype=float)
    out = np.eye(4, dtype=float)
    out[:3, :3] = aff[:3, :3].T
    out[:3, 3] = aff[:3, 3]
    return out


def _affinetform_from_internal(matrix: np.ndarray) -> np.ndarray:
    """Convert internal transform storage to affinetform3d row convention."""
    internal = np.asarray(matrix, dtype=float)
    out = np.eye(4, dtype=float)
    out[:3, :3] = internal[:3, :3].T
    out[:3, 3] = internal[:3, 3]
    return out


def _build_global_matrix(
    transform_type: TransformType,
    *,
    scale: np.ndarray | None = None,
    rotation: np.ndarray | None = None,
    translation: np.ndarray | None = None,
    affine: np.ndarray | None = None,
) -> np.ndarray:
    """Build moving->fixed affinetform3d matrix (MATLAB row convention)."""
    matrix = np.eye(4, dtype=float)
    translation_vec = np.asarray(translation, dtype=float).ravel()[:3]
    if transform_type in {"similarity", "similarity_nonrigid"}:
        rotation_clean = _clean_rotation(rotation)
        uniform_scale = float(np.asarray(scale, dtype=float).ravel()[0])
        matrix[:3, :3] = rotation_clean * uniform_scale
        matrix[:3, 3] = translation_vec
    elif transform_type == "rigid":
        rotation_clean = _clean_rotation(rotation)
        matrix[:3, :3] = rotation_clean
        matrix[:3, 3] = translation_vec
    elif transform_type in {"affine", "affine_nonrigid"}:
        affine_matrix = np.asarray(affine, dtype=float)
        if affine_matrix.shape == (4, 4):
            matrix[:3, :3] = affine_matrix[:3, :3]
            matrix[:3, 3] = affine_matrix[:3, 3]
        else:
            matrix[:3, :3] = affine_matrix
            matrix[:3, 3] = translation_vec
    return matrix


def _fit_similarity_from_registered(
    moving: np.ndarray,
    registered: np.ndarray,
    *,
    scale: float,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> np.ndarray:
    """Pick the BCPD similarity convention that best reproduces ``registered``."""
    candidates: list[tuple[float, np.ndarray]] = []
    rotation_clean = _clean_rotation(rotation)
    translation_vec = np.asarray(translation, dtype=float).ravel()[:3]
    variants = (
        ("matlab", rotation_clean * scale, translation_vec),
        ("transpose", rotation_clean.T * scale, translation_vec),
        ("scale_left", scale * rotation_clean, translation_vec),
        ("scale_left_transpose", scale * rotation_clean.T, translation_vec),
    )
    for _name, linear, trans in variants:
        predicted = _transform_rows(moving, linear, trans)
        err = float(np.linalg.norm(predicted - registered, axis=1).mean())
        aff = np.eye(4, dtype=float)
        aff[:3, :3] = linear
        aff[:3, 3] = trans
        candidates.append((err, aff))
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1]


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

        try:
            result = subprocess.run(command, capture_output=True, text=True, check=False)
        except OSError as exc:
            if exc.errno == 8:
                hint = _exec_format_hint(executable)
                msg = f"BCPD is not executable on this platform: {executable}"
                if hint:
                    msg = f"{msg}\n{hint}"
                raise RuntimeError(msg) from exc
            raise
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
            uniform_scale = float(np.asarray(scale, dtype=float).ravel()[0])
            moving_to_fixed = _fit_similarity_from_registered(
                moving,
                registered,
                scale=uniform_scale,
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
            predicted = _transform_rows(
                moving,
                moving_to_fixed[:3, :3],
                moving_to_fixed[:3, 3],
            )
            err = float(np.linalg.norm(predicted - registered, axis=1).mean())
            if err > 1.0:
                alt = np.eye(4, dtype=float)
                alt[:3, :3] = _clean_rotation(rotation).T
                alt[:3, 3] = np.asarray(translation, dtype=float).ravel()[:3]
                alt_err = float(
                    np.linalg.norm(
                        _transform_rows(moving, alt[:3, :3], alt[:3, 3]) - registered,
                        axis=1,
                    ).mean()
                )
                if alt_err < err:
                    moving_to_fixed = alt
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
