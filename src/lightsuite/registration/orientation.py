"""Brain volume axis permutation (getBrainOrientation.m / permuteBrainVolume.m)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from lightsuite.config.models import BrainPipelineConfig
from lightsuite.registration.volume import permute_brain_volume

# Map each atlas output axis to a signed source dimension (1-based, MATLAB convention).
PERMUTATION_OPTIONS: list[tuple[int, str]] = [
    (1, "Dim 1 (+)"),
    (-1, "Dim 1 (flipped)"),
    (2, "Dim 2 (+)"),
    (-2, "Dim 2 (flipped)"),
    (3, "Dim 3 (+)"),
    (-3, "Dim 3 (flipped)"),
]

DEFAULT_PERMVEC: list[int] = [1, 2, 3]


def orientation_path(save_path: Path) -> Path:
    return save_path.expanduser() / "brain_orientation.txt"


def validate_permvec(permvec: list[int]) -> None:
    if len(permvec) != 3:
        msg = f"Orientation must have 3 values, got {permvec}"
        raise ValueError(msg)
    if len(set(abs(v) for v in permvec)) != 3:
        msg = f"Orientation dimensions must be unique (ignoring sign): {permvec}"
        raise ValueError(msg)
    for value in permvec:
        if abs(value) not in {1, 2, 3}:
            msg = f"Orientation values must be +/-1, +/-2, or +/-3: {permvec}"
            raise ValueError(msg)


def save_orientation(save_path: Path, permvec: list[int]) -> Path:
    validate_permvec(permvec)
    path = orientation_path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.array(permvec, dtype=int), fmt="%d")
    return path


def load_orientation_file(path: Path) -> list[int]:
    values = np.loadtxt(path.expanduser(), dtype=np.int64)
    permvec = values.flatten().astype(int).tolist()
    validate_permvec(permvec)
    return permvec


def resolve_orientation(config: BrainPipelineConfig, save_path: Path) -> list[int]:
    """Load orientation from config, brain_orientation.txt, or default."""
    if config.registration.orientation is not None:
        permvec = list(config.registration.orientation)
    else:
        orient_file = orientation_path(save_path)
        if orient_file.is_file():
            permvec = load_orientation_file(orient_file)
        else:
            permvec = DEFAULT_PERMVEC.copy()
    validate_permvec(permvec)
    return permvec


def permute_for_atlas(volume: np.ndarray, permvec: list[int]) -> np.ndarray:
    validate_permvec(permvec)
    return permute_brain_volume(volume, permvec)


def permvec_from_indices(indices: tuple[int, int, int]) -> list[int]:
    """Convert 0-based dropdown indices to a permvec."""
    values = [PERMUTATION_OPTIONS[i][0] for i in indices]
    validate_permvec(values)
    return values


def indices_from_permvec(permvec: list[int]) -> tuple[int, int, int]:
    lookup = {value: idx for idx, (value, _label) in enumerate(PERMUTATION_OPTIONS)}
    return lookup[permvec[0]], lookup[permvec[1]], lookup[permvec[2]]
