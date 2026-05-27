"""Tests for BCPD registration wrapper."""

from __future__ import annotations

import numpy as np
import pytest

from lightsuite.registration.bcpd import find_bcpd_executable, register_bcpd


@pytest.fixture(scope="module")
def bcpd_executable() -> str | None:
    found = find_bcpd_executable()
    return str(found) if found is not None else None


@pytest.mark.skipif(find_bcpd_executable() is None, reason="bcpd not installed")
def test_register_bcpd_similarity_roundtrip(bcpd_executable: str) -> None:
    rng = np.random.default_rng(0)
    fixed = rng.random((80, 3)) * 40.0
    scale = 1.08
    rotation = np.array(
        [
            [0.98, -0.17, 0.0],
            [0.17, 0.98, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.array([5.0, -3.0, 2.0])
    moving = scale * (fixed @ rotation.T) + translation

    registered, moving_to_fixed = register_bcpd(
        moving,
        fixed,
        "similarity",
        bcpd_path=bcpd_executable,
        gamma=1.0,
        convergence_tolerance=1e-6,
    )

    assert registered.shape == moving.shape
    err = np.linalg.norm(registered - fixed, axis=1).mean()
    assert err < 3.0

    recovered = _transform_points(moving, moving_to_fixed)
    err_map = np.linalg.norm(recovered - fixed, axis=1).mean()
    assert err_map < 3.0


def _transform_points(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    hom = np.column_stack([points, np.ones(points.shape[0])])
    return (hom @ matrix.T)[:, :3]
