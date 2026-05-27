"""Tests for BCPD registration wrapper."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from lightsuite.registration.bcpd import (
    _fit_similarity_from_registered,
    _native_bcpd_candidates,
    _transform_rows,
    find_bcpd_executable,
    register_bcpd,
    to_matlab_voxel_points,
)


def test_native_bcpd_candidates_from_win_exe() -> None:
    configured = Path("/opt/bcpd-master/win/bcpd.exe")
    candidates = _native_bcpd_candidates(configured)
    assert candidates[0] == Path("/opt/bcpd-master/win/bcpd")
    assert candidates[1] == Path("/opt/bcpd-master/bcpd")


def test_find_bcpd_prefers_native_binary_over_win_exe(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source_root = tmp_path / "bcpd-master"
    win_dir = source_root / "win"
    win_dir.mkdir(parents=True)
    (win_dir / "bcpd.exe").write_bytes(b"MZ")
    native = source_root / "bcpd"
    native.write_text("#!/bin/sh\necho bcpd\n", encoding="utf-8")
    native.chmod(0o755)

    monkeypatch.delenv("PATH", raising=False)
    found = find_bcpd_executable(win_dir / "bcpd.exe")
    assert found == native.resolve()


def test_fit_similarity_prefers_matlab_row_convention() -> None:
    rng = np.random.default_rng(0)
    moving = rng.random((40, 3)) * 30.0
    rotation = np.array(
        [
            [0.96, -0.28, 0.0],
            [0.28, 0.96, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    scale = 1.12
    translation = np.array([3.0, -5.0, 2.0])
    registered = _transform_rows(moving, rotation * scale, translation)
    matrix = _fit_similarity_from_registered(
        moving,
        registered,
        scale=scale,
        rotation=rotation,
        translation=translation,
    )
    predicted = _transform_rows(moving, matrix[:3, :3], matrix[:3, 3])
    assert np.linalg.norm(predicted - registered, axis=1).mean() < 1e-6


def test_to_matlab_voxel_points() -> None:
    points = np.array([[0.0, 1.0, 2.0]])
    assert np.allclose(to_matlab_voxel_points(points), [[1.0, 2.0, 3.0]])


@pytest.fixture(scope="module")
def bcpd_executable() -> str | None:
    from lightsuite.registration.bcpd import find_bcpd_executable

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
    moving = _transform_rows(fixed, rotation * scale, translation)

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

    recovered = _transform_rows(moving, moving_to_fixed[:3, :3], moving_to_fixed[:3, 3])
    err_map = np.linalg.norm(recovered - fixed, axis=1).mean()
    assert err_map < 3.0
