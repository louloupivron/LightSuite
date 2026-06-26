"""Tests for spinal cord straightening optimizer."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from lightsuite.gui.straighten_cord import (
    clear_all_alignment_points,
    fit_overlay_geometry,
    pop_last_alignment_edit,
)
from lightsuite.registration.straightening import (
    compute_straightening_transforms,
    transform_cord_images_slices,
)
from lightsuite.registration.straightening_optimizer import run_straightening_optimizer, solve_spline, unwrap_sparse
from lightsuite.validation.spinal_parity import compare_alignment_optimizer


def test_unwrap_sparse() -> None:
    theta = np.array([0.0, np.pi - 0.1, -np.pi + 0.1])
    out = unwrap_sparse(theta)
    assert out[1] > out[0]
    assert out[2] > out[1]


def test_solve_spline_returns_finite_curve() -> None:
    n = 10
    obs_idx = np.array([0, 4, 9])
    obs_val = np.array([1.0, 5.0, 2.0])
    fit = solve_spline(n, obs_idx, obs_val, lam=100.0)
    assert fit.shape == (n,)
    assert np.all(np.isfinite(fit))


def test_optimizer_parity_fixture() -> None:
    fixture_root = Path(__file__).resolve().parent / "fixtures" / "spinal_cord"
    ref_path = fixture_root / "parity" / "align_out_reference.json"
    if not ref_path.is_file():
        from tests.fixtures.spinal_cord.build_fixtures import build_fixtures

        build_fixtures(fixture_root)
    report = compare_alignment_optimizer(ref_path)
    assert report.passed, report.messages

    ref = json.loads(ref_path.read_text(encoding="utf-8"))
    fit = run_straightening_optimizer(
        np.asarray(ref["user_cen"]),
        np.asarray(ref["user_ant"]),
        np.asarray(ref["user_pos"]),
    )
    assert fit["fit_x"].shape[0] == len(ref["fit_x"])


def test_fit_overlay_geometry() -> None:
    n = 5
    user_cen = np.full((n, 2), np.nan)
    user_ant = np.full((n, 2), np.nan)
    user_pos = np.full((n, 2), np.nan)
    for z in range(n):
        user_cen[z] = [10.0 + z, 20.0 + z]
        user_ant[z] = [12.0 + z, 18.0 + z]
        user_pos[z] = [8.0 + z, 22.0 + z]
    fit = run_straightening_optimizer(user_cen, user_ant, user_pos)
    geom = fit_overlay_geometry(fit, 2)
    assert geom is not None
    assert geom["center"].shape == (1, 2)
    assert geom["axis"].shape == (2, 2)
    assert np.all(np.isfinite(geom["center"]))
    assert fit_overlay_geometry(fit, 0) is not None
    empty = run_straightening_optimizer(
        np.full((3, 2), np.nan),
        np.full((3, 2), np.nan),
        np.full((3, 2), np.nan),
    )
    assert fit_overlay_geometry(empty, 0) is None


def test_straightening_warp_preserves_intensity() -> None:
    """Regression: row/col vs x/y swap must not flatten straightvol to fill value."""
    n_slices = 5
    user_cen = np.full((n_slices, 2), np.nan)
    user_ant = np.full((n_slices, 2), np.nan)
    user_pos = np.full((n_slices, 2), np.nan)
    for z in range(n_slices):
        user_cen[z] = [20.0, 30.0 + z * 0.5]
        user_ant[z] = [24.0, 28.0 + z * 0.5]
        user_pos[z] = [16.0, 32.0 + z * 0.5]
    fit = run_straightening_optimizer(user_cen, user_ant, user_pos)
    tforms = compute_straightening_transforms(
        fit["fit_x"],
        fit["fit_y"],
        fit["fit_theta"],
        target_center=(40.0, 60.0),
        target_orientation_deg=90.0,
    )
    cordvol = np.zeros((80, 80, n_slices), dtype=np.uint16)
    for z in range(n_slices):
        cy = int(round(fit["fit_y"][z]))
        cx = int(round(fit["fit_x"][z]))
        cordvol[cy - 3 : cy + 4, cx - 3 : cx + 4, z] = 100 + z
    straight = transform_cord_images_slices(cordvol, tforms, output_size=(120, 80), fill_value=1.0)
    assert straight.max() > 1
    assert np.any(straight[:, :, 2] > 50)


def test_clear_all_and_undo_last_point() -> None:
    user_cen = np.array([[1.0, 2.0], [np.nan, np.nan]])
    user_ant = np.array([[3.0, 4.0], [5.0, 6.0]])
    user_pos = np.full((2, 2), np.nan)
    history: list[tuple[int, str, tuple[float, float] | None]] = [
        (0, "cen", None),
        (0, "ant", None),
    ]
    assert pop_last_alignment_edit(history, user_cen, user_ant, user_pos) == 0
    assert np.isnan(user_ant[0, 0])
    assert history == [(0, "cen", None)]
    clear_all_alignment_points(user_cen, user_ant, user_pos)
    assert np.all(np.isnan(user_cen))
    assert np.all(np.isnan(user_ant))
