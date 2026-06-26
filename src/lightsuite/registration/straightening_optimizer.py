"""Spinal cord straightening spline optimizer (spinal_cord_aligner.m port)."""

from __future__ import annotations

import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def unwrap_sparse(theta: np.ndarray) -> np.ndarray:
    """Port of unwrap_sparse from spinal_cord_aligner.m."""
    out = theta.astype(float).copy()
    for k in range(1, len(out)):
        diff = out[k] - out[k - 1]
        shift = round(diff / (2 * math.pi))
        out[k] = out[k] - shift * 2 * math.pi
    return out


def solve_spline(n_slices: int, obs_idx: np.ndarray, obs_val: np.ndarray, lam: float) -> np.ndarray:
    """Regularized 1D spline fit over slice index."""
    n = int(n_slices)
    e = np.ones(n)
    d2 = sparse.spdiags([e, -2 * e, e], [-1, 0, 1], n, n, format="csr")
    d2 = d2.tolil()
    d2[0, :] = 0
    d2[0, 0:2] = [1, -1]
    d2[-1, :] = 0
    d2[-1, -2:] = [-1, 1]
    d2 = d2.tocsr()
    lap = d2.T @ d2
    obs_idx = np.asarray(obs_idx, dtype=int)
    obs_val = np.asarray(obs_val, dtype=float)
    w = sparse.csr_matrix(
        (np.ones(obs_idx.shape[0]), (obs_idx, obs_idx)),
        shape=(n, n),
    )
    b = np.zeros(n, dtype=float)
    b[obs_idx] = obs_val
    return spsolve(w + lam * lap, b)


def run_straightening_optimizer(
    user_cen: np.ndarray,
    user_ant: np.ndarray,
    user_pos: np.ndarray,
    *,
    lambda_pos: float = 5000.0,
    lambda_ang: float = 5000.0,
) -> dict[str, np.ndarray]:
    """Compute fit_x, fit_y, fit_theta from user clicks (run_optimizer)."""
    n = user_cen.shape[0]
    obs_x = np.full(n, np.nan)
    obs_y = np.full(n, np.nan)
    obs_th = np.full(n, np.nan)
    obs_rad = np.full(n, np.nan)

    for z in range(n):
        c = user_cen[z]
        a = user_ant[z]
        p = user_pos[z]
        has_c = not np.isnan(c[0])
        has_a = not np.isnan(a[0])
        has_p = not np.isnan(p[0])

        if has_c:
            obs_x[z], obs_y[z] = c
        elif has_a and has_p:
            obs_x[z] = (a[0] + p[0]) / 2
            obs_y[z] = (a[1] + p[1]) / 2

        vec = np.array([np.nan, np.nan])
        curr_rad = np.nan
        if has_a and has_p:
            vec = a - p
            curr_rad = float(np.linalg.norm(vec) / 2)
        elif has_c and has_a:
            vec = a - c
            curr_rad = float(np.linalg.norm(vec))
        elif has_c and has_p:
            vec = c - p
            curr_rad = float(np.linalg.norm(vec))

        if not np.isnan(vec[0]):
            obs_th[z] = math.atan2(vec[1], vec[0])
            obs_rad[z] = curr_rad

    idx_x = np.flatnonzero(~np.isnan(obs_x))
    if idx_x.size > 1:
        fit_x = solve_spline(n, idx_x, obs_x[idx_x], lambda_pos)
        fit_y = solve_spline(n, idx_x, obs_y[idx_x], lambda_pos)
    else:
        fit_x = np.full(n, np.nan)
        fit_y = np.full(n, np.nan)

    idx_th = np.flatnonzero(~np.isnan(obs_th))
    if idx_th.size > 1:
        unwrapped = unwrap_sparse(obs_th[idx_th])
        fit_unwrapped = solve_spline(n, idx_th, unwrapped, lambda_ang)
        fit_theta = np.angle(np.exp(1j * fit_unwrapped))
    else:
        fit_theta = np.full(n, np.nan)

    idx_r = np.flatnonzero(~np.isnan(obs_rad))
    if idx_r.size > 1:
        fit_rad = solve_spline(n, idx_r, obs_rad[idx_r], lambda_pos)
    elif idx_r.size == 1:
        fit_rad = np.full(n, float(obs_rad[idx_r[0]]))
    else:
        fit_rad = np.full(n, 1.0)

    return {
        "obs_x": obs_x,
        "obs_y": obs_y,
        "obs_th": obs_th,
        "fit_x": fit_x,
        "fit_y": fit_y,
        "fit_theta": fit_theta,
        "fit_rad": fit_rad,
    }
