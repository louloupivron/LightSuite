"""Generate control-point slice list (generate_cp_list_alt.m)."""

from __future__ import annotations

import numpy as np


def generate_control_point_list(volume_shape: tuple[int, int, int]) -> np.ndarray:
    """Return chooselist array with columns [slice_index, axis, flag_a, flag_b]."""
    ny, nx, nz = volume_shape
    nmin = min(ny, nx, nz)
    minstart = int(np.ceil(nmin / 20))
    n_per_side = 20
    indmat = np.array([[1, 1], [1, 2], [2, 1], [2, 2]], dtype=int)
    n_types = indmat.shape[0]

    data_all: list[np.ndarray] = []
    for idim in range(3):
        axis_size = volume_shape[idim]
        sids = np.round(np.linspace(minstart, axis_size - minstart, n_per_side * n_types)).astype(int)
        sids = sids.reshape(n_types, n_per_side)
        for ii in range(n_types):
            block = np.column_stack(
                [
                    sids[ii, :],
                    np.full(n_per_side, idim + 1, dtype=int),
                    np.full(n_per_side, indmat[ii, 0], dtype=int),
                    np.full(n_per_side, indmat[ii, 1], dtype=int),
                ]
            )
            data_all.append(block)

    minidx = min(arr.shape[0] for arr in data_all)
    data_all = [arr[:minidx, :] for arr in data_all]
    cplist = np.vstack(data_all)

    rng = np.random.default_rng(1)
    iperm = rng.permutation(n_per_side)
    reshaped = cplist.reshape(n_per_side, 3, n_types, n_types)
    cplist = reshaped[iperm, :, :, :].transpose(3, 1, 2, 0).reshape(n_types, -1).T
    return cplist.astype(int)
