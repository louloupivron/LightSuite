"""Generate control-point slice list (generate_cp_list_alt.m)."""

from __future__ import annotations

import numpy as np


def default_ap_cut_axis(volume_shape: tuple[int, int, int]) -> int:
    """Heuristic AP axis: longest volume extent (1-based axis index)."""
    return int(np.argmax(volume_shape) + 1)


def generate_ap_alignment_list(
    volume_shape: tuple[int, int, int],
    *,
    cut_axis: int | None = None,
    n_slices: int = 20,
) -> np.ndarray:
    """Return chooselist rows for AP-only slice alignment: [index, cut_axis, 1, 1]."""
    axis = cut_axis if cut_axis is not None else default_ap_cut_axis(volume_shape)
    if axis not in {1, 2, 3}:
        msg = f"cut_axis must be 1, 2, or 3, got {axis}"
        raise ValueError(msg)
    axis_size = volume_shape[axis - 1]
    nmin = min(volume_shape)
    minstart = max(1, int(np.ceil(nmin / 20)))
    if axis_size <= 2 * minstart:
        minstart = 1
    sids = np.round(np.linspace(minstart, axis_size - minstart, n_slices)).astype(int)
    sids = np.clip(sids, 1, axis_size)
    return np.column_stack(
        [
            sids,
            np.full(n_slices, axis, dtype=int),
            np.ones(n_slices, dtype=int),
            np.ones(n_slices, dtype=int),
        ]
    ).astype(int)


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


def generate_cord_longitudinal_list(
    length_axis_size: int,
    n_slices: int = 20,
    *,
    cut_axis: int = 3,
) -> np.ndarray:
    """Return ordered chooselist rows for longitudinal align (rostrocaudal anchors)."""
    if cut_axis not in {1, 2, 3}:
        msg = f"cut_axis must be 1, 2, or 3, got {cut_axis}"
        raise ValueError(msg)
    if length_axis_size < 1:
        msg = f"length_axis_size must be >= 1, got {length_axis_size}"
        raise ValueError(msg)
    minstart = max(1, int(np.ceil(length_axis_size / 20)))
    if length_axis_size <= 2 * minstart:
        minstart = 1
    sids = np.round(np.linspace(minstart, length_axis_size - minstart, n_slices)).astype(int)
    sids = np.clip(sids, 1, length_axis_size)
    return np.column_stack(
        [
            sids,
            np.full(n_slices, cut_axis, dtype=int),
            np.ones(n_slices, dtype=int),
            np.ones(n_slices, dtype=int),
        ]
    ).astype(int)


def generate_cord_control_point_list(
    length_axis_size: int,
    n_slices: int = 100,
    *,
    cut_axis: int = 3,
) -> np.ndarray:
    """Return chooselist rows for spinal cord match-points (matchControlPointsSpine.m).

    Cord registration always cuts transverse slices along the rostrocaudal axis (volume
    dimension 3 after straightening). Rows are ``[slice_index, cut_axis, 1, 1]`` and
    shuffled with seed 1 to match MATLAB ``rng(1); randperm(...)``.
    """
    if cut_axis not in {1, 2, 3}:
        msg = f"cut_axis must be 1, 2, or 3, got {cut_axis}"
        raise ValueError(msg)
    if length_axis_size < 1:
        msg = f"length_axis_size must be >= 1, got {length_axis_size}"
        raise ValueError(msg)
    sids = np.round(np.linspace(1, length_axis_size, n_slices)).astype(int)
    sids = np.clip(sids, 1, length_axis_size)
    chooselist = np.column_stack(
        [
            sids,
            np.full(n_slices, cut_axis, dtype=int),
            np.ones(n_slices, dtype=int),
            np.ones(n_slices, dtype=int),
        ]
    )
    rng = np.random.default_rng(1)
    return chooselist[rng.permutation(n_slices)].astype(int)
