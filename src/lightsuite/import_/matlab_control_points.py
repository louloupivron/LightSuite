"""Load MATLAB ``atlas2histology_tform.mat`` control-point sessions."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat, savemat

from lightsuite.gui.control_points import (
    COORD_SCHEMA_VERSION,
    POINT_COORD_SOURCE_MATLAB,
    ControlPointSession,
)


def find_matlab_control_point_session(save_path: Path) -> Path | None:
    """Return a MATLAB control-point file in ``save_path`` (MATLAB ``*tform.mat`` glob)."""
    save_path = save_path.expanduser()
    preferred = save_path / "atlas2histology_tform.mat"
    if preferred.is_file():
        return preferred
    legacy = save_path / "atlas2histology.mat"
    if legacy.is_file():
        return legacy
    matches = sorted(save_path.glob("*tform.mat"))
    return matches[0] if matches else None


def _as_matrix(value: Any, *, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (4, 4):
        msg = f"Expected {name} to be 4x4, got shape {arr.shape}"
        raise ValueError(msg)
    return arr


def _swap_xy_points(points: np.ndarray) -> np.ndarray:
    """MATLAB multiobjRegistration ``(:, [2 1 3])`` on stored GUI points."""
    pts = np.asarray(points, dtype=float)
    if pts.size == 0:
        return pts.reshape(0, 3)
    return pts[:, [1, 0, 2]]


def _cell_slices_to_lists(cell: Any, *, swap_xy: bool = False) -> list[list[list[float]]]:
    """Convert a MATLAB cell array of Nx3/Nx4 point tables to nested Python lists.

    Points are stored in the same layout as the MATLAB GUI (volume dims 1–3 in
    columns 1–3). The ``[2 1 3]`` axis swap for registration is applied later in
    :meth:`ControlPointSession.paired_points_xyz`, not at import time.
    """
    raw = np.asarray(cell, dtype=object)
    if raw.size == 0:
        return []
    flat = raw.ravel(order="C")
    slices: list[list[list[float]]] = []
    for item in flat:
        if item is None:
            slices.append([])
            continue
        pts = np.asarray(item, dtype=float)
        if pts.size == 0:
            slices.append([])
            continue
        if pts.ndim == 1:
            pts = pts.reshape(1, -1)
        row = pts[:, :3]
        if swap_xy:
            row = _swap_xy_points(row)
        slices.append(row.tolist())
    return slices


def _infer_chooserow_from_points(points: list[list[float]]) -> list[int] | None:
    """Recover ``[slice_index, cut_axis, 0, 0]`` from one cell of MATLAB sample points.

    ``matchControlPoints_unified.m`` writes ``cpt(idim) = chooselist(curr_slice, 1)``,
    so the cut axis is the single column that is a constant integer across the cell.
    The two blanking flags are not recoverable and are left as 0 ("unknown").
    """
    if not points:
        return None
    pts = np.asarray(points, dtype=float)[:, :3]
    matches = [
        (axis + 1, int(round(float(pts[0, axis]))))
        for axis in range(3)
        if np.allclose(pts[:, axis], np.round(pts[:, axis]))
        and np.allclose(pts[:, axis], pts[0, axis])
    ]
    if len(matches) != 1:
        return None
    axis, index = matches[0]
    return [index, axis, 0, 0]


def infer_chooselist_from_points(
    histology: list[list[list[float]]],
    *,
    fallback: np.ndarray | list[list[int]] | None = None,
) -> list[list[int]] | None:
    """Rebuild the MATLAB chooselist for cells that carry sample control points.

    Cells without points keep ``fallback`` (the generated list) when supplied,
    otherwise they are marked unknown with a zero row.
    """
    if not histology:
        return None
    base = np.asarray(fallback, dtype=int) if fallback is not None else None
    rows: list[list[int]] = []
    recovered = 0
    for idx, points in enumerate(histology):
        row = _infer_chooserow_from_points(points)
        if row is not None:
            recovered += 1
        elif base is not None and idx < base.shape[0]:
            row = [int(v) for v in base[idx]]
        else:
            row = [0, 0, 0, 0]
        rows.append(row)
    return rows if recovered else None


def load_control_point_session_from_mat(
    path: Path,
    *,
    original_trans: np.ndarray | list[list[float]] | None = None,
    fallback_chooselist: np.ndarray | list[list[int]] | None = None,
) -> ControlPointSession:
    """Load ``ControlPointSession`` from a MATLAB ``atlas2histology_tform.mat`` file.

    Point tables are copied as stored in MATLAB (1-based volume indices, optional
    timestamp column). ``ori_trans`` is taken from the file when present; otherwise
    ``original_trans`` is used so register matches MATLAB ``multiobjRegistration.m``.

    The slice geometry each point was placed on is recovered from the points
    themselves unless the file stores ``chooselist``; ``fallback_chooselist`` fills
    the cells that hold no points.
    """
    path = path.expanduser()
    if not path.is_file():
        msg = f"MATLAB control-point file not found: {path}"
        raise FileNotFoundError(msg)

    data = loadmat(path, squeeze_me=True, struct_as_record=False)
    keys = {key for key in data if not key.startswith("__")}

    if "histology_control_points" not in keys or "atlas_control_points" not in keys:
        msg = (
            f"{path.name} must contain histology_control_points and atlas_control_points "
            f"(found: {sorted(keys)})."
        )
        raise ValueError(msg)

    histology = _cell_slices_to_lists(data["histology_control_points"])
    atlas = _cell_slices_to_lists(data["atlas_control_points"])
    n_slices = max(len(histology), len(atlas))
    if len(histology) < n_slices:
        histology.extend([[] for _ in range(n_slices - len(histology))])
    if len(atlas) < n_slices:
        atlas.extend([[] for _ in range(n_slices - len(atlas))])

    if "atlas2histology_tform" in keys:
        manual_tform = _as_matrix(data["atlas2histology_tform"], name="atlas2histology_tform")
    elif "histology_ccf_manual_alignment" in keys:
        manual_tform = _as_matrix(
            data["histology_ccf_manual_alignment"],
            name="histology_ccf_manual_alignment",
        )
    else:
        manual_tform = np.eye(4)

    if "ori_trans" in keys:
        ori_trans = _as_matrix(data["ori_trans"], name="ori_trans")
    elif original_trans is not None:
        ori_trans = np.asarray(original_trans, dtype=float)
    else:
        ori_trans = np.eye(4)

    if "chooselist" in keys:
        chooselist = np.asarray(data["chooselist"], dtype=int).reshape(-1, 4).tolist()
    else:
        chooselist = infer_chooselist_from_points(histology, fallback=fallback_chooselist)

    return ControlPointSession(
        atlas2histology_tform=manual_tform.tolist(),
        histology_control_points=histology,
        atlas_control_points=atlas,
        ori_trans=ori_trans.tolist(),
        chooselist=chooselist,
        coord_schema_version=COORD_SCHEMA_VERSION,
        point_coord_source=POINT_COORD_SOURCE_MATLAB,
    )


def import_matlab_control_points(
    mat_path: Path,
    *,
    output_json: Path | None = None,
    original_trans: np.ndarray | list[list[float]] | None = None,
) -> Path:
    """Convert a MATLAB session file to ``atlas2histology_tform.json``."""
    session = load_control_point_session_from_mat(mat_path, original_trans=original_trans)
    out = output_json.expanduser() if output_json is not None else mat_path.with_suffix(".json")
    session.save(out)
    return out


def _lists_to_cell_slices(
    slices: list[list[list[float]]],
    *,
    four_columns: bool,
) -> np.ndarray:
    """Convert nested Python lists to an ``Nx1`` MATLAB cell array of point tables."""
    cell = np.empty((len(slices), 1), dtype=object)
    for row_idx, pts_list in enumerate(slices):
        if not pts_list:
            cell[row_idx, 0] = np.zeros((0, 4 if four_columns else 3), dtype=float)
            continue
        pts = np.asarray(pts_list, dtype=float)
        if pts.ndim == 1:
            pts = pts.reshape(1, -1)
        n_cols = 4 if four_columns or pts.shape[1] >= 4 else 3
        out = np.zeros((pts.shape[0], n_cols), dtype=float)
        xyz = pts[:, :3]
        out[:, :3] = xyz
        if n_cols == 4:
            out[:, 3] = pts[:, 3] if pts.shape[1] >= 4 else 0.0
        cell[row_idx, 0] = out
    return cell


def export_control_point_session_to_mat(
    session: ControlPointSession,
    path: Path,
    *,
    include_ori_trans: bool = True,
) -> Path:
    """Write ``atlas2histology_tform.mat`` for MATLAB ``multiobjRegistration.m``.

    Point tables are stored in the same layout as ``matchControlPoints_unified.m``
    (inverse of :func:`load_control_point_session_from_mat`).
    """
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)

    histology = _lists_to_cell_slices(
        session.histology_control_points,
        four_columns=True,
    )
    atlas_has_ts = any(
        len(pts) >= 4 for slice_pts in session.atlas_control_points for pts in slice_pts
    )
    atlas = _lists_to_cell_slices(
        session.atlas_control_points,
        four_columns=atlas_has_ts,
    )

    payload: dict[str, Any] = {
        "atlas2histology_tform": np.asarray(session.atlas2histology_tform, dtype=float),
        "histology_control_points": histology,
        "atlas_control_points": atlas,
    }
    if include_ori_trans:
        payload["ori_trans"] = np.asarray(session.ori_trans, dtype=float)
    if session.chooselist:
        payload["chooselist"] = np.asarray(session.chooselist, dtype=float)

    savemat(path, payload, do_compression=False)
    return path


def export_matlab_control_points(
    session_path: Path,
    *,
    output_mat: Path | None = None,
) -> Path:
    """Convert ``atlas2histology_tform.json`` to MATLAB ``atlas2histology_tform.mat``."""
    session_path = session_path.expanduser()
    if not session_path.is_file():
        msg = f"Control-point session not found: {session_path}"
        raise FileNotFoundError(msg)
    session = (
        ControlPointSession.load(session_path)
        if session_path.suffix.lower() == ".json"
        else load_control_point_session_from_mat(session_path)
    )
    out = output_mat.expanduser() if output_mat is not None else session_path.with_name(
        "atlas2histology_tform.mat"
    )
    return export_control_point_session_to_mat(session, out)
