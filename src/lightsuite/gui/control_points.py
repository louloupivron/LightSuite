"""Control point session persistence (atlas2histology_tform.mat equivalent)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

COORD_SCHEMA_VERSION = 2
POINT_COORD_SOURCE_MATLAB = "matlab"
POINT_COORD_SOURCE_NAPARI = "napari"


def in_plane_volume_axes(cut_axis_1based: int) -> tuple[int, int]:
    """0-based volume axes in the slice plane for a 1-based cut axis."""
    cut = int(cut_axis_1based) - 1
    axes = [d for d in range(3) if d != cut]
    return axes[0], axes[1]


def swap_in_plane_volume_coords(point: list[float], cut_axis_1based: int) -> list[float]:
    """Swap the two in-plane components (undo legacy napari transpose storage)."""
    out = list(point)
    axis_a, axis_b = in_plane_volume_axes(cut_axis_1based)
    out[axis_a], out[axis_b] = out[axis_b], out[axis_a]
    return out


def migrate_napari_transposed_control_points(
    session: "ControlPointSession",
    chooselist: np.ndarray | list[list[int]],
) -> int:
    """Fix in-plane axis swap from pre-v2 napari match-points sessions.

    Returns the number of point rows updated. Skips sessions imported from MATLAB
    (``point_coord_source == "matlab"``) or already at ``coord_schema_version >= 2``.
    """
    if session.coord_schema_version >= COORD_SCHEMA_VERSION:
        return 0
    if session.point_coord_source == POINT_COORD_SOURCE_MATLAB:
        session.coord_schema_version = COORD_SCHEMA_VERSION
        return 0

    rows = np.asarray(chooselist, dtype=int).reshape(-1, 4)
    updated = 0
    for slice_idx, chooserow in enumerate(rows):
        cut_axis = int(chooserow[1])
        if cut_axis not in (1, 2, 3):
            continue
        for store in (session.histology_control_points, session.atlas_control_points):
            pts = store[slice_idx]
            if not pts:
                continue
            for pt_idx, point in enumerate(pts):
                store[slice_idx][pt_idx] = swap_in_plane_volume_coords(point, cut_axis)
                updated += 1
    session.coord_schema_version = COORD_SCHEMA_VERSION
    if session.point_coord_source is None:
        session.point_coord_source = POINT_COORD_SOURCE_NAPARI
    return updated


def session_needs_napari_transpose_migration(session: ControlPointSession) -> bool:
    """True when a legacy Napari session likely has transposed in-plane storage."""
    return (
        session.coord_schema_version < COORD_SCHEMA_VERSION
        and session.point_coord_source != POINT_COORD_SOURCE_MATLAB
    )


def through_axis_column_for_cut(cut_axis_1based: int) -> int:
    """Column in ``paired_points_xyz`` ``[X, Y, Z]`` for the slice-normal volume axis."""
    mapping = {1: 1, 2: 0, 3: 2}
    cut = int(cut_axis_1based)
    if cut not in mapping:
        msg = f"cut_axis must be 1, 2, or 3, got {cut_axis_1based}"
        raise ValueError(msg)
    return mapping[cut]


def mark_session_saved_from_napari(session: ControlPointSession) -> None:
    """Record that coordinates were written with the v2 napari (row, col) schema."""
    session.coord_schema_version = COORD_SCHEMA_VERSION
    if session.point_coord_source != POINT_COORD_SOURCE_MATLAB:
        session.point_coord_source = POINT_COORD_SOURCE_NAPARI


@dataclass
class ControlPointSession:
    """User-defined control points and manual atlas-to-sample alignment."""

    atlas2histology_tform: list[list[float]]
    histology_control_points: list[list[list[float]]]
    atlas_control_points: list[list[list[float]]]
    ori_trans: list[list[float]]
    chooselist: list[list[int]] | None = None
    # Per chooselist entry: atlas plane index along the cut axis (1-based). None = auto on load.
    atlas_slice_indices: list[int] | None = None
    coord_schema_version: int = COORD_SCHEMA_VERSION
    point_coord_source: str | None = None

    @classmethod
    def empty(cls, ori_trans: np.ndarray, n_slices: int) -> ControlPointSession:
        return cls(
            atlas2histology_tform=np.eye(4).tolist(),
            histology_control_points=[[] for _ in range(n_slices)],
            atlas_control_points=[[] for _ in range(n_slices)],
            ori_trans=ori_trans.tolist(),
            coord_schema_version=COORD_SCHEMA_VERSION,
            point_coord_source=POINT_COORD_SOURCE_NAPARI,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "atlas2histology_tform": self.atlas2histology_tform,
            "histology_control_points": self.histology_control_points,
            "atlas_control_points": self.atlas_control_points,
            "ori_trans": self.ori_trans,
            "chooselist": self.chooselist,
            "atlas_slice_indices": self.atlas_slice_indices,
            "coord_schema_version": self.coord_schema_version,
            "point_coord_source": self.point_coord_source,
        }

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> ControlPointSession:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        raw.setdefault("atlas_slice_indices", None)
        raw.setdefault("chooselist", None)
        raw.setdefault("ori_trans", np.eye(4).tolist())
        raw.setdefault("coord_schema_version", 1)
        raw.setdefault("point_coord_source", None)
        if "atlas2histology_tform" not in raw:
            raw["atlas2histology_tform"] = np.eye(4).tolist()
        return cls(**raw)

    def point_counts(self) -> tuple[int, int, int]:
        """Return (matched_pair_count, total_sample_points, total_atlas_points)."""
        total_sample = sum(len(s) for s in self.histology_control_points)
        total_atlas = sum(len(s) for s in self.atlas_control_points)
        matched = sum(
            len(a_list)
            for a_list, s_list in zip(
                self.atlas_control_points,
                self.histology_control_points,
                strict=True,
            )
            if len(a_list) == len(s_list) and len(a_list) > 0
        )
        return matched, total_sample, total_atlas

    def paired_points_xyz(self) -> tuple[np.ndarray, np.ndarray]:
        """Return matched atlas/sample points as Nx3 (MATLAB dim order [2,1,3])."""
        atlas_pts: list[np.ndarray] = []
        sample_pts: list[np.ndarray] = []
        for a_list, s_list in zip(self.atlas_control_points, self.histology_control_points, strict=True):
            if len(a_list) != len(s_list) or len(a_list) == 0:
                continue
            a = np.asarray(a_list, dtype=float)[:, :3]
            s = np.asarray(s_list, dtype=float)[:, :3]
            atlas_pts.append(a[:, [1, 0, 2]])
            sample_pts.append(s[:, [1, 0, 2]])
        if not atlas_pts:
            return np.zeros((0, 3)), np.zeros((0, 3))
        return np.vstack(atlas_pts), np.vstack(sample_pts)

    def paired_points_volume_yxz(self) -> tuple[np.ndarray, np.ndarray]:
        """Return matched points in 0-based volume array order (Y, X, Z).

        Cord match-points store in-plane pixel indices as 0-based volume indices and
        the slice coordinate on the cut axis as 1-based (MATLAB-style).
        """
        atlas_pts: list[np.ndarray] = []
        sample_pts: list[np.ndarray] = []
        for a_list, s_list in zip(self.atlas_control_points, self.histology_control_points, strict=True):
            if len(a_list) != len(s_list) or len(a_list) == 0:
                continue
            a = np.asarray(a_list, dtype=float)[:, :3].copy()
            s = np.asarray(s_list, dtype=float)[:, :3].copy()
            a[:, 2] -= 1.0
            s[:, 2] -= 1.0
            atlas_pts.append(a)
            sample_pts.append(s)
        if not atlas_pts:
            return np.zeros((0, 3)), np.zeros((0, 3))
        return np.vstack(atlas_pts), np.vstack(sample_pts)

    def update_manual_alignment(
        self,
        min_pairs: int = 16,
        *,
        fallback_tform: np.ndarray | list[list[float]] | None = None,
        constrain_cut_axis: int | None = None,
    ) -> float | None:
        """Recompute manual alignment from paired slices; return MSE if fit.

        Below ``min_pairs`` (MATLAB ``Nmin``), keeps or restores ``fallback_tform``
        instead of fitting an under-constrained affine.

        When ``constrain_cut_axis`` is set (spinal cord match-points), the through-plane
        coordinate is copied from sample to atlas before fitting so a global affine only
        refines in-plane alignment. Longitudinal correspondence and elastix already
        place atlas and sample in the same straightened grid; atlas-plane scrolling only
        picks the visible slice and should not drive a 3D Z shear in the overlay.
        """
        from lightsuite.gui.affine import fit_affine_transform

        atlas_pts, sample_pts = self.paired_points_xyz()
        if atlas_pts.shape[0] < min_pairs:
            if fallback_tform is not None:
                self.atlas2histology_tform = np.asarray(fallback_tform, dtype=float).tolist()
            return None
        if constrain_cut_axis is not None:
            atlas_pts = atlas_pts.copy()
            col = through_axis_column_for_cut(constrain_cut_axis)
            atlas_pts[:, col] = sample_pts[:, col]
        matrix, mse = fit_affine_transform(atlas_pts, sample_pts)
        self.atlas2histology_tform = matrix.tolist()
        return mse


def default_session_path(save_path: Path) -> Path:
    return save_path / "atlas2histology_tform.json"


def load_registration_control_point_session(
    save_path: Path,
    *,
    original_trans: list[list[float]] | np.ndarray,
    prefer_matlab: bool = False,
) -> ControlPointSession:
    """Load manual control points if present; otherwise an empty session (MATLAB optional *tform.mat)."""
    from lightsuite.import_.matlab_control_points import (
        find_matlab_control_point_session,
        load_control_point_session_from_mat,
    )

    save_path = save_path.expanduser()
    json_path = default_session_path(save_path)
    mat_path = find_matlab_control_point_session(save_path)
    matrix = np.asarray(original_trans, dtype=float)

    if prefer_matlab and mat_path is not None:
        return load_control_point_session_from_mat(mat_path, original_trans=matrix)
    if json_path.is_file():
        return ControlPointSession.load(json_path)
    if mat_path is not None:
        return load_control_point_session_from_mat(mat_path, original_trans=matrix)
    return ControlPointSession.empty(matrix, n_slices=1)
