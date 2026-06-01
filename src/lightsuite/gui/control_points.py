"""Control point session persistence (atlas2histology_tform.mat equivalent)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


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

    @classmethod
    def empty(cls, ori_trans: np.ndarray, n_slices: int) -> ControlPointSession:
        return cls(
            atlas2histology_tform=np.eye(4).tolist(),
            histology_control_points=[[] for _ in range(n_slices)],
            atlas_control_points=[[] for _ in range(n_slices)],
            ori_trans=ori_trans.tolist(),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "atlas2histology_tform": self.atlas2histology_tform,
            "histology_control_points": self.histology_control_points,
            "atlas_control_points": self.atlas_control_points,
            "ori_trans": self.ori_trans,
            "chooselist": self.chooselist,
            "atlas_slice_indices": self.atlas_slice_indices,
        }

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> ControlPointSession:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        raw.setdefault("atlas_slice_indices", None)
        return cls(**raw)

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

    def update_manual_alignment(self, min_pairs: int = 16) -> float | None:  # MATLAB Nmin
        """Recompute manual alignment from paired slices; return MSE if fit."""
        from lightsuite.gui.affine import fit_affine_transform

        atlas_pts, sample_pts = self.paired_points_xyz()
        if atlas_pts.shape[0] < min_pairs:
            return None
        matrix, mse = fit_affine_transform(atlas_pts, sample_pts)
        self.atlas2histology_tform = matrix.tolist()
        return mse


def default_session_path(save_path: Path) -> Path:
    return save_path / "atlas2histology_tform.json"
