"""Landmark session persistence for multiresolution overview / ROI placement."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

import numpy as np

LandmarkFitMode = Literal["similarity", "affine", "rigid"]


@dataclass
class MultiresLandmarkSession:
    """Paired landmark points linking ROI voxels to overview voxels."""

    overview_points_zyx: list[list[float]] = field(default_factory=list)
    roi_points_zyx: list[list[float]] = field(default_factory=list)
    fit_mode: LandmarkFitMode = "similarity"
    roi_to_overview_tform: list[list[float]] | None = None
    rms_error_um: float | None = None
    fit_point_errors_um: list[float] | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> MultiresLandmarkSession:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls(**raw)

    def paired_points_zyx(self) -> tuple[np.ndarray, np.ndarray]:
        """Return matched overview and ROI points as Nx3 arrays in ZYX order."""
        if len(self.overview_points_zyx) != len(self.roi_points_zyx):
            msg = (
                "Overview and ROI landmark counts must match: "
                f"{len(self.overview_points_zyx)} vs {len(self.roi_points_zyx)}"
            )
            raise ValueError(msg)
        if not self.overview_points_zyx:
            return np.zeros((0, 3)), np.zeros((0, 3))
        overview = np.asarray(self.overview_points_zyx, dtype=float)
        roi = np.asarray(self.roi_points_zyx, dtype=float)
        if overview.shape[1] != 3 or roi.shape[1] != 3:
            msg = "Landmark points must be 3D coordinates in [Z, Y, X] order"
            raise ValueError(msg)
        return overview, roi

    def point_counts(self) -> tuple[int, int, int]:
        """Return (matched_pairs, n_overview, n_roi)."""
        n_overview = len(self.overview_points_zyx)
        n_roi = len(self.roi_points_zyx)
        return min(n_overview, n_roi), n_overview, n_roi


def default_landmark_session_path(save_path: Path, pair_label: str | None = None) -> Path:
    """Default landmark JSON path; include pair label when available to avoid collisions."""
    root = save_path.expanduser()
    if pair_label:
        from lightsuite.multires.registration import sanitize_experiment_name

        slug = sanitize_experiment_name(str(pair_label))
        return root / f"multires_landmarks_{slug}.json"
    return root / "multires_landmarks.json"
