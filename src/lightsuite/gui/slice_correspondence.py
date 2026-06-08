"""Sample-to-atlas slice correspondence along the primary anatomical axis."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class SliceAnchor:
    """One sample slice index paired with an atlas plane along the cut axis."""

    sample_index: int
    atlas_plane: int
    confirmed: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_index": int(self.sample_index),
            "atlas_plane": int(self.atlas_plane),
            "confirmed": bool(self.confirmed),
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> SliceAnchor:
        return cls(
            sample_index=int(raw["sample_index"]),
            atlas_plane=int(raw["atlas_plane"]),
            confirmed=bool(raw.get("confirmed", False)),
        )


@dataclass
class SliceCorrespondence:
    """AP (or primary-axis) slice map between sample and atlas volumes."""

    cut_axis: int
    original_trans: list[list[float]]
    anchors: list[SliceAnchor] = field(default_factory=list)
    version: int = 1
    source: str = "manual"

    def to_dict(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "source": self.source,
            "cut_axis": int(self.cut_axis),
            "original_trans": self.original_trans,
            "anchors": [anchor.to_dict() for anchor in self.anchors],
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> SliceCorrespondence:
        return cls(
            version=int(raw.get("version", 1)),
            source=str(raw.get("source", "manual")),
            cut_axis=int(raw["cut_axis"]),
            original_trans=raw["original_trans"],
            anchors=[SliceAnchor.from_dict(item) for item in raw.get("anchors", [])],
        )

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> SliceCorrespondence:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls.from_dict(raw)

    def confirmed_anchors(self) -> list[SliceAnchor]:
        return [anchor for anchor in self.anchors if anchor.confirmed]

    def interpolate_atlas_plane(
        self,
        sample_index: int,
        cut_axis: int,
        atlas_axis_size: int,
    ) -> int | None:
        """Return atlas plane for a sample index; None if cut axis differs or no anchors."""
        if int(cut_axis) != int(self.cut_axis):
            return None
        confirmed = self.confirmed_anchors()
        if not confirmed:
            return None
        sample_vals = np.array([anchor.sample_index for anchor in confirmed], dtype=float)
        atlas_vals = np.array([anchor.atlas_plane for anchor in confirmed], dtype=float)
        if sample_vals.size == 1:
            plane = float(atlas_vals[0])
        else:
            order = np.argsort(sample_vals)
            sample_vals = sample_vals[order]
            atlas_vals = atlas_vals[order]
            plane = float(np.interp(float(sample_index), sample_vals, atlas_vals))
        return int(np.clip(int(np.round(plane)), 1, atlas_axis_size))


def default_correspondence_path(save_path: Path) -> Path:
    return save_path.expanduser() / "slice_correspondence.json"
