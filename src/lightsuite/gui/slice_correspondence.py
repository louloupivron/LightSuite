"""Sample-to-atlas slice correspondence along volume axes (Y, X, Z)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

VOLUME_AXES = (1, 2, 3)


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


def _interpolate_anchor_curve(
    anchors: list[SliceAnchor],
    sample_index: int,
    atlas_axis_size: int,
    *,
    confirmed_only: bool = True,
    allow_extrapolation: bool = True,
) -> int | None:
    working = [anchor for anchor in anchors if anchor.confirmed] if confirmed_only else list(anchors)
    if not working:
        return None
    sample_vals = np.array([anchor.sample_index for anchor in working], dtype=float)
    atlas_vals = np.array([anchor.atlas_plane for anchor in working], dtype=float)
    if sample_vals.size == 1:
        if not allow_extrapolation and float(sample_index) != sample_vals[0]:
            return None
        plane = float(atlas_vals[0])
    else:
        order = np.argsort(sample_vals)
        sample_vals = sample_vals[order]
        atlas_vals = atlas_vals[order]
        if not allow_extrapolation and (
            float(sample_index) < sample_vals[0] or float(sample_index) > sample_vals[-1]
        ):
            return None
        plane = float(np.interp(float(sample_index), sample_vals, atlas_vals))
    return int(np.clip(int(np.round(plane)), 1, atlas_axis_size))


def resolve_correspondence_atlas_plane(
    correspondence: SliceCorrespondence,
    sample_index: int,
    cut_axis: int,
    atlas_axis_size: int,
) -> int | None:
    """Atlas plane for match-points from saved align-slices correspondence.

    Interpolates within the confirmed-anchor sample-index range. Outside that
    range, falls back to the full anchor curve (including auto-estimated
    unconfirmed anchors) so partial align-slices sessions still cover the volume.
    """
    anchors = correspondence.anchors_for_axis(cut_axis)
    plane = _interpolate_anchor_curve(
        anchors,
        sample_index,
        atlas_axis_size,
        confirmed_only=True,
        allow_extrapolation=False,
    )
    if plane is not None:
        return plane
    return _interpolate_anchor_curve(
        anchors,
        sample_index,
        atlas_axis_size,
        confirmed_only=False,
        allow_extrapolation=True,
    )


@dataclass
class SliceCorrespondence:
    """Per-axis sample index ↔ atlas plane maps between sample and atlas volumes."""

    original_trans: list[list[float]]
    axes: dict[int, list[SliceAnchor]] = field(default_factory=dict)
    version: int = 2
    source: str = "manual"

    @classmethod
    def single_axis(
        cls,
        cut_axis: int,
        original_trans: list[list[float]],
        anchors: list[SliceAnchor],
        *,
        version: int = 2,
        source: str = "manual",
    ) -> SliceCorrespondence:
        """Convenience builder for one axis (tests and v1 migration)."""
        return cls(
            original_trans=original_trans,
            axes={int(cut_axis): list(anchors)},
            version=version,
            source=source,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "source": self.source,
            "original_trans": self.original_trans,
            "axes": {
                str(axis): [anchor.to_dict() for anchor in self.axes[axis]]
                for axis in sorted(self.axes)
            },
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> SliceCorrespondence:
        version = int(raw.get("version", 1))
        source = str(raw.get("source", "manual"))
        original_trans = raw["original_trans"]
        if "axes" in raw:
            axes = {
                int(key): [SliceAnchor.from_dict(item) for item in value]
                for key, value in raw["axes"].items()
            }
        else:
            axes = {
                int(raw["cut_axis"]): [
                    SliceAnchor.from_dict(item) for item in raw.get("anchors", [])
                ]
            }
            version = max(version, 1)
        return cls(
            version=version if "axes" in raw else 2,
            source=source,
            original_trans=original_trans,
            axes=axes,
        )

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        self.version = 2
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> SliceCorrespondence:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls.from_dict(raw)

    @property
    def cut_axis(self) -> int:
        """First axis with anchors (backward compatibility for single-axis callers)."""
        if not self.axes:
            return 1
        return int(min(self.axes))

    @property
    def anchors(self) -> list[SliceAnchor]:
        """Anchors for :pyattr:`cut_axis` only (v1 compatibility)."""
        return list(self.axes.get(self.cut_axis, []))

    def anchors_for_axis(self, cut_axis: int) -> list[SliceAnchor]:
        return self.axes.setdefault(int(cut_axis), [])

    def confirmed_anchors(self, cut_axis: int | None = None) -> list[SliceAnchor]:
        if cut_axis is not None:
            return [anchor for anchor in self.anchors_for_axis(cut_axis) if anchor.confirmed]
        return [
            anchor
            for axis in sorted(self.axes)
            for anchor in self.axes[axis]
            if anchor.confirmed
        ]

    def has_confirmed_anchors(self, cut_axis: int | None = None) -> bool:
        return bool(self.confirmed_anchors(cut_axis))

    def confirmed_axis_count(self) -> int:
        return sum(1 for axis in VOLUME_AXES if self.has_confirmed_anchors(axis))

    def interpolate_atlas_plane(
        self,
        sample_index: int,
        cut_axis: int,
        atlas_axis_size: int,
        *,
        confirmed_only: bool = True,
        allow_extrapolation: bool = True,
    ) -> int | None:
        """Return atlas plane for a sample index along ``cut_axis``; None if no anchors."""
        return _interpolate_anchor_curve(
            self.anchors_for_axis(cut_axis),
            sample_index,
            atlas_axis_size,
            confirmed_only=confirmed_only,
            allow_extrapolation=allow_extrapolation,
        )


def default_correspondence_path(save_path: Path) -> Path:
    return save_path.expanduser() / "slice_correspondence.json"
