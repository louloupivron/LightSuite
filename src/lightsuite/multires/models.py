"""Data models for multiresolution pair manifests (LightSuite Multires Pair v1)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

MANIFEST_FORMAT = "lightsuite_multires_pair_v1"


@dataclass
class ManifestVolumeSpec:
    """Physical frame for one overview or ROI volume."""

    volume_path: str
    shape_zyx: list[int]
    spacing_um: list[float]
    origin_um: list[float]
    direction: list[float] = field(
        default_factory=lambda: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    )

    def to_dict(self) -> dict[str, Any]:
        return {
            "volume_path": self.volume_path,
            "shape_zyx": self.shape_zyx,
            "spacing_um": self.spacing_um,
            "origin_um": self.origin_um,
            "direction": self.direction,
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> ManifestVolumeSpec:
        direction = raw.get("direction")
        if direction is None:
            direction = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
        return cls(
            volume_path=str(raw["volume_path"]),
            shape_zyx=[int(v) for v in raw["shape_zyx"]],
            spacing_um=[float(v) for v in raw["spacing_um"]],
            origin_um=[float(v) for v in raw["origin_um"]],
            direction=[float(v) for v in direction],
        )


@dataclass
class MultiresPairManifest:
    """Canonical contract between conversion notebooks and the multires pipeline."""

    format: str
    sample_name: str
    pair_label: str
    overview: ManifestVolumeSpec
    roi: ManifestVolumeSpec
    provenance: dict[str, str] = field(default_factory=dict)
    landmarks_path: str | None = None

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "format": self.format,
            "sample_name": self.sample_name,
            "pair_label": self.pair_label,
            "overview": self.overview.to_dict(),
            "roi": self.roi.to_dict(),
        }
        if self.provenance:
            out["provenance"] = self.provenance
        if self.landmarks_path is not None:
            out["landmarks_path"] = self.landmarks_path
        return out

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> MultiresPairManifest:
        fmt = raw.get("format")
        if fmt != MANIFEST_FORMAT:
            msg = f"Unsupported manifest format {fmt!r}; expected {MANIFEST_FORMAT!r}"
            raise ValueError(msg)
        landmarks = raw.get("landmarks_path")
        return cls(
            format=MANIFEST_FORMAT,
            sample_name=str(raw["sample_name"]),
            pair_label=str(raw["pair_label"]),
            overview=ManifestVolumeSpec.from_dict(raw["overview"]),
            roi=ManifestVolumeSpec.from_dict(raw["roi"]),
            provenance={str(k): str(v) for k, v in (raw.get("provenance") or {}).items()},
            landmarks_path=str(landmarks) if landmarks is not None else None,
        )


def serialize_report(report: dict[str, object]) -> dict[str, object]:
    out: dict[str, object] = {}
    for key, value in report.items():
        if isinstance(value, np.ndarray):
            out[key] = value.tolist()
        else:
            out[key] = value
    return out
