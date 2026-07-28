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
class MultiresChannelSpecs:
    """Overview and ROI volume specs for one imaging channel."""

    overview: ManifestVolumeSpec
    roi: ManifestVolumeSpec

    def to_dict(self) -> dict[str, Any]:
        return {
            "overview": self.overview.to_dict(),
            "roi": self.roi.to_dict(),
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> MultiresChannelSpecs:
        return cls(
            overview=ManifestVolumeSpec.from_dict(raw["overview"]),
            roi=ManifestVolumeSpec.from_dict(raw["roi"]),
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
    reference_channel: str | None = None
    channels: dict[str, MultiresChannelSpecs] | None = None

    def channel_names(self) -> list[str]:
        if self.channels:
            return list(self.channels.keys())
        return []

    def resolved_reference_channel(self) -> str | None:
        if self.reference_channel is not None:
            return self.reference_channel
        if self.channels:
            return next(iter(self.channels))
        return None

    def channel_specs(self, channel: str) -> MultiresChannelSpecs:
        if self.channels is None:
            msg = "Manifest has no multichannel entries"
            raise KeyError(msg)
        if channel not in self.channels:
            msg = f"Unknown channel {channel!r}; available: {sorted(self.channels)}"
            raise KeyError(msg)
        return self.channels[channel]

    def non_reference_channels(
        self,
        *,
        reference_channel: str | None = None,
        apply_transform_to: list[str] | None = None,
    ) -> list[str]:
        if not self.channels:
            return []
        ref = reference_channel or self.resolved_reference_channel()
        if ref is None:
            return []
        if apply_transform_to is not None:
            return [name for name in apply_transform_to if name != ref]
        return [name for name in self.channels if name != ref]

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
        if self.reference_channel is not None:
            out["reference_channel"] = self.reference_channel
        if self.channels:
            out["channels"] = {name: specs.to_dict() for name, specs in self.channels.items()}
        return out

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> MultiresPairManifest:
        fmt = raw.get("format")
        if fmt != MANIFEST_FORMAT:
            msg = f"Unsupported manifest format {fmt!r}; expected {MANIFEST_FORMAT!r}"
            raise ValueError(msg)
        landmarks = raw.get("landmarks_path")
        channels_raw = raw.get("channels")
        channels = (
            {str(name): MultiresChannelSpecs.from_dict(specs) for name, specs in channels_raw.items()}
            if channels_raw
            else None
        )
        reference_channel = raw.get("reference_channel")
        if reference_channel is not None:
            reference_channel = str(reference_channel)
        overview = ManifestVolumeSpec.from_dict(raw["overview"])
        roi = ManifestVolumeSpec.from_dict(raw["roi"])
        if channels and reference_channel and reference_channel in channels:
            overview = channels[reference_channel].overview
            roi = channels[reference_channel].roi
        return cls(
            format=MANIFEST_FORMAT,
            sample_name=str(raw["sample_name"]),
            pair_label=str(raw["pair_label"]),
            overview=overview,
            roi=roi,
            provenance={str(k): str(v) for k, v in (raw.get("provenance") or {}).items()},
            landmarks_path=str(landmarks) if landmarks is not None else None,
            reference_channel=reference_channel,
            channels=channels,
        )


def serialize_report(report: dict[str, object]) -> dict[str, object]:
    out: dict[str, object] = {}
    for key, value in report.items():
        if isinstance(value, np.ndarray):
            out[key] = value.tolist()
        else:
            out[key] = value
    return out
