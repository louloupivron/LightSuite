"""Pydantic models for manifest-driven multiresolution registration."""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from lightsuite.config.models import ImportConfig
from lightsuite.multires.landmark_session import LandmarkFitMode, default_landmark_session_path


class MesospimGeometryOverride(BaseModel):
    """Partial mesoSPIM geometry fields merged onto defaults when building manifests."""

    stage_xy_is_center: bool | None = None
    itk_lateral_dim0_motor: Literal["x", "y"] | None = None
    lateral_flip: Annotated[list[int], Field(min_length=2, max_length=2)] | None = None

    @field_validator("lateral_flip")
    @classmethod
    def validate_lateral_flip(cls, value: list[int] | None) -> list[int] | None:
        if value is None:
            return None
        if any(v not in (-1, 1) for v in value):
            msg = "lateral_flip values must be +1 or -1"
            raise ValueError(msg)
        return value


class MultiresMesospimGeometryConfig(BaseModel):
    """Optional per-volume mesoSPIM geometry when building manifests from ``channels``."""

    overview: MesospimGeometryOverride | None = None
    roi: MesospimGeometryOverride | None = None


class MultiresGeometryMode(StrEnum):
    METADATA = "metadata"
    HYBRID = "hybrid"


class MultiresGeometryCheckLevel(StrEnum):
    """How much voxel data ``check-geometry`` loads."""

    METADATA_ONLY = "metadata-only"
    SLICE_QC = "slice-qc"
    FULL = "full"


class MultiresSampleConfig(BaseModel):
    name: str = Field(min_length=1)
    save_path: Path
    scratch: Path | None = None

    @field_validator("save_path", "scratch")
    @classmethod
    def expand_paths(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser()


class MultiresLandmarkConfig(BaseModel):
    session_path: Path | None = None
    fit_mode: LandmarkFitMode = "similarity"
    min_pairs: int = Field(default=3, ge=3)

    @field_validator("session_path")
    @classmethod
    def expand_session_path(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser()


class MultiresChannelPathConfig(BaseModel):
    """Overview / ROI volume paths for one imaging channel."""

    overview: Path
    roi: Path
    overview_meta_path: Path | None = None
    roi_meta_path: Path | None = None

    @field_validator("overview", "roi", "overview_meta_path", "roi_meta_path")
    @classmethod
    def expand_paths(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser().resolve()

    @field_validator("overview", "roi")
    @classmethod
    def path_must_exist(cls, value: Path) -> Path:
        if not value.exists():
            msg = f"Channel volume path does not exist: {value}"
            raise ValueError(msg)
        return value


class MultiresRegistrationSettings(BaseModel):
    overlap_margin_um: float = 0.0
    registration_bin: int = Field(default=1, ge=1)
    max_slab_bytes: int = Field(default=500_000_000, ge=50_000_000)
    experiment_name: str = "default"
    elastix_stages: Annotated[list[str], Field(min_length=1)] = ["translation", "rigid"]
    write_full_overview_canvas: bool = True
    reference_channel: str | None = None
    apply_transform_to: list[str] | None = None


class MultiresConfig(BaseModel):
    pair_manifest: Path | None = None
    pair_label: str | None = None
    overview_meta_path: Path | None = None
    channels: dict[str, MultiresChannelPathConfig] | None = None
    mesospim_geometry: MultiresMesospimGeometryConfig | None = None
    geometry_mode: MultiresGeometryMode = MultiresGeometryMode.METADATA
    landmarks: MultiresLandmarkConfig = Field(default_factory=MultiresLandmarkConfig)
    registration: MultiresRegistrationSettings = Field(default_factory=MultiresRegistrationSettings)

    @field_validator("pair_manifest", "overview_meta_path")
    @classmethod
    def expand_optional_paths(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser().resolve()

    @field_validator("geometry_mode", mode="before")
    @classmethod
    def normalize_legacy_landmarks_mode(cls, value: object) -> object:
        if value == "landmarks":
            return MultiresGeometryMode.HYBRID
        return value

    @model_validator(mode="after")
    def validate_manifest_or_channels(self) -> MultiresConfig:
        has_channels = bool(self.channels)
        has_manifest = self.pair_manifest is not None
        if not has_channels and not has_manifest:
            msg = "Provide multires.pair_manifest and/or multires.channels"
            raise ValueError(msg)
        if has_manifest and self.pair_manifest is not None and not has_channels:
            if not self.pair_manifest.is_file():
                msg = f"Pair manifest does not exist: {self.pair_manifest}"
                raise ValueError(msg)
        if has_channels:
            if self.registration.reference_channel is None:
                # Default to first declared channel when building from config paths.
                self.registration.reference_channel = next(iter(self.channels))
            ref = self.registration.reference_channel
            if ref not in self.channels:
                msg = (
                    f"registration.reference_channel {ref!r} missing from multires.channels "
                    f"(available: {sorted(self.channels)})"
                )
                raise ValueError(msg)
            if self.apply_targets_missing():
                missing = [
                    name
                    for name in (self.registration.apply_transform_to or [])
                    if name not in self.channels
                ]
                msg = f"registration.apply_transform_to channels missing from multires.channels: {missing}"
                raise ValueError(msg)
        return self

    def apply_targets_missing(self) -> bool:
        if not self.channels or not self.registration.apply_transform_to:
            return False
        return any(name not in self.channels for name in self.registration.apply_transform_to)

    def resolved_pair_label(self, sample_name: str) -> str:
        if self.pair_label:
            return self.pair_label
        if self.registration.experiment_name and self.registration.experiment_name != "default":
            return self.registration.experiment_name
        return f"{sample_name}_pair"

    def resolved_pair_manifest_path(self, save_path: Path, sample_name: str) -> Path:
        if self.pair_manifest is not None:
            return self.pair_manifest.expanduser().resolve()
        label = self.resolved_pair_label(sample_name)
        return (save_path.expanduser() / "converted" / f"{label}_pair.json").resolve()

    def resolved_landmark_session_path(self, save_path: Path, manifest: object) -> Path:
        if self.landmarks.session_path is not None:
            return self.landmarks.session_path.expanduser().resolve()
        landmarks_path = getattr(manifest, "landmarks_path", None)
        if landmarks_path:
            path = Path(str(landmarks_path)).expanduser()
            if not path.is_absolute():
                manifest_path = self.resolved_pair_manifest_path(
                    save_path,
                    getattr(manifest, "sample_name", "sample"),
                )
                path = manifest_path.parent / path
            return path.resolve()
        pair_label = getattr(manifest, "pair_label", None)
        return default_landmark_session_path(
            save_path,
            pair_label=str(pair_label) if pair_label else None,
        )


class MultiresPipelineConfig(BaseModel):
    """Top-level manifest-driven multiresolution registration configuration."""

    sample: MultiresSampleConfig
    multires: MultiresConfig
    import_config: ImportConfig | None = Field(default=None, alias="import")

    model_config = {"populate_by_name": True}
