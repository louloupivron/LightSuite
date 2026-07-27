"""Pydantic models for manifest-driven multiresolution registration."""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Annotated

from pydantic import BaseModel, Field, field_validator, model_validator

from lightsuite.multires.landmark_session import LandmarkFitMode, default_landmark_session_path


class MultiresGeometryMode(StrEnum):
    METADATA = "metadata"
    LANDMARKS = "landmarks"
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


class MultiresRegistrationSettings(BaseModel):
    overlap_margin_um: float = 0.0
    registration_bin: int = Field(default=1, ge=1)
    experiment_name: str = "default"
    elastix_stages: Annotated[list[str], Field(min_length=1)] = ["translation", "rigid"]
    write_full_overview_canvas: bool = True


class MultiresConfig(BaseModel):
    pair_manifest: Path
    geometry_mode: MultiresGeometryMode = MultiresGeometryMode.METADATA
    landmarks: MultiresLandmarkConfig = Field(default_factory=MultiresLandmarkConfig)
    registration: MultiresRegistrationSettings = Field(default_factory=MultiresRegistrationSettings)

    @field_validator("pair_manifest")
    @classmethod
    def expand_manifest_path(cls, value: Path) -> Path:
        return value.expanduser().resolve()

    @field_validator("pair_manifest")
    @classmethod
    def manifest_must_exist(cls, value: Path) -> Path:
        if not value.is_file():
            msg = f"Pair manifest does not exist: {value}"
            raise ValueError(msg)
        return value

    @model_validator(mode="after")
    def validate_landmark_mode(self) -> MultiresConfig:
        if self.geometry_mode in (MultiresGeometryMode.LANDMARKS, MultiresGeometryMode.HYBRID):
            return self
        return self

    def resolved_landmark_session_path(self, save_path: Path, manifest: object) -> Path:
        if self.landmarks.session_path is not None:
            return self.landmarks.session_path.expanduser().resolve()
        landmarks_path = getattr(manifest, "landmarks_path", None)
        if landmarks_path:
            path = Path(str(landmarks_path)).expanduser()
            if not path.is_absolute():
                path = self.pair_manifest.parent / path
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
