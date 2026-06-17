"""Pydantic models for mesoSPIM overview / ROI registration."""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from lightsuite.mesospim.landmark_session import LandmarkFitMode, default_landmark_session_path
from lightsuite.mesospim.meta import meta_path_for_tiff


class MesospimGeometryMode(StrEnum):
    METADATA = "metadata"
    LANDMARKS = "landmarks"
    HYBRID = "hybrid"


class MesospimSampleConfig(BaseModel):
    name: str = Field(min_length=1)
    save_path: Path
    scratch: Path | None = None

    @field_validator("save_path", "scratch")
    @classmethod
    def expand_paths(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser()


class MesospimVolumeConfig(BaseModel):
    path: Path
    meta_path: Path | None = None
    voxel_um: Annotated[list[float], Field(min_length=3, max_length=3)] | None = None

    @field_validator("path", "meta_path")
    @classmethod
    def expand_paths(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser().resolve()

    @field_validator("path")
    @classmethod
    def path_must_exist(cls, value: Path) -> Path:
        if not value.is_file():
            msg = f"Volume TIFF does not exist: {value}"
            raise ValueError(msg)
        return value

    @field_validator("voxel_um")
    @classmethod
    def voxel_um_positive(cls, value: list[float] | None) -> list[float] | None:
        if value is None:
            return None
        if any(v <= 0 for v in value):
            msg = "voxel_um values must be positive"
            raise ValueError(msg)
        return value

    def resolved_meta_path(self) -> Path:
        if self.meta_path is not None:
            return self.meta_path
        return meta_path_for_tiff(self.path).resolve()

    def has_meta_sidecar(self) -> bool:
        return self.resolved_meta_path().is_file()


class MesospimGeometryConfig(BaseModel):
    stage_xy_is_center: bool = True
    itk_lateral_dim0_motor: Literal["x", "y"] = "x"
    lateral_flip: Annotated[list[int], Field(min_length=2, max_length=2)] = [1, -1]

    @field_validator("lateral_flip")
    @classmethod
    def validate_lateral_flip(cls, value: list[int]) -> list[int]:
        if any(v not in (-1, 1) for v in value):
            msg = "lateral_flip values must be +1 or -1"
            raise ValueError(msg)
        return value


class MesospimLandmarkConfig(BaseModel):
    session_path: Path | None = None
    fit_mode: LandmarkFitMode = "similarity"
    min_pairs: int = Field(default=3, ge=3)

    @field_validator("session_path")
    @classmethod
    def expand_session_path(cls, value: Path | None) -> Path | None:
        if value is None:
            return None
        return value.expanduser()


class MesospimTiffRemapConfig(BaseModel):
    reverse_z: bool = False
    rot90_k_overview: int = 0
    rot90_k_roi: int = 0
    flip_row: bool = False
    flip_col: bool = False
    swap_xy: bool = False
    extra_flip_row_overview: bool = False
    extra_flip_col_overview: bool = False
    extra_flip_row_roi: bool = False
    extra_flip_col_roi: bool = False

    @field_validator("rot90_k_overview", "rot90_k_roi")
    @classmethod
    def normalize_rot90(cls, value: int) -> int:
        return int(value) % 4


class MesospimRegistrationSettings(BaseModel):
    overlap_margin_um: float = 0.0
    registration_bin: int = Field(default=1, ge=1)
    experiment_name: str = "default"
    elastix_stages: Annotated[list[str], Field(min_length=1)] = ["translation", "rigid"]
    write_full_overview_canvas: bool = True


class MesospimConfig(BaseModel):
    overview: MesospimVolumeConfig
    roi: MesospimVolumeConfig
    geometry_mode: MesospimGeometryMode = MesospimGeometryMode.METADATA
    geometry: MesospimGeometryConfig = Field(default_factory=MesospimGeometryConfig)
    landmarks: MesospimLandmarkConfig = Field(default_factory=MesospimLandmarkConfig)
    tiff_remap: MesospimTiffRemapConfig = Field(default_factory=MesospimTiffRemapConfig)
    registration: MesospimRegistrationSettings = Field(default_factory=MesospimRegistrationSettings)

    @model_validator(mode="after")
    def validate_geometry_requirements(self) -> MesospimConfig:
        mode = self.geometry_mode
        if mode == MesospimGeometryMode.METADATA:
            for label, volume in (("overview", self.overview), ("roi", self.roi)):
                if not volume.has_meta_sidecar():
                    msg = (
                        f"geometry_mode=metadata requires a meta sidecar for {label}: "
                        f"{volume.resolved_meta_path()}"
                    )
                    raise ValueError(msg)
        elif mode == MesospimGeometryMode.LANDMARKS:
            for label, volume in (("overview", self.overview), ("roi", self.roi)):
                if volume.voxel_um is None:
                    msg = (
                        f"geometry_mode=landmarks requires mesospim.{label}.voxel_um "
                        f"when metadata is unavailable"
                    )
                    raise ValueError(msg)
        elif mode == MesospimGeometryMode.HYBRID:
            for label, volume in (("overview", self.overview), ("roi", self.roi)):
                if not volume.has_meta_sidecar():
                    msg = (
                        f"geometry_mode=hybrid requires a meta sidecar for {label}: "
                        f"{volume.resolved_meta_path()}"
                    )
                    raise ValueError(msg)
        return self

    def resolved_landmark_session_path(self, save_path: Path) -> Path:
        if self.landmarks.session_path is not None:
            return self.landmarks.session_path.expanduser().resolve()
        return default_landmark_session_path(save_path)


class MesospimPipelineConfig(BaseModel):
    """Top-level mesoSPIM overview ↔ ROI registration configuration."""

    sample: MesospimSampleConfig
    mesospim: MesospimConfig
