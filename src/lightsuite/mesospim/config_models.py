"""Pydantic models for mesoSPIM overview / ROI registration."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from lightsuite.mesospim.meta import meta_path_for_tiff


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

    def resolved_meta_path(self) -> Path:
        if self.meta_path is not None:
            return self.meta_path
        return meta_path_for_tiff(self.path).resolve()

    @model_validator(mode="after")
    def meta_must_exist(self) -> MesospimVolumeConfig:
        meta = self.resolved_meta_path()
        if not meta.is_file():
            msg = f"Meta sidecar does not exist: {meta}"
            raise ValueError(msg)
        return self


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


class MesospimConfig(BaseModel):
    overview: MesospimVolumeConfig
    roi: MesospimVolumeConfig
    geometry: MesospimGeometryConfig = Field(default_factory=MesospimGeometryConfig)
    tiff_remap: MesospimTiffRemapConfig = Field(default_factory=MesospimTiffRemapConfig)
    registration: MesospimRegistrationSettings = Field(default_factory=MesospimRegistrationSettings)


class MesospimPipelineConfig(BaseModel):
    """Top-level mesoSPIM overview ↔ ROI registration configuration."""

    sample: MesospimSampleConfig
    mesospim: MesospimConfig
