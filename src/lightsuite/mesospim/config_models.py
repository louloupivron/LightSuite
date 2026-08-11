"""Pydantic models for mesoSPIM geometry (used by multires vendor helpers)."""

from __future__ import annotations

from typing import Annotated, Literal

from pydantic import BaseModel, Field, field_validator


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
