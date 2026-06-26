"""Spinal cord pipeline checkpoint I/O (replaces regopts.mat / spinal_alignment_opt.mat)."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any


@dataclass
class CordRegOptsCheckpoint:
    """Registration/preprocess state for spinal cord pipeline stages."""

    sample_name: str
    data_folder: str
    lsfolder: str
    orisize: list[int]
    nchans: int
    sampleres_um: list[float]
    registrationres_um: list[float]
    reg_channel: int
    sample_perm: list[int]
    tofliprc: bool
    ikeeprange: list[int]
    xrange: list[int]
    yrange: list[int]
    regvol_path: str
    tv_path: str
    av_path: str
    smpts_path: str
    tvpts_path: str
    atlas_res_um: list[float]
    segments_path: str
    regions_path: str
    tiff_type: str
    regvolpaths: dict[str, str] | None = None
    straightvol_path: str | None = None
    slicetforms_path: str | None = None
    affine_atlas_to_samp: list[list[float]] | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> CordRegOptsCheckpoint:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        known = {field.name for field in fields(cls)}
        filtered = {key: value for key, value in raw.items() if key in known}
        return cls(**filtered)


@dataclass
class SpinalAlignmentCheckpoint:
    """Straightening GUI output (spinal_alignment_opt.mat)."""

    user_cen: list[list[float | None]]
    user_ant: list[list[float | None]]
    user_pos: list[list[float | None]]
    fit_x: list[float]
    fit_y: list[float]
    fit_theta: list[float]
    lambda_pos: float
    lambda_ang: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> SpinalAlignmentCheckpoint:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        known = {field.name for field in fields(cls)}
        filtered = {key: value for key, value in raw.items() if key in known}
        return cls(**filtered)


@dataclass
class CordTransformParamsCheckpoint:
    """Post-registration transforms (transform_params.mat)."""

    tform_bspline_samp20um_to_atlas_20um_px: str
    tform_affine_samp20um_to_atlas_20um_px: list[list[float]]
    control_point_weight: float
    samp_ikeeplong: list[int]
    samp_ikeepx: list[int]
    samp_ikeepy: list[int]
    how_to_perm: list[int]
    slicetforms_path: str
    sampleres_um: list[float]
    registrationres_um: list[float]
    tofliprc: bool
    atlassize: list[int]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> CordTransformParamsCheckpoint:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        known = {field.name for field in fields(cls)}
        filtered = {key: value for key, value in raw.items() if key in known}
        return cls(**filtered)
