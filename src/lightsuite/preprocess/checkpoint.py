"""Pipeline checkpoint I/O (replaces regopts.mat)."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any

_DOWNSTREAM_CHECKPOINT_FIELDS = (
    "permute_sample_to_atlas",
    "original_trans",
    "downfac_reg",
    "autocpsample",
    "autocpatlas",
    "brain_atlas",
    "auto_points_source",
    "auto_points_refined",
    "auto_points_mode",
    "auto_points_correspondence_path",
)


def compute_preprocess_fingerprint(
    *,
    ny: int,
    nx: int,
    nz: int,
    nchans: int,
    voxel_um: list[float],
    registres_um: float,
    tiff_type: str,
    sample_content_crop: str = "off",
    sample_content_box: list[int] | None = None,
) -> dict[str, Any]:
    """Inputs that determine registration TIFF downsampling; YAML-only changes do not affect this."""
    fp = {
        "ny": int(ny),
        "nx": int(nx),
        "nz": int(nz),
        "nchans": int(nchans),
        "voxel_um": [float(v) for v in voxel_um],
        "registres_um": float(registres_um),
        "tiff_type": str(tiff_type),
        "sample_content_crop": str(sample_content_crop),
    }
    if sample_content_box is not None:
        fp["sample_content_box"] = [int(v) for v in sample_content_box]
    return fp


@dataclass
class RegOptsCheckpoint:
    """Registration/preprocess state persisted between pipeline stages."""

    sample_name: str
    ny: int
    nx: int
    nz: int
    nchans: int
    voxel_um: list[float]
    registres_um: float
    regvolpath: str
    regvolpath_secondary: str | None
    regvolpaths: dict[str, str]
    tiff_type: str
    channel_primary: int
    channel_secondary: int | None
    permute_sample_to_atlas: list[int] | None = None
    original_trans: list[list[float]] | None = None
    downfac_reg: float | None = None
    autocpsample: list[list[float]] | None = None
    autocpatlas: list[list[float]] | None = None
    brain_atlas: str | None = None
    auto_points_source: str | None = None
    auto_points_refined: bool = False
    auto_points_mode: str | None = None
    auto_points_correspondence_path: str | None = None
    preprocess_fingerprint: dict[str, Any] | None = None
    content_crop_start: list[int] | None = None
    content_crop_size: list[int] | None = None
    native_crop_offset_yxz: list[float] | None = None
    atlas_crop_start_native: list[int] | None = None
    atlas_native_shape: list[int] | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def merge_downstream_from(self, previous: RegOptsCheckpoint | None) -> RegOptsCheckpoint:
        """Keep init-registration / align-slices state when only refreshing preprocess metadata."""
        if previous is None:
            return self
        for name in _DOWNSTREAM_CHECKPOINT_FIELDS:
            setattr(self, name, getattr(previous, name))
        return self

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> RegOptsCheckpoint:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        known = {field.name for field in fields(cls)}
        filtered = {key: value for key, value in raw.items() if key in known}
        return cls(**filtered)
