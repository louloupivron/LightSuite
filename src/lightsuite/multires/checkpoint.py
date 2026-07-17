"""Checkpoint I/O for manifest-driven multiresolution registration."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass
class MultiresRegOptsCheckpoint:
    """State persisted between multires pipeline stages."""

    sample_name: str
    pair_label: str
    pair_manifest_path: str
    experiment_slug: str
    overview_volume_path: str
    roi_volume_path: str
    overlap_box_um: list[list[float]] | None = None
    elastix_output_dir: str | None = None
    transform_paths: list[str] | None = None
    cropped_overview_path: str | None = None
    registered_roi_path: str | None = None
    registered_roi_full_overview_path: str | None = None
    crop_start_index: list[int] | None = None
    geometry_report_paths: dict[str, str] | None = None
    geometry_mode: str | None = None
    landmark_session_path: str | None = None
    roi_to_overview_tform: list[list[float]] | None = None
    landmark_rms_error_um: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> MultiresRegOptsCheckpoint:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls(**raw)


def multires_checkpoint_path(save_path: Path) -> Path:
    return save_path.expanduser() / "multires_regopts.json"
