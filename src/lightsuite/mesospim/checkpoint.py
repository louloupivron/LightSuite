"""Checkpoint I/O for mesoSPIM overview / ROI registration."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass
class MesospimRegOptsCheckpoint:
    """State persisted between mesoSPIM pipeline stages."""

    sample_name: str
    overview_path: str
    roi_path: str
    overview_meta_path: str
    roi_meta_path: str
    experiment_slug: str
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
    def load(cls, path: Path) -> MesospimRegOptsCheckpoint:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls(**raw)


def mesospim_checkpoint_path(save_path: Path) -> Path:
    return save_path.expanduser() / "mesospim_regopts.json"
