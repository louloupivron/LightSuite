"""Pipeline checkpoint I/O (replaces regopts.mat)."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


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

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> RegOptsCheckpoint:
        raw = json.loads(path.expanduser().read_text(encoding="utf-8"))
        return cls(**raw)
