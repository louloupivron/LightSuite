"""VolumeReader protocol for lazy microscopy volume access."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable


@dataclass(frozen=True)
class VolumeMetadata:
    """Shape and spacing for a microscopy volume."""

    shape_yxzt: tuple[int, int, int, int]  # Ny, Nx, Nz, Nchans
    voxel_um: tuple[float, float, float] | None = None
    axes: str = "zyxc"


@runtime_checkable
class VolumeReader(Protocol):
    """Lazy, chunked access to microscopy volumes."""

    @classmethod
    def supports(cls, path: Path) -> bool: ...

    @property
    def metadata(self) -> VolumeMetadata: ...

    def get_slice(self, z: int, channel: int = 0) -> object: ...

    def close(self) -> None: ...
