"""VolumeReader auto-detection and registration."""

from __future__ import annotations

from pathlib import Path
from typing import Type

from lightsuite.config.models import SourceFormat, TiffLayout
from lightsuite.io.readers.tiff_stack import TiffStackReader

_READERS: list[Type] = [TiffStackReader]


def register_reader(reader_cls: Type) -> None:
    if reader_cls not in _READERS:
        _READERS.insert(0, reader_cls)


def open_volume(
    path: Path,
    source_format: SourceFormat = SourceFormat.AUTO,
    tiff_type: TiffLayout = TiffLayout.CHANNEL_PER_FILE,
    voxel_um: tuple[float, float, float] | None = None,
) -> TiffStackReader:
    """Open a volume using the first matching reader."""
    path = path.expanduser().resolve()

    if source_format in (SourceFormat.AUTO, SourceFormat.TIFF_STACK):
        if TiffStackReader.supports(path):
            return TiffStackReader.open(path, tiff_type=tiff_type, voxel_um=voxel_um)

    if source_format == SourceFormat.OME_ZARR:
        msg = "OME-Zarr reader not yet installed. Use: uv sync --extra formats"
        raise NotImplementedError(msg)

    if source_format == SourceFormat.IMARIS:
        msg = "Imaris reader not yet installed. Use: uv sync --extra formats"
        raise NotImplementedError(msg)

    if source_format == SourceFormat.CZI:
        msg = "CZI reader not yet implemented."
        raise NotImplementedError(msg)

    msg = f"No VolumeReader supports: {path} (format={source_format.value})"
    raise ValueError(msg)
