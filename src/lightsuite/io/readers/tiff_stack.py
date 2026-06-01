"""TIFF stack VolumeReader implementation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from lightsuite.config.models import TiffLayout
from lightsuite.io.discover import TiffStackDiscovery, discover_tiff_stack
from lightsuite.io.readers.base import VolumeMetadata
from lightsuite.io.tiff_volume import read_volumetric_slice


class TiffStackReader:
    """Read lightsheet TIFF stacks (planeperfile or channelperfile)."""

    def __init__(self, discovery: TiffStackDiscovery, voxel_um: tuple[float, float, float] | None = None):
        self._discovery = discovery
        self._voxel_um = voxel_um
        self._page_infos: dict[Path, tifffile.TiffFile] = {}

    @classmethod
    def supports(cls, path: Path) -> bool:
        path = path.expanduser()
        if path.is_file():
            return path.suffix.lower() in {".tif", ".tiff"}
        if path.is_dir():
            return bool(list(path.glob("*.tif")) or list(path.glob("*.tiff")))
        return False

    @classmethod
    def open(
        cls,
        path: Path,
        tiff_type: TiffLayout = TiffLayout.CHANNEL_PER_FILE,
        voxel_um: tuple[float, float, float] | None = None,
    ) -> TiffStackReader:
        folder = path if path.is_dir() else path.parent
        discovery = discover_tiff_stack(folder, tiff_type=tiff_type)
        return cls(discovery, voxel_um=voxel_um)

    @property
    def metadata(self) -> VolumeMetadata:
        d = self._discovery
        return VolumeMetadata(
            shape_yxzt=(d.ny, d.nx, d.nz, d.nchans),
            voxel_um=self._voxel_um,
        )

    @property
    def discovery(self) -> TiffStackDiscovery:
        return self._discovery

    def _open_tiff(self, path: Path) -> tifffile.TiffFile:
        handle = self._page_infos.get(path)
        if handle is None:
            handle = tifffile.TiffFile(path)
            self._page_infos[path] = handle
        return handle

    def _read_planeperfile(self, z: int, channel: int) -> np.ndarray:
        if channel != 0:
            msg = "planeperfile layout has a single channel."
            raise IndexError(msg)
        index = z - 1
        if index < 0 or index >= len(self._discovery.tfiles):
            raise IndexError(f"Slice index {z} out of range 1..{self._discovery.nz}")
        path = self._discovery.tfiles[index]
        with tifffile.TiffFile(path) as tif:
            return tif.pages[0].asarray()

    def _read_channelperfile(self, z: int, channel: int) -> np.ndarray:
        d = self._discovery
        chan_index = channel
        if chan_index < 0 or chan_index >= d.nchans:
            raise IndexError(f"Channel {channel} out of range 0..{d.nchans - 1}")
        z_index = z
        if z_index < 1 or z_index > d.nz:
            raise IndexError(f"Slice index {z} out of range 1..{d.nz}")

        if d.multitiffs:
            path = d.tfiles[chan_index]
            if d.use_native_tiff_pages:
                tif = self._open_tiff(path)
                return tif.pages[z_index - 1].asarray()
            return tifffile.imread(path)
        path = d.tfiles[0]
        if d.stack_read_mode != "pages":
            return read_volumetric_slice(path, z_index - 1, d.stack_read_mode)
        if d.use_native_tiff_pages:
            tif = self._open_tiff(path)
            return tif.pages[z_index - 1].asarray()
        return tifffile.imread(path, key=z_index - 1)

    def get_slice(self, z: int, channel: int = 0) -> np.ndarray:
        if self._discovery.tiff_type == TiffLayout.PLANE_PER_FILE:
            return self._read_planeperfile(z, channel)
        return self._read_channelperfile(z, channel)

    def close(self) -> None:
        for handle in self._page_infos.values():
            handle.close()
        self._page_infos.clear()
