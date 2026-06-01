"""TIFF stack discovery (port of readLightsheetOpts.m)."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import tifffile

from lightsuite.config.models import TiffLayout


@dataclass(frozen=True)
class TiffStackDiscovery:
    """Discovered TIFF stack layout and dimensions."""

    tiff_type: TiffLayout
    tfiles: tuple[Path, ...]
    ny: int
    nx: int
    nz: int
    nchans: int
    multitiffs: bool
    planes_in_time: bool
    use_native_tiff_pages: bool
    stack_read_mode: str = "pages"


def _list_tiffs(folder: Path) -> list[Path]:
    paths: list[Path] = []
    for pattern in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
        paths.extend(sorted(folder.glob(pattern)))
    return paths


def _ome_size(pattern: str, ome_metadata: str) -> int | None:
    match = re.search(pattern, ome_metadata)
    if match:
        return int(match.group(1))
    return None


def _volume_layout_from_shape(shape: tuple[int, ...]) -> tuple[int, int, int, str]:
    """Infer ny, nx, nz from a 3D stack shape (Z often the smallest axis)."""
    if len(shape) != 3:
        msg = f"Expected 3D stack shape, got {shape}"
        raise ValueError(msg)
    z_dim, y_dim, x_dim = int(shape[0]), int(shape[1]), int(shape[2])
    if shape[0] <= shape[1] and shape[0] <= shape[2]:
        return y_dim, x_dim, z_dim, "memmap_zyx"
    if shape[2] <= shape[0] and shape[2] <= shape[1]:
        return z_dim, y_dim, x_dim, "memmap_yxz"
    return y_dim, x_dim, z_dim, "memmap_zyx"


def _single_tiff_stack_info(path: Path) -> tuple[int, int, int, bool, str]:
    """Return ny, nx, nz, use_native_pages, stack_read_mode for one stack file."""
    with tifffile.TiffFile(path) as tif:
        n_pages = len(tif.pages)
        page0 = tif.pages[0]
        arr0 = page0.asarray()

        if tif.series:
            shape = tif.series[0].shape
            if len(shape) == 3:
                ny, nx, nz, mode = _volume_layout_from_shape(shape)
                return ny, nx, nz, False, mode

        if arr0.ndim == 3:
            ny, nx, nz, mode = _volume_layout_from_shape(arr0.shape)
            return ny, nx, nz, False, mode

        ny, nx = int(arr0.shape[0]), int(arr0.shape[1])

        if n_pages > 1:
            return ny, nx, n_pages, True, "pages"

        ij = tif.imagej_metadata or {}
        for key in ("images", "slices"):
            if key in ij:
                nz = int(ij[key])
                if nz > 1:
                    return ny, nx, nz, False, "memmap_zyx"

        if tif.ome_metadata:
            nz = _ome_size(r'SizeZ="(\d+)"', tif.ome_metadata)
            if nz is not None and nz > 1:
                return ny, nx, nz, False, "memmap_zyx"

    return ny, nx, 1, False, "pages"


def _tiff_channel_stack_dims(path: Path) -> tuple[int, int, int, bool, bool, str]:
    """Return ny, nx, nz, planes_in_time, use_native_pages, stack_read_mode."""
    ny, nx, nz, use_native, mode = _single_tiff_stack_info(path)
    planes_in_time = use_native and nz > 1
    return ny, nx, nz, planes_in_time, use_native, mode


def discover_tiff_stack(
    data_folder: Path,
    tiff_type: TiffLayout = TiffLayout.CHANNEL_PER_FILE,
) -> TiffStackDiscovery:
    """Discover TIFF files and volume dimensions under data_folder."""
    folder = data_folder.expanduser().resolve()
    if not folder.is_dir():
        msg = f"Data folder not found: {folder}"
        raise FileNotFoundError(msg)

    tfiles = _list_tiffs(folder)
    if not tfiles:
        msg = f"No TIFF files found in {folder}"
        raise FileNotFoundError(msg)

    if tiff_type == TiffLayout.PLANE_PER_FILE:
        nz = len(tfiles)
        with tifffile.TiffFile(tfiles[0]) as tif:
            ny, nx = tif.pages[0].shape[:2]
            n_pages = len(tif.pages)
        if nz == 1 and n_pages > 1:
            msg = (
                f"planeperfile found one TIFF with {n_pages} pages in {folder}. "
                "Use tiff_type: channelperfile for a multi-page stack, or point "
                "source.path at a folder with one TIFF file per Z plane."
            )
            raise ValueError(msg)
        if nz == 1:
            msg = (
                f"planeperfile found only one TIFF in {folder}. "
                "Expected many plane files (e.g. z_0001.tif, z_0002.tif, ...). "
                "Check source.path, or use channelperfile for a single stack file."
            )
            raise ValueError(msg)
        return TiffStackDiscovery(
            tiff_type=tiff_type,
            tfiles=tuple(tfiles),
            ny=ny,
            nx=nx,
            nz=nz,
            nchans=1,
            multitiffs=False,
            planes_in_time=False,
            use_native_tiff_pages=True,
            stack_read_mode="pages",
        )

    # channelperfile
    if len(tfiles) == 1:
        path = tfiles[0]
        ny, nx, nz, planes_in_time, use_native, stack_read_mode = _tiff_channel_stack_dims(path)
        nchans = 1
        with tifffile.TiffFile(path) as tif:
            if tif.ome_metadata:
                nc = _ome_size(r'SizeC="(\d+)"', tif.ome_metadata)
                if nc is not None and nc >= 1:
                    nchans = nc
        return TiffStackDiscovery(
            tiff_type=tiff_type,
            tfiles=(path,),
            ny=ny,
            nx=nx,
            nz=nz,
            nchans=nchans,
            multitiffs=False,
            planes_in_time=planes_in_time,
            use_native_tiff_pages=use_native,
            stack_read_mode=stack_read_mode,
        )

    nchans = len(tfiles)
    dims = [_tiff_channel_stack_dims(p) for p in tfiles]
    first = dims[0]
    if not all(d[:3] == first[:3] for d in dims[1:]):
        msg = "Channel TIFF files have mismatched dimensions."
        raise ValueError(msg)

    planes_in_time = any(d[3] for d in dims)
    use_native = all(d[4] for d in dims)
    stack_read_mode = first[5]
    if not all(d[5] == stack_read_mode for d in dims):
        msg = "Channel TIFF files use different stack layouts."
        raise ValueError(msg)
    ny, nx, nz = first[0], first[1], first[2]

    return TiffStackDiscovery(
        tiff_type=tiff_type,
        tfiles=tuple(tfiles),
        ny=ny,
        nx=nx,
        nz=nz,
        nchans=nchans,
        multitiffs=True,
        planes_in_time=planes_in_time,
        use_native_tiff_pages=use_native,
        stack_read_mode=stack_read_mode,
    )
