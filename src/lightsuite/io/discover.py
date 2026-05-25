"""TIFF stack discovery (port of readLightsheetOpts.m)."""

from __future__ import annotations

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


def _list_tiffs(folder: Path) -> list[Path]:
    paths: list[Path] = []
    for pattern in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
        paths.extend(sorted(folder.glob(pattern)))
    return paths


def _tiff_channel_stack_dims(path: Path) -> tuple[int, int, int, bool, bool]:
    """Return ny, nx, nz, planes_in_time, prefer_native_pages."""
    with tifffile.TiffFile(path) as tif:
        page = tif.pages[0]
        ny, nx = page.shape[:2]
        n_pages = len(tif.pages)
        if n_pages > 1:
            return ny, nx, n_pages, False, True
    # Single IFD — depth unknown without OME; caller may use Bioformats later.
    return ny, nx, 1, False, False


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
        )

    # channelperfile
    if len(tfiles) == 1:
        path = tfiles[0]
        ny, nx, nz, planes_in_time, use_native = _tiff_channel_stack_dims(path)
        with tifffile.TiffFile(path) as tif:
            n_pages = len(tif.pages)
            if nz == 1 and n_pages > 1:
                nz = n_pages
                planes_in_time = True
                use_native = True
        # Estimate channel count from OME if present
        nchans = 1
        with tifffile.TiffFile(path) as tif:
            if tif.ome_metadata:
                import re

                match = re.search(r'SizeC="(\d+)"', tif.ome_metadata)
                if match:
                    nchans = int(match.group(1))
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
        )

    nchans = len(tfiles)
    dims = [_tiff_channel_stack_dims(p) for p in tfiles]
    first = dims[0]
    if not all(d[:3] == first[:3] for d in dims[1:]):
        msg = "Channel TIFF files have mismatched dimensions."
        raise ValueError(msg)

    planes_in_time = any(d[3] for d in dims)
    use_native = all(d[4] for d in dims)
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
    )
