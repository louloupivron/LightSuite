"""TIFF stack discovery (port of readLightsheetOpts.m)."""

from __future__ import annotations

import os
import re
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import tifffile

from lightsuite.config.models import TiffLayout
from lightsuite.reporter import emit_pipeline_message, format_duration


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
    channel_plane_files: tuple[tuple[Path, ...], ...] | None = None


_TIFF_SUFFIXES = {".tif", ".tiff"}
_SLOW_STORAGE_MARKERS = ("gvfs", "smb-share", "fuse.gvfs")
_LIST_HEARTBEAT_S = 2.0
_SLOW_LOG_S = 2.0


def _is_tiff_filename(name: str) -> bool:
    return Path(name).suffix.lower() in _TIFF_SUFFIXES


def looks_like_slow_storage(path: Path) -> bool:
    """True when *path* is (or points at) GVFS/SMB, which makes listing/opens slow."""
    try:
        target = os.readlink(path) if path.is_symlink() else str(path)
    except OSError:
        target = str(path)
    lowered = target.lower()
    return any(marker in lowered for marker in _SLOW_STORAGE_MARKERS)


def downsample_duration_hint(nz: int, nchans: int, sample_path: Path) -> str:
    """Coarse per-run downsample estimate from plane count and storage kind."""
    nplanes = max(nz, 1)
    nch = max(nchans, 1)
    if looks_like_slow_storage(sample_path):
        lo = format_duration(nplanes * 0.5)
        hi = format_duration(nplanes * 3.0)
        return (
            f"Downsampling {nplanes} planes × {nch} channel(s) over GVFS/SMB "
            f"(typically {lo}–{hi} per channel)."
        )
    return f"Downsampling {nplanes} planes × {nch} channel(s)."


def _list_tiffs(folder: Path, *, log_progress: bool = True) -> list[Path]:
    """List TIFF paths in *folder* with a single directory scan (no 4× glob)."""
    t0 = time.perf_counter()
    last_log = t0
    n_entries = 0
    paths: list[Path] = []
    announced = False
    with os.scandir(folder) as iterator:
        for entry in iterator:
            n_entries += 1
            if _is_tiff_filename(entry.name):
                paths.append(Path(entry.path))
            now = time.perf_counter()
            if log_progress and now - last_log >= _LIST_HEARTBEAT_S:
                if not announced:
                    emit_pipeline_message(f"Listing TIFFs in {folder} (slow directory)…")
                    announced = True
                elapsed = now - t0
                rate = n_entries / elapsed if elapsed > 0 else 0.0
                emit_pipeline_message(
                    f"  {n_entries} entries, {len(paths)} TIFFs, "
                    f"{format_duration(elapsed)} (~{rate:.0f}/s)"
                )
                last_log = now
    paths.sort()
    elapsed = time.perf_counter() - t0
    if log_progress and elapsed >= _SLOW_LOG_S:
        emit_pipeline_message(
            f"Listed {len(paths)} TIFF files in {format_duration(elapsed)}"
        )
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


def _volume_layout_from_labeled_axes(
    shape: tuple[int, ...],
    axes: str,
) -> tuple[int, int, int, str] | None:
    """Map a tifffile series ``axes`` label (e.g. ``ZYX``) to LightSuite Y/X/Z layout."""
    if len(shape) != 3 or len(axes) != 3:
        return None
    labels = axes.upper()
    if set(labels) != {"X", "Y", "Z"}:
        return None
    sizes = {label: int(shape[index]) for index, label in enumerate(labels)}
    ny, nx, nz = sizes["Y"], sizes["X"], sizes["Z"]
    z_axis = labels.index("Z")
    if z_axis == 0:
        mode = "memmap_zyx"
    elif z_axis == 1:
        mode = "memmap_xzy"
    elif z_axis == 2:
        mode = "memmap_yxz"
    else:
        return None
    return ny, nx, nz, mode


def _volume_layout_from_series(shape: tuple[int, ...], axes: str | None) -> tuple[int, int, int, str]:
    if axes:
        labeled = _volume_layout_from_labeled_axes(shape, axes)
        if labeled is not None:
            return labeled
    return _volume_layout_from_shape(shape)


def _single_tiff_stack_info(path: Path) -> tuple[int, int, int, bool, str]:
    """Return ny, nx, nz, use_native_pages, stack_read_mode for one stack file."""
    with tifffile.TiffFile(path) as tif:
        n_pages = len(tif.pages)
        page0 = tif.pages[0]
        arr0 = page0.asarray()

        # Prefer native TIFF pages when each IFD is one 2D plane (e.g. compressed
        # BigTIFF stacks). tifffile series may still report a 3D shape, but
        # memmap only works on uncompressed contiguous IFDs.
        if n_pages > 1 and arr0.ndim == 2:
            ny, nx = int(arr0.shape[0]), int(arr0.shape[1])
            return ny, nx, n_pages, True, "pages"

        if arr0.ndim == 3:
            ny, nx, nz, mode = _volume_layout_from_shape(arr0.shape)
            return ny, nx, nz, False, mode

        if tif.series:
            series = tif.series[0]
            shape = series.shape
            if len(shape) == 3:
                ny, nx, nz, mode = _volume_layout_from_series(shape, series.axes)
                return ny, nx, nz, False, mode

        ny, nx = int(arr0.shape[0]), int(arr0.shape[1])

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


_SMARTSPIM_CHANNEL_PLANE_RE = re.compile(
    r"^.+_Ch(?P<channel>\d+)\.(?:tif|tiff)$",
    re.IGNORECASE,
)


def _smartspim_interleaved_channel_indices(tfiles: Sequence[Path]) -> set[int]:
    channels: set[int] = set()
    for path in tfiles:
        match = _SMARTSPIM_CHANNEL_PLANE_RE.match(path.name)
        if match is not None:
            channels.add(int(match.group("channel")))
    return channels


def _raise_if_smartspim_interleaved_planes(folder: Path, tfiles: list[Path]) -> None:
    channels = _smartspim_interleaved_channel_indices(tfiles)
    if len(channels) < 2:
        return
    msg = (
        f"SmartSPIM interleaved planes detected in {folder} "
        f"(channels {sorted(channels)} in one folder). "
        "Run: lightsuite brain split-smartspim-channels -s <All_Channels> "
        "then set source.channels to the Ch0/Ch1 folders (planeperfile)."
    )
    raise ValueError(msg)


def _discover_planeperfile_folder(folder: Path) -> tuple[tuple[Path, ...], int, int, int]:
    """Return sorted plane TIFFs and ny, nx, nz for one planeperfile root."""
    tfiles = _list_tiffs(folder)
    if not tfiles:
        msg = f"No TIFF files found in {folder}"
        raise FileNotFoundError(msg)

    _raise_if_smartspim_interleaved_planes(folder, tfiles)

    nz = len(tfiles)
    first = tfiles[0]
    t_open = time.perf_counter()
    with tifffile.TiffFile(first) as tif:
        ny, nx = tif.pages[0].shape[:2]
        n_pages = len(tif.pages)
    open_s = time.perf_counter() - t_open
    if open_s >= _SLOW_LOG_S:
        emit_pipeline_message(
            f"Read plane header from {first.name} in {format_duration(open_s)}"
            + (" (GVFS/SMB)" if looks_like_slow_storage(first) else "")
        )
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
    return tuple(tfiles), ny, nx, nz


def discover_tiff_stack(
    data_folder: Path,
    tiff_type: TiffLayout = TiffLayout.CHANNEL_PER_FILE,
    *,
    channel_folders: Sequence[Path] | None = None,
) -> TiffStackDiscovery:
    """Discover TIFF files and volume dimensions under data_folder."""
    folder = data_folder.expanduser().resolve()
    if not folder.is_dir():
        msg = f"Data folder not found: {folder}"
        raise FileNotFoundError(msg)

    if tiff_type == TiffLayout.PLANE_PER_FILE and channel_folders:
        roots = tuple(p.expanduser().resolve() for p in channel_folders)
        channel_planes: list[tuple[Path, ...]] = []
        ny = nx = nz = 0
        for i, root in enumerate(roots):
            if not root.is_dir():
                msg = f"Channel folder not found: {root}"
                raise FileNotFoundError(msg)
            if len(roots) > 1:
                emit_pipeline_message(
                    f"Channel folder {i + 1}/{len(roots)}: {root}"
                )
            planes, cny, cnx, cnz = _discover_planeperfile_folder(root)
            channel_planes.append(planes)
            if i == 0:
                ny, nx, nz = cny, cnx, cnz
            elif (cny, cnx, cnz) != (ny, nx, nz):
                msg = (
                    f"Channel folders have mismatched dimensions: {root} is "
                    f"{cny}x{cnx}x{cnz}, expected {ny}x{nx}x{nz}."
                )
                raise ValueError(msg)
        first_planes = channel_planes[0]
        return TiffStackDiscovery(
            tiff_type=tiff_type,
            tfiles=first_planes,
            ny=ny,
            nx=nx,
            nz=nz,
            nchans=len(channel_planes),
            multitiffs=False,
            planes_in_time=False,
            use_native_tiff_pages=True,
            stack_read_mode="pages",
            channel_plane_files=tuple(channel_planes),
        )

    if tiff_type == TiffLayout.PLANE_PER_FILE:
        plane_files, ny, nx, nz = _discover_planeperfile_folder(folder)
        return TiffStackDiscovery(
            tiff_type=tiff_type,
            tfiles=plane_files,
            ny=ny,
            nx=nx,
            nz=nz,
            nchans=1,
            multitiffs=False,
            planes_in_time=False,
            use_native_tiff_pages=True,
            stack_read_mode="pages",
        )

    tfiles = _list_tiffs(folder)
    if not tfiles:
        msg = f"No TIFF files found in {folder}"
        raise FileNotFoundError(msg)

    if len(tfiles) > 1:
        emit_pipeline_message(f"Inspecting {len(tfiles)} TIFF(s) as channelperfile…")

    # channelperfile: one multi-page (or volumetric) TIFF per imaging channel.
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
    first_probe = _tiff_channel_stack_dims(tfiles[0])
    if first_probe[2] == 1 and nchans > 4:
        msg = (
            f"channelperfile found {nchans} TIFF files with 1 plane each in {folder}. "
            "Each file would be treated as a separate imaging channel "
            "(one slice per 'channel'). This folder is one TIFF per Z plane — "
            "set tiff_type: planeperfile."
        )
        raise ValueError(msg)
    dims = [first_probe, *(_tiff_channel_stack_dims(p) for p in tfiles[1:])]
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
