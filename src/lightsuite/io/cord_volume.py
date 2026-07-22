"""Spinal cord volume I/O, layout detection, and streaming downsampling."""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
from skimage.transform import resize

from lightsuite.config.models import CordTiffLayout
from lightsuite.io.discover import _list_tiffs, _single_tiff_stack_info, discover_tiff_stack
from lightsuite.config.models import TiffLayout
from lightsuite.preprocess.slice_ops import read_source_plane, resize_xy_fast, SliceLoadJob

from rich.console import Console

console = Console()


def _print_slice_progress(current: int, total: int) -> None:
    if current == 1 or current % 100 == 0 or current == total:
        console.print(f"  read slice {current} / {total}")


@dataclass(frozen=True)
class SkippedSlice:
    slice_index: int
    file_name: str
    file_path: str
    message: str


def normalize_res_um(value: float | list[float] | np.ndarray) -> np.ndarray:
    """Expand scalar voxel size to [x, y, z] micrometers."""
    arr = np.atleast_1d(np.asarray(value, dtype=float)).ravel()
    if arr.size == 1:
        return np.full(3, arr[0], dtype=float)
    return arr[:3].astype(float)


def cord_volume_downsample_spec(
    native_size: tuple[int, int, int],
    sampleres_um: np.ndarray,
    registrationres_um: np.ndarray,
) -> tuple[np.ndarray, tuple[int, int, int]]:
    """Pixel scale factors and target shape for the registration grid."""
    native = np.asarray(native_size, dtype=float)
    resfac = normalize_res_um(sampleres_um) / normalize_res_um(registrationres_um)
    target = tuple(max(1, int(np.ceil(n * f))) for n, f in zip(native, resfac, strict=True))
    return resfac, target


def cord_downsample_volume(
    cordvol: np.ndarray,
    sampleres_um: np.ndarray,
    registrationres_um: np.ndarray,
    *,
    native_orisize: tuple[int, int, int] | None = None,
) -> np.ndarray:
    """Resample to registration grid; no-op when already at target size."""
    if native_orisize is None:
        native_orisize = tuple(int(v) for v in cordvol.shape[:3])
    _, target = cord_volume_downsample_spec(native_orisize, sampleres_um, registrationres_um)
    if tuple(cordvol.shape[:3]) == target:
        return cordvol
    nchan = cordvol.shape[3] if cordvol.ndim == 4 else 1
    if cordvol.ndim == 3:
        cordvol = cordvol[:, :, :, np.newaxis]
    out = np.zeros((*target, nchan), dtype=cordvol.dtype)
    for ich in range(nchan):
        out[:, :, :, ich] = resize(
            cordvol[:, :, :, ich],
            target,
            order=1,
            preserve_range=True,
            anti_aliasing=True,
        ).astype(cordvol.dtype)
    return out


def read_plane_tiff(path: Path) -> np.ndarray:
    """Read a single 2D plane TIFF (one page)."""
    path = path.expanduser()
    with tifffile.TiffFile(path) as tif:
        if len(tif.pages) != 1:
            msg = f"Expected single-page plane TIFF, got {len(tif.pages)} pages: {path}"
            raise ValueError(msg)
        img = tif.pages[0].asarray()
    if img.ndim > 2:
        img = img[:, :, 0]
    return np.asarray(img, dtype=np.uint16)


def _natural_sort_key(name: str) -> list:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", name)]


def _sorted_tiff_files(folder: Path) -> list[Path]:
    files = _list_tiffs(folder)
    return sorted(files, key=lambda p: _natural_sort_key(p.name))


def looks_like_plane_per_file(files: list[Path]) -> bool:
    """Heuristic: sample up to 8 files; each must be single-page."""
    if len(files) < 2:
        return False
    indices = np.unique(np.round(np.linspace(0, len(files) - 1, min(8, len(files)))).astype(int))
    for idx in indices:
        with tifffile.TiffFile(files[int(idx)]) as tif:
            if len(tif.pages) != 1:
                return False
    return True


def resolve_cord_tiff_layout(
    folder: Path,
    layout: CordTiffLayout,
    *,
    channel_folders: Sequence[Path] | None = None,
) -> CordTiffLayout:
    if channel_folders:
        if layout not in (CordTiffLayout.AUTO, CordTiffLayout.PLANE_PER_FILE):
            msg = "source.channels requires tiff_type planeperfile (or auto)"
            raise ValueError(msg)
        for root in channel_folders:
            files = _sorted_tiff_files(Path(root))
            if not files:
                msg = f"No .tif/.tiff files found in channel folder {root}"
                raise FileNotFoundError(msg)
        return CordTiffLayout.PLANE_PER_FILE

    files = _sorted_tiff_files(folder)
    if not files:
        msg = f"No .tif/.tiff files found in {folder}"
        raise FileNotFoundError(msg)
    if layout != CordTiffLayout.AUTO:
        return layout
    if len(files) == 1:
        return CordTiffLayout.CHANNEL_PER_FILE
    if looks_like_plane_per_file(files):
        return CordTiffLayout.PLANE_PER_FILE
    return CordTiffLayout.CHANNEL_PER_FILE


def filter_spinal_cord_slices(
    files: list[Path],
    *,
    skip_corrupt: bool = False,
) -> tuple[list[Path], list[SkippedSlice]]:
    """Validate plane TIFFs; optionally drop unreadable slices."""
    if not files:
        return [], []
    ny0, nx0 = read_plane_tiff(files[0]).shape
    bad: list[SkippedSlice] = []
    good: list[Path] = []
    for iz, path in enumerate(files, start=1):
        try:
            plane = read_plane_tiff(path)
            if plane.shape != (ny0, nx0):
                msg = f"size {plane.shape[0]}x{plane.shape[1]} px, expected {ny0}x{nx0} px"
                raise ValueError(msg)
            good.append(path)
        except (OSError, ValueError, tifffile.TiffFileError) as exc:
            bad.append(
                SkippedSlice(
                    slice_index=iz,
                    file_name=path.name,
                    file_path=str(path),
                    message=str(exc),
                )
            )
            if not skip_corrupt:
                lines = "\n".join(f"  [{s.slice_index}] {s.file_path}: {s.message}" for s in bad)
                msg = f"Corrupt or unreadable spinal cord slice TIFF(s):\n{lines}"
                raise RuntimeError(msg) from exc
    if not good:
        msg = "No readable slice TIFFs remain after filtering."
        raise RuntimeError(msg)
    return good, bad


def _assert_plane_stack_loadable(
    native_size: tuple[int, int, int],
    sampleres_um: np.ndarray,
    registrationres_um: np.ndarray,
) -> None:
    resfac, target = cord_volume_downsample_spec(native_size, sampleres_um, registrationres_um)
    native_bytes = float(np.prod(native_size)) * 2.0
    if native_bytes <= 200e9 and np.allclose(resfac, 1.0):
        return
    if np.allclose(resfac, 1.0):
        msg = (
            f"Native stack is {native_bytes / 1e9:.1f} GB ({native_size[0]} x {native_size[1]} x "
            f"{native_size[2]} px) but sample.voxel_um matches registration.resolution_um "
            f"(no downsampling). Set sample.voxel_um to your native microscope voxel size."
        )
        raise ValueError(msg)
    target_bytes = float(np.prod(target)) * 2.0
    if target_bytes > 200e9:
        msg = (
            f"Downsampled stack would still be {target_bytes / 1e9:.1f} GB ({target} px). "
            "Check sample.voxel_um and registration.resolution_um."
        )
        raise ValueError(msg)


def _plane_z_key(path: Path) -> int:
    """Parse a sortable plane index from common lightsheet / SmartSPIM TIFF names."""
    name = path.name
    for pattern in (
        r"_(\d+)_Ch\d+\.tif$",  # SmartSPIM: ..._000040_Ch0.tif
        r"_(\d+)\.tif$",  # Terastitcher / plane_003.tif
    ):
        match = re.search(pattern, name, re.IGNORECASE)
        if match:
            return int(match.group(1))
    stem_match = re.search(r"(\d+)(?:\D*)$", path.stem)
    if stem_match:
        return int(stem_match.group(1))
    msg = f"Could not parse plane index from TIFF filename: {name}"
    raise ValueError(msg)


def _channel_plane_key_map(
    folder: Path,
    *,
    skip_corrupt_slices: bool,
) -> tuple[dict[int, Path], list[SkippedSlice]]:
    files = _sorted_tiff_files(folder)
    if skip_corrupt_slices:
        files, skipped = filter_spinal_cord_slices(files, skip_corrupt=True)
    else:
        files, skipped = filter_spinal_cord_slices(files, skip_corrupt=False)
    key_map: dict[int, Path] = {}
    for path in files:
        key = _plane_z_key(path)
        if key in key_map:
            msg = (
                f"Duplicate plane index {key} in {folder}: "
                f"{key_map[key].name} and {path.name}"
            )
            raise ValueError(msg)
        key_map[key] = path
    return key_map, skipped


def align_multi_channel_plane_files(
    roots: Sequence[Path],
    *,
    skip_corrupt_slices: bool = False,
) -> tuple[list[list[Path]], list[SkippedSlice], dict[str, int]]:
    """Intersect plane indices across channel folders and return aligned file lists."""
    if len(roots) < 2:
        msg = "align_multi_channel_plane_files expects at least two channel folders."
        raise ValueError(msg)

    key_maps: list[dict[int, Path]] = []
    all_skipped: list[SkippedSlice] = []
    for root in roots:
        key_map, skipped = _channel_plane_key_map(root, skip_corrupt_slices=skip_corrupt_slices)
        key_maps.append(key_map)
        all_skipped.extend(skipped)

    common_keys = set.intersection(*(set(key_map) for key_map in key_maps))
    if not common_keys:
        msg = "No common plane indices across channel folders."
        raise ValueError(msg)

    sorted_keys = sorted(common_keys)
    aligned = [[key_map[key] for key in sorted_keys] for key_map in key_maps]
    dropped = {
        root.name: len(key_map) - len(common_keys)
        for root, key_map in zip(roots, key_maps, strict=True)
        if len(key_map) != len(common_keys)
    }
    return aligned, all_skipped, dropped


def _load_one_plane_per_file_channel(
    folder: Path,
    *,
    sampleres_um: np.ndarray,
    registrationres_um: np.ndarray,
    skip_corrupt_slices: bool = False,
    channel_label: str | None = None,
    max_slices: int | None = None,
    files: Sequence[Path] | None = None,
) -> tuple[np.ndarray, tuple[int, int, int], list[SkippedSlice]]:
    """Load one Terastitcher-style folder into a (Y, X, Z) uint16 volume."""
    skipped: list[SkippedSlice] = []
    if files is None:
        files = _sorted_tiff_files(folder)
        if not skip_corrupt_slices:
            files, skipped = filter_spinal_cord_slices(files, skip_corrupt=False)
        else:
            files, skipped = filter_spinal_cord_slices(files, skip_corrupt=True)
    else:
        files = list(files)

    if max_slices is not None and max_slices > 0:
        files = files[:max_slices]

    ny0, nx0 = read_plane_tiff(files[0]).shape
    nz_listed = len(files)
    native_listed = (ny0, nx0, nz_listed)
    resfac, target_size = cord_volume_downsample_spec(native_listed, sampleres_um, registrationres_um)
    _assert_plane_stack_loadable(native_listed, sampleres_um, registrationres_um)

    target_ny, target_nx, _target_nz = target_size
    scale_xy = float(resfac[0])
    label = channel_label or folder.name

    console.print(f"Using planeperfile loading for {label} ({len(files)} slice TIFFs)")
    if skipped:
        console.print(f"  filtered out {len(skipped)} corrupt slice(s) before loading")
        for row in skipped[:5]:
            console.print(f"  skipping corrupt slice {row.slice_index}: {row.file_name}")
        if len(skipped) > 5:
            console.print(f"    ... and {len(skipped) - 5} more")

    if target_size != native_listed:
        console.print(
            f"  downsampling while loading: {native_listed[0]} x {native_listed[1]} x {native_listed[2]} "
            f"-> {target_size[0]} x {target_size[1]} x {target_size[2]} px "
            f"(sampleres {sampleres_um.tolist()} um -> registration {registrationres_um.tolist()} um)"
        )

    if target_size == native_listed:
        vol = np.zeros((ny0, nx0, nz_listed), dtype=np.uint16)
        for iz, path in enumerate(files, start=1):
            vol[:, :, iz - 1] = read_plane_tiff(path)
            _print_slice_progress(iz, nz_listed)
        return vol, (ny0, nx0, nz_listed), skipped

    backvol = np.zeros((target_ny, target_nx, 0), dtype=np.uint16)
    for iz, path in enumerate(files, start=1):
        plane = read_plane_tiff(path)
        if np.isclose(scale_xy, 1.0):
            xy = plane
        else:
            xy = resize_xy_fast(plane, scale_xy)
        backvol = np.concatenate([backvol, xy[:, :, np.newaxis]], axis=2)
        _print_slice_progress(iz, len(files))

    if np.isclose(resfac[2], 1.0):
        finvol = backvol
    else:
        _, target_loaded = cord_volume_downsample_spec(
            (ny0, nx0, backvol.shape[2]),
            sampleres_um,
            registrationres_um,
        )
        finvol = resize(
            backvol,
            target_loaded,
            order=1,
            preserve_range=True,
            anti_aliasing=True,
        ).astype(np.uint16)

    return finvol, native_listed, skipped


def load_plane_per_file_stack(
    folder: Path,
    *,
    sampleres_um: np.ndarray,
    registrationres_um: np.ndarray,
    skip_corrupt_slices: bool = False,
    channel_folders: Sequence[Path] | None = None,
) -> tuple[np.ndarray, tuple[int, int, int], list[SkippedSlice], CordTiffLayout]:
    """Load Terastitcher-style slice series with optional streaming XY/Z downsampling.

    When ``channel_folders`` is set (brain-style multi-channel planeperfile), each
    folder is loaded and stacked along the channel axis.
    """
    roots = [Path(p) for p in channel_folders] if channel_folders else [folder]
    channel_vols: list[np.ndarray] = []
    native_orisize: tuple[int, int, int] | None = None
    all_skipped: list[SkippedSlice] = []

    aligned_files: list[list[Path]] | None = None
    if len(roots) > 1:
        aligned_files, all_skipped, dropped = align_multi_channel_plane_files(
            roots,
            skip_corrupt_slices=skip_corrupt_slices,
        )
        if dropped:
            common_nz = len(aligned_files[0])
            console.print(
                f"Multi-channel plane alignment: using [bold]{common_nz}[/bold] common plane(s)"
            )
            for root in roots:
                if dropped.get(root.name, 0):
                    console.print(
                        f"  {root.name}: dropped {dropped[root.name]} slice(s) without a match "
                        "in all channels"
                    )

    for ich, root in enumerate(roots, start=1):
        label = f"channel {ich}/{len(roots)} ({root.name})" if len(roots) > 1 else None
        channel_files = aligned_files[ich - 1] if aligned_files is not None else None
        vol, native, skipped = _load_one_plane_per_file_channel(
            root,
            sampleres_um=sampleres_um,
            registrationres_um=registrationres_um,
            skip_corrupt_slices=skip_corrupt_slices,
            channel_label=label,
            files=channel_files,
        )
        if native_orisize is None:
            native_orisize = native
        elif native != native_orisize:
            msg = (
                f"Channel folders have mismatched native dimensions: {root} is "
                f"{native[0]}x{native[1]}x{native[2]}, expected "
                f"{native_orisize[0]}x{native_orisize[1]}x{native_orisize[2]}."
            )
            raise ValueError(msg)
        if channel_vols and vol.shape != channel_vols[0].shape:
            msg = (
                f"Channel folders have mismatched loaded shapes: {root} is "
                f"{vol.shape}, expected {channel_vols[0].shape}."
            )
            raise ValueError(msg)
        channel_vols.append(vol)
        all_skipped.extend(skipped)

    assert native_orisize is not None
    volume = np.stack(channel_vols, axis=-1)
    return volume, native_orisize, all_skipped, CordTiffLayout.PLANE_PER_FILE


def load_channel_per_file_stack(
    folder: Path,
    *,
    sampleres_um: np.ndarray,
    registrationres_um: np.ndarray,
) -> tuple[np.ndarray, tuple[int, int, int], CordTiffLayout]:
    """Load channel-per-file or single multi-page stack."""
    discovery = discover_tiff_stack(folder, TiffLayout.CHANNEL_PER_FILE)
    if discovery.tiff_type == TiffLayout.PLANE_PER_FILE:
        msg = "Internal error: channel loader received plane-per-file discovery."
        raise RuntimeError(msg)

    channels: list[np.ndarray] = []
    native_sizes: list[tuple[int, int, int]] = []
    for chan_idx, path in enumerate(discovery.tfiles, start=1):
        console.print(f"  channel {chan_idx}/{len(discovery.tfiles)}: {path.name}")
        if discovery.multitiffs:
            with tifffile.TiffFile(path) as tif:
                if discovery.use_native_tiff_pages and discovery.stack_read_mode == "pages":
                    planes = [
                        read_source_plane(
                            SliceLoadJob(
                                source_path=str(path),
                                z_page=page,
                                scale_xy=1.0,
                                fill_background=False,
                                capture_binary=False,
                                stack_read_mode="pages",
                            )
                        )
                        for page in range(len(tif.pages))
                    ]
                    stack = np.stack(planes, axis=2).astype(np.uint16)
                else:
                    stack = tifffile.imread(path).astype(np.uint16)
                    if stack.ndim == 2:
                        stack = stack[:, :, np.newaxis]
        else:
            stack = tifffile.imread(path).astype(np.uint16)
            if stack.ndim == 2:
                stack = stack[:, :, np.newaxis]
        if stack.ndim != 3:
            msg = f"Expected 3D channel stack at {path}, got shape {stack.shape}"
            raise ValueError(msg)
        native_sizes.append(stack.shape)
        channels.append(stack)

    if len(set(native_sizes)) != 1:
        msg = "Channel TIFF stacks have mismatched dimensions."
        raise ValueError(msg)
    native_orisize = native_sizes[0]
    volume = np.stack(channels, axis=-1)
    _, target = cord_volume_downsample_spec(native_orisize, sampleres_um, registrationres_um)
    if target != native_orisize:
        console.print(
            f"  downsampling {native_orisize[0]} x {native_orisize[1]} x {native_orisize[2]} "
            f"-> {target[0]} x {target[1]} x {target[2]} px"
        )
    volume = cord_downsample_volume(volume, sampleres_um, registrationres_um, native_orisize=native_orisize)
    return volume, native_orisize, CordTiffLayout.CHANNEL_PER_FILE


def load_multichannel_single_stack(
    path: Path,
    *,
    sampleres_um: np.ndarray,
    registrationres_um: np.ndarray,
) -> tuple[np.ndarray, tuple[int, int, int], CordTiffLayout]:
    """Load a single TIFF using stack metadata (Bioformats-style layouts)."""
    ny, nx, nz, _, stack_read_mode = _single_tiff_stack_info(path)
    if stack_read_mode == "pages" and nz > 1:
        with tifffile.TiffFile(path) as tif:
            planes = [tif.pages[i].asarray().astype(np.uint16) for i in range(len(tif.pages))]
        stack = np.stack(planes, axis=2)
    else:
        stack = tifffile.imread(path).astype(np.uint16)
        if stack.ndim == 2:
            stack = stack[:, :, np.newaxis]
    if stack.ndim != 3:
        msg = f"Expected 3D stack in {path}, got {stack.shape}"
        raise ValueError(msg)
    native_orisize = stack.shape
    volume = stack[:, :, :, np.newaxis]
    _, target = cord_volume_downsample_spec(native_orisize, sampleres_um, registrationres_um)
    if target != native_orisize:
        console.print(
            f"  downsampling {native_orisize[0]} x {native_orisize[1]} x {native_orisize[2]} "
            f"-> {target[0]} x {target[1]} x {target[2]} px"
        )
    volume = cord_downsample_volume(volume, sampleres_um, registrationres_um, native_orisize=native_orisize)
    return volume, native_orisize, CordTiffLayout.MULTICHANNEL_SINGLE


def registration_volume_path(save_path: Path, channel: int, resolution_um: float) -> Path:
    label = int(round(resolution_um))
    return save_path / f"chan_{channel}_sample_register_{label}um.tif"


def load_registration_volumes(checkpoint) -> np.ndarray:
    """Load cached per-channel registration TIFFs from regopts."""
    paths = checkpoint.regvolpaths
    if not paths:
        msg = "regopts.json missing regvolpaths; re-run spinal preprocess."
        raise RuntimeError(msg)
    channels = sorted((int(k), Path(v)) for k, v in paths.items())
    first = tifffile.imread(channels[0][1])
    volume = np.zeros((*first.shape, len(channels)), dtype=np.uint16)
    volume[:, :, :, 0] = first
    for idx, (_, path) in enumerate(channels[1:], start=1):
        volume[:, :, :, idx] = tifffile.imread(path)
    return volume
