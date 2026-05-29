"""Fast slice IO and downsampling for brain preprocess."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile
from scipy.ndimage import median_filter, zoom


@dataclass(frozen=True)
class SliceLoadJob:
    """Pickleable work item for one source Z plane."""

    source_path: str
    z_page: int | None  # None = whole file (planeperfile); else 0-based page index
    scale_xy: float
    fill_background: bool
    capture_binary: bool


@dataclass(frozen=True)
class SliceProcessResult:
    plane_xy: np.ndarray
    binary_bytes: bytes | None


def output_xy_shape(ny: int, nx: int, scale_xy: float) -> tuple[int, int]:
    out_h = max(1, int(np.ceil(ny * scale_xy)))
    out_w = max(1, int(np.ceil(nx * scale_xy)))
    return out_h, out_w


def background_fill_fast(slice_2d: np.ndarray) -> np.ndarray:
    """Replace zero pixels with mode of positive samples (MATLAB preprocess step)."""
    data = np.asarray(slice_2d, dtype=np.uint16)
    flat = data.ravel()
    probe = 50_000 if flat.size > 50_000 else flat.size
    if probe > 0:
        rng = np.random.default_rng(0)
        idx = rng.choice(flat.size, size=probe, replace=False) if flat.size > probe else np.arange(flat.size)
        if not np.any(flat[idx] == 0):
            return data
    zeros = data == 0
    if not np.any(zeros):
        return data
    sample_size = min(flat.size, 20_000)
    rng = np.random.default_rng(0)
    if flat.size > sample_size:
        idx = rng.choice(flat.size, size=sample_size, replace=False)
        sample = flat[idx]
    else:
        sample = flat
    positive = sample[sample > 0]
    if positive.size == 0:
        return data
    values, counts = np.unique(positive.astype(np.int64), return_counts=True)
    background = int(values[int(counts.argmax())])
    out = data.copy()
    out[zeros] = background
    return out


def resize_xy_fast(slice_2d: np.ndarray, scale_xy: float) -> np.ndarray:
    """Downsample one XY plane using striding (crop) with a small zoom fallback."""
    if np.isclose(scale_xy, 1.0):
        return np.asarray(slice_2d, dtype=np.uint16)

    h, w = slice_2d.shape
    new_h = max(1, int(np.ceil(h * scale_xy)))
    new_w = max(1, int(np.ceil(w * scale_xy)))

    source = np.asarray(slice_2d, dtype=np.uint16)
    stride = max(1, min(int(h // new_h), int(w // new_w)))
    if stride > 1:
        source = source[::stride, ::stride]

    sh, sw = source.shape
    if sh >= new_h and sw >= new_w:
        return source[:new_h, :new_w].copy()

    if sh == new_h and sw == new_w:
        return source

    zy = new_h / sh
    zx = new_w / sw
    out = zoom(source, (zy, zx), order=1, prefilter=False)
    return np.clip(out, 0, np.iinfo(np.uint16).max).astype(np.uint16)


def read_source_plane(job: SliceLoadJob) -> np.ndarray:
    if job.z_page is None:
        return tifffile.imread(job.source_path)
    return tifffile.imread(job.source_path, key=job.z_page)


def process_slice_job(job: SliceLoadJob) -> SliceProcessResult:
    currim = read_source_plane(job)
    if job.fill_background:
        currim = background_fill_fast(currim)
    plane = resize_xy_fast(currim, job.scale_xy)
    binary_bytes = None
    if job.capture_binary:
        filtered = median_filter(np.asarray(currim, dtype=np.uint16), size=3, mode="mirror")
        binary_bytes = filtered.astype("<u2", copy=False).tobytes()
    return SliceProcessResult(plane_xy=plane, binary_bytes=binary_bytes)


def output_z_count(nz: int, scale_z: float) -> int:
    return max(1, int(np.ceil(nz * scale_z)))


def z_downsample_plane(
    vol: np.ndarray | np.memmap,
    out_index: int,
    scale_z: float,
) -> np.ndarray:
    """Linear Z interpolation between nearest native Z planes."""
    _, _, z = vol.shape
    if z == 1:
        return np.asarray(vol[:, :, 0], dtype=np.uint16)

    new_z = output_z_count(z, scale_z)
    if new_z == 1:
        return np.asarray(vol[:, :, 0], dtype=np.uint16)

    src = out_index / scale_z
    src = min(max(src, 0.0), float(z - 1))
    z0 = int(np.floor(src))
    z1 = min(z0 + 1, z - 1)
    t = float(src - z0)
    if z0 == z1 or t <= 0.0:
        return np.asarray(vol[:, :, z0], dtype=np.uint16)
    a = vol[:, :, z0].astype(np.float32, copy=False)
    b = vol[:, :, z1].astype(np.float32, copy=False)
    blended = (1.0 - t) * a + t * b
    return np.clip(blended, 0, np.iinfo(np.uint16).max).astype(np.uint16)


def write_z_downsampled_volume(
    vol: np.ndarray | np.memmap,
    path: Path,
    scale_z: float,
    *,
    compression: str = "lzw",
) -> int:
    """Stream Z-downsampled pages to disk without holding the full 3D output in RAM."""
    _, _, z = vol.shape
    new_z = output_z_count(z, scale_z)
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    with tifffile.TiffWriter(path) as writer:
        for k in range(new_z):
            page = z_downsample_plane(vol, k, scale_z)
            writer.write(
                page,
                compression=compression,
                photometric="minisblack",
            )
    return new_z
