"""Brain lightsheet preprocessing (preprocessLightSheetVolume.m port)."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rich.console import Console
from scipy.ndimage import median_filter
from skimage.transform import resize

from lightsuite.config.models import BrainPipelineConfig
from lightsuite.io.discover import discover_tiff_stack
from lightsuite.io.readers.tiff_stack import TiffStackReader
from lightsuite.io.tiff_write import save_registration_volume
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint

console = Console()


@dataclass(frozen=True)
class PreprocessResult:
    checkpoint: RegOptsCheckpoint
    regvolpaths: dict[int, Path]


def _require_voxel_um(config: BrainPipelineConfig) -> tuple[float, float, float]:
    if config.sample.voxel_um is None:
        msg = "sample.voxel_um is required for preprocessing (e.g. [5.26, 5.26, 5.0])."
        raise ValueError(msg)
    vx, vy, vz = config.sample.voxel_um
    return float(vx), float(vy), float(vz)


def _background_fill(slice_2d: np.ndarray) -> np.ndarray:
    """Replace zero pixels with mode of positive samples (MATLAB preprocess step)."""
    data = slice_2d.astype(np.float64, copy=True)
    flat = data.ravel()
    sample_size = min(flat.size, 20_000)
    rng = np.random.default_rng(0)
    sample = rng.choice(flat, size=sample_size, replace=False) if flat.size > sample_size else flat
    positive = sample[sample > 0]
    if positive.size == 0:
        return slice_2d
    values, counts = np.unique(positive.astype(np.int64), return_counts=True)
    background = values[int(counts.argmax())]
    data[data == 0] = background
    return data.astype(slice_2d.dtype, copy=False)


def _resize_xy(slice_2d: np.ndarray, scale_xy: float) -> np.ndarray:
    if np.isclose(scale_xy, 1.0):
        return slice_2d.astype(np.uint16, copy=False)
    h, w = slice_2d.shape
    new_h = max(1, int(np.ceil(h * scale_xy)))
    new_w = max(1, int(np.ceil(w * scale_xy)))
    out = resize(
        slice_2d,
        (new_h, new_w),
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    )
    return np.clip(out, 0, np.iinfo(np.uint16).max).astype(np.uint16)


def _resize_z(volume_yxz: np.ndarray, scale_z: float) -> np.ndarray:
    if np.isclose(scale_z, 1.0):
        return volume_yxz.astype(np.uint16, copy=False)
    h, w, z = volume_yxz.shape
    new_z = max(1, int(np.ceil(z * scale_z)))
    out = resize(
        volume_yxz,
        (h, w, new_z),
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    )
    return np.clip(out, 0, np.iinfo(np.uint16).max).astype(np.uint16)


def _channel_for_cells(config: BrainPipelineConfig) -> int | None:
    if not config.detection.enabled:
        return None
    return config.detection.channel


def preprocess_lightsheet_volume(config: BrainPipelineConfig) -> PreprocessResult:
    """Downsample channels and write registration TIFFs under save_path."""
    vx, vy, vz = _require_voxel_um(config)
    registres = config.registration.resolution_um
    scale_xy = vx / registres
    scale_z = vz / registres

    config.sample.scratch.mkdir(parents=True, exist_ok=True)
    config.sample.save_path.mkdir(parents=True, exist_ok=True)

    discovery = discover_tiff_stack(
        config.sample.source.path,
        tiff_type=config.sample.source.tiff_type,
    )
    reader = TiffStackReader(discovery, voxel_um=(vx, vy, vz))

    ny, nx, nz, nchans = discovery.ny, discovery.nx, discovery.nz, discovery.nchans
    regvolpaths: dict[int, Path] = {}
    cell_channel = _channel_for_cells(config)

    out_h = max(1, int(np.ceil(ny * scale_xy)))
    out_w = max(1, int(np.ceil(nx * scale_xy)))

    for ichannel in range(1, nchans + 1):
        chan0 = ichannel - 1
        has_cells = cell_channel is not None and ichannel == cell_channel
        backvol = np.zeros((out_h, out_w, nz), dtype=np.uint16)

        binary_path: Path | None = None
        binary_handle = None
        if has_cells:
            binary_path = config.sample.scratch / f"chan_{ichannel}_binary_{config.sample.name}.dat"
            if binary_path.exists():
                binary_path.unlink()
            binary_handle = binary_path.open("wb")

        t0 = time.perf_counter()
        for islice in range(1, nz + 1):
            currim = reader.get_slice(islice, channel=chan0)
            currim = _background_fill(currim)
            backvol[:, :, islice - 1] = _resize_xy(currim, scale_xy)

            if binary_handle is not None:
                filtered = median_filter(currim.astype(np.uint16), size=3, mode="mirror")
                binary_handle.write(filtered.astype("<u2", copy=False).tobytes())

            if islice == 1 or islice % 20 == 0 or islice == nz:
                elapsed = time.perf_counter() - t0
                per_slice = elapsed / islice
                console.print(
                    f"Channel {ichannel}/{nchans}. Slice {islice}/{nz}. "
                    f"Time per slice {per_slice:.2f} s. Elapsed {elapsed:.2f} s..."
                )

        if binary_handle is not None:
            binary_handle.close()
            console.print(
                "[yellow]Cell detection not yet ported;[/yellow] "
                f"binary scratch written to {binary_path}"
            )
            if binary_path is not None and binary_path.exists():
                binary_path.unlink()

        console.print(f"Saving volume for registration (channel {ichannel})...", end=" ")
        voldown = _resize_z(backvol, scale_z)
        sample_path = (
            config.sample.save_path
            / f"chan_{ichannel}_sample_register_{int(registres)}um.tif"
        )
        if sample_path.exists():
            sample_path.unlink()
        save_registration_volume(voldown, sample_path)
        regvolpaths[ichannel] = sample_path
        console.print("Done.")

    reader.close()

    primary = config.registration.channel_primary
    secondary = config.registration.channel_secondary
    if primary < 1 or primary > nchans:
        msg = f"registration.channel_primary={primary} out of range 1..{nchans}"
        raise ValueError(msg)
    if secondary is not None:
        if secondary < 1 or secondary > nchans:
            msg = f"registration.channel_secondary={secondary} out of range 1..{nchans}"
            raise ValueError(msg)
        if secondary == primary:
            msg = "registration.channel_secondary must differ from channel_primary"
            raise ValueError(msg)

    checkpoint = RegOptsCheckpoint(
        sample_name=config.sample.name,
        ny=ny,
        nx=nx,
        nz=nz,
        nchans=nchans,
        voxel_um=[vx, vy, vz],
        registres_um=registres,
        regvolpath=str(regvolpaths[primary]),
        regvolpath_secondary=str(regvolpaths[secondary]) if secondary else None,
        regvolpaths={str(k): str(v) for k, v in regvolpaths.items()},
        tiff_type=config.sample.source.tiff_type.value,
        channel_primary=primary,
        channel_secondary=secondary,
    )
    regopts_path = config.sample.save_path / "regopts.json"
    checkpoint.save(regopts_path)
    console.print(f"Wrote checkpoint [bold]{regopts_path}[/bold]")

    return PreprocessResult(checkpoint=checkpoint, regvolpaths=regvolpaths)
