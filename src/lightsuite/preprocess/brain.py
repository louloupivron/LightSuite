"""Brain lightsheet preprocessing (preprocessLightSheetVolume.m port)."""

from __future__ import annotations

import time
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.config.models import BrainPipelineConfig, TiffLayout
from lightsuite.io.discover import TiffStackDiscovery, discover_tiff_stack
from lightsuite.io.readers.tiff_stack import TiffStackReader
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.preprocess.slice_ops import (
    SliceLoadJob,
    SliceProcessResult,
    output_xy_shape,
    process_slice_job,
    write_z_downsampled_volume,
)

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


def _channel_for_cells(config: BrainPipelineConfig) -> int | None:
    if not config.detection.enabled:
        return None
    return config.detection.channel


def _slice_jobs_for_channel(
    discovery: TiffStackDiscovery,
    chan0: int,
    scale_xy: float,
    *,
    fill_background: bool,
    capture_binary: bool,
) -> list[SliceLoadJob]:
    jobs: list[SliceLoadJob] = []
    if discovery.tiff_type == TiffLayout.PLANE_PER_FILE:
        for path in discovery.tfiles:
            jobs.append(
                SliceLoadJob(
                    source_path=str(path),
                    z_page=None,
                    scale_xy=scale_xy,
                    fill_background=fill_background,
                    capture_binary=capture_binary,
                    stack_read_mode=discovery.stack_read_mode,
                )
            )
        return jobs

    if discovery.multitiffs:
        path = str(discovery.tfiles[chan0])
        for z_page in range(discovery.nz):
            jobs.append(
                SliceLoadJob(
                    source_path=path,
                    z_page=z_page,
                    scale_xy=scale_xy,
                    fill_background=fill_background,
                    capture_binary=capture_binary,
                    stack_read_mode=discovery.stack_read_mode,
                )
            )
        return jobs

    path = str(discovery.tfiles[0])
    for z_page in range(discovery.nz):
        jobs.append(
            SliceLoadJob(
                source_path=path,
                z_page=z_page,
                scale_xy=scale_xy,
                fill_background=fill_background,
                capture_binary=capture_binary,
                stack_read_mode=discovery.stack_read_mode,
            )
        )
    return jobs


def _effective_preprocess_workers(
    discovery: TiffStackDiscovery,
    requested: int,
) -> int:
    """Many small TIFFs (planeperfile) are faster sequentially — parallel reads thrash disk."""
    if discovery.tiff_type != TiffLayout.PLANE_PER_FILE:
        return max(1, requested)
    if requested <= 1:
        return 1
    console.print(
        "[yellow]planeperfile: using 1 worker "
        f"(requested {requested}; parallel reads slow on many TIFF files).[/yellow]"
    )
    return 1


def _iter_processed_slices(
    jobs: list[SliceLoadJob],
    workers: int,
) -> Iterator[SliceProcessResult]:
    if workers <= 1:
        for job in jobs:
            yield process_slice_job(job)
        return

    with ProcessPoolExecutor(max_workers=workers) as pool:
        yield from pool.map(process_slice_job, jobs, chunksize=1)


def _allocate_xy_stack(
    *,
    scratch_dir: Path,
    sample_name: str,
    channel: int,
    shape: tuple[int, int, int],
    max_in_memory_bytes: int,
) -> tuple[np.ndarray | np.memmap, Path | None, bool]:
    """Return (volume, optional memmap path, uses_disk)."""
    nbytes = int(np.prod(shape)) * np.dtype(np.uint16).itemsize
    if nbytes <= max_in_memory_bytes:
        return np.zeros(shape, dtype=np.uint16), None, False

    mmap_path = scratch_dir / f"chan_{channel}_xy_{sample_name}.dat"
    mmap_path.parent.mkdir(parents=True, exist_ok=True)
    if mmap_path.exists():
        mmap_path.unlink()
    return (
        np.memmap(mmap_path, dtype=np.uint16, mode="w+", shape=shape),
        mmap_path,
        True,
    )


def _process_channel_to_registration_tiff(
    *,
    jobs: list[SliceLoadJob],
    scratch_dir: Path,
    sample_name: str,
    channel: int,
    out_h: int,
    out_w: int,
    nz: int,
    scale_z: float,
    workers: int,
    output_path: Path,
    binary_path: Path | None,
    max_in_memory_bytes: int,
) -> None:
    xy_stack, mmap_path, on_disk = _allocate_xy_stack(
        scratch_dir=scratch_dir,
        sample_name=sample_name,
        channel=channel,
        shape=(out_h, out_w, nz),
        max_in_memory_bytes=max_in_memory_bytes,
    )

    binary_handle = None
    if binary_path is not None:
        if binary_path.exists():
            binary_path.unlink()
        binary_handle = binary_path.open("wb")

    t0 = time.perf_counter()
    try:
        for islice, result in enumerate(_iter_processed_slices(jobs, workers), start=1):
            xy_stack[:, :, islice - 1] = result.plane_xy
            if binary_handle is not None and result.binary_bytes is not None:
                binary_handle.write(result.binary_bytes)
            if islice == 1 or islice % 20 == 0 or islice == nz:
                elapsed = time.perf_counter() - t0
                console.print(
                    f"Slice {islice}/{nz}. "
                    f"Time per slice {elapsed / islice:.2f} s. Elapsed {elapsed:.2f} s..."
                )
    finally:
        if on_disk:
            del xy_stack
        if binary_handle is not None:
            binary_handle.close()

    if binary_path is not None:
        console.print(
            "[yellow]Cell detection not yet ported;[/yellow] "
            f"binary scratch written to {binary_path}"
        )
        if binary_path.exists():
            binary_path.unlink()

    console.print(f"Saving volume for registration (channel {channel})...", end=" ")
    try:
        if on_disk and mmap_path is not None:
            xy_read = np.memmap(mmap_path, dtype=np.uint16, mode="r", shape=(out_h, out_w, nz))
            write_z_downsampled_volume(xy_read, output_path, scale_z)
            del xy_read
        else:
            write_z_downsampled_volume(xy_stack, output_path, scale_z)
    finally:
        if on_disk and mmap_path is not None and mmap_path.exists():
            mmap_path.unlink()
    console.print("Done.")


def preprocess_lightsheet_volume(config: BrainPipelineConfig) -> PreprocessResult:
    """Downsample channels and write registration TIFFs under save_path."""
    vx, vy, vz = _require_voxel_um(config)
    registres = config.registration.resolution_um
    scale_xy = vx / registres
    scale_z = vz / registres
    requested_workers = config.compute.workers

    config.sample.scratch.mkdir(parents=True, exist_ok=True)
    config.sample.save_path.mkdir(parents=True, exist_ok=True)

    discovery = discover_tiff_stack(
        config.sample.source.path,
        tiff_type=config.sample.source.tiff_type,
    )
    reader = TiffStackReader(discovery, voxel_um=(vx, vy, vz))
    workers = _effective_preprocess_workers(discovery, requested_workers)

    ny, nx, nz, nchans = discovery.ny, discovery.nx, discovery.nz, discovery.nchans
    regvolpaths: dict[int, Path] = {}
    cell_channel = _channel_for_cells(config)
    out_h, out_w = output_xy_shape(ny, nx, scale_xy)

    scratch_bytes = out_h * out_w * nz * 2
    max_ram_bytes = int(config.compute.max_in_memory_scratch_gb * (1024**3))
    scratch_in_ram = scratch_bytes <= max_ram_bytes
    if discovery.stack_read_mode != "pages":
        console.print(
            f"Single-file volumetric TIFF: {nz} Z planes "
            f"({discovery.stack_read_mode}, channelperfile)."
        )
    elif discovery.tiff_type == TiffLayout.PLANE_PER_FILE and nchans == 1:
        est_xy_gb = scratch_bytes / (1024**3)
        where = "RAM" if scratch_in_ram else f"disk memmap on {config.sample.scratch}"
        console.print(
            f"planeperfile: {nz} planes, XY scratch ~{est_xy_gb:.1f} GB in {where} "
            f"(workers={workers})."
        )

    for ichannel in range(1, nchans + 1):
        chan0 = ichannel - 1
        has_cells = cell_channel is not None and ichannel == cell_channel
        console.print(f"Channel {ichannel}/{nchans}.")

        jobs = _slice_jobs_for_channel(
            discovery,
            chan0,
            scale_xy,
            fill_background=True,
            capture_binary=has_cells,
        )

        binary_path = None
        if has_cells:
            binary_path = config.sample.scratch / f"chan_{ichannel}_binary_{config.sample.name}.dat"

        sample_path = (
            config.sample.save_path / f"chan_{ichannel}_sample_register_{int(registres)}um.tif"
        )
        if sample_path.exists():
            sample_path.unlink()

        _process_channel_to_registration_tiff(
            jobs=jobs,
            scratch_dir=config.sample.scratch,
            sample_name=config.sample.name,
            channel=ichannel,
            out_h=out_h,
            out_w=out_w,
            nz=nz,
            scale_z=scale_z,
            workers=workers,
            output_path=sample_path,
            binary_path=binary_path,
            max_in_memory_bytes=max_ram_bytes,
        )
        regvolpaths[ichannel] = sample_path

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
