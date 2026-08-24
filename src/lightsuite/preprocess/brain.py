"""Brain lightsheet preprocessing (preprocessLightSheetVolume.m port)."""

from __future__ import annotations

import time
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
from rich.console import Console

from lightsuite.config.models import BrainPipelineConfig, TiffLayout
from lightsuite.reporter import (
    check_stage_cancelled,
    emit_pipeline_message,
    format_duration,
    iter_cancellable_process_map,
    report_step_progress,
)
from lightsuite.io.discover import (
    TiffStackDiscovery,
    discover_tiff_stack,
    downsample_duration_hint,
)
from lightsuite.io.readers.tiff_stack import TiffStackReader
from lightsuite.preprocess.checkpoint import (
    RegOptsCheckpoint,
    compute_preprocess_fingerprint,
)
from lightsuite.preprocess.sample_crop import apply_sample_content_crop
from lightsuite.preprocess.slice_ops import (
    SliceLoadJob,
    SliceProcessResult,
    output_xy_shape,
    output_z_count,
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
        if discovery.channel_plane_files is not None:
            plane_paths = discovery.channel_plane_files[chan0]
        else:
            plane_paths = discovery.tfiles
        for path in plane_paths:
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
    return 1


def _iter_processed_slices(
    jobs: list[SliceLoadJob],
    workers: int,
) -> Iterator[SliceProcessResult]:
    yield from iter_cancellable_process_map(
        process_slice_job,
        jobs,
        max_workers=max(1, workers),
        chunksize=1,
    )


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

    label = f"channel {channel}"
    t0 = time.perf_counter()
    try:
        for islice, result in enumerate(_iter_processed_slices(jobs, workers), start=1):
            xy_stack[:, :, islice - 1] = result.plane_xy
            if binary_handle is not None and result.binary_bytes is not None:
                binary_handle.write(result.binary_bytes)
            report_step_progress(
                islice,
                nz,
                label=label,
                t0=t0,
                workers=workers,
                every=max(20, nz // 10) if nz >= 20 else 1,
                unit="slice",
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

    emit_pipeline_message(f"[{label}] writing registration TIFF ({output_path.name})…")
    t_save = time.perf_counter()
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
    emit_pipeline_message(
        f"[{label}] saved registration TIFF in {time.perf_counter() - t_save:.1f}s"
    )


def _registration_volume_status(
    path: Path,
    *,
    expected_shape: tuple[int, int],
    expected_pages: int,
) -> tuple[bool, str]:
    want = f"{expected_pages} pages at {expected_shape[0]}×{expected_shape[1]}"
    if not path.is_file():
        return False, f"{path.name} not found (need {want})"
    try:
        with tifffile.TiffFile(path) as tif:
            n_pages = len(tif.pages)
            shape = tuple(int(v) for v in tif.pages[0].shape)
    except (OSError, ValueError, tifffile.TiffFileError) as exc:
        return False, f"{path.name} could not be read ({exc})"
    if n_pages != expected_pages or shape != expected_shape:
        return False, f"{path.name} is {n_pages} pages at {shape[0]}×{shape[1]}, expected {want}"
    return True, f"{path.name} matches cache ({want})"


def _registration_tiff_path(save_path: Path, channel: int, registres_um: float) -> Path:
    return save_path / f"chan_{channel}_sample_register_{int(registres_um)}um.tif"


def _fingerprint_matches(
    stored: dict[str, Any] | None,
    current: dict[str, Any],
) -> bool:
    return stored is not None and stored == current


def _load_existing_checkpoint(save_path: Path) -> RegOptsCheckpoint | None:
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        return None
    return RegOptsCheckpoint.load(regopts_path)


def _build_preprocess_checkpoint(
    config: BrainPipelineConfig,
    *,
    ny: int,
    nx: int,
    nz: int,
    nchans: int,
    voxel_um: list[float],
    registres_um: float,
    regvolpaths: dict[int, Path],
    fingerprint: dict[str, Any],
) -> RegOptsCheckpoint:
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

    return RegOptsCheckpoint(
        sample_name=config.sample.name,
        ny=ny,
        nx=nx,
        nz=nz,
        nchans=nchans,
        voxel_um=voxel_um,
        registres_um=registres_um,
        regvolpath=str(regvolpaths[primary]),
        regvolpath_secondary=str(regvolpaths[secondary]) if secondary else None,
        regvolpaths={str(k): str(v) for k, v in regvolpaths.items()},
        tiff_type=config.sample.source.tiff_type.value,
        channel_primary=primary,
        channel_secondary=secondary,
        preprocess_fingerprint=fingerprint,
    )


def preprocess_lightsheet_volume(
    config: BrainPipelineConfig,
    *,
    force: bool = False,
) -> PreprocessResult:
    """Downsample channels and write registration TIFFs under save_path."""
    vx, vy, vz = _require_voxel_um(config)
    registres = config.registration.resolution_um
    scale_xy = vx / registres
    scale_z = vz / registres
    requested_workers = config.compute.workers

    config.sample.scratch.mkdir(parents=True, exist_ok=True)
    config.sample.save_path.mkdir(parents=True, exist_ok=True)

    emit_pipeline_message(
        f"Discovering TIFF stack at {config.sample.source.path} "
        f"({config.sample.source.tiff_type.value})…"
    )
    discovery = discover_tiff_stack(
        config.sample.source.path,
        tiff_type=config.sample.source.tiff_type,
        channel_folders=config.sample.source.channel_roots,
    )
    reader = TiffStackReader(discovery, voxel_um=(vx, vy, vz))
    workers = _effective_preprocess_workers(discovery, requested_workers)

    ny, nx, nz, nchans = discovery.ny, discovery.nx, discovery.nz, discovery.nchans
    emit_pipeline_message(
        f"Stack: {ny}×{nx}, {nz} planes, {nchans} channel(s) "
        f"({discovery.tiff_type.value}) → {registres:g} µm"
    )
    regvolpaths: dict[int, Path] = {}
    cell_channel = _channel_for_cells(config)
    out_h, out_w = output_xy_shape(ny, nx, scale_xy)
    expected_pages = output_z_count(nz, scale_z)
    expected_shape = (out_h, out_w)
    fingerprint = compute_preprocess_fingerprint(
        ny=ny,
        nx=nx,
        nz=nz,
        nchans=nchans,
        voxel_um=[vx, vy, vz],
        registres_um=registres,
        tiff_type=config.sample.source.tiff_type.value,
        sample_content_crop=config.registration.sample_content_crop.value,
        sample_content_box=config.registration.sample_content_box,
    )
    existing = _load_existing_checkpoint(config.sample.save_path)
    fingerprint_unchanged = _fingerprint_matches(
        existing.preprocess_fingerprint if existing is not None else None,
        fingerprint,
    )
    cache_statuses = [
        _registration_volume_status(
            _registration_tiff_path(config.sample.save_path, ichannel, registres),
            expected_shape=expected_shape,
            expected_pages=expected_pages,
        )
        for ichannel in range(1, nchans + 1)
    ]
    cached_tiffs_valid = all(ok for ok, _reason in cache_statuses)
    skip_downsample = not force and cached_tiffs_valid
    if force and cached_tiffs_valid:
        emit_pipeline_message(
            "force=True: ignoring matching registration TIFFs and downsampling again."
        )
    elif skip_downsample:
        console.print("[green]Using cached registration TIFFs.[/green]")
        if existing is not None and (
            existing.channel_primary != config.registration.channel_primary
            or existing.channel_secondary != config.registration.channel_secondary
        ):
            console.print(
                "[yellow]channel_primary / channel_secondary changed — "
                "re-run init-registration if you switch the primary channel.[/yellow]"
            )
    else:
        for _ok, reason in cache_statuses:
            if not _ok:
                emit_pipeline_message(f"Not reusing cache: {reason}")
        sample_file = discovery.tfiles[0] if discovery.tfiles else config.sample.source.path
        emit_pipeline_message(downsample_duration_hint(nz, nchans, sample_file))

    scratch_bytes = out_h * out_w * nz * 2
    max_ram_bytes = int(config.compute.max_in_memory_scratch_gb * (1024**3))
    scratch_in_ram = scratch_bytes <= max_ram_bytes
    if not skip_downsample:
        if discovery.stack_read_mode != "pages":
            console.print(
                f"Single-file volumetric TIFF: {nz} Z planes "
                f"({discovery.stack_read_mode}, channelperfile)."
            )
        elif discovery.tiff_type == TiffLayout.PLANE_PER_FILE and nchans >= 1:
            est_xy_gb = scratch_bytes / (1024**3)
            where = "RAM" if scratch_in_ram else f"disk memmap on {config.sample.scratch}"
            channel_note = f", {nchans} channels" if nchans > 1 else ""
            console.print(
                f"Downsampling {nz} planes{channel_note}; XY scratch ~{est_xy_gb:.1f} GB in {where}."
            )

    for ichannel in range(1, nchans + 1):
        check_stage_cancelled()
        chan0 = ichannel - 1
        sample_path = _registration_tiff_path(config.sample.save_path, ichannel, registres)
        regvolpaths[ichannel] = sample_path

        if skip_downsample:
            continue

        has_cells = cell_channel is not None and ichannel == cell_channel
        if nchans > 1:
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

    crop_meta: tuple[list[int], list[int], list[float]] | None = None
    needs_crop = config.registration.sample_content_crop.value != "off"
    if needs_crop and (
        not skip_downsample
        or existing is None
        or existing.content_crop_start is None
    ):
        crop_meta = apply_sample_content_crop(
            config,
            regvolpaths,
            primary_channel=config.registration.channel_primary,
        )
    elif existing is not None and existing.content_crop_start is not None:
        crop_meta = (
            existing.content_crop_start,
            existing.content_crop_size or [],
            existing.native_crop_offset_yxz or [],
        )

    reader.close()

    checkpoint = _build_preprocess_checkpoint(
        config,
        ny=ny,
        nx=nx,
        nz=nz,
        nchans=nchans,
        voxel_um=[vx, vy, vz],
        registres_um=registres,
        regvolpaths=regvolpaths,
        fingerprint=fingerprint,
    )
    if existing is not None and (skip_downsample or fingerprint_unchanged):
        checkpoint = checkpoint.merge_downstream_from(existing)
    if crop_meta is not None:
        crop_start, crop_size, native_off = crop_meta
        checkpoint.content_crop_start = crop_start
        checkpoint.content_crop_size = crop_size
        checkpoint.native_crop_offset_yxz = native_off
    regopts_path = config.sample.save_path / "regopts.json"
    checkpoint.save(regopts_path)

    from lightsuite.import_.sample_reference import write_sample_reference

    ref_path = write_sample_reference(
        config.sample.save_path,
        sample_name=config.sample.name,
        ny=ny,
        nx=nx,
        nz=nz,
        voxel_um=[vx, vy, vz],
    )
    console.print(f"Wrote {regopts_path.name} and {ref_path.name}")

    return PreprocessResult(checkpoint=checkpoint, regvolpaths=regvolpaths)
