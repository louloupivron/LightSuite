"""Read spinal cord TIFF samples (readSpinalCordSample.m port)."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.config.models import CordTiffLayout, SpinalCordPipelineConfig
from lightsuite.io.cord_volume import (
    SkippedSlice,
    load_channel_per_file_stack,
    load_multichannel_single_stack,
    load_plane_per_file_stack,
    normalize_res_um,
    resolve_cord_tiff_layout,
    _sorted_tiff_files,
)

console = Console()


@dataclass(frozen=True)
class CordSampleVolume:
    volume: np.ndarray
    native_orisize: tuple[int, int, int]
    n_channels: int
    layout: CordTiffLayout
    skipped_slices: tuple[SkippedSlice, ...] = ()

    @property
    def n_slices(self) -> int:
        return int(self.volume.shape[0])

    @property
    def ny(self) -> int:
        return int(self.volume.shape[1])

    @property
    def nx(self) -> int:
        return int(self.volume.shape[2])


def read_spinal_cord_sample(config: SpinalCordPipelineConfig) -> CordSampleVolume:
    """Load a spinal cord volume, downsampling to the registration grid when needed."""
    folder = config.sample.source.path
    requested = config.sample.source.tiff_type
    layout = resolve_cord_tiff_layout(folder, requested)
    if requested == CordTiffLayout.AUTO and layout != CordTiffLayout.AUTO:
        console.print(f"Auto-detected TIFF layout: {layout.value}")
    sampleres = normalize_res_um(config.sample.voxel_um)
    regres = normalize_res_um([config.registration.resolution_um] * 3)
    skip_corrupt = config.sample.source.skip_corrupt_slices

    t0 = time.perf_counter()
    if layout == CordTiffLayout.PLANE_PER_FILE:
        volume, native_orisize, skipped, layout = load_plane_per_file_stack(
            folder,
            sampleres_um=sampleres,
            registrationres_um=regres,
            skip_corrupt_slices=skip_corrupt,
        )
    elif layout == CordTiffLayout.MULTICHANNEL_SINGLE:
        files = _sorted_tiff_files(folder)
        if len(files) != 1:
            msg = "multichannel_single layout expects exactly one TIFF file."
            raise ValueError(msg)
        console.print(f"Loading multichannel_single stack: {files[0].name}")
        volume, native_orisize, layout = load_multichannel_single_stack(
            files[0],
            sampleres_um=sampleres,
            registrationres_um=regres,
        )
        skipped = []
    else:
        console.print(f"Loading channel-per-file stack from {folder}")
        volume, native_orisize, layout = load_channel_per_file_stack(
            folder,
            sampleres_um=sampleres,
            registrationres_um=regres,
        )
        skipped = []

    if volume.ndim != 4:
        msg = f"Expected 4D cord volume (Y, X, Z, C), got {volume.shape}"
        raise ValueError(msg)
    if volume.shape[3] < 1 or min(volume.shape[:3]) < 2:
        msg = f"Degenerate cord stack shape {volume.shape}; check TIFF layout and corrupt slices."
        raise ValueError(msg)

    elapsed = time.perf_counter() - t0
    y, x, z, nchan = volume.shape
    console.print(
        f"Parsed spinal cord sample in {elapsed:.1f} s. "
        f"Size {y} x {x} x {z} with {nchan} channel(s)"
    )
    if native_orisize != volume.shape[:3]:
        ny, nx, nz = native_orisize
        console.print(
            f"  (native size was {ny} x {nx} x {nz} px; downsampled to registration grid)"
        )
    if skipped:
        console.print(f"  {len(skipped)} slice(s) skipped as corrupt")

    return CordSampleVolume(
        volume=volume.astype(np.uint16),
        native_orisize=native_orisize,
        n_channels=int(volume.shape[3]),
        layout=layout,
        skipped_slices=tuple(skipped),
    )
