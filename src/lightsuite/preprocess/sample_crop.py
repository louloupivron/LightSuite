"""Crop registration TIFFs to sample foreground after brain preprocess."""

from __future__ import annotations

from pathlib import Path

import tifffile
from rich.console import Console

from lightsuite.config.models import BrainPipelineConfig, SampleContentCropMode
from lightsuite.registration.content_bbox import (
    ContentBox,
    crop_volume_yxz,
    sample_foreground_bbox,
)
from lightsuite.registration.coordinates import content_box_to_lists, native_crop_offset_yxz
from lightsuite.registration.volume import load_registration_volume

console = Console()


def resolve_sample_content_box(
    volume: object,
    config: BrainPipelineConfig,
) -> ContentBox | None:
    mode = config.registration.sample_content_crop
    if mode == SampleContentCropMode.OFF:
        return None
    if mode == SampleContentCropMode.MANUAL:
        if config.registration.sample_content_box is None:
            msg = "registration.sample_content_box required when sample_content_crop=manual"
            raise ValueError(msg)
        return ContentBox.from_manual_box(config.registration.sample_content_box)
    return sample_foreground_bbox(
        volume,
        margin_vox=config.registration.sample_content_margin_vox,
        trim_z=config.registration.sample_content_trim_z,
    )


def crop_registration_tiff(path: Path, box: ContentBox) -> tuple[int, int, int]:
    """Crop a plane-per-page registration TIFF in place; return new (H, W, Z)."""
    volume = load_registration_volume(path)
    cropped = crop_volume_yxz(volume, box)
    if path.exists():
        path.unlink()
    with tifffile.TiffWriter(path) as writer:
        for z in range(cropped.shape[2]):
            writer.write(
                cropped[:, :, z],
                compression="lzw",
                photometric="minisblack",
            )
    return tuple(int(v) for v in cropped.shape)


def apply_sample_content_crop(
    config: BrainPipelineConfig,
    regvolpaths: dict[int, Path],
    *,
    primary_channel: int,
) -> tuple[list[int], list[int], list[float]] | None:
    """Crop all registration channels to a shared box; return crop metadata or None."""
    primary_path = regvolpaths[primary_channel]
    primary_vol = load_registration_volume(primary_path)
    box = resolve_sample_content_box(primary_vol, config)
    if box is None:
        return None
    if box.size_yxz == primary_vol.shape and box.start_yxz == (0, 0, 0):
        console.print("[dim]Sample content crop: foreground spans full registration grid.[/dim]")
        return None

    console.print(
        f"[green]Sample content crop:[/green] {primary_vol.shape} → {box.size_yxz} "
        f"(origin {box.start_yxz})"
    )
    for path in regvolpaths.values():
        crop_registration_tiff(path, box)

    crop_start, crop_size = content_box_to_lists(box)
    native_off = native_crop_offset_yxz(
        box.start_yxz,
        voxel_um=config.sample.voxel_um or [1.0, 1.0, 1.0],
        registres_um=config.registration.resolution_um,
    )
    return crop_start, crop_size, list(native_off)
