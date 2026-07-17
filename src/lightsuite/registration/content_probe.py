"""Load volumes and auto-detect foreground boxes for atlas/sample content cropping."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np
from skimage.transform import resize

from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import resolve_brain_atlas_from_config
from lightsuite.atlas.trim import resolve_atlas_content_box
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.content_bbox import ContentBox, sample_foreground_bbox
from lightsuite.registration.volume import load_registration_volume

CONTENT_PROBE_PREVIEW_MAX_BYTES = 256 * 1024 * 1024


class ContentProbeTarget(str, Enum):
    ATLAS = "atlas"
    SAMPLE = "sample"


@dataclass
class ContentProbeData:
    """Volume and initial crop box for interactive or headless probing."""

    target: ContentProbeTarget
    volume: np.ndarray
    full_shape: tuple[int, int, int]
    box: ContentBox
    source_label: str
    preview_downsampled: bool = False

    @property
    def yaml_section(self) -> str:
        return "atlas" if self.target == ContentProbeTarget.ATLAS else "registration"

    @property
    def yaml_mode_key(self) -> str:
        if self.target == ContentProbeTarget.ATLAS:
            return "content_trim"
        return "sample_content_crop"

    @property
    def yaml_box_key(self) -> str:
        if self.target == ContentProbeTarget.ATLAS:
            return "content_box"
        return "sample_content_box"


def _downsample_for_preview(
    volume: np.ndarray,
    *,
    max_bytes: int = CONTENT_PROBE_PREVIEW_MAX_BYTES,
) -> tuple[np.ndarray, bool]:
    vol = np.asarray(volume, dtype=np.float32)
    if vol.nbytes <= max_bytes:
        return vol, False
    scale = (max_bytes / vol.nbytes) ** (1.0 / 3.0)
    new_shape = tuple(max(1, int(dim * scale)) for dim in vol.shape)
    preview = resize(
        vol,
        new_shape,
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    ).astype(np.float32)
    return preview, True


def _box_for_preview(box: ContentBox, full_shape: tuple[int, int, int], preview_shape: tuple[int, int, int]) -> ContentBox:
    def map_index(value: int, full_len: int, preview_len: int) -> int:
        if full_len <= 1:
            return 0
        return int(round(value * (preview_len - 1) / (full_len - 1)))

    return ContentBox(
        y0=map_index(box.y0, full_shape[0], preview_shape[0]),
        y1=map_index(box.y1, full_shape[0], preview_shape[0]),
        x0=map_index(box.x0, full_shape[1], preview_shape[1]),
        x1=map_index(box.x1, full_shape[1], preview_shape[1]),
        z0=map_index(box.z0, full_shape[2], preview_shape[2]),
        z1=map_index(box.z1, full_shape[2], preview_shape[2]),
    )


def apply_box_mask(volume: np.ndarray, box: ContentBox) -> np.ndarray:
    """Dim voxels outside ``box`` for Napari preview."""
    vol = np.asarray(volume, dtype=np.float32)
    masked = vol * 0.25
    ys, xs, zs = box.slices()
    masked[ys, xs, zs] = vol[ys, xs, zs]
    return masked


def auto_detect_atlas_box(config: BrainPipelineConfig) -> ContentBox:
    atlas = resolve_brain_atlas_from_config(config.atlas)
    annotation = load_atlas_volume(atlas.annotation_path)
    template = load_atlas_volume(atlas.template_path)
    return resolve_atlas_content_box(annotation, template, config.atlas)


def auto_detect_sample_box(config: BrainPipelineConfig) -> ContentBox:
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run 'lightsuite brain preprocess' first."
        raise FileNotFoundError(msg)
    checkpoint = RegOptsCheckpoint.load(regopts_path)
    volume = load_registration_volume(Path(checkpoint.regvolpath))
    box = sample_foreground_bbox(
        volume,
        margin_vox=config.registration.sample_content_margin_vox,
        trim_z=config.registration.sample_content_trim_z,
    )
    if box is None:
        return ContentBox.full_volume(tuple(int(v) for v in volume.shape))
    return box


def load_content_probe_data(
    config: BrainPipelineConfig,
    *,
    target: ContentProbeTarget,
    box: ContentBox | None = None,
) -> ContentProbeData:
    if target == ContentProbeTarget.ATLAS:
        atlas = resolve_brain_atlas_from_config(config.atlas)
        volume = load_atlas_volume(atlas.template_path).astype(np.float32)
        source_label = str(atlas.template_path)
        initial_box = box or auto_detect_atlas_box(config)
    elif target == ContentProbeTarget.SAMPLE:
        save_path = config.sample.save_path.expanduser()
        regopts_path = save_path / "regopts.json"
        if not regopts_path.is_file():
            msg = f"Missing {regopts_path}. Run 'lightsuite brain preprocess' first."
            raise FileNotFoundError(msg)
        checkpoint = RegOptsCheckpoint.load(regopts_path)
        volume = load_registration_volume(Path(checkpoint.regvolpath)).astype(np.float32)
        source_label = str(checkpoint.regvolpath)
        initial_box = box or auto_detect_sample_box(config)
    else:
        msg = f"Unknown content probe target: {target!r}"
        raise ValueError(msg)

    full_shape = tuple(int(v) for v in volume.shape)
    preview, downsampled = _downsample_for_preview(volume)
    return ContentProbeData(
        target=target,
        volume=preview,
        full_shape=full_shape,
        box=initial_box,
        source_label=source_label,
        preview_downsampled=downsampled,
    )


def format_content_box_report(data: ContentProbeData) -> str:
    box = data.box
    lines = [
        f"target: {data.target.value}",
        f"source: {data.source_label}",
        f"full shape (Y, X, Z): {data.full_shape}",
        f"crop box [y0, y1, x0, x1, z0, z1]: {box.to_manual_list()}",
        f"crop start (Y, X, Z): {box.start_yxz}",
        f"crop size (Y, X, Z): {box.size_yxz}",
        "",
        "YAML snippet:",
        f"{data.yaml_section}:",
        f"  {data.yaml_mode_key}: manual",
        f"  {data.yaml_box_key}: {box.to_manual_list()}",
    ]
    return "\n".join(lines)


def preview_box_for_display(data: ContentProbeData, box: ContentBox) -> ContentBox:
    if not data.preview_downsampled:
        return box
    return _box_for_preview(box, data.full_shape, tuple(int(v) for v in data.volume.shape))
