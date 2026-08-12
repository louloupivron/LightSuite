"""Link multires ROI registration outputs into brain view-registration (20 µm sample grid)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from skimage.transform import resize as sk_resize

from lightsuite.config.models import BrainMultiresLinkConfig, BrainPipelineConfig
from lightsuite.import_.brain_import import _registration_shape_native_yxz
from lightsuite.import_.sample_reference import SampleReference, load_sample_reference
from lightsuite.multires.checkpoint import MultiresRegOptsCheckpoint, multires_checkpoint_path
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.volume import (
    _all_pages_are_2d,
    _open_tiff_file,
    load_registration_volume,
    permute_brain_volume,
)


def resolve_multires_checkpoint_path(link: BrainMultiresLinkConfig) -> Path | None:
    """Return multires_regopts.json from an explicit path or linked multires YAML."""
    if link.checkpoint is not None:
        return link.checkpoint.expanduser()
    if link.config is not None:
        from lightsuite.config.loader import load_multires_config

        multires_cfg = load_multires_config(link.config)
        return multires_checkpoint_path(multires_cfg.sample.save_path)
    return None


def discover_multires_registered_roi_paths(
    checkpoint: MultiresRegOptsCheckpoint,
    *,
    use_full_overview: bool = True,
) -> dict[str, Path]:
    """Map channel name to registered ROI TIFF on the overview grid."""
    channels: dict[str, Path] = {}
    reference = checkpoint.reference_channel or "reference"

    if use_full_overview and checkpoint.registered_roi_full_overview_path:
        full_path = Path(checkpoint.registered_roi_full_overview_path).expanduser()
        if full_path.is_file():
            channels[str(reference)] = full_path.resolve()

    if str(reference) not in channels and checkpoint.registered_roi_path:
        crop_path = Path(checkpoint.registered_roi_path).expanduser()
        if crop_path.is_file():
            channels[str(reference)] = crop_path.resolve()

    for channel, path_str in sorted((checkpoint.additional_channel_paths or {}).items()):
        crop_path = Path(path_str).expanduser()
        if not crop_path.is_file():
            continue
        if use_full_overview:
            full_path = crop_path.parent / f"{crop_path.stem}_in_full_overview.tif"
            if full_path.is_file():
                channels[str(channel)] = full_path.resolve()
                continue
        channels[str(channel)] = crop_path.resolve()

    return channels


def load_overview_native_volume_yxz(path: Path) -> np.ndarray:
    """Load an overview-grid TIFF as (Y, X, Z)."""
    return load_registration_volume(path).astype(np.float32, copy=False)


def _tiff_shape_yxz(path: Path) -> tuple[int, int, int]:
    with _open_tiff_file(path) as tif:
        if _all_pages_are_2d(tif):
            ny, nx = (int(v) for v in tif.pages[0].shape)
            return ny, nx, len(tif.pages)
        shape = tif.series[0].shape
        if len(shape) == 2:
            return int(shape[0]), int(shape[1]), 1
        if len(shape) == 3:
            axes = tif.series[0].axes
            if axes in {"ZYX", "IYX"}:
                return int(shape[1]), int(shape[2]), int(shape[0])
            if axes == "XYZ":
                return int(shape[0]), int(shape[1]), int(shape[2])
    msg = f"Could not infer YXZ shape for {path}"
    raise ValueError(msg)


def _overview_start_to_reg_start_yxz(
    crop_start_index: list[int],
    voxel_um: list[float],
    registres_um: float,
) -> tuple[int, int, int]:
    """Map overview XYZ voxel indices to registration-grid YXZ indices (before permute)."""
    ix0, iy0, iz0 = (int(v) for v in crop_start_index)
    vx, vy, vz = (float(v) for v in voxel_um)
    return (
        int(round(iy0 * vy / registres_um)),
        int(round(ix0 * vx / registres_um)),
        int(round(iz0 * vz / registres_um)),
    )


def resample_crop_to_registration_grid(
    crop_yxz: np.ndarray,
    crop_start_index: list[int],
    overview_shape_yxz: tuple[int, int, int],
    *,
    voxel_um: list[float],
    registres_um: float,
    permute: list[int],
    expected_shape: tuple[int, int, int],
) -> np.ndarray:
    """Resample an overlap crop and paste it on the 20 µm grid (no full-overview canvas)."""
    crop_reg = _resample_intensity_to_shape(
        np.asarray(crop_yxz, dtype=np.float32),
        _registration_shape_native_yxz(crop_yxz.shape, voxel_um, registres_um),
    )
    full_reg_native = _registration_shape_native_yxz(overview_shape_yxz, voxel_um, registres_um)
    canvas_native = np.zeros(full_reg_native, dtype=np.float32)
    ry0, rx0, rz0 = _overview_start_to_reg_start_yxz(crop_start_index, voxel_um, registres_um)
    cy, cx, cz = crop_reg.shape
    y1 = min(ry0 + cy, full_reg_native[0])
    x1 = min(rx0 + cx, full_reg_native[1])
    z1 = min(rz0 + cz, full_reg_native[2])
    yy0, xx0, zz0 = max(ry0, 0), max(rx0, 0), max(rz0, 0)
    sy0, sx0, sz0 = yy0 - ry0, xx0 - rx0, zz0 - rz0
    canvas_native[yy0:y1, xx0:x1, zz0:z1] = crop_reg[
        sy0 : sy0 + (y1 - yy0),
        sx0 : sx0 + (x1 - xx0),
        sz0 : sz0 + (z1 - zz0),
    ]
    sample = permute_brain_volume(canvas_native, permute)
    if tuple(sample.shape) != expected_shape:
        sample = _resample_intensity_to_shape(sample, expected_shape)
    return sample


def embed_crop_in_overview_yxz(
    crop_yxz: np.ndarray,
    crop_start_index: list[int],
    overview_shape_yxz: tuple[int, int, int],
) -> np.ndarray:
    """Paste an overlap-crop ROI volume into a zero-filled full-overview canvas."""
    ix0, iy0, iz0 = (int(v) for v in crop_start_index)
    oy, ox, oz = overview_shape_yxz
    canvas = np.zeros(overview_shape_yxz, dtype=crop_yxz.dtype)
    cy, cx, cz = crop_yxz.shape
    y1 = min(iy0 + cy, oy)
    x1 = min(ix0 + cx, ox)
    z1 = min(iz0 + cz, oz)
    yy0, xx0, zz0 = max(iy0, 0), max(ix0, 0), max(iz0, 0)
    sy0, sx0, sz0 = yy0 - iy0, xx0 - ix0, zz0 - iz0
    canvas[yy0:y1, xx0:x1, zz0:z1] = crop_yxz[
        sy0 : sy0 + (y1 - yy0),
        sx0 : sx0 + (x1 - xx0),
        sz0 : sz0 + (z1 - zz0),
    ]
    return canvas


def _resample_intensity_to_shape(
    volume: np.ndarray,
    target_shape_yxz: tuple[int, int, int],
) -> np.ndarray:
    vol = np.asarray(volume, dtype=np.float32)
    if tuple(vol.shape) == target_shape_yxz:
        return vol
    return sk_resize(
        vol,
        target_shape_yxz,
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    ).astype(np.float32, copy=False)


def resample_overview_volume_to_registration_grid(
    volume_yxz: np.ndarray,
    *,
    voxel_um: list[float],
    registres_um: float,
    permute: list[int],
    expected_shape: tuple[int, int, int],
) -> np.ndarray:
    """Downsample an overview-native volume onto the permuted 20 µm registration grid."""
    vol = np.asarray(volume_yxz, dtype=np.float32)
    reg_shape = _registration_shape_native_yxz(vol.shape, voxel_um, registres_um)
    if tuple(vol.shape) != reg_shape:
        vol = _resample_intensity_to_shape(vol, reg_shape)
    sample = permute_brain_volume(vol, permute)
    if tuple(sample.shape) != expected_shape:
        sample = _resample_intensity_to_shape(sample, expected_shape)
    return sample


def _channel_on_registration_grid(
    path: Path,
    multires_checkpoint: MultiresRegOptsCheckpoint,
    *,
    sample_reference: SampleReference,
    voxel_um: list[float],
    registres_um: float,
    permute: list[int],
    expected_shape: tuple[int, int, int],
) -> np.ndarray | None:
    path = path.expanduser()
    if not path.is_file():
        return None

    overview_shape = sample_reference.shape_tuple
    on_disk_shape = _tiff_shape_yxz(path)
    crop_start = multires_checkpoint.crop_start_index

    if on_disk_shape == overview_shape:
        volume = load_overview_native_volume_yxz(path)
        return resample_overview_volume_to_registration_grid(
            volume,
            voxel_um=voxel_um,
            registres_um=registres_um,
            permute=permute,
            expected_shape=expected_shape,
        )

    if crop_start is None:
        return None

    crop = load_overview_native_volume_yxz(path)
    return resample_crop_to_registration_grid(
        crop,
        crop_start,
        overview_shape,
        voxel_um=voxel_um,
        registres_um=registres_um,
        permute=permute,
        expected_shape=expected_shape,
    )


def _overview_volume_for_channel(
    path: Path,
    checkpoint: MultiresRegOptsCheckpoint,
    *,
    sample_reference: SampleReference,
) -> np.ndarray | None:
    """Deprecated path kept for tests; prefer ``_channel_on_registration_grid``."""
    path = path.expanduser()
    if not path.is_file():
        return None
    volume = load_overview_native_volume_yxz(path)
    overview_shape = sample_reference.shape_tuple
    if tuple(volume.shape) == overview_shape:
        return volume
    if checkpoint.crop_start_index is None:
        return None
    return embed_crop_in_overview_yxz(volume, checkpoint.crop_start_index, overview_shape)


def load_multires_roi_channels_on_registration_grid(
    config: BrainPipelineConfig,
    *,
    checkpoint: RegOptsCheckpoint,
    transform_params: object,
    expected_shape: tuple[int, int, int],
    sample_reference: SampleReference | None = None,
) -> dict[str, np.ndarray]:
    """Load multires registered ROI channel(s) resampled to the brain 20 µm grid."""
    link = config.multires_link
    if link is None:
        return {}

    checkpoint_path = resolve_multires_checkpoint_path(link)
    if checkpoint_path is None or not checkpoint_path.is_file():
        return {}

    multires_checkpoint = MultiresRegOptsCheckpoint.load(checkpoint_path)
    roi_paths = discover_multires_registered_roi_paths(
        multires_checkpoint,
        use_full_overview=link.use_full_overview,
    )
    if not roi_paths:
        return {}

    save_path = config.sample.save_path.expanduser()
    reference = sample_reference or load_sample_reference(save_path)
    voxel_um = [float(v) for v in checkpoint.voxel_um]
    registres_um = float(checkpoint.registres_um)
    permute = getattr(transform_params, "permute_sample_to_atlas", None) or [1, 2, 3]

    channels: dict[str, np.ndarray] = {}
    for channel, path in roi_paths.items():
        grid = _channel_on_registration_grid(
            path,
            multires_checkpoint,
            sample_reference=reference,
            voxel_um=voxel_um,
            registres_um=registres_um,
            permute=permute,
            expected_shape=expected_shape,
        )
        if grid is not None:
            channels[channel] = grid
    return channels
