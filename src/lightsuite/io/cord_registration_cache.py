"""Cache registration-grid spinal cord volumes between pipeline stages."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
from rich.console import Console

from lightsuite.config.models import CordTiffLayout, SpinalCordPipelineConfig, TiffLayout
from lightsuite.io.cord_reader import CordSampleVolume, read_spinal_cord_sample
from lightsuite.io.cord_volume import (
    SkippedSlice,
    _single_tiff_stack_info,
    _sorted_tiff_files,
    align_multi_channel_plane_files,
    cord_volume_downsample_spec,
    filter_spinal_cord_slices,
    normalize_res_um,
    read_plane_tiff,
    registration_volume_path,
    resolve_cord_tiff_layout,
)
from lightsuite.io.discover import discover_tiff_stack
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint
from lightsuite.registration.cord_paths import cord_cache_dir, cord_save_path

console = Console()

REGISTER_CACHE_MANIFEST = "register_cache.json"


@dataclass(frozen=True)
class CordSourceProbe:
    native_orisize: tuple[int, int, int]
    n_channels: int
    layout: CordTiffLayout
    n_skipped_slices: int = 0


@dataclass(frozen=True)
class CordRegistrationVolume:
    """Sample volume on the registration grid plus metadata for preprocess."""

    volume: np.ndarray
    native_orisize: tuple[int, int, int]
    n_channels: int
    layout: CordTiffLayout
    regvolpaths: dict[str, str]
    skipped_slices: tuple[SkippedSlice, ...] = ()
    from_cache: bool = False


def _source_roots(config: SpinalCordPipelineConfig) -> list[str]:
    roots = config.sample.source.channel_roots
    if roots is not None:
        return [str(path) for path in roots]
    path = config.sample.source.path
    if path is None:
        msg = "sample.source.path is required"
        raise ValueError(msg)
    return [str(path)]


def probe_cord_source(config: SpinalCordPipelineConfig) -> CordSourceProbe:
    """Infer native stack shape and channel count without loading the full volume."""
    folder = config.sample.source.path
    if folder is None:
        msg = "sample.source.path is required"
        raise ValueError(msg)
    channel_folders = config.sample.source.channel_roots
    layout = resolve_cord_tiff_layout(
        folder,
        config.sample.source.tiff_type,
        channel_folders=channel_folders,
    )
    skip_corrupt = config.sample.source.skip_corrupt_slices

    if layout == CordTiffLayout.PLANE_PER_FILE:
        roots = list(channel_folders) if channel_folders else [folder]
        if len(roots) > 1:
            aligned, skipped, _dropped = align_multi_channel_plane_files(
                roots,
                skip_corrupt_slices=skip_corrupt,
            )
            files = aligned[0]
            ny, nx = read_plane_tiff(files[0]).shape
            return CordSourceProbe(
                native_orisize=(ny, nx, len(files)),
                n_channels=len(roots),
                layout=layout,
                n_skipped_slices=len(skipped),
            )
        files, skipped = filter_spinal_cord_slices(
            _sorted_tiff_files(roots[0]),
            skip_corrupt=skip_corrupt,
        )
        ny, nx = read_plane_tiff(files[0]).shape
        return CordSourceProbe(
            native_orisize=(ny, nx, len(files)),
            n_channels=len(roots),
            layout=layout,
            n_skipped_slices=len(skipped),
        )

    if layout == CordTiffLayout.MULTICHANNEL_SINGLE:
        files = _sorted_tiff_files(folder)
        ny, nx, nz, _, _ = _single_tiff_stack_info(files[0])
        return CordSourceProbe(
            native_orisize=(ny, nx, nz),
            n_channels=1,
            layout=layout,
        )

    discovery = discover_tiff_stack(folder, TiffLayout.CHANNEL_PER_FILE)
    path = discovery.tfiles[0]
    if discovery.multitiffs:
        with tifffile.TiffFile(path) as tif:
            if discovery.use_native_tiff_pages and discovery.stack_read_mode == "pages":
                plane = tif.pages[0].asarray()
                nz = len(tif.pages)
            else:
                stack = tifffile.imread(path)
                if stack.ndim == 2:
                    nz = 1
                    plane = stack
                else:
                    nz = int(stack.shape[2])
                    plane = stack[:, :, 0]
    else:
        stack = tifffile.imread(path)
        if stack.ndim == 2:
            nz = 1
            plane = stack
        else:
            nz = int(stack.shape[2])
            plane = stack[:, :, 0]
    if plane.ndim > 2:
        plane = plane[:, :, 0]
    ny, nx = plane.shape
    return CordSourceProbe(
        native_orisize=(int(ny), int(nx), int(nz)),
        n_channels=len(discovery.tfiles),
        layout=layout,
    )


def compute_cord_register_fingerprint(
    config: SpinalCordPipelineConfig,
    probe: CordSourceProbe,
) -> dict[str, Any]:
    """Inputs that determine registration-grid downsampling."""
    sampleres = normalize_res_um(config.sample.voxel_um)
    regres = normalize_res_um([config.registration.resolution_um] * 3)
    return {
        "source_roots": _source_roots(config),
        "native_orisize": [int(v) for v in probe.native_orisize],
        "nchans": int(probe.n_channels),
        "voxel_um": sampleres.tolist(),
        "registrationres_um": regres.tolist(),
        "tiff_type": probe.layout.value,
        "skip_corrupt_slices": bool(config.sample.source.skip_corrupt_slices),
        "n_skipped_slices": int(probe.n_skipped_slices),
    }


def _expected_registration_shape(
    probe: CordSourceProbe,
    config: SpinalCordPipelineConfig,
) -> tuple[int, int, int]:
    sampleres = normalize_res_um(config.sample.voxel_um)
    regres = normalize_res_um([config.registration.resolution_um] * 3)
    _, target = cord_volume_downsample_spec(probe.native_orisize, sampleres, regres)
    return target


def _probe_from_manifest(manifest: dict[str, Any]) -> CordSourceProbe | None:
    """Rebuild a source probe from a register-cache manifest (no raw TIFF reads)."""
    native = manifest.get("native_orisize")
    tiff_type = manifest.get("tiff_type")
    fingerprint = manifest.get("fingerprint")
    if not isinstance(native, list) or len(native) != 3:
        return None
    if not isinstance(tiff_type, str):
        return None
    if not isinstance(fingerprint, dict):
        return None
    try:
        layout = CordTiffLayout(tiff_type)
    except ValueError:
        return None
    nchans = fingerprint.get("nchans", 1)
    n_skipped = fingerprint.get("n_skipped_slices", 0)
    return CordSourceProbe(
        native_orisize=(int(native[0]), int(native[1]), int(native[2])),
        n_channels=int(nchans),
        layout=layout,
        n_skipped_slices=int(n_skipped),
    )


def _probe_from_checkpoint(checkpoint: CordRegOptsCheckpoint) -> CordSourceProbe:
    return CordSourceProbe(
        native_orisize=(int(checkpoint.orisize[0]), int(checkpoint.orisize[1]), int(checkpoint.orisize[2])),
        n_channels=int(checkpoint.nchans),
        layout=CordTiffLayout(checkpoint.tiff_type),
    )


def _plane_per_file_slice_count(config: SpinalCordPipelineConfig) -> int | None:
    """Return the number of plane TIFFs via directory listing only (single-folder stacks)."""
    folder = config.sample.source.path
    if folder is None:
        return None
    channel_folders = config.sample.source.channel_roots
    try:
        layout = resolve_cord_tiff_layout(
            folder,
            config.sample.source.tiff_type,
            channel_folders=channel_folders,
        )
    except (FileNotFoundError, ValueError):
        return None
    if layout != CordTiffLayout.PLANE_PER_FILE:
        return None
    roots = list(channel_folders) if channel_folders else [folder]
    if len(roots) != 1:
        return None
    return len(_sorted_tiff_files(roots[0]))


def _manifest_fingerprint_matches(
    config: SpinalCordPipelineConfig,
    manifest: dict[str, Any],
) -> bool:
    """True when manifest fingerprint matches config and plane count is unchanged."""
    probe = _probe_from_manifest(manifest)
    if probe is None:
        return False
    if compute_cord_register_fingerprint(config, probe) != manifest.get("fingerprint"):
        return False
    if probe.layout == CordTiffLayout.PLANE_PER_FILE:
        listed = _plane_per_file_slice_count(config)
        if listed is not None:
            expected = probe.native_orisize[2] + probe.n_skipped_slices
            if listed != expected:
                return False
    return True


def _tiff_stack_shape(path: Path) -> tuple[int, int, int] | None:
    """Return (Y, X, Z) for a volume TIFF without decoding pixel data."""
    if not path.is_file():
        return None
    try:
        with tifffile.TiffFile(path) as tif:
            if tif.series:
                shape = tif.series[0].shape
                if len(shape) == 2:
                    return (int(shape[0]), int(shape[1]), 1)
                if len(shape) == 3:
                    return (int(shape[0]), int(shape[1]), int(shape[2]))
            if not tif.pages:
                return None
            page = tif.pages[0]
            if len(tif.pages) == 1:
                return (int(page.shape[0]), int(page.shape[1]), 1)
            return (int(page.shape[0]), int(page.shape[1]), len(tif.pages))
    except (OSError, ValueError, tifffile.TiffFileError):
        return None
    return None


def _registration_tiff_matches(path: Path, expected_shape: tuple[int, int, int]) -> bool:
    shape = _tiff_stack_shape(path)
    return shape == expected_shape


def _regvolpaths_valid(
    regvolpaths: dict[str, str],
    *,
    n_channels: int,
    expected_shape: tuple[int, int, int],
) -> bool:
    if len(regvolpaths) != n_channels:
        return False
    for ich in range(1, n_channels + 1):
        path = regvolpaths.get(str(ich))
        if path is None or not _registration_tiff_matches(Path(path), expected_shape):
            return False
    return True


def _manifest_path(cache_dir: Path) -> Path:
    return cache_dir / REGISTER_CACHE_MANIFEST


def _load_manifest(cache_dir: Path) -> dict[str, Any] | None:
    path = _manifest_path(cache_dir)
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def register_cache_valid(
    cache_dir: Path,
    config: SpinalCordPipelineConfig,
    probe: CordSourceProbe | None = None,
) -> bool:
    """Return True when cached registration TIFFs match the current sample inputs."""
    manifest = _load_manifest(cache_dir)
    if manifest is None:
        return False
    if probe is None:
        if not _manifest_fingerprint_matches(config, manifest):
            return False
        probe = _probe_from_manifest(manifest)
        if probe is None:
            return False
    else:
        fingerprint = compute_cord_register_fingerprint(config, probe)
        if manifest.get("fingerprint") != fingerprint:
            return False
    regvolpaths = manifest.get("regvolpaths")
    if not isinstance(regvolpaths, dict):
        return False
    expected_shape = _expected_registration_shape(probe, config)
    return _regvolpaths_valid(
        {str(k): str(v) for k, v in regvolpaths.items()},
        n_channels=probe.n_channels,
        expected_shape=expected_shape,
    )


def _checkpoint_matches_config(
    checkpoint: CordRegOptsCheckpoint,
    config: SpinalCordPipelineConfig,
) -> bool:
    sampleres = normalize_res_um(config.sample.voxel_um)
    regres = normalize_res_um([config.registration.resolution_um] * 3)
    if checkpoint.data_folder != str(config.sample.source.path):
        return False
    if not np.allclose(checkpoint.sampleres_um, sampleres.tolist()):
        return False
    if not np.allclose(checkpoint.registrationres_um, regres.tolist()):
        return False
    return True


def _checkpoint_matches_inputs(
    checkpoint: CordRegOptsCheckpoint,
    config: SpinalCordPipelineConfig,
    probe: CordSourceProbe,
) -> bool:
    if not _checkpoint_matches_config(checkpoint, config):
        return False
    if checkpoint.nchans != probe.n_channels:
        return False
    if list(checkpoint.orisize) != list(probe.native_orisize):
        return False
    if checkpoint.tiff_type != probe.layout.value:
        return False
    return True


def _try_load_register_cache(
    config: SpinalCordPipelineConfig,
    cache_dir: Path,
) -> CordRegistrationVolume | None:
    """Load cached registration TIFFs without probing raw source planes."""
    manifest = _load_manifest(cache_dir)
    if manifest is None or not _manifest_fingerprint_matches(config, manifest):
        return None
    probe = _probe_from_manifest(manifest)
    if probe is None:
        return None
    regvolpaths = manifest.get("regvolpaths")
    if not isinstance(regvolpaths, dict):
        return None
    expected_shape = _expected_registration_shape(probe, config)
    paths = {str(k): str(v) for k, v in regvolpaths.items()}
    if not _regvolpaths_valid(
        paths,
        n_channels=probe.n_channels,
        expected_shape=expected_shape,
    ):
        return None
    console.print(
        "[green]Using cached registration-grid sample TIFFs[/green] "
        "(inputs unchanged; written by check-orientation or a prior run)."
    )
    return _load_cached_volume(paths, probe=probe, layout=probe.layout)


def _try_load_regopts_cache(
    config: SpinalCordPipelineConfig,
    save_path: Path,
) -> CordRegistrationVolume | None:
    """Load registration TIFFs referenced by regopts.json without probing raw planes."""
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        return None
    checkpoint = CordRegOptsCheckpoint.load(regopts_path)
    if not _checkpoint_matches_config(checkpoint, config):
        return None
    probe = _probe_from_checkpoint(checkpoint)
    if probe.layout == CordTiffLayout.PLANE_PER_FILE:
        listed = _plane_per_file_slice_count(config)
        if listed is not None and listed != probe.native_orisize[2]:
            return None
    regvolpaths = checkpoint.regvolpaths or {}
    expected_shape = _expected_registration_shape(probe, config)
    if not _regvolpaths_valid(
        regvolpaths,
        n_channels=probe.n_channels,
        expected_shape=expected_shape,
    ):
        return None
    console.print(
        "[green]Using cached registration-grid sample TIFFs[/green] "
        "(from regopts.json; inputs unchanged)."
    )
    return _load_cached_volume(regvolpaths, probe=probe, layout=probe.layout)


def _write_register_cache(
    cache_dir: Path,
    *,
    config: SpinalCordPipelineConfig,
    sample: CordSampleVolume,
    fingerprint: dict[str, Any],
) -> dict[str, str]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    regvolpaths: dict[str, str] = {}
    for ich in range(sample.n_channels):
        out_path = registration_volume_path(
            cache_dir,
            ich + 1,
            config.registration.resolution_um,
        )
        tifffile.imwrite(out_path, sample.volume[:, :, :, ich].astype(np.uint16))
        regvolpaths[str(ich + 1)] = str(out_path)
    manifest = {
        "fingerprint": fingerprint,
        "regvolpaths": regvolpaths,
        "native_orisize": list(sample.native_orisize),
        "tiff_type": sample.layout.value,
    }
    _manifest_path(cache_dir).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    console.print(f"Cached {len(regvolpaths)} channel registration TIFF(s) under {cache_dir}")
    return regvolpaths


def _load_cached_volume(
    regvolpaths: dict[str, str],
    *,
    probe: CordSourceProbe,
    layout: CordTiffLayout,
    skipped_slices: tuple[SkippedSlice, ...] = (),
) -> CordRegistrationVolume:
    channels = sorted((int(k), Path(v)) for k, v in regvolpaths.items())
    first = tifffile.imread(channels[0][1])
    volume = np.zeros((*first.shape, len(channels)), dtype=np.uint16)
    volume[:, :, :, 0] = first
    for idx, (_, path) in enumerate(channels[1:], start=1):
        volume[:, :, :, idx] = tifffile.imread(path)
    return CordRegistrationVolume(
        volume=volume,
        native_orisize=probe.native_orisize,
        n_channels=probe.n_channels,
        layout=layout,
        regvolpaths={str(k): str(v) for k, v in regvolpaths.items()},
        skipped_slices=skipped_slices,
        from_cache=True,
    )


def manifest_register_fingerprint(
    config: SpinalCordPipelineConfig,
) -> dict[str, Any] | None:
    """Return the register-cache fingerprint when manifest matches config (no raw TIFF reads)."""
    manifest = _load_manifest(cord_cache_dir(config))
    if manifest is None or not _manifest_fingerprint_matches(config, manifest):
        return None
    probe = _probe_from_manifest(manifest)
    if probe is None:
        return None
    return compute_cord_register_fingerprint(config, probe)


def load_or_cache_cord_registration(
    config: SpinalCordPipelineConfig,
    *,
    force: bool = False,
) -> CordRegistrationVolume:
    """Load the registration-grid sample volume, reusing cache when inputs are unchanged."""
    cache_dir = cord_cache_dir(config)
    save_path = cord_save_path(config)

    if not force:
        cached = _try_load_regopts_cache(config, save_path)
        if cached is not None:
            return cached
        cached = _try_load_register_cache(config, cache_dir)
        if cached is not None:
            return cached

    probe = probe_cord_source(config)
    fingerprint = compute_cord_register_fingerprint(config, probe)

    sample = read_spinal_cord_sample(config)
    regvolpaths = _write_register_cache(
        cache_dir,
        config=config,
        sample=sample,
        fingerprint=fingerprint,
    )
    return CordRegistrationVolume(
        volume=sample.volume,
        native_orisize=sample.native_orisize,
        n_channels=sample.n_channels,
        layout=sample.layout,
        regvolpaths=regvolpaths,
        skipped_slices=sample.skipped_slices,
        from_cache=False,
    )
