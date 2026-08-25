"""Load and build registered spinal cord volumes for export and Napari QC."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import tifffile
from rich.console import Console

from lightsuite.analysis.cord_hemisphere import (
    REGISTERED_HEMISPHERE_FILENAME,
    load_fiederling_hemisphere_native,
)
from lightsuite.atlas.fiederling import load_fiederling_atlas_volumes
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.preprocess.cord_checkpoint import (
    CordRegOptsCheckpoint,
    CordTransformParamsCheckpoint,
)

console = Console()

REGISTERED_ANNOTATION_FILENAME = "annotation_registered.tiff"
REGISTERED_TEMPLATE_FILENAME = "template_registered.tiff"
_CHANNEL_PATTERN = re.compile(r"chan(\d+)_channel\d+\.tiff$", re.IGNORECASE)


@dataclass(frozen=True)
class CordRegisteredInspectPaths:
    volume_registered_dir: Path
    annotation_path: Path
    registered_channels: dict[int, Path] = field(default_factory=dict)


@dataclass
class CordRegisteredInspectVolumes:
    annotation: np.ndarray
    registered_channels: dict[int, np.ndarray]


def volume_registered_dir(config: SpinalCordPipelineConfig) -> Path:
    return config.sample.save_path.expanduser() / "volume_registered"


def cord_stats_dir(config: SpinalCordPipelineConfig | Path) -> Path:
    """Return ``<save_path>/stats`` (created)."""
    save_path = (
        config.sample.save_path.expanduser()
        if isinstance(config, SpinalCordPipelineConfig)
        else Path(config).expanduser()
    )
    path = save_path / "stats"
    path.mkdir(parents=True, exist_ok=True)
    return path


def export_layout_from_native(volume_native: np.ndarray) -> np.ndarray:
    """Map native atlas layout (Z, Y, X) to slice TIFF layout (Y, X, Z)."""
    return np.transpose(volume_native, (1, 2, 0))


def native_layout_from_export(volume_export: np.ndarray) -> np.ndarray:
    """Invert :func:`export_layout_from_native`."""
    return np.transpose(volume_export, (2, 0, 1))


def warp_output_to_uint16(volume: np.ndarray) -> np.ndarray:
    """Clip warped float volumes before uint16 export (MATLAB imwarp clips at 0)."""
    return np.clip(np.rint(volume), 0, np.iinfo(np.uint16).max).astype(np.uint16)


def load_native_template_export_layout(config: SpinalCordPipelineConfig) -> np.ndarray:
    """Native Fiederling template in the same (Y, X, Z) layout as exported channels."""
    atlas_volumes = load_fiederling_atlas_volumes(config.atlas)
    return export_layout_from_native(atlas_volumes.template.astype(np.float32))


def discover_registered_cord_paths(config: SpinalCordPipelineConfig) -> CordRegisteredInspectPaths:
    """Discover exported registered channel TIFFs and annotation under volume_registered/."""
    vr = volume_registered_dir(config)
    if not vr.is_dir():
        msg = f"Missing {vr}. Run 'lightsuite spinal export' first."
        raise FileNotFoundError(msg)

    registered_channels: dict[int, Path] = {}
    for path in sorted(vr.glob("chan*_channel*.tiff")):
        match = _CHANNEL_PATTERN.fullmatch(path.name)
        if match is None:
            continue
        registered_channels[int(match.group(1))] = path.resolve()

    annotation_path = vr / REGISTERED_ANNOTATION_FILENAME
    if not registered_channels:
        msg = (
            f"No registered channel TIFFs found in {vr}. "
            f"Expected files like chan01_channel1.tiff."
        )
        raise FileNotFoundError(msg)
    return CordRegisteredInspectPaths(
        volume_registered_dir=vr.resolve(),
        annotation_path=annotation_path.resolve(),
        registered_channels=registered_channels,
    )


def load_registered_stack(path: Path) -> np.ndarray:
    """Load one exported registered volume (Y, X, Z)."""
    volume = tifffile.imread(path)
    if volume.ndim != 3:
        msg = f"Expected 3D registered TIFF at {path}, got shape {volume.shape}"
        raise ValueError(msg)
    return np.asarray(volume)


def _load_transform_params(save_path: Path) -> CordTransformParamsCheckpoint:
    json_path = save_path / "transform_params.json"
    if not json_path.is_file():
        msg = f"Missing {json_path}. Run 'lightsuite spinal register' first."
        raise FileNotFoundError(msg)
    return CordTransformParamsCheckpoint.load(json_path)


def compute_registered_hemisphere_volume(config: SpinalCordPipelineConfig) -> np.ndarray:
    """Return the Fiederling hemisphere mask in registered export layout (Y, X, Z)."""
    return export_layout_from_native(
        load_fiederling_hemisphere_native(config.atlas.atlas_dir).astype(np.uint8)
    )


def ensure_registered_hemisphere_volume(
    config: SpinalCordPipelineConfig,
    *,
    paths: CordRegisteredInspectPaths | None = None,
    force_recompute: bool = False,
) -> Path:
    """Write hemisphere_registered.tiff if missing (or when forced)."""
    paths = paths or discover_registered_cord_paths(config)
    out_path = paths.volume_registered_dir / REGISTERED_HEMISPHERE_FILENAME
    if out_path.is_file() and not force_recompute:
        return out_path

    console.print("Computing registered hemisphere volume...")
    hemisphere = compute_registered_hemisphere_volume(config)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(out_path, hemisphere, imagej=True)
    console.print(f"Wrote {out_path}")
    return out_path


def load_registered_hemisphere_volume(
    config: SpinalCordPipelineConfig,
    register_path: Path,
) -> np.ndarray:
    """Load or create the registered hemisphere mask aligned to annotation_registered.tiff."""
    out_path = register_path / REGISTERED_HEMISPHERE_FILENAME
    if not out_path.is_file():
        ensure_registered_hemisphere_volume(config)
    return load_registered_stack(out_path).astype(np.uint8, copy=False)


def compute_registered_annotation_volume(
    config: SpinalCordPipelineConfig,
    *,
    regopts: CordRegOptsCheckpoint | None = None,
    transform_params: CordTransformParamsCheckpoint | None = None,
    save_path: Path | None = None,
) -> np.ndarray:
    """Return the atlas annotation in the registered export grid (Y, X, Z).

    Registered sample channels are warped *into* the native atlas space, so the matching
    annotation is simply the native atlas annotation (which already lives in that space and
    aligns with ``template_registered.tiff``). Round-tripping the annotation through the
    sample-space transforms only re-introduces registration error and is unnecessary.
    """
    atlas_volumes = load_fiederling_atlas_volumes(config.atlas)
    return export_layout_from_native(atlas_volumes.annotation.astype(np.uint16))


def ensure_registered_annotation_volume(
    config: SpinalCordPipelineConfig,
    *,
    paths: CordRegisteredInspectPaths | None = None,
    force_recompute: bool = False,
) -> Path:
    """Write annotation_registered.tiff if missing (or when forced)."""
    paths = paths or discover_registered_cord_paths(config)
    if paths.annotation_path.is_file() and not force_recompute:
        return paths.annotation_path

    console.print("Computing registered annotation volume...")
    annotation = compute_registered_annotation_volume(config)
    paths.annotation_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(paths.annotation_path, annotation, imagej=True)
    console.print(f"Wrote {paths.annotation_path}")
    return paths.annotation_path


def load_registered_cord_volumes(
    config: SpinalCordPipelineConfig,
    *,
    paths: CordRegisteredInspectPaths | None = None,
    recompute_annotation: bool = False,
) -> CordRegisteredInspectVolumes:
    """Load exported channels and registered annotation for Napari."""
    paths = paths or discover_registered_cord_paths(config)
    ensure_registered_annotation_volume(config, paths=paths, force_recompute=recompute_annotation)

    registered_channels: dict[int, np.ndarray] = {}
    expected_shape: tuple[int, ...] | None = None
    for ichan, path in paths.registered_channels.items():
        vol = load_registered_stack(path).astype(np.float32, copy=False)
        if expected_shape is None:
            expected_shape = vol.shape
        elif vol.shape != expected_shape:
            msg = f"{path.name} shape {vol.shape} != expected {expected_shape}"
            raise ValueError(msg)
        registered_channels[ichan] = vol

    annotation = load_registered_stack(paths.annotation_path).astype(np.int32, copy=False)
    if expected_shape is not None and annotation.shape != expected_shape:
        msg = (
            f"Annotation shape {annotation.shape} != registered channel shape {expected_shape}. "
            "Re-run 'lightsuite spinal export'."
        )
        raise ValueError(msg)
    return CordRegisteredInspectVolumes(
        annotation=annotation,
        registered_channels=registered_channels,
    )
