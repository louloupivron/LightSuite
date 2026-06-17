"""Napari viewer for atlas-space export + annotation import QA."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import nibabel as nib
import numpy as np

from lightsuite.atlas.registry import resolve_brain_atlas
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.registration.volume import load_registration_volume


@dataclass(frozen=True)
class BrainImportInspectPaths:
    """Resolved inputs under ``<save_path>/volume_registered/``."""

    volume_registered_dir: Path
    template_path: Path
    annotation_path: Path
    registered_channels: dict[int, Path] = field(default_factory=dict)
    point_npz_paths: dict[str, Path] = field(default_factory=dict)
    mask_paths: dict[str, Path] = field(default_factory=dict)


@dataclass
class BrainImportInspectVolumes:
    template: np.ndarray
    annotation: np.ndarray
    registered_channels: dict[int, np.ndarray]
    point_layers: dict[str, np.ndarray]
    mask_layers: dict[str, np.ndarray]


def _volume_registered_dir(config: BrainPipelineConfig) -> Path:
    return config.sample.save_path.expanduser() / "volume_registered"


def _label_from_stem(stem: str, suffix: str) -> str:
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def discover_brain_import_inspect_paths(config: BrainPipelineConfig) -> BrainImportInspectPaths:
    """Discover atlas TIFFs, registered channels, and import annotation outputs."""
    save_path = config.sample.save_path.expanduser()
    vr = _volume_registered_dir(config)
    if not vr.is_dir():
        msg = (
            f"Missing {vr}. Run 'lightsuite brain export --save-volume' and "
            "'lightsuite brain import-annotations' first."
        )
        raise FileNotFoundError(msg)

    transform_params_path = save_path / "transform_params.json"
    if not transform_params_path.is_file():
        msg = f"Missing {transform_params_path}. Run 'lightsuite brain register' first."
        raise FileNotFoundError(msg)

    atlas = resolve_brain_atlas(config.atlas.provider.value, config.atlas.atlas_dir)

    registered_channels: dict[int, Path] = {}
    for path in sorted(vr.glob("chan_*_registered_atlas.tif")):
        match = re.fullmatch(r"chan_(\d+)_registered_atlas\.tif", path.name)
        if match is None:
            continue
        registered_channels[int(match.group(1))] = path.resolve()

    point_npz_paths: dict[str, Path] = {}
    for path in sorted(vr.glob("*_atlas_coords.npz")):
        label = _label_from_stem(path.stem, "_atlas_coords")
        point_npz_paths[label] = path.resolve()

    mask_paths: dict[str, Path] = {}
    for path in sorted(vr.glob("*_registered_atlas.tif")):
        if path.name.startswith("chan_"):
            continue
        label = _label_from_stem(path.stem, "_registered_atlas")
        mask_paths[label] = path.resolve()

    if not registered_channels and not point_npz_paths and not mask_paths:
        msg = (
            f"No inspectable layers found in {vr}. "
            "Expected chan_*_registered_atlas.tif and/or import outputs "
            "(*_atlas_coords.npz, *_registered_atlas.tif)."
        )
        raise FileNotFoundError(msg)

    return BrainImportInspectPaths(
        volume_registered_dir=vr.resolve(),
        template_path=atlas.template_path.resolve(),
        annotation_path=atlas.annotation_path.resolve(),
        registered_channels=registered_channels,
        point_npz_paths=point_npz_paths,
        mask_paths=mask_paths,
    )


def _load_atlas_volume_yxz(path: Path) -> np.ndarray:
    data = np.asanyarray(nib.load(str(path)).dataobj)
    if data.ndim != 3:
        msg = f"Expected 3D atlas volume in {path}, got shape {data.shape}"
        raise ValueError(msg)
    return data.astype(np.float32, copy=False)


def _load_points_npz(path: Path) -> np.ndarray:
    with np.load(path) as archive:
        if "atlasptcoords" not in archive:
            msg = f"NPZ missing atlasptcoords array: {path}"
            raise KeyError(msg)
        coords = np.asarray(archive["atlasptcoords"], dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] < 3:
        msg = f"Expected Nx3+ atlasptcoords in {path}, got shape {coords.shape}"
        raise ValueError(msg)
    return coords


def _load_points_csv(path: Path) -> np.ndarray:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            msg = f"CSV missing header row: {path}"
            raise ValueError(msg)
        field_map = {name.strip().lower(): name for name in reader.fieldnames}
        for required in ("atlas_x", "atlas_y", "atlas_z"):
            if required not in field_map:
                msg = f"CSV missing column {required!r} in {path}"
                raise KeyError(msg)
        rows = list(reader)
    if not rows:
        msg = f"No point rows in {path}"
        raise ValueError(msg)
    keys = field_map["atlas_x"], field_map["atlas_y"], field_map["atlas_z"]
    return np.column_stack([[float(row[k]) for row in rows] for k in keys])


def atlas_points_to_napari_zyx(coords_xyz: np.ndarray) -> np.ndarray:
    """Map 1-based atlas voxel indices (x, y, z) to Napari point coordinates (z, y, x)."""
    pts = np.asarray(coords_xyz, dtype=np.float64)
    if pts.size == 0:
        return np.zeros((0, 3), dtype=np.float64)
    return np.column_stack([pts[:, 2] - 1.0, pts[:, 1] - 1.0, pts[:, 0] - 1.0])


def volume_yxz_to_napari_zyx(volume: np.ndarray) -> np.ndarray:
    """LightSuite volumes are (Y, X, Z); Napari 3D images use (Z, Y, X)."""
    return np.transpose(volume, (2, 0, 1))


def _contrast_limits(volume: np.ndarray) -> tuple[float, float]:
    positive = volume[volume > 0]
    if positive.size:
        lo, hi = np.percentile(positive, (1.0, 99.5))
    else:
        lo, hi = float(np.min(volume)), float(np.max(volume))
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def load_brain_import_inspect_volumes(
    config: BrainPipelineConfig,
    *,
    paths: BrainImportInspectPaths | None = None,
) -> BrainImportInspectVolumes:
    """Load atlas, registered channels, masks, and point layers for Napari."""
    paths = paths or discover_brain_import_inspect_paths(config)
    save_path = config.sample.save_path.expanduser()
    transform_params = _load_transform_params(save_path)
    expected_shape = tuple(int(v) for v in transform_params.atlassize)

    template = _load_atlas_volume_yxz(paths.template_path)
    annotation = _load_atlas_volume_yxz(paths.annotation_path)
    if template.shape != expected_shape:
        msg = (
            f"Allen template shape {template.shape} != transform atlassize {expected_shape}. "
            "Re-run register/export on this sample."
        )
        raise ValueError(msg)
    if annotation.shape != expected_shape:
        msg = (
            f"Allen annotation shape {annotation.shape} != transform atlassize {expected_shape}."
        )
        raise ValueError(msg)

    registered_channels: dict[int, np.ndarray] = {}
    for ichan, path in paths.registered_channels.items():
        vol = load_registration_volume(path)
        if tuple(vol.shape) != expected_shape:
            msg = f"{path.name} shape {vol.shape} != expected atlas shape {expected_shape}"
            raise ValueError(msg)
        registered_channels[ichan] = vol.astype(np.float32, copy=False)

    mask_layers: dict[str, np.ndarray] = {}
    for label, path in paths.mask_paths.items():
        vol = load_registration_volume(path)
        if tuple(vol.shape) != expected_shape:
            msg = f"{path.name} shape {vol.shape} != expected atlas shape {expected_shape}"
            raise ValueError(msg)
        mask_layers[label] = vol.astype(np.float32, copy=False)

    point_layers: dict[str, np.ndarray] = {}
    for label, npz_path in paths.point_npz_paths.items():
        point_layers[label] = _load_points_npz(npz_path)
    # Optional CSV fallback when NPZ is absent
    for path in sorted(paths.volume_registered_dir.glob("*_atlas_coords.csv")):
        label = _label_from_stem(path.stem, "_atlas_coords")
        if label in point_layers:
            continue
        point_layers[label] = _load_points_csv(path)

    return BrainImportInspectVolumes(
        template=template,
        annotation=annotation,
        registered_channels=registered_channels,
        point_layers=point_layers,
        mask_layers=mask_layers,
    )


def run_brain_inspect_imports(
    config: BrainPipelineConfig,
    *,
    headless: bool = False,
) -> BrainImportInspectPaths:
    """Open Napari to QC registered channels and imported annotations in atlas space."""
    paths = discover_brain_import_inspect_paths(config)
    if headless:
        load_brain_import_inspect_volumes(config, paths=paths)
        return paths

    try:
        import napari
        from napari.utils.notifications import show_info
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    volumes = load_brain_import_inspect_volumes(config, paths=paths)
    viewer = napari.Viewer(title=f"LightSuite brain import QC — {config.sample.name}")

    template_limits = _contrast_limits(volumes.template)
    viewer.add_image(
        volume_yxz_to_napari_zyx(volumes.template),
        name="Allen template",
        colormap="gray",
        blending="opaque",
        contrast_limits=template_limits,
    )
    viewer.add_image(
        volume_yxz_to_napari_zyx(volumes.annotation),
        name="Allen annotation",
        colormap="green",
        blending="additive",
        opacity=0.2,
        contrast_limits=(0.0, float(np.max(volumes.annotation)) or 1.0),
    )

    channel_cmaps = ["magenta", "cyan", "yellow", "red"]
    for idx, (ichan, vol) in enumerate(sorted(volumes.registered_channels.items())):
        cmap = channel_cmaps[idx % len(channel_cmaps)]
        viewer.add_image(
            volume_yxz_to_napari_zyx(vol),
            name=f"channel {ichan} registered",
            colormap=cmap,
            blending="additive",
            opacity=0.55,
            contrast_limits=_contrast_limits(vol),
        )

    for label, mask in volumes.mask_layers.items():
        viewer.add_image(
            volume_yxz_to_napari_zyx(mask),
            name=f"mask: {label}",
            colormap="red",
            blending="additive",
            opacity=0.35,
            contrast_limits=(0.0, 1.0),
        )

    for label, coords in volumes.point_layers.items():
        napari_pts = atlas_points_to_napari_zyx(coords)
        viewer.add_points(
            napari_pts,
            name=f"points: {label}",
            size=3,
            face_color="red",
            border_color="white",
        )

    summary_path = paths.volume_registered_dir / "import_annotations_summary.json"
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        show_info(
            f"Loaded {len(volumes.registered_channels)} channel(s), "
            f"{len(volumes.point_layers)} point layer(s), "
            f"{len(volumes.mask_layers)} mask layer(s). "
            f"Summary: {summary_path.name}"
        )
    else:
        show_info(
            f"Loaded {len(volumes.registered_channels)} channel(s), "
            f"{len(volumes.point_layers)} point layer(s), "
            f"{len(volumes.mask_layers)} mask layer(s)."
        )
    napari.run()
    return paths
