"""Napari viewer for atlas-space export + annotation import QA."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np

from lightsuite.atlas.registry import (
    atlas_display_provider_from_config,
    resolve_brain_atlas_from_config,
)
from lightsuite.atlas.display import (
    atlas_points_xyz_to_napari_zyx,
    atlas_volume_yxz_to_napari_zyx,
    registration_points_xyz_to_napari_zyx,
    registration_volume_yxz_to_napari_zyx,
)
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.registration.volume import load_registration_volume

ViewSpace = Literal["atlas", "sample"]


@dataclass(frozen=True)
class BrainImportInspectPaths:
    """Resolved inputs for brain import Napari QC."""

    volume_registered_dir: Path
    template_path: Path | None
    annotation_path: Path | None
    registered_channels: dict[int, Path] = field(default_factory=dict)
    point_npz_paths: dict[str, Path] = field(default_factory=dict)
    mask_paths: dict[str, Path] = field(default_factory=dict)
    space: ViewSpace = "atlas"
    sample_space_dir: Path | None = None


@dataclass
class BrainImportInspectVolumes:
    template: np.ndarray | None
    annotation: np.ndarray | None
    registered_channels: dict[int, np.ndarray]
    point_layers: dict[str, np.ndarray]
    mask_layers: dict[str, np.ndarray]
    channels_warped_on_the_fly: bool = False


def _volume_registered_dir(config: BrainPipelineConfig) -> Path:
    return config.sample.save_path.expanduser() / "volume_registered"


def _label_from_stem(stem: str, suffix: str) -> str:
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def discover_brain_import_inspect_paths(
    config: BrainPipelineConfig,
    *,
    space: ViewSpace = "atlas",
) -> BrainImportInspectPaths:
    """Discover registered volumes and imported annotation outputs for Napari QC."""
    if space == "sample":
        return _discover_brain_import_inspect_paths_sample(config)
    return _discover_brain_import_inspect_paths_atlas(config)


def _discover_brain_import_inspect_paths_atlas(
    config: BrainPipelineConfig,
) -> BrainImportInspectPaths:
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

    atlas = resolve_brain_atlas_from_config(config.atlas)

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
        space="atlas",
    )


def _discover_brain_import_inspect_paths_sample(
    config: BrainPipelineConfig,
) -> BrainImportInspectPaths:
    from lightsuite.export.brain_sample_space import discover_brain_sample_space_inspect_paths

    sample_paths = discover_brain_sample_space_inspect_paths(config)
    return BrainImportInspectPaths(
        volume_registered_dir=sample_paths.volume_registered_dir,
        template_path=sample_paths.template_path,
        annotation_path=sample_paths.annotation_path,
        registered_channels=sample_paths.registered_channels,
        point_npz_paths=sample_paths.point_npz_paths,
        mask_paths=sample_paths.mask_paths,
        space="sample",
        sample_space_dir=sample_paths.sample_space_dir,
    )


def _load_atlas_volume_yxz(path: Path) -> np.ndarray:
    from lightsuite.atlas.io import load_atlas_volume

    return load_atlas_volume(path).astype(np.float32, copy=False)


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
    """Map 1-based voxel indices (x, y, z) to Napari point coordinates (z, y, x)."""
    pts = np.asarray(coords_xyz, dtype=np.float64)
    if pts.size == 0:
        return np.zeros((0, 3), dtype=np.float64)
    return np.column_stack([pts[:, 2] - 1.0, pts[:, 1] - 1.0, pts[:, 0] - 1.0])


def volume_yxz_to_napari_zyx(volume: np.ndarray) -> np.ndarray:
    """LightSuite volumes are (Y, X, Z); Napari 3D images use (Z, Y, X)."""
    return np.transpose(volume, (2, 0, 1))


def brain_volume_to_napari_zyx(
    volume: np.ndarray,
    *,
    space: ViewSpace,
    atlas_provider: str | None = None,
    permute_sample_to_atlas: list[int] | None = None,
) -> np.ndarray:
    """Map a brain volume to Napari ZYX using atlas QC layout when in atlas space."""
    if space == "atlas":
        if atlas_provider is None:
            msg = "atlas_provider is required for atlas-space Napari display"
            raise ValueError(msg)
        return atlas_volume_yxz_to_napari_zyx(volume, atlas_provider=atlas_provider)
    if atlas_provider is None or permute_sample_to_atlas is None:
        return volume_yxz_to_napari_zyx(volume)
    return registration_volume_yxz_to_napari_zyx(
        volume,
        atlas_provider=atlas_provider,
        permute_sample_to_atlas=permute_sample_to_atlas,
    )


def brain_points_to_napari_zyx(
    coords_xyz: np.ndarray,
    *,
    space: ViewSpace,
    volume_shape_yxz: tuple[int, int, int],
    atlas_provider: str | None = None,
    permute_sample_to_atlas: list[int] | None = None,
) -> np.ndarray:
    """Map 1-based point coordinates to Napari ZYX."""
    if space == "atlas":
        if atlas_provider is None:
            msg = "atlas_provider is required for atlas-space Napari display"
            raise ValueError(msg)
        return atlas_points_xyz_to_napari_zyx(
            coords_xyz,
            atlas_provider=atlas_provider,
            volume_shape_yxz=volume_shape_yxz,
        )
    if atlas_provider is None or permute_sample_to_atlas is None:
        return atlas_points_to_napari_zyx(coords_xyz)
    return registration_points_xyz_to_napari_zyx(
        coords_xyz,
        atlas_provider=atlas_provider,
        volume_shape_yxz=volume_shape_yxz,
        permute_sample_to_atlas=permute_sample_to_atlas,
    )


def _contrast_limits(volume: np.ndarray) -> tuple[float, float]:
    positive = volume[volume > 0]
    if positive.size:
        lo, hi = np.percentile(positive, (1.0, 99.5))
    else:
        lo, hi = float(np.min(volume)), float(np.max(volume))
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def _warp_registration_channels_for_atlas_inspect(
    config: BrainPipelineConfig,
    transform_params,
    *,
    expected_shape: tuple[int, int, int],
) -> dict[int, np.ndarray]:
    """Warp registration-grid channels to atlas space when export TIFFs are absent."""
    import shutil

    from lightsuite.export.atlas_space import transform_volume_to_atlas
    from lightsuite.preprocess.checkpoint import RegOptsCheckpoint

    if shutil.which("transformix") is None:
        return {}

    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        return {}

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    channel_paths = {
        int(k): Path(v).expanduser()
        for k, v in (checkpoint.regvolpaths or {}).items()
    }
    if not channel_paths:
        return {}

    permute = transform_params.permute_sample_to_atlas or [1, 2, 3]
    spacing_mm = checkpoint.registres_um * 1e-3
    transformix_root = save_path / "transformix_inspect_temp"
    transformix_root.mkdir(parents=True, exist_ok=True)

    registered_channels: dict[int, np.ndarray] = {}
    for ichan, volpath in sorted(channel_paths.items()):
        if not volpath.is_file():
            continue
        volume = load_registration_volume(volpath)
        registered = transform_volume_to_atlas(
            volume,
            transform_params,
            permute=permute,
            spacing_mm=spacing_mm,
            temp_dir=transformix_root / f"chan_{ichan:02d}",
        )
        if tuple(registered.shape) != expected_shape:
            msg = (
                f"On-the-fly atlas warp for {volpath.name} produced shape "
                f"{registered.shape}, expected {expected_shape}."
            )
            raise ValueError(msg)
        registered_channels[ichan] = registered.astype(np.float32, copy=False)
    return registered_channels


def load_brain_import_inspect_volumes(
    config: BrainPipelineConfig,
    *,
    paths: BrainImportInspectPaths | None = None,
    space: ViewSpace = "atlas",
) -> BrainImportInspectVolumes:
    """Load registered channels, overlays, masks, and point layers for Napari."""
    paths = paths or discover_brain_import_inspect_paths(config, space=space)
    view_space = paths.space if paths is not None else space
    if view_space == "sample":
        from lightsuite.export.brain_sample_space import (
            BrainSampleSpaceInspectPaths,
            load_brain_sample_space_inspect_volumes,
        )

        sample_paths = BrainSampleSpaceInspectPaths(
            volume_registered_dir=paths.volume_registered_dir,
            sample_space_dir=paths.sample_space_dir,
            template_path=paths.template_path,
            annotation_path=paths.annotation_path,
            registered_channels=paths.registered_channels,
            point_npz_paths=paths.point_npz_paths,
            mask_paths=paths.mask_paths,
        )
        return load_brain_sample_space_inspect_volumes(config, paths=sample_paths)
    return _load_brain_import_inspect_volumes_atlas(config, paths=paths)


def _load_brain_import_inspect_volumes_atlas(
    config: BrainPipelineConfig,
    *,
    paths: BrainImportInspectPaths,
) -> BrainImportInspectVolumes:
    save_path = config.sample.save_path.expanduser()
    transform_params = _load_transform_params(save_path)
    expected_shape = tuple(int(v) for v in transform_params.atlassize)

    if paths.template_path is None or paths.annotation_path is None:
        msg = "Atlas inspect paths missing template or annotation."
        raise ValueError(msg)

    template = _load_atlas_volume_yxz(paths.template_path)
    annotation = _load_atlas_volume_yxz(paths.annotation_path)
    if template.shape != expected_shape:
        msg = (
            f"Atlas template shape {template.shape} != transform atlassize {expected_shape}. "
            "Re-run register/export on this sample."
        )
        raise ValueError(msg)
    if annotation.shape != expected_shape:
        msg = (
            f"Atlas annotation shape {annotation.shape} != transform atlassize {expected_shape}."
        )
        raise ValueError(msg)

    registered_channels: dict[int, np.ndarray] = {}
    for ichan, path in paths.registered_channels.items():
        vol = load_registration_volume(path)
        if tuple(vol.shape) != expected_shape:
            msg = f"{path.name} shape {vol.shape} != expected atlas shape {expected_shape}"
            raise ValueError(msg)
        registered_channels[ichan] = vol.astype(np.float32, copy=False)

    channels_warped_on_the_fly = False
    if not registered_channels:
        registered_channels = _warp_registration_channels_for_atlas_inspect(
            config,
            transform_params,
            expected_shape=expected_shape,
        )
        channels_warped_on_the_fly = bool(registered_channels)

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
        channels_warped_on_the_fly=channels_warped_on_the_fly,
    )


def run_brain_inspect_imports(
    config: BrainPipelineConfig,
    *,
    space: ViewSpace = "atlas",
    headless: bool = False,
) -> BrainImportInspectPaths:
    """Open Napari to QC registered channels and imported annotations."""
    paths = discover_brain_import_inspect_paths(config, space=space)
    if headless:
        load_brain_import_inspect_volumes(config, paths=paths, space=space)
        return paths

    try:
        import napari
        from napari.utils.notifications import show_info
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    volumes = load_brain_import_inspect_volumes(config, paths=paths, space=space)
    atlas_provider = atlas_display_provider_from_config(config.atlas)
    permute_sample_to_atlas: list[int] | None = None
    if space == "sample":
        transform_params = _load_transform_params(config.sample.save_path.expanduser())
        permute_sample_to_atlas = transform_params.permute_sample_to_atlas or [1, 2, 3]
    reference_shape = (
        volumes.template.shape
        if volumes.template is not None
        else next(iter(volumes.registered_channels.values())).shape
        if volumes.registered_channels
        else next(iter(volumes.mask_layers.values())).shape
    )

    if space == "sample":
        title = f"LightSuite brain import QC — {config.sample.name} (sample, 20 µm registration grid)"
        template_name = "atlas template (warped)"
        annotation_name = "atlas annotation (warped)"
        channel_suffix = "(sample warped)"
    else:
        title = f"LightSuite brain import QC — {config.sample.name} (atlas)"
        template_name = "atlas template"
        annotation_name = "atlas annotation"
        channel_suffix = "(sample warped)"

    viewer = napari.Viewer(title=title)

    if volumes.template is not None:
        viewer.add_image(
            brain_volume_to_napari_zyx(
                volumes.template,
                space=space,
                atlas_provider=atlas_provider,
                permute_sample_to_atlas=permute_sample_to_atlas,
            ),
            name=template_name,
            colormap="gray",
            blending="opaque",
            contrast_limits=_contrast_limits(volumes.template),
        )

    if volumes.annotation is not None:
        viewer.add_labels(
            brain_volume_to_napari_zyx(
                volumes.annotation,
                space=space,
                atlas_provider=atlas_provider,
                permute_sample_to_atlas=permute_sample_to_atlas,
            ).astype(np.int64, copy=False),
            name=annotation_name,
            opacity=0.45,
        )

    channel_cmaps = ["magenta", "cyan", "yellow", "red"]
    for idx, (ichan, vol) in enumerate(sorted(volumes.registered_channels.items())):
        cmap = channel_cmaps[idx % len(channel_cmaps)]
        viewer.add_image(
            brain_volume_to_napari_zyx(
                vol,
                space=space,
                atlas_provider=atlas_provider,
                permute_sample_to_atlas=permute_sample_to_atlas,
            ),
            name=f"channel {ichan} {channel_suffix}",
            colormap=cmap,
            blending="additive",
            opacity=0.55,
            contrast_limits=_contrast_limits(vol),
        )

    for label, mask in volumes.mask_layers.items():
        viewer.add_image(
            brain_volume_to_napari_zyx(
                mask,
                space=space,
                atlas_provider=atlas_provider,
                permute_sample_to_atlas=permute_sample_to_atlas,
            ),
            name=f"mask: {label}",
            colormap="red",
            blending="additive",
            opacity=0.35,
            contrast_limits=(0.0, 1.0),
        )

    for label, coords in volumes.point_layers.items():
        napari_pts = brain_points_to_napari_zyx(
            coords,
            space=space,
            atlas_provider=atlas_provider,
            volume_shape_yxz=reference_shape,
            permute_sample_to_atlas=permute_sample_to_atlas,
        )
        viewer.add_points(
            napari_pts,
            name=f"points: {label}",
            size=3,
            face_color="red",
            border_color="white",
        )

    summary_path = paths.volume_registered_dir / "import_annotations_summary.json"
    space_note = "sample-space " if space == "sample" else ""
    atlas_overlay_note = ""
    if space == "sample" and volumes.annotation is None:
        atlas_overlay_note = (
            " No warped atlas annotation — run "
            "'lightsuite brain export -c <config> --space sample' first."
        )
    warp_note = ""
    if space == "atlas" and volumes.channels_warped_on_the_fly:
        warp_note = (
            " Sample channels warped on the fly from regopts.json "
            "(run 'brain export --space atlas --save-volume' to cache)."
        )
    channel_info = (
        f"Loaded {len(volumes.registered_channels)} {space_note}channel(s), "
        f"{len(volumes.point_layers)} point layer(s), "
        f"{len(volumes.mask_layers)} mask layer(s)."
        f"{atlas_overlay_note}{warp_note}"
    )
    if summary_path.is_file():
        show_info(f"{channel_info} Summary: {summary_path.name}")
    else:
        show_info(channel_info)
    napari.run()
    return paths
