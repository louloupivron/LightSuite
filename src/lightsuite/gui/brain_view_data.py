"""Load volumes and paths for brain registration review (view-registration stage)."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from lightsuite.analysis.counts import SAMPLE_POINTS_KEY, load_atlas_points
from lightsuite.analysis.division_map import ensure_division_map
from lightsuite.atlas.display import (
    atlas_points_xyz_to_napari_zyx,
    atlas_volume_yxz_to_napari_zyx,
    registration_points_xyz_to_napari_zyx,
    registration_volume_yxz_to_napari_zyx,
)
from lightsuite.atlas.registry import (
    resolve_brain_atlas_from_config,
    resolve_brain_atlas_with_config,
)
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.export.brain_sample_space import (
    BOUNDARY_IN_SAMPLE,
    DIVISION_IN_SAMPLE,
    atlas_volumes_permuted_on_disk,
    discover_brain_sample_space_inspect_paths,
    load_sample_space_atlas_volume,
)
from lightsuite.import_.adapters import load_annotation
from lightsuite.import_.brain_import import (
    _registration_shape_native_yxz,
    _resample_mask_native_to_registration,
)
from lightsuite.import_.models import ImportedMask, ImportedPoints
from lightsuite.import_.orchestrator import slug_for_label
from lightsuite.import_.sample_reference import load_sample_reference
from lightsuite.import_.transform import sample_points_to_registration_voxels
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.brain_paths import (
    REGISTRATION_DIAGNOSTICS_FILENAME,
    TRANSFORMIX_VIEW_TEMP,
    brain_volume_registered_dir,
    brain_work_dir,
    cleanup_brain_work,
    iter_brain_import_paths,
    resolve_brain_qc_file,
)
from lightsuite.registration.points import volume_indices_to_cloud_xyz
from lightsuite.registration.register_diagnostics import RegistrationDiagnostics
from lightsuite.registration.volume import (
    load_permuted_registration_volume,
    load_registration_volume,
    permute_brain_volume,
)

ViewSpace = Literal["atlas", "sample"]


@dataclass(frozen=True)
class BrainViewPaths:
    """Resolved inputs for brain registration review."""

    volume_registered_dir: Path
    template_path: Path | None
    annotation_path: Path | None
    boundary_path: Path | None
    division_labels_path: Path | None
    division_legend_path: Path | None
    registered_channels: dict[int, Path] = field(default_factory=dict)
    point_npz_paths: dict[str, Path] = field(default_factory=dict)
    mask_paths: dict[str, Path] = field(default_factory=dict)
    space: ViewSpace = "sample"
    sample_space_dir: Path | None = None
    diagnostics_path: Path | None = None


@dataclass
class BrainViewVolumes:
    template: np.ndarray | None
    annotation: np.ndarray | None
    boundary: np.ndarray | None
    division_labels: np.ndarray | None
    division_legend: pd.DataFrame | None
    registered_channels: dict[int, np.ndarray]
    point_layers: dict[str, np.ndarray]
    mask_layers: dict[str, np.ndarray]
    resampled_roi_masks: dict[str, np.ndarray] = field(default_factory=dict)
    resampled_roi_points: dict[str, np.ndarray] = field(default_factory=dict)
    multires_roi_channels: dict[str, np.ndarray] = field(default_factory=dict)
    channels_warped_on_the_fly: bool = False
    diagnostics: RegistrationDiagnostics | None = None


# Backward-compatible aliases used by sample-space export loaders.
BrainImportInspectPaths = BrainViewPaths
BrainImportInspectVolumes = BrainViewVolumes


def _volume_registered_dir(config: BrainPipelineConfig) -> Path:
    return brain_volume_registered_dir(config.sample.save_path.expanduser())


def _label_from_stem(stem: str, suffix: str) -> str:
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def discover_brain_view_paths(
    config: BrainPipelineConfig,
    *,
    space: ViewSpace = "sample",
) -> BrainViewPaths:
    """Discover layers for registration review in sample or atlas space."""
    if space == "sample":
        return _discover_brain_view_paths_sample(config)
    return _discover_brain_view_paths_atlas(config)


discover_brain_import_inspect_paths = discover_brain_view_paths


def _discover_brain_view_paths_sample(config: BrainPipelineConfig) -> BrainViewPaths:
    sample_paths = discover_brain_sample_space_inspect_paths(config)
    save_path = config.sample.save_path.expanduser()
    division_labels_path: Path | None = None
    division_legend_path: Path | None = None
    boundary_path: Path | None = None
    if sample_paths.sample_space_dir is not None:
        ss = sample_paths.sample_space_dir
        candidate = ss / DIVISION_IN_SAMPLE
        if candidate.is_file():
            division_labels_path = candidate.resolve()
        candidate = ss / BOUNDARY_IN_SAMPLE
        if candidate.is_file():
            boundary_path = candidate.resolve()
    try:
        transform_params = _load_transform_params(save_path)
        atlas = resolve_brain_atlas_with_config(transform_params.brain_atlas, config.atlas)
        division = ensure_division_map(atlas)
        if division_labels_path is None:
            division_labels_path = division.paths.labels_tiff.resolve()
        division_legend_path = division.paths.legend_csv.resolve()
    except (FileNotFoundError, OSError, ValueError):
        pass

    diagnostics_path = resolve_brain_qc_file(save_path, REGISTRATION_DIAGNOSTICS_FILENAME)
    if not diagnostics_path.is_file():
        diagnostics_path = None

    return BrainViewPaths(
        volume_registered_dir=sample_paths.volume_registered_dir,
        template_path=sample_paths.template_path,
        annotation_path=sample_paths.annotation_path,
        boundary_path=boundary_path,
        division_labels_path=division_labels_path,
        division_legend_path=division_legend_path,
        registered_channels=sample_paths.registered_channels,
        point_npz_paths=sample_paths.point_npz_paths,
        mask_paths=sample_paths.mask_paths,
        space="sample",
        sample_space_dir=sample_paths.sample_space_dir,
        diagnostics_path=diagnostics_path,
    )


def _can_warp_registration_channels_to_atlas(config: BrainPipelineConfig) -> bool:
    """True when registration channels can be warped to atlas space for QC on demand."""
    import shutil

    if shutil.which("transformix") is None:
        return False
    save_path = config.sample.save_path.expanduser()
    if not (save_path / "transform_params.json").is_file():
        return False
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        return False
    from lightsuite.preprocess.checkpoint import RegOptsCheckpoint

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    if not checkpoint.regvolpaths:
        return False
    return any(Path(v).expanduser().is_file() for v in checkpoint.regvolpaths.values())


def _discover_brain_view_paths_atlas(config: BrainPipelineConfig) -> BrainViewPaths:
    save_path = config.sample.save_path.expanduser()
    vr = _volume_registered_dir(config)

    transform_params_path = save_path / "transform_params.json"
    if not transform_params_path.is_file():
        msg = f"Missing {transform_params_path}. Run 'lightsuite brain register' first."
        raise FileNotFoundError(msg)

    atlas = resolve_brain_atlas_from_config(config.atlas)
    transform_params = _load_transform_params(save_path)
    atlas_resolved = resolve_brain_atlas_with_config(transform_params.brain_atlas, config.atlas)
    division = ensure_division_map(atlas_resolved)

    registered_channels: dict[int, Path] = {}
    if vr.is_dir():
        for path in sorted(vr.glob("chan_*_registered_atlas.tif")):
            match = re.fullmatch(r"chan_(\d+)_registered_atlas\.tif", path.name)
            if match is None:
                continue
            registered_channels[int(match.group(1))] = path.resolve()

    point_npz_paths: dict[str, Path] = {}
    mask_paths: dict[str, Path] = {}
    for path in iter_brain_import_paths(save_path, "*_atlas_coords.npz"):
        label = _label_from_stem(path.stem, "_atlas_coords")
        point_npz_paths[label] = path

    for path in iter_brain_import_paths(save_path, "*_registered_atlas.tif"):
        if path.name.startswith("chan_"):
            continue
        label = _label_from_stem(path.stem, "_registered_atlas")
        mask_paths[label] = path

    has_on_disk_exports = bool(registered_channels or point_npz_paths or mask_paths)
    if not has_on_disk_exports and not _can_warp_registration_channels_to_atlas(config):
        msg = (
            f"No atlas-space exports found under {vr}. "
            "Run export with Atlas enabled, or keep registration channel TIFFs "
            "and transformix available for on-the-fly warping."
        )
        raise FileNotFoundError(msg)

    diagnostics_path = resolve_brain_qc_file(save_path, REGISTRATION_DIAGNOSTICS_FILENAME)
    return BrainViewPaths(
        volume_registered_dir=vr.resolve() if vr.is_dir() else save_path,
        template_path=atlas.template_path.resolve(),
        annotation_path=atlas.annotation_path.resolve(),
        boundary_path=None,
        division_labels_path=division.paths.labels_tiff.resolve(),
        division_legend_path=division.paths.legend_csv.resolve(),
        registered_channels=registered_channels,
        point_npz_paths=point_npz_paths,
        mask_paths=mask_paths,
        space="atlas",
        diagnostics_path=diagnostics_path if diagnostics_path.is_file() else None,
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


def label_from_stem(stem: str, suffix: str) -> str:
    return _label_from_stem(stem, suffix)


load_points_npz = _load_points_npz
load_points_csv = _load_points_csv


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


def contrast_limits(volume: np.ndarray) -> tuple[float, float]:
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
    transform_params: Any,
    *,
    expected_shape: tuple[int, int, int],
) -> dict[int, np.ndarray]:
    """Warp registration-grid channels to atlas space when export TIFFs are absent."""
    import shutil

    from lightsuite.export.atlas_space import transform_volume_to_atlas

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
    transformix_root = brain_work_dir(save_path, TRANSFORMIX_VIEW_TEMP)

    registered_channels: dict[int, np.ndarray] = {}
    try:
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
    finally:
        cleanup_brain_work(save_path, TRANSFORMIX_VIEW_TEMP)
    return registered_channels


def _already_imported(output_dir: Path, slug: str) -> bool:
    return (output_dir / f"{slug}_in_sample_20um.tif").is_file() or (
        output_dir / f"{slug}_sample_coords.npz"
    ).is_file()


def load_resampled_config_annotations(
    config: BrainPipelineConfig,
    *,
    expected_shape: tuple[int, int, int],
    output_dir: Path,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Preview native / overview-grid ROIs from import.annotations on the registration grid."""
    import_cfg = config.import_config
    if import_cfg is None or not import_cfg.annotations:
        return {}, {}

    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        return {}, {}

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    transform_params = _load_transform_params(save_path)
    permute = transform_params.permute_sample_to_atlas or [1, 2, 3]
    voxel_um = [float(v) for v in checkpoint.voxel_um]
    registres_um = float(checkpoint.registres_um)
    sample_reference = load_sample_reference(save_path)

    mask_layers: dict[str, np.ndarray] = {}
    point_layers: dict[str, np.ndarray] = {}

    for spec in import_cfg.annotations:
        label = spec.label or Path(spec.path).stem
        slug = slug_for_label(label)
        if _already_imported(output_dir, slug):
            continue
        try:
            loaded = load_annotation(spec)
        except (OSError, ValueError, KeyError):
            continue

        if isinstance(loaded, ImportedMask):
            reg_shape = _registration_shape_native_yxz(loaded.volume.shape, voxel_um, registres_um)
            mask_reg = _resample_mask_native_to_registration(
                loaded.volume,
                target_shape_yxz=reg_shape,
                native_voxel_um=voxel_um,
                registres_um=registres_um,
            )
            mask_sample = permute_brain_volume(mask_reg, permute)
            if mask_sample.shape != expected_shape:
                mask_sample = _resample_mask_native_to_registration(
                    mask_sample,
                    target_shape_yxz=expected_shape,
                    native_voxel_um=[registres_um, registres_um, registres_um],
                    registres_um=registres_um,
                )
            mask_layers[f"ROI: {label}"] = mask_sample.astype(np.float32, copy=False)
        elif isinstance(loaded, ImportedPoints):
            prepared = loaded
            if loaded.coordinates.size:
                from lightsuite.import_.adapters import prepare_points_for_sample

                prepared = prepare_points_for_sample(loaded, sample_reference)
            reg_yxz = sample_points_to_registration_voxels(
                prepared.coordinates,
                transform_params,
                registres_um=registres_um,
                content_crop_start=checkpoint.content_crop_start,
            )
            point_layers[f"ROI: {label}"] = volume_indices_to_cloud_xyz(reg_yxz)

    return mask_layers, point_layers


def load_brain_view_volumes(
    config: BrainPipelineConfig,
    *,
    paths: BrainViewPaths | None = None,
    space: ViewSpace = "sample",
    load_multires_roi: bool = True,
) -> BrainViewVolumes:
    """Load registration review layers for Napari."""
    paths = paths or discover_brain_view_paths(config, space=space)
    if paths.space == "sample":
        return _load_brain_view_volumes_sample(
            config,
            paths=paths,
            load_multires_roi=load_multires_roi,
        )
    return _load_brain_view_volumes_atlas(config, paths=paths)


load_brain_import_inspect_volumes = load_brain_view_volumes


def _load_brain_view_volumes_sample(
    config: BrainPipelineConfig,
    *,
    paths: BrainViewPaths,
    load_multires_roi: bool = True,
) -> BrainViewVolumes:
    save_path = config.sample.save_path.expanduser()
    transform_params = _load_transform_params(save_path)
    expected_shape = tuple(int(v) for v in transform_params.regvolsize)
    permute = transform_params.permute_sample_to_atlas or [1, 2, 3]
    atlas_permuted = atlas_volumes_permuted_on_disk(paths.sample_space_dir)

    registered_channels: dict[int, np.ndarray] = {}
    for ichan, path in paths.registered_channels.items():
        vol = load_permuted_registration_volume(path, permute)
        if tuple(vol.shape) != expected_shape:
            msg = f"{path.name} shape {vol.shape} != expected registration shape {expected_shape}"
            raise ValueError(msg)
        registered_channels[ichan] = vol

    template: np.ndarray | None = None
    if paths.template_path is not None and paths.template_path.is_file():
        template = load_sample_space_atlas_volume(
            paths.template_path,
            permute,
            permuted_on_disk=atlas_permuted,
        )
        if tuple(template.shape) != expected_shape:
            msg = (
                f"Template shape {template.shape} != expected registration shape {expected_shape}. "
                "Re-run sample-space export."
            )
            raise ValueError(msg)

    annotation: np.ndarray | None = None
    if paths.annotation_path is not None and paths.annotation_path.is_file():
        annotation = load_sample_space_atlas_volume(
            paths.annotation_path,
            permute,
            permuted_on_disk=atlas_permuted,
        )

    boundary: np.ndarray | None = None
    if paths.boundary_path is not None and paths.boundary_path.is_file():
        boundary = load_sample_space_atlas_volume(
            paths.boundary_path,
            permute,
            permuted_on_disk=atlas_permuted,
        )

    division_labels: np.ndarray | None = None
    division_legend: pd.DataFrame | None = None
    if paths.division_labels_path is not None and paths.division_labels_path.is_file():
        if paths.sample_space_dir is not None and paths.division_labels_path.parent == paths.sample_space_dir:
            division_labels = load_sample_space_atlas_volume(
                paths.division_labels_path,
                permute,
                permuted_on_disk=atlas_permuted,
            ).astype(np.int32, copy=False)
        else:
            division_labels = load_registration_volume(paths.division_labels_path).astype(
                np.int32, copy=False
            )
        if paths.division_legend_path is not None and paths.division_legend_path.is_file():
            division_legend = pd.read_csv(paths.division_legend_path)

    mask_layers: dict[str, np.ndarray] = {}
    for label, path in paths.mask_paths.items():
        vol = load_registration_volume(path)
        mask_layers[label] = vol.astype(np.float32, copy=False)

    point_layers: dict[str, np.ndarray] = {}
    for label, npz_path in paths.point_npz_paths.items():
        point_layers[label] = load_atlas_points(npz_path, key=SAMPLE_POINTS_KEY)

    resampled_masks, resampled_points = load_resampled_config_annotations(
        config,
        expected_shape=expected_shape,
        output_dir=paths.volume_registered_dir,
    )

    from lightsuite.gui.brain_multires_link import load_multires_roi_channels_on_registration_grid

    multires_roi_channels: dict[str, np.ndarray] = {}
    if load_multires_roi and config.multires_link is not None:
        regopts_path = save_path / "regopts.json"
        if regopts_path.is_file():
            regopts = RegOptsCheckpoint.load(regopts_path)
            multires_roi_channels = load_multires_roi_channels_on_registration_grid(
                config,
                checkpoint=regopts,
                transform_params=transform_params,
                expected_shape=expected_shape,
            )

    diagnostics: RegistrationDiagnostics | None = None
    if paths.diagnostics_path is not None and paths.diagnostics_path.is_file():
        diagnostics = RegistrationDiagnostics.load(paths.diagnostics_path)

    return BrainViewVolumes(
        template=template,
        annotation=annotation,
        boundary=boundary,
        division_labels=division_labels,
        division_legend=division_legend,
        registered_channels=registered_channels,
        point_layers=point_layers,
        mask_layers=mask_layers,
        resampled_roi_masks=resampled_masks,
        resampled_roi_points=resampled_points,
        multires_roi_channels=multires_roi_channels,
        diagnostics=diagnostics,
    )


def _load_brain_view_volumes_atlas(
    config: BrainPipelineConfig,
    *,
    paths: BrainViewPaths,
) -> BrainViewVolumes:
    save_path = config.sample.save_path.expanduser()
    transform_params = _load_transform_params(save_path)
    expected_shape = tuple(int(v) for v in transform_params.atlassize)

    if paths.template_path is None or paths.annotation_path is None:
        msg = "Atlas view paths missing template or annotation."
        raise ValueError(msg)

    template = _load_atlas_volume_yxz(paths.template_path)
    annotation = _load_atlas_volume_yxz(paths.annotation_path)

    registered_channels: dict[int, np.ndarray] = {}
    for ichan, path in paths.registered_channels.items():
        vol = load_registration_volume(path)
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
        mask_layers[label] = load_registration_volume(path).astype(np.float32, copy=False)

    point_layers: dict[str, np.ndarray] = {}
    for label, npz_path in paths.point_npz_paths.items():
        point_layers[label] = _load_points_npz(npz_path)
    save_path = paths.volume_registered_dir.parent
    for path in iter_brain_import_paths(save_path, "*_atlas_coords.csv"):
        label = _label_from_stem(path.stem, "_atlas_coords")
        if label in point_layers:
            continue
        point_layers[label] = _load_points_csv(path)

    division_labels: np.ndarray | None = None
    division_legend: pd.DataFrame | None = None
    if paths.division_labels_path is not None and paths.division_labels_path.is_file():
        import tifffile

        division_labels = np.asarray(tifffile.imread(paths.division_labels_path), dtype=np.int32)
        if paths.division_legend_path is not None and paths.division_legend_path.is_file():
            division_legend = pd.read_csv(paths.division_legend_path)

    diagnostics: RegistrationDiagnostics | None = None
    if paths.diagnostics_path is not None and paths.diagnostics_path.is_file():
        diagnostics = RegistrationDiagnostics.load(paths.diagnostics_path)

    return BrainViewVolumes(
        template=template,
        annotation=annotation,
        boundary=None,
        division_labels=division_labels,
        division_legend=division_legend,
        registered_channels=registered_channels,
        point_layers=point_layers,
        mask_layers=mask_layers,
        channels_warped_on_the_fly=channels_warped_on_the_fly,
        diagnostics=diagnostics,
    )


def brain_view_spaces_available(config: BrainPipelineConfig) -> dict[ViewSpace, bool]:
    """Return which coordinate spaces can be opened in view-registration."""
    available: dict[ViewSpace, bool] = {"sample": False, "atlas": False}
    save_path = config.sample.save_path.expanduser()
    vr = save_path / "volume_registered"

    if (save_path / "transform_params.json").is_file():
        try:
            resolve_brain_atlas_from_config(config.atlas)
        except (FileNotFoundError, OSError, ValueError):
            pass
        else:
            if vr.is_dir() and any(vr.glob("chan_*_registered_atlas.tif")):
                available["atlas"] = True
            elif _can_warp_registration_channels_to_atlas(config):
                available["atlas"] = True

    try:
        sample_paths = discover_brain_sample_space_inspect_paths(config)
        available["sample"] = bool(sample_paths.registered_channels)
        if sample_paths.sample_space_dir is not None and sample_paths.sample_space_dir.is_dir():
            available["sample"] = True
    except FileNotFoundError:
        pass

    return available


def resolve_brain_view_space(
    config: BrainPipelineConfig,
    *,
    preferred: ViewSpace = "sample",
) -> ViewSpace:
    """Pick sample- or atlas-space review when the preferred export is missing."""
    available = brain_view_spaces_available(config)
    if not any(available.values()):
        msg = (
            f"Missing view-registration exports under {config.sample.save_path.expanduser() / 'volume_registered'}. "
            "Run 'lightsuite brain export' (atlas and/or sample space) first."
        )
        raise FileNotFoundError(msg)
    if available.get(preferred, False):
        return preferred
    for space in ("sample", "atlas"):
        if available.get(space, False):
            return space
    msg = "No view-registration exports are available."
    raise FileNotFoundError(msg)


def brain_view_load_summary(
    paths: BrainViewPaths,
    volumes: BrainViewVolumes,
    *,
    space: ViewSpace,
    config: BrainPipelineConfig,
) -> str:
    """Short status line after loading a brain view-registration space."""
    del paths, config
    space_note = "sample-space " if space == "sample" else "atlas (template) space "
    extras: list[str] = []
    if volumes.template is not None:
        extras.append("template")
    if volumes.annotation is not None:
        extras.append("annotation")
    if volumes.channels_warped_on_the_fly:
        extras.append("channels warped on the fly")
    extra_text = f" ({', '.join(extras)})" if extras else ""
    return (
        f"Loaded {space_note}with {len(volumes.registered_channels)} channel(s), "
        f"{len(volumes.point_layers)} import point layer(s), "
        f"{len(volumes.mask_layers)} import mask layer(s){extra_text}."
    )


def format_diagnostics_summary(diagnostics: RegistrationDiagnostics | None) -> str:
    if diagnostics is None:
        return "No registration_diagnostics.json found."
    lines = [
        f"Status: {diagnostics.status.upper()} — {diagnostics.status_message}",
        (
            f"Affine landmarks: median {diagnostics.affine_median_error_vox:.1f} vox, "
            f"p95 {diagnostics.affine_p95_error_vox:.1f} vox"
        ),
        f"Annotation overlap: {diagnostics.annotation_label_voxels:,} label voxels",
    ]
    for warning in diagnostics.warnings:
        lines.append(f"Warning: {warning}")
    return "\n".join(lines)
