"""Brain pipeline sample-space export (atlas warped onto registration grid)."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from rich.console import Console

from lightsuite.analysis.brain_runner import maybe_write_brain_sample_region_stats
from lightsuite.analysis.division_map import ensure_division_map
from lightsuite.analysis.ontology import RegionTable, load_region_table
from lightsuite.analysis.region_stats import (
    concat_tidy,
    parcellation_result_to_tidy,
    write_region_stats_csv,
)
from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import (
    resolve_brain_atlas_from_config,
    uses_ccf_id_parcellation,
)
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.parcellation import (
    compute_allen_parcellation,
    compute_perens_parcellation,
)
from lightsuite.export.sample_space import transform_atlas_volume_to_sample
from lightsuite.io.tiff_write import save_registration_volume
from lightsuite.preprocess.checkpoint import RegOptsCheckpoint
from lightsuite.registration.brain_register import TransformParamsCheckpoint
from lightsuite.registration.plots import boundary_volume_from_annotation
from lightsuite.registration.volume import (
    load_permuted_registration_volume,
    load_registration_volume,
    permute_brain_volume,
)

console = Console()

SAMPLE_SPACE_SUBDIR = "sample_space"
ANNOTATION_IN_SAMPLE = "annotation_in_sample_20um.tif"
TEMPLATE_IN_SAMPLE = "template_in_sample_20um.tif"
BOUNDARY_IN_SAMPLE = "annotation_boundary_in_sample_20um.tif"
DIVISION_IN_SAMPLE = "division_labels_in_sample_20um.tif"
MANIFEST_NAME = "sample_space_manifest.json"


def sample_space_dir(save_path: Path) -> Path:
    return save_path / "volume_registered" / SAMPLE_SPACE_SUBDIR


def atlas_volumes_permuted_on_disk(sample_dir: Path | None) -> bool:
    """Whether exported atlas TIFFs are already on the permuted registration grid."""
    if sample_dir is None:
        return False
    manifest_path = sample_dir / MANIFEST_NAME
    if not manifest_path.is_file():
        return False
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    return bool(manifest.get("atlas_volumes_permuted", False))


def load_sample_space_atlas_volume(
    path: Path,
    permute_sample_to_atlas: list[int],
    *,
    permuted_on_disk: bool,
) -> np.ndarray:
    """Load a warped atlas TIFF for overlay with permuted registration channels."""
    vol = load_registration_volume(path).astype(np.float32, copy=False)
    if permuted_on_disk:
        return vol
    return permute_brain_volume(vol, permute_sample_to_atlas)


def export_brain_sample_space(
    config: BrainPipelineConfig,
    *,
    transform_params: TransformParamsCheckpoint,
    checkpoint: RegOptsCheckpoint,
    write_csv: bool,
    save_volume: bool,
    region_table: RegionTable | None,
) -> tuple[Path, dict[int, Path] | None]:
    """Warp atlas labels onto the registration grid and write sample-space stats."""
    save_path = config.sample.save_path.expanduser()
    out_dir = sample_space_dir(save_path)
    if save_volume or write_csv:
        out_dir.mkdir(parents=True, exist_ok=True)

    atlas = resolve_brain_atlas_from_config(config.atlas)
    spacing_mm = checkpoint.registres_um * 1e-3
    registres_um = float(checkpoint.registres_um)
    transformix_root = save_path / "transformix_sample_export_temp"
    transformix_root.mkdir(parents=True, exist_ok=True)

    console.print("Warping atlas volumes to sample registration grid...")
    t0 = time.perf_counter()

    av_native = load_atlas_volume(atlas.annotation_path)
    annotation_sample = transform_atlas_volume_to_sample(
        av_native,
        transform_params,
        save_path=save_path,
        spacing_mm=spacing_mm,
        temp_dir=transformix_root / "annotation",
        nearest=True,
    )

    tv_native = load_atlas_volume(atlas.template_path)
    template_sample = transform_atlas_volume_to_sample(
        tv_native,
        transform_params,
        save_path=save_path,
        spacing_mm=spacing_mm,
        temp_dir=transformix_root / "template",
        nearest=False,
    )

    if atlas.boundary_path is not None and atlas.boundary_path.is_file():
        boundary_native = load_atlas_volume(atlas.boundary_path)
        boundary_sample = transform_atlas_volume_to_sample(
            boundary_native,
            transform_params,
            save_path=save_path,
            spacing_mm=spacing_mm,
            temp_dir=transformix_root / "boundary",
            nearest=True,
        )
    else:
        boundary_sample = boundary_volume_from_annotation(annotation_sample)

    division_result = ensure_division_map(atlas)
    division_sample = transform_atlas_volume_to_sample(
        division_result.labels.astype(np.float32),
        transform_params,
        save_path=save_path,
        spacing_mm=spacing_mm,
        temp_dir=transformix_root / "division",
        nearest=True,
    )

    if save_volume:
        save_registration_volume(
            np.rint(annotation_sample).astype(np.uint16),
            out_dir / ANNOTATION_IN_SAMPLE,
        )
        save_registration_volume(
            np.clip(template_sample, 0, np.iinfo(np.uint16).max).astype(np.uint16),
            out_dir / TEMPLATE_IN_SAMPLE,
        )
        save_registration_volume(
            (boundary_sample > 0).astype(np.uint16),
            out_dir / BOUNDARY_IN_SAMPLE,
        )
        save_registration_volume(
            np.rint(division_sample).astype(np.uint16),
            out_dir / DIVISION_IN_SAMPLE,
        )

    channel_paths = {int(k): Path(v) for k, v in (checkpoint.regvolpaths or {}).items()}
    manifest = {
        "space": "sample",
        "grid": "registration",
        "shape_yxz": list(transform_params.regvolsize),
        "voxel_um": [registres_um, registres_um, registres_um],
        "straightened": False,
        "permute_sample_to_atlas": transform_params.permute_sample_to_atlas,
        "registration_canvas": transform_params.registration_canvas,
        "channel_paths": {str(k): str(v) for k, v in channel_paths.items()},
        "annotation_path": str(out_dir / ANNOTATION_IN_SAMPLE),
        "template_path": str(out_dir / TEMPLATE_IN_SAMPLE),
        "boundary_path": str(out_dir / BOUNDARY_IN_SAMPLE),
        "division_labels_path": str(out_dir / DIVISION_IN_SAMPLE),
        "atlas_volumes_permuted": True,
    }
    manifest_path = out_dir / MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    region_stats_paths: dict[int, Path] | None = None
    tidy_frames: list[pd.DataFrame] = []
    if write_csv and atlas.supports_parcellation:
        ann_int = np.rint(annotation_sample).astype(np.int32)
        for ichan, volpath in sorted(channel_paths.items()):
            vol = load_permuted_registration_volume(
                volpath,
                transform_params.permute_sample_to_atlas,
            )
            if uses_ccf_id_parcellation(atlas):
                result = compute_perens_parcellation(
                    vol,
                    atlas,
                    transform_params.atlas_resolution_um,
                    annotation=ann_int,
                    voxel_um=registres_um,
                )
            else:
                result = compute_allen_parcellation(
                    vol,
                    atlas,
                    transform_params.atlas_resolution_um,
                    annotation=ann_int,
                    voxel_um=registres_um,
                )
            tidy = parcellation_result_to_tidy(
                result,
                region_table,
                sample=config.sample.name,
                channel=ichan,
                atlas=atlas.brain_atlas,
                intensity_metrics=config.analysis.intensity_metrics,
            )
            tidy_path = out_dir / f"chan{ichan:02d}_region_stats_sample.csv"
            write_region_stats_csv(tidy_path, tidy)
            if region_stats_paths is None:
                region_stats_paths = {}
            region_stats_paths[ichan] = tidy_path
            tidy_frames.append(tidy)

        if tidy_frames:
            register_path = save_path / "volume_registered"
            maybe_write_brain_sample_region_stats(
                config,
                tidy_frames=tidy_frames,
                annotation_sample=np.rint(annotation_sample).astype(np.int32),
                region_table=region_table,
                atlas=atlas,
                registres_um=registres_um,
                transform_params=transform_params,
                register_path=register_path,
            )

    console.print(
        f"Sample-space export done in {time.perf_counter() - t0:.1f}s under {out_dir}"
    )
    return manifest_path, region_stats_paths


@dataclass(frozen=True)
class BrainSampleSpaceInspectPaths:
    """Resolved sample-space inputs for Napari import QC."""

    volume_registered_dir: Path
    sample_space_dir: Path | None
    template_path: Path | None
    annotation_path: Path | None
    registered_channels: dict[int, Path] = field(default_factory=dict)
    point_npz_paths: dict[str, Path] = field(default_factory=dict)
    mask_paths: dict[str, Path] = field(default_factory=dict)


def _label_from_stem(stem: str, suffix: str) -> str:
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def discover_brain_sample_space_inspect_paths(
    config: BrainPipelineConfig,
) -> BrainSampleSpaceInspectPaths:
    """Discover registration-grid channels and imported sample-space annotations."""
    from lightsuite.preprocess.checkpoint import RegOptsCheckpoint

    save_path = config.sample.save_path.expanduser()
    vr = save_path / "volume_registered"
    out_dir = sample_space_dir(save_path)

    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run 'lightsuite brain preprocess' first."
        raise FileNotFoundError(msg)

    checkpoint = RegOptsCheckpoint.load(regopts_path)
    registered_channels: dict[int, Path] = {}
    for key, path_str in (checkpoint.regvolpaths or {}).items():
        path = Path(path_str).expanduser()
        if path.is_file():
            registered_channels[int(key)] = path.resolve()

    point_npz_paths: dict[str, Path] = {}
    mask_paths: dict[str, Path] = {}
    if vr.is_dir():
        for path in sorted(vr.glob("*_sample_coords.npz")):
            label = _label_from_stem(path.stem, "_sample_coords")
            point_npz_paths[label] = path.resolve()
        for path in sorted(vr.glob("*_in_sample_20um.tif")):
            label = _label_from_stem(path.stem, "_in_sample_20um")
            mask_paths[label] = path.resolve()

    template_path: Path | None = None
    annotation_path: Path | None = None
    sample_space_dir_resolved: Path | None = out_dir.resolve() if out_dir.is_dir() else None
    if sample_space_dir_resolved is not None:
        manifest_path = sample_space_dir_resolved / MANIFEST_NAME
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if manifest.get("template_path"):
                template_path = Path(manifest["template_path"]).expanduser().resolve()
            if manifest.get("annotation_path"):
                annotation_path = Path(manifest["annotation_path"]).expanduser().resolve()
        if template_path is None:
            candidate = sample_space_dir_resolved / TEMPLATE_IN_SAMPLE
            if candidate.is_file():
                template_path = candidate.resolve()
        if annotation_path is None:
            candidate = sample_space_dir_resolved / ANNOTATION_IN_SAMPLE
            if candidate.is_file():
                annotation_path = candidate.resolve()

    if not registered_channels and not point_npz_paths and not mask_paths:
        msg = (
            "No inspectable sample-space layers found. Run "
            "'lightsuite brain import-annotations' and/or ensure regopts.json "
            "regvolpaths point at chan_*_sample_register_*um.tif volumes."
        )
        raise FileNotFoundError(msg)

    return BrainSampleSpaceInspectPaths(
        volume_registered_dir=vr.resolve() if vr.is_dir() else save_path / "volume_registered",
        sample_space_dir=sample_space_dir_resolved,
        template_path=template_path,
        annotation_path=annotation_path,
        registered_channels=registered_channels,
        point_npz_paths=point_npz_paths,
        mask_paths=mask_paths,
    )


def load_brain_sample_space_inspect_volumes(
    config: BrainPipelineConfig,
    *,
    paths: BrainSampleSpaceInspectPaths | None = None,
) -> "BrainViewVolumes":
    """Load registration-grid channels, warped atlas overlays, and sample imports."""
    from lightsuite.analysis.counts import SAMPLE_POINTS_KEY, load_atlas_points
    from lightsuite.export.brain_export import _load_transform_params
    from lightsuite.gui.brain_view_data import BrainViewVolumes

    paths = paths or discover_brain_sample_space_inspect_paths(config)
    save_path = config.sample.save_path.expanduser()
    transform_params = _load_transform_params(save_path)
    expected_shape = tuple(int(v) for v in transform_params.regvolsize)
    permute = transform_params.permute_sample_to_atlas or [1, 2, 3]

    registered_channels: dict[int, np.ndarray] = {}
    for ichan, path in paths.registered_channels.items():
        vol = load_permuted_registration_volume(path, permute)
        if tuple(vol.shape) != expected_shape:
            msg = f"{path.name} shape {vol.shape} != expected registration shape {expected_shape}"
            raise ValueError(msg)
        registered_channels[ichan] = vol

    atlas_permuted = atlas_volumes_permuted_on_disk(paths.sample_space_dir)

    template: np.ndarray | None = None
    if paths.template_path is not None and paths.template_path.is_file():
        template = load_sample_space_atlas_volume(
            paths.template_path,
            permute,
            permuted_on_disk=atlas_permuted,
        )
        if tuple(template.shape) != expected_shape:
            msg = (
                f"Template shape {template.shape} != registration shape {expected_shape}. "
                "Re-run 'lightsuite brain export --space sample'."
            )
            raise ValueError(msg)

    annotation: np.ndarray | None = None
    if paths.annotation_path is not None and paths.annotation_path.is_file():
        annotation = load_sample_space_atlas_volume(
            paths.annotation_path,
            permute,
            permuted_on_disk=atlas_permuted,
        )
        if tuple(annotation.shape) != expected_shape:
            msg = (
                f"Annotation shape {annotation.shape} != registration shape {expected_shape}. "
                "Re-run 'lightsuite brain export --space sample'."
            )
            raise ValueError(msg)

    mask_layers: dict[str, np.ndarray] = {}
    for label, path in paths.mask_paths.items():
        vol = load_registration_volume(path)
        if tuple(vol.shape) != expected_shape:
            msg = f"{path.name} shape {vol.shape} != expected registration shape {expected_shape}"
            raise ValueError(msg)
        mask_layers[label] = vol.astype(np.float32, copy=False)

    point_layers: dict[str, np.ndarray] = {}
    for label, npz_path in paths.point_npz_paths.items():
        point_layers[label] = load_atlas_points(npz_path, key=SAMPLE_POINTS_KEY)

    return BrainViewVolumes(
        template=template,
        annotation=annotation,
        boundary=None,
        division_labels=None,
        division_legend=None,
        registered_channels=registered_channels,
        point_layers=point_layers,
        mask_layers=mask_layers,
    )
