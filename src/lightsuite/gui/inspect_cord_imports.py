"""Napari viewer for spinal cord registration export + annotation import QA."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import numpy as np

from lightsuite.analysis.counts import SAMPLE_POINTS_KEY, load_atlas_points
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import (
    REGISTERED_ANNOTATION_FILENAME,
    REGISTERED_TEMPLATE_FILENAME,
    CordRegisteredInspectPaths,
    discover_registered_cord_paths,
    load_native_template_export_layout,
    load_registered_cord_volumes,
    load_registered_stack,
)
from lightsuite.export.cord_sample_space import (
    discover_cord_sample_space_paths,
    load_cord_sample_space_volumes,
)
from lightsuite.gui.cord_napari_display import align_sample_space_for_atlas_qc, load_cord_tofliprc
from lightsuite.gui.brain_view_data import (
    atlas_points_to_napari_zyx,
    label_from_stem,
    load_points_csv,
    load_points_npz,
    contrast_limits,
    volume_yxz_to_napari_zyx,
)
from lightsuite.gui.stage_controller import DockStageController

ViewSpace = Literal["atlas", "sample"]


def _add_channel_layers(viewer, channels: dict[int, np.ndarray], *, name_suffix: str = "") -> None:
    channel_cmaps = ["gray", "magenta", "cyan", "yellow", "green", "red"]
    suffix = f" {name_suffix}".rstrip()
    for idx, (ichan, vol) in enumerate(sorted(channels.items())):
        cmap = channel_cmaps[idx % len(channel_cmaps)]
        opacity = 1.0 if len(channels) == 1 else 0.65
        viewer.add_image(
            volume_yxz_to_napari_zyx(vol),
            name=f"channel {ichan}{suffix}",
            colormap=cmap,
            blending="additive" if len(channels) > 1 else "opaque",
            opacity=opacity,
            contrast_limits=contrast_limits(vol),
        )


def _add_point_layers(viewer, point_layers: dict[str, np.ndarray]) -> None:
    point_colors = ["red", "yellow", "cyan", "magenta", "orange", "lime"]
    for idx, (label, coords) in enumerate(sorted(point_layers.items())):
        napari_pts = atlas_points_to_napari_zyx(coords)
        color = point_colors[idx % len(point_colors)]
        viewer.add_points(
            napari_pts,
            name=f"points: {label}",
            size=4,
            face_color=color,
            border_color="white",
        )


@dataclass(frozen=True)
class CordImportInspectPaths:
    """Resolved inputs under ``<save_path>/volume_registered/``."""

    volume_registered_dir: Path
    template_path: Path
    annotation_path: Path
    registered_channels: dict[int, Path] = field(default_factory=dict)
    point_npz_paths: dict[str, Path] = field(default_factory=dict)
    space: ViewSpace = "atlas"
    sample_space_dir: Path | None = None


@dataclass
class CordImportInspectVolumes:
    template: np.ndarray
    annotation: np.ndarray
    registered_channels: dict[int, np.ndarray]
    point_layers: dict[str, np.ndarray]
    hemisphere: np.ndarray | None = None


def discover_cord_import_inspect_paths(
    config: SpinalCordPipelineConfig,
    *,
    space: ViewSpace = "atlas",
) -> CordImportInspectPaths:
    """Discover registered cord volumes and imported annotation outputs."""
    if space == "sample":
        return _discover_cord_import_inspect_paths_sample(config)
    return _discover_cord_import_inspect_paths_atlas(config)


def resolve_cord_view_space(
    config: SpinalCordPipelineConfig,
    *,
    preferred: ViewSpace = "sample",
) -> ViewSpace:
    """Pick sample- or atlas-space review when the preferred export is missing."""
    available = cord_view_spaces_available(config)
    if preferred == "sample" and available["sample"]:
        return "sample"
    if available["atlas"]:
        return "atlas"
    if available["sample"]:
        return "sample"

    save_path = config.sample.save_path.expanduser()
    from lightsuite.export.cord_sample_space import sample_space_dir

    msg = (
        f"Missing view-registration exports under {save_path / 'volume_registered'}. "
        "Run 'lightsuite spinal export --space both' first."
    )
    raise FileNotFoundError(msg)


def cord_view_spaces_available(config: SpinalCordPipelineConfig) -> dict[ViewSpace, bool]:
    """Return which coordinate spaces can be opened in view-registration."""
    available: dict[ViewSpace, bool] = {"sample": False, "atlas": False}
    try:
        discover_cord_import_inspect_paths(config, space="sample")
    except FileNotFoundError:
        pass
    else:
        available["sample"] = True
    try:
        discover_cord_import_inspect_paths(config, space="atlas")
    except FileNotFoundError:
        pass
    else:
        available["atlas"] = True
    return available


def cord_view_load_summary(
    paths: CordImportInspectPaths,
    volumes: CordImportInspectVolumes,
    *,
    space: ViewSpace,
    config: SpinalCordPipelineConfig,
) -> str:
    """Short status line after loading a spinal view-registration space."""
    summary_path = paths.volume_registered_dir / "import_annotations_summary.json"
    space_note = "sample-space " if space == "sample" else "atlas-space "
    parts = [
        f"Loaded {space_note}{len(volumes.registered_channels)} channel(s) and "
        f"{len(volumes.point_layers)} point layer(s)."
    ]
    if summary_path.is_file():
        parts.append(f"Summary: {summary_path.name}.")
    if space == "sample" and load_cord_tofliprc(config.sample.save_path):
        parts.append("Sample Z flipped to match atlas rostrocaudal orientation (tofliprc).")
    return " ".join(parts)


def add_cord_view_layers(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    paths: CordImportInspectPaths,
    volumes: CordImportInspectVolumes,
    space: ViewSpace,
) -> None:
    """Add Napari layers for one spinal view-registration coordinate space."""
    if space == "sample":
        template_name = "atlas template (warped)"
        annotation_name = "atlas annotation (warped)"
        channel_suffix = "straightened"
    else:
        template_name = "atlas template"
        annotation_name = "atlas annotation"
        channel_suffix = "registered"

    viewer.add_image(
        volume_yxz_to_napari_zyx(volumes.template),
        name=template_name,
        colormap="green",
        blending="additive",
        opacity=0.35,
        contrast_limits=contrast_limits(volumes.template),
    )
    _add_channel_layers(viewer, volumes.registered_channels, name_suffix=channel_suffix)
    viewer.add_labels(
        volume_yxz_to_napari_zyx(volumes.annotation),
        name=annotation_name,
        opacity=0.45,
    )
    if space == "sample" and volumes.hemisphere is not None:
        viewer.add_labels(
            volume_yxz_to_napari_zyx(volumes.hemisphere),
            name="hemisphere (warped)",
            opacity=0.2,
            visible=False,
        )
    if volumes.point_layers:
        _add_point_layers(viewer, volumes.point_layers)


def iter_cord_import_paths(save_path: Path, pattern: str) -> list[Path]:
    paths: list[Path] = []
    for sub in ("volume_registered", "imports", ""):
        folder = save_path / sub if sub else save_path
        if folder.is_dir():
            paths.extend(sorted(folder.glob(pattern)))
    return sorted(dict.fromkeys(p.resolve() for p in paths))


def _discover_cord_import_inspect_paths_atlas(
    config: SpinalCordPipelineConfig,
) -> CordImportInspectPaths:
    save_path = config.sample.save_path.expanduser()
    vr = save_path / "volume_registered"
    if not vr.is_dir():
        msg = (
            f"Missing {vr}. Run 'lightsuite spinal export' and "
            "'lightsuite spinal import-annotations' first."
        )
        raise FileNotFoundError(msg)

    transform_params_path = save_path / "transform_params.json"
    if not transform_params_path.is_file():
        msg = f"Missing {transform_params_path}. Run 'lightsuite spinal register' first."
        raise FileNotFoundError(msg)

    registered = discover_registered_cord_paths(config)

    point_npz_paths: dict[str, Path] = {}
    for path in iter_cord_import_paths(save_path, "*_atlas_coords.npz"):
        label = label_from_stem(path.stem, "_atlas_coords")
        point_npz_paths.setdefault(label, path.resolve())
    for path in iter_cord_import_paths(save_path, "*_atlas_coords.csv"):
        label = label_from_stem(path.stem, "_atlas_coords")
        point_npz_paths.setdefault(label, path.resolve())

    if not registered.registered_channels and not point_npz_paths:
        msg = (
            f"No inspectable layers found in {vr}. "
            "Expected chan*_channel*.tiff and/or import outputs (*_atlas_coords.npz)."
        )
        raise FileNotFoundError(msg)

    template_path = vr / REGISTERED_TEMPLATE_FILENAME
    if not template_path.is_file():
        template_path = vr / "template_registered.tiff"

    return CordImportInspectPaths(
        volume_registered_dir=vr.resolve(),
        template_path=template_path.resolve(),
        annotation_path=registered.annotation_path.resolve(),
        registered_channels=registered.registered_channels,
        point_npz_paths=point_npz_paths,
        space="atlas",
    )


def _discover_cord_import_inspect_paths_sample(
    config: SpinalCordPipelineConfig,
) -> CordImportInspectPaths:
    save_path = config.sample.save_path.expanduser()
    vr = save_path / "volume_registered"
    sample_paths = discover_cord_sample_space_paths(config)

    point_npz_paths: dict[str, Path] = {}
    for path in iter_cord_import_paths(save_path, "*_sample_coords.npz"):
        label = label_from_stem(path.stem, "_sample_coords")
        point_npz_paths.setdefault(label, path.resolve())
    for path in iter_cord_import_paths(save_path, "*_sample_coords.csv"):
        label = label_from_stem(path.stem, "_sample_coords")
        point_npz_paths.setdefault(label, path.resolve())

    if not sample_paths.channel_paths and not point_npz_paths:
        msg = (
            f"No inspectable sample-space layers found. Run "
            "'lightsuite spinal export --space sample' and/or "
            "'lightsuite spinal import-annotations' first."
        )
        raise FileNotFoundError(msg)

    return CordImportInspectPaths(
        volume_registered_dir=vr.resolve() if vr.is_dir() else sample_paths.sample_space_dir,
        template_path=sample_paths.template_path,
        annotation_path=sample_paths.annotation_path,
        registered_channels=sample_paths.channel_paths,
        point_npz_paths=point_npz_paths or sample_paths.point_npz_paths,
        space="sample",
        sample_space_dir=sample_paths.sample_space_dir,
    )


def load_cord_import_inspect_volumes(
    config: SpinalCordPipelineConfig,
    *,
    paths: CordImportInspectPaths | None = None,
    space: ViewSpace = "atlas",
    recompute_annotation: bool = False,
) -> CordImportInspectVolumes:
    """Load registered channels, atlas labels, and imported point layers for Napari."""
    paths = paths or discover_cord_import_inspect_paths(config, space=space)
    view_space = paths.space if paths is not None else space
    if view_space == "sample":
        return _load_cord_import_inspect_volumes_sample(config, paths=paths)
    return _load_cord_import_inspect_volumes_atlas(
        config,
        paths=paths,
        recompute_annotation=recompute_annotation,
    )


def _load_cord_import_inspect_volumes_atlas(
    config: SpinalCordPipelineConfig,
    *,
    paths: CordImportInspectPaths,
    recompute_annotation: bool,
) -> CordImportInspectVolumes:
    registered = load_registered_cord_volumes(
        config,
        paths=CordRegisteredInspectPaths(
            volume_registered_dir=paths.volume_registered_dir,
            annotation_path=paths.annotation_path,
            registered_channels=paths.registered_channels,
        ),
        recompute_annotation=recompute_annotation,
    )

    if paths.template_path.is_file():
        template = load_registered_stack(paths.template_path).astype(np.float32, copy=False)
    else:
        template = load_native_template_export_layout(config).astype(np.float32, copy=False)

    if template.shape != registered.annotation.shape:
        msg = (
            f"Template shape {template.shape} != annotation shape {registered.annotation.shape}. "
            "Re-run 'lightsuite spinal export'."
        )
        raise ValueError(msg)

    point_layers: dict[str, np.ndarray] = {}
    for label, pt_path in paths.point_npz_paths.items():
        if pt_path.suffix.lower() == ".npz":
            try:
                point_layers[label] = load_atlas_points(pt_path, key="atlasptcoords")
            except (KeyError, ValueError):
                point_layers[label] = load_points_npz(pt_path)
        elif pt_path.suffix.lower() == ".csv":
            point_layers[label] = load_points_csv(pt_path)

    return CordImportInspectVolumes(
        template=template,
        annotation=registered.annotation.astype(np.int32, copy=False),
        registered_channels=registered.registered_channels,
        point_layers=point_layers,
    )


def _load_cord_import_inspect_volumes_sample(
    config: SpinalCordPipelineConfig,
    *,
    paths: CordImportInspectPaths,
) -> CordImportInspectVolumes:
    sample_paths = discover_cord_sample_space_paths(config)
    volumes = load_cord_sample_space_volumes(config, paths=sample_paths)

    point_layers: dict[str, np.ndarray] = dict(volumes.point_layers)
    for label, pt_path in paths.point_npz_paths.items():
        if label in point_layers:
            continue
        if pt_path.suffix.lower() == ".npz":
            point_layers[label] = load_atlas_points(pt_path, key=SAMPLE_POINTS_KEY)
        elif pt_path.suffix.lower() == ".csv":
            point_layers[label] = load_points_csv(pt_path)

    tofliprc = load_cord_tofliprc(config.sample.save_path)
    template, annotation, channels, point_layers, hemisphere = align_sample_space_for_atlas_qc(
        template=volumes.template,
        annotation=volumes.annotation,
        channels=volumes.channels,
        point_layers=point_layers,
        hemisphere=volumes.hemisphere,
        tofliprc=tofliprc,
    )

    return CordImportInspectVolumes(
        template=template,
        annotation=annotation,
        registered_channels=channels,
        point_layers=point_layers,
        hemisphere=hemisphere,
    )


def attach_cord_inspect_imports(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    paths: CordImportInspectPaths,
    volumes: CordImportInspectVolumes,
    space: ViewSpace = "atlas",
) -> DockStageController:
    """Attach spinal import QC layers to an existing napari viewer."""
    from napari.utils.notifications import show_info

    add_cord_view_layers(
        viewer,
        config,
        paths=paths,
        volumes=volumes,
        space=space,
    )

    def _notify() -> None:
        show_info(cord_view_load_summary(paths, volumes, space=space, config=config))

    return DockStageController(
        dock_widgets=[],
        _refresh_fn=_notify,
        result=paths,
    )


def run_cord_inspect_imports(
    config: SpinalCordPipelineConfig,
    *,
    space: ViewSpace = "atlas",
    headless: bool = False,
    recompute_annotation: bool = False,
) -> CordImportInspectPaths:
    """Open Napari to QC registered channels and imported annotations.

    Deprecated alias — use :func:`lightsuite.gui.view_registered_cord.run_spinal_view_registration`.
    """
    from lightsuite.gui.view_registered_cord import run_spinal_view_registration

    return run_spinal_view_registration(
        config,
        space=space,
        headless=headless,
        recompute_annotation=recompute_annotation,
    )
