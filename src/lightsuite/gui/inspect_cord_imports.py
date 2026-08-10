"""Napari viewer for spinal cord registration export + annotation import QA."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

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
from lightsuite.gui.inspect_brain_imports import (
    _contrast_limits,
    _label_from_stem,
    _load_points_csv,
    _load_points_npz,
    volume_yxz_to_napari_zyx,
)
from lightsuite.gui.view_registered_cord import _add_channel_layers, _add_point_layers

ViewSpace = Literal["atlas", "sample"]


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
    for path in sorted(vr.glob("*_atlas_coords.npz")):
        label = _label_from_stem(path.stem, "_atlas_coords")
        point_npz_paths[label] = path.resolve()

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
    if vr.is_dir():
        for path in sorted(vr.glob("*_sample_coords.npz")):
            label = _label_from_stem(path.stem, "_sample_coords")
            point_npz_paths[label] = path.resolve()

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
    for label, npz_path in paths.point_npz_paths.items():
        point_layers[label] = _load_points_npz(npz_path)
    for path in sorted(paths.volume_registered_dir.glob("*_atlas_coords.csv")):
        label = _label_from_stem(path.stem, "_atlas_coords")
        if label in point_layers:
            continue
        point_layers[label] = _load_points_csv(path)

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
    for label, npz_path in paths.point_npz_paths.items():
        if label in point_layers:
            continue
        point_layers[label] = load_atlas_points(npz_path, key=SAMPLE_POINTS_KEY)

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


def run_cord_inspect_imports(
    config: SpinalCordPipelineConfig,
    *,
    space: ViewSpace = "atlas",
    headless: bool = False,
    recompute_annotation: bool = False,
) -> CordImportInspectPaths:
    """Open Napari to QC registered channels and imported annotations."""
    paths = discover_cord_import_inspect_paths(config, space=space)
    if headless:
        load_cord_import_inspect_volumes(
            config,
            paths=paths,
            space=space,
            recompute_annotation=recompute_annotation,
        )
        return paths

    try:
        import napari
        from napari.utils.notifications import show_info
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    volumes = load_cord_import_inspect_volumes(
        config,
        paths=paths,
        space=space,
        recompute_annotation=recompute_annotation,
    )

    if space == "sample":
        title = f"LightSuite spinal import QC — {config.sample.name} (sample, 20 µm straightened)"
        template_name = "atlas template (warped)"
        annotation_name = "atlas annotation (warped)"
        channel_suffix = "straightened"
    else:
        title = f"LightSuite spinal import QC — {config.sample.name} (atlas)"
        template_name = "atlas template"
        annotation_name = "atlas annotation"
        channel_suffix = "registered"

    viewer = napari.Viewer(title=title)

    viewer.add_image(
        volume_yxz_to_napari_zyx(volumes.template),
        name=template_name,
        colormap="green",
        blending="additive",
        opacity=0.35,
        contrast_limits=_contrast_limits(volumes.template),
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

    summary_path = paths.volume_registered_dir / "import_annotations_summary.json"
    space_note = "sample-space " if space == "sample" else ""
    if summary_path.is_file():
        show_info(
            f"Loaded {len(volumes.registered_channels)} {space_note}channel(s) and "
            f"{len(volumes.point_layers)} point layer(s). "
            f"Summary: {summary_path.name}"
            + (
                " Sample Z flipped to match atlas rostrocaudal orientation (tofliprc)."
                if space == "sample" and load_cord_tofliprc(config.sample.save_path)
                else ""
            )
        )
    else:
        show_info(
            f"Loaded {len(volumes.registered_channels)} {space_note}channel(s) and "
            f"{len(volumes.point_layers)} point layer(s)."
        )
    napari.run()
    return paths
