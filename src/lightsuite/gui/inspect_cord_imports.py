"""Napari viewer for spinal cord registration export + annotation import QA."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

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
from lightsuite.gui.inspect_brain_imports import (
    _contrast_limits,
    _label_from_stem,
    _load_points_csv,
    _load_points_npz,
    atlas_points_to_napari_zyx,
    volume_yxz_to_napari_zyx,
)


@dataclass(frozen=True)
class CordImportInspectPaths:
    """Resolved inputs under ``<save_path>/volume_registered/``."""

    volume_registered_dir: Path
    template_path: Path
    annotation_path: Path
    registered_channels: dict[int, Path] = field(default_factory=dict)
    point_npz_paths: dict[str, Path] = field(default_factory=dict)


@dataclass
class CordImportInspectVolumes:
    template: np.ndarray
    annotation: np.ndarray
    registered_channels: dict[int, np.ndarray]
    point_layers: dict[str, np.ndarray]


def discover_cord_import_inspect_paths(config: SpinalCordPipelineConfig) -> CordImportInspectPaths:
    """Discover registered cord volumes and imported annotation outputs."""
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
    )


def load_cord_import_inspect_volumes(
    config: SpinalCordPipelineConfig,
    *,
    paths: CordImportInspectPaths | None = None,
    recompute_annotation: bool = False,
) -> CordImportInspectVolumes:
    """Load registered channels, atlas labels, and imported point layers for Napari."""
    paths = paths or discover_cord_import_inspect_paths(config)
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


def run_cord_inspect_imports(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
    recompute_annotation: bool = False,
) -> CordImportInspectPaths:
    """Open Napari to QC registered channels and imported annotations in atlas space."""
    paths = discover_cord_import_inspect_paths(config)
    if headless:
        load_cord_import_inspect_volumes(
            config,
            paths=paths,
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
        recompute_annotation=recompute_annotation,
    )
    viewer = napari.Viewer(title=f"LightSuite spinal import QC — {config.sample.name}")

    viewer.add_image(
        volume_yxz_to_napari_zyx(volumes.template),
        name="atlas template",
        colormap="green",
        blending="additive",
        opacity=0.35,
        contrast_limits=_contrast_limits(volumes.template),
    )

    channel_cmaps = ["gray", "magenta", "cyan", "yellow", "green", "red"]
    for idx, (ichan, vol) in enumerate(sorted(volumes.registered_channels.items())):
        cmap = channel_cmaps[idx % len(channel_cmaps)]
        opacity = 1.0 if len(volumes.registered_channels) == 1 else 0.65
        viewer.add_image(
            volume_yxz_to_napari_zyx(vol),
            name=f"channel {ichan} registered",
            colormap=cmap,
            blending="additive" if len(volumes.registered_channels) > 1 else "opaque",
            opacity=opacity,
            contrast_limits=_contrast_limits(vol),
        )

    viewer.add_labels(
        volume_yxz_to_napari_zyx(volumes.annotation),
        name="atlas annotation",
        opacity=0.45,
    )

    point_colors = ["red", "yellow", "cyan", "magenta", "orange", "lime"]
    for idx, (label, coords) in enumerate(sorted(volumes.point_layers.items())):
        napari_pts = atlas_points_to_napari_zyx(coords)
        color = point_colors[idx % len(point_colors)]
        viewer.add_points(
            napari_pts,
            name=f"points: {label}",
            size=4,
            face_color=color,
            border_color="white",
        )

    summary_path = paths.volume_registered_dir / "import_annotations_summary.json"
    if summary_path.is_file():
        show_info(
            f"Loaded {len(volumes.registered_channels)} channel(s) and "
            f"{len(volumes.point_layers)} point layer(s). "
            f"Summary: {summary_path.name}"
        )
    else:
        show_info(
            f"Loaded {len(volumes.registered_channels)} channel(s) and "
            f"{len(volumes.point_layers)} point layer(s)."
        )
    napari.run()
    return paths
