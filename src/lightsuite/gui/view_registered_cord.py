"""Napari viewer for registered spinal cord channels and atlas annotations."""

from __future__ import annotations

from typing import Any, Literal

import numpy as np

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import (
    CordRegisteredInspectPaths,
    discover_registered_cord_paths,
    load_native_template_export_layout,
    load_registered_cord_volumes,
)
from lightsuite.export.cord_sample_space import (
    CordSampleSpaceInspectPaths,
    discover_cord_sample_space_paths,
    load_cord_sample_space_volumes,
)
from lightsuite.gui.cord_napari_display import align_sample_space_for_atlas_qc, load_cord_tofliprc
from lightsuite.gui.inspect_brain_imports import (
    _contrast_limits,
    atlas_points_to_napari_zyx,
    volume_yxz_to_napari_zyx,
)
from lightsuite.gui.stage_controller import DockStageController, run_attached_stage

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
            contrast_limits=_contrast_limits(vol),
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


def attach_spinal_registered_atlas_view(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    paths: CordRegisteredInspectPaths,
    volumes,
) -> DockStageController:
    """Attach atlas-space registered cord layers to an existing napari viewer."""
    from napari.utils.notifications import show_info

    template = load_native_template_export_layout(config)
    viewer.add_image(
        volume_yxz_to_napari_zyx(template),
        name="atlas template",
        colormap="green",
        blending="additive",
        opacity=0.35,
        contrast_limits=_contrast_limits(template),
    )

    _add_channel_layers(viewer, volumes.registered_channels, name_suffix="registered")

    viewer.add_labels(
        volume_yxz_to_napari_zyx(volumes.annotation),
        name="atlas annotation",
        opacity=0.45,
    )

    def _notify() -> None:
        show_info(
            f"Loaded {len(volumes.registered_channels)} registered channel(s) "
            f"and annotation labels from {paths.volume_registered_dir.name}/."
        )

    return DockStageController(
        dock_widgets=[],
        _refresh_fn=_notify,
        result=paths,
    )


def attach_spinal_registered_sample_view(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    paths: CordSampleSpaceInspectPaths,
    template: np.ndarray,
    annotation: np.ndarray,
    channels: dict[int, np.ndarray],
    point_layers: dict[str, np.ndarray],
    tofliprc: bool,
) -> DockStageController:
    """Attach sample-space registered cord layers to an existing napari viewer."""
    from napari.utils.notifications import show_info

    viewer.add_image(
        volume_yxz_to_napari_zyx(template),
        name="atlas template (warped)",
        colormap="green",
        blending="additive",
        opacity=0.35,
        contrast_limits=_contrast_limits(template),
    )

    _add_channel_layers(viewer, channels, name_suffix="straightened")

    viewer.add_labels(
        volume_yxz_to_napari_zyx(annotation),
        name="atlas annotation (warped)",
        opacity=0.45,
    )

    if point_layers:
        _add_point_layers(viewer, point_layers)

    def _notify() -> None:
        summary = (
            f"Loaded {len(channels)} straightened channel(s) and warped annotation "
            f"from {paths.sample_space_dir.name}/ on the 20 µm registration grid."
        )
        if point_layers:
            summary += f" Point layers: {len(point_layers)}."
        if tofliprc:
            summary += " Sample Z flipped for atlas-aligned rostrocaudal QC."
        show_info(summary)

    return DockStageController(
        dock_widgets=[],
        _refresh_fn=_notify,
        result=paths,
    )


def run_spinal_registered_view(
    config: SpinalCordPipelineConfig,
    *,
    space: ViewSpace = "atlas",
    headless: bool = False,
    recompute_annotation: bool = False,
) -> CordRegisteredInspectPaths | CordSampleSpaceInspectPaths:
    """Open Napari showing registered channels and warped atlas annotations."""
    if space == "sample":
        return _run_spinal_sample_space_view(config, headless=headless)
    return _run_spinal_atlas_space_view(
        config,
        headless=headless,
        recompute_annotation=recompute_annotation,
    )


def _run_spinal_atlas_space_view(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool,
    recompute_annotation: bool,
) -> CordRegisteredInspectPaths:
    paths = discover_registered_cord_paths(config)
    if headless:
        load_registered_cord_volumes(
            config,
            paths=paths,
            recompute_annotation=recompute_annotation,
        )
        return paths

    volumes = load_registered_cord_volumes(
        config,
        paths=paths,
        recompute_annotation=recompute_annotation,
    )
    title = f"LightSuite spinal registration — {config.sample.name} (atlas)"

    def _attach(viewer: Any) -> DockStageController:
        return attach_spinal_registered_atlas_view(
            viewer,
            config,
            paths=paths,
            volumes=volumes,
        )

    final = run_attached_stage(title, _attach)
    return final if final is not None else paths


def _run_spinal_sample_space_view(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool,
) -> CordSampleSpaceInspectPaths:
    paths = discover_cord_sample_space_paths(config)
    if headless:
        load_cord_sample_space_volumes(config, paths=paths)
        return paths

    volumes = load_cord_sample_space_volumes(config, paths=paths)
    tofliprc = load_cord_tofliprc(config.sample.save_path)
    template, annotation, channels, point_layers, _hemisphere = align_sample_space_for_atlas_qc(
        template=volumes.template,
        annotation=volumes.annotation,
        channels=volumes.channels,
        point_layers=volumes.point_layers,
        hemisphere=volumes.hemisphere,
        tofliprc=tofliprc,
    )
    title = (
        f"LightSuite spinal registration — {config.sample.name} (sample, 20 µm straightened)"
    )

    def _attach(viewer: Any) -> DockStageController:
        return attach_spinal_registered_sample_view(
            viewer,
            config,
            paths=paths,
            template=template,
            annotation=annotation,
            channels=channels,
            point_layers=point_layers,
            tofliprc=tofliprc,
        )

    final = run_attached_stage(title, _attach)
    return final if final is not None else paths
