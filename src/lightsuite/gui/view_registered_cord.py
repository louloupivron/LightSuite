"""Napari viewer for registered spinal cord channels and atlas annotations."""

from __future__ import annotations

from typing import Literal

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
from lightsuite.gui.inspect_brain_imports import (
    _contrast_limits,
    atlas_points_to_napari_zyx,
    volume_yxz_to_napari_zyx,
)

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

    try:
        import napari
        from napari.utils.notifications import show_info
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    volumes = load_registered_cord_volumes(
        config,
        paths=paths,
        recompute_annotation=recompute_annotation,
    )
    viewer = napari.Viewer(title=f"LightSuite spinal registration — {config.sample.name} (atlas)")

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

    show_info(
        f"Loaded {len(volumes.registered_channels)} registered channel(s) "
        f"and annotation labels from {paths.volume_registered_dir.name}/."
    )
    napari.run()
    return paths


def _run_spinal_sample_space_view(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool,
) -> CordSampleSpaceInspectPaths:
    paths = discover_cord_sample_space_paths(config)
    if headless:
        load_cord_sample_space_volumes(config, paths=paths)
        return paths

    try:
        import napari
        from napari.utils.notifications import show_info
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    volumes = load_cord_sample_space_volumes(config, paths=paths)
    viewer = napari.Viewer(
        title=f"LightSuite spinal registration — {config.sample.name} (sample, 20 µm straightened)"
    )

    viewer.add_image(
        volume_yxz_to_napari_zyx(volumes.template),
        name="atlas template (warped)",
        colormap="green",
        blending="additive",
        opacity=0.35,
        contrast_limits=_contrast_limits(volumes.template),
    )

    _add_channel_layers(viewer, volumes.channels, name_suffix="straightened")

    viewer.add_labels(
        volume_yxz_to_napari_zyx(volumes.annotation),
        name="atlas annotation (warped)",
        opacity=0.45,
    )

    if volumes.point_layers:
        _add_point_layers(viewer, volumes.point_layers)

    summary = (
        f"Loaded {len(volumes.channels)} straightened channel(s) and warped annotation "
        f"from {paths.sample_space_dir.name}/ on the 20 µm registration grid."
    )
    if volumes.point_layers:
        summary += f" Point layers: {len(volumes.point_layers)}."
    show_info(summary)
    napari.run()
    return paths
