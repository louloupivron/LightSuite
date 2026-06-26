"""Napari viewer for registered spinal cord channels and atlas annotations."""

from __future__ import annotations

import numpy as np

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.export.cord_registered import (
    CordRegisteredInspectPaths,
    discover_registered_cord_paths,
    load_native_template_export_layout,
    load_registered_cord_volumes,
)
from lightsuite.gui.inspect_brain_imports import _contrast_limits, volume_yxz_to_napari_zyx


def run_spinal_registered_view(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
    recompute_annotation: bool = False,
) -> CordRegisteredInspectPaths:
    """Open Napari showing registered sample channel(s) and warped annotation labels."""
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
    viewer = napari.Viewer(title=f"LightSuite spinal registration — {config.sample.name}")

    template = load_native_template_export_layout(config)
    viewer.add_image(
        volume_yxz_to_napari_zyx(template),
        name="atlas template",
        colormap="green",
        blending="additive",
        opacity=0.35,
        contrast_limits=_contrast_limits(template),
    )

    channel_cmaps = ["gray", "magenta", "cyan", "yellow", "green", "red"]
    for idx, (ichan, vol) in enumerate(sorted(volumes.registered_channels.items())):
        cmap = channel_cmaps[idx % len(channel_cmaps)]
        opacity = 1.0 if len(volumes.registered_channels) == 1 else 0.65
        viewer.add_image(
            volume_yxz_to_napari_zyx(vol),
            name=f"channel {ichan}",
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

    show_info(
        f"Loaded {len(volumes.registered_channels)} registered channel(s) "
        f"and annotation labels from {paths.volume_registered_dir.name}/."
    )
    napari.run()
    return paths
