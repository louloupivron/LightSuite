"""Napari GUI for atlas/sample manual content-box picking."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.config.loader import save_content_box_to_config
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.brain_data import _normalize_display
from lightsuite.registration.content_bbox import ContentBox
from lightsuite.registration.content_probe import (
    ContentProbeData,
    ContentProbeTarget,
    apply_box_mask,
    format_content_box_report,
    load_content_probe_data,
    preview_box_for_display,
)

console = Console()


def run_content_box_picker(
    config: BrainPipelineConfig,
    config_path: Path,
    *,
    target: ContentProbeTarget,
    headless: bool = False,
    write_config: bool = False,
) -> ContentBox:
    """Probe or interactively pick a content crop box for atlas or sample."""
    data = load_content_probe_data(config, target=target)
    if headless:
        if write_config:
            save_content_box_to_config(config_path, target=target, box=data.box.to_manual_list())
        return data.box

    try:
        import napari
        from magicgui import magicgui
        from napari.utils.notifications import show_info
        from qtpy.QtCore import QTimer
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    viewer = napari.Viewer(title=f"LightSuite content box — {config.sample.name} ({target.value})")
    display = _normalize_display(data.volume)
    viewer.add_image(display, name=f"{target.value} volume", colormap="gray", blending="opaque")
    mask_layer = viewer.add_image(
        display,
        name="crop preview (outside dimmed)",
        colormap="green",
        blending="additive",
        opacity=0.35,
    )

    state = {"box": data.box}

    def _refresh_layers(box: ContentBox) -> None:
        preview_box = preview_box_for_display(data, box)
        masked = apply_box_mask(data.volume, preview_box)
        mask_layer.data = _normalize_display(masked)
        viewer.status = (
            f"box={box.to_manual_list()} | crop size {box.size_yxz} | "
            f"full shape {data.full_shape}"
        )

    @magicgui(
        y0={"min": 0, "max": data.full_shape[0] - 1, "label": "Y start"},
        y1={"min": 0, "max": data.full_shape[0] - 1, "label": "Y end"},
        x0={"min": 0, "max": data.full_shape[1] - 1, "label": "X start"},
        x1={"min": 0, "max": data.full_shape[1] - 1, "label": "X end"},
        z0={"min": 0, "max": data.full_shape[2] - 1, "label": "Z start"},
        z1={"min": 0, "max": data.full_shape[2] - 1, "label": "Z end"},
        call_button="Update preview",
    )
    def controls(
        y0: int = data.box.y0,
        y1: int = data.box.y1,
        x0: int = data.box.x0,
        x1: int = data.box.x1,
        z0: int = data.box.z0,
        z1: int = data.box.z1,
    ) -> None:
        try:
            box = ContentBox(y0=y0, y1=y1, x0=x0, x1=x1, z0=z0, z1=z1)
        except ValueError as exc:
            show_info(str(exc))
            return
        state["box"] = box
        _refresh_layers(box)

    @magicgui(call_button="Auto-detect foreground")
    def auto_controls() -> None:
        redetected = load_content_probe_data(config, target=target).box
        state["box"] = redetected
        controls.y0.value = redetected.y0
        controls.y1.value = redetected.y1
        controls.x0.value = redetected.x0
        controls.x1.value = redetected.x1
        controls.z0.value = redetected.z0
        controls.z1.value = redetected.z1
        _refresh_layers(redetected)
        show_info(f"Auto box: {redetected.to_manual_list()}")

    @magicgui(call_button="Reset to full volume")
    def reset_controls() -> None:
        full = ContentBox.full_volume(data.full_shape)
        state["box"] = full
        controls.y0.value = full.y0
        controls.y1.value = full.y1
        controls.x0.value = full.x0
        controls.x1.value = full.x1
        controls.z0.value = full.z0
        controls.z1.value = full.z1
        _refresh_layers(full)

    @magicgui(call_button="Save box to config && close")
    def save_controls() -> None:
        box = state["box"]
        path = save_content_box_to_config(
            config_path,
            target=target.value,
            box=box.to_manual_list(),
        )
        show_info(f"Saved {target.value} manual crop to {path}")
        QTimer.singleShot(0, viewer.close)

    viewer.window.add_dock_widget(controls, area="right", name="Crop box")
    viewer.window.add_dock_widget(auto_controls, area="right", name="Auto")
    viewer.window.add_dock_widget(reset_controls, area="right", name="Reset")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")

    if data.preview_downsampled:
        console.print(
            "[yellow]Preview downsampled for Napari "
            f"({data.full_shape} → {tuple(int(v) for v in data.volume.shape)}). "
            "Slider indices remain in full-volume voxel coordinates.[/yellow]"
        )
    console.print(
        f"[bold]Content box picker[/bold] ({target.value}) — "
        "adjust Y/X/Z start/end, then Update preview. "
        "Outside the box is dimmed in green."
    )
    console.print(format_content_box_report(data))
    _refresh_layers(data.box)
    napari.run()
    return state["box"]
