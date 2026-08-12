"""Napari GUI for rostrocaudal sample-to-atlas alignment before cord init-registration."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from rich.console import Console

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.cord_align_data import (
    align_longitudinal_pair,
    estimate_cord_atlas_plane,
    load_cord_align_longitudinal_data,
    prepare_cord_align_longitudinal_session,
)
from lightsuite.gui.stage_controller import (
    DockStageController,
    close_stage_or_viewer,
    require_magicgui,
    run_attached_stage,
)
from lightsuite.registration.cord_longitudinal import CORD_LONGITUDINAL_AXIS

console = Console()
PANEL_GAP_X = 24


def _anchor_label(chooselist: np.ndarray, slice_idx: int) -> str:
    row = np.asarray(chooselist[slice_idx - 1], dtype=int)
    return f"sample z={row[0]} / axis {row[1]} cut"


def attach_spinal_align_longitudinal(
    viewer: Any,
    config: SpinalCordPipelineConfig,
) -> DockStageController:
    """Attach cord longitudinal alignment controls to an existing napari viewer."""
    from napari.utils.notifications import show_info
    from qtpy.QtCore import QTimer

    magicgui = require_magicgui()
    data = load_cord_align_longitudinal_data(config)
    state = {
        "slice": 1,
        "_nav_syncing": False,
        "_view_shape": None,
        "_atlas_plane": None,
    }

    def _n_slices() -> int:
        return int(data.chooselist.shape[0])

    viewer.dims.ndisplay = 2
    sample_layer = viewer.add_image(np.zeros((10, 10)), name="sample", colormap="gray")
    atlas_layer = viewer.add_image(np.zeros((10, 10)), name="atlas", colormap="gray")

    def _layout_panels() -> None:
        _h, w = sample_layer.data.shape
        sample_layer.translate = (0.0, 0.0)
        atlas_layer.translate = (0.0, float(w + PANEL_GAP_X))

    def _current_anchor():
        return data.correspondence.anchors_for_axis(CORD_LONGITUDINAL_AXIS)[state["slice"] - 1]

    def _atlas_plane_limits() -> tuple[int, int]:
        return 1, data.atlas_depth

    def _default_plane(slice_idx: int) -> int:
        anchor = data.correspondence.anchors_for_axis(CORD_LONGITUDINAL_AXIS)[slice_idx - 1]
        if anchor.atlas_plane > 0:
            return int(anchor.atlas_plane)
        row = np.asarray(data.chooselist[slice_idx - 1], dtype=int)
        return estimate_cord_atlas_plane(
            int(row[0]),
            nslices=data.nslices,
            atlas_depth=data.atlas_depth,
        )

    def _current_atlas_plane() -> int:
        plane = state.get("_atlas_plane")
        if plane is None:
            plane = _default_plane(state["slice"])
            state["_atlas_plane"] = plane
        return int(plane)

    def _set_atlas_plane(plane: int) -> None:
        plane = int(np.clip(plane, *_atlas_plane_limits()))
        state["_atlas_plane"] = plane
        _current_anchor().atlas_plane = plane
        _refresh()

    def _update_status() -> None:
        idx = state["slice"]
        plane = _current_atlas_plane()
        anchors = data.correspondence.anchors_for_axis(CORD_LONGITUDINAL_AXIS)
        confirmed = sum(anchor.confirmed for anchor in anchors)
        anchor = anchors[idx - 1]
        flag = "confirmed" if anchor.confirmed else "pending"
        viewer.status = (
            f"Anchor {idx}/{_n_slices()} ({_anchor_label(data.chooselist, idx)}) | {flag} "
            f"| atlas z {plane}/{data.atlas_depth} "
            f"| confirmed {confirmed}/{_n_slices()}"
        )

    def _refresh() -> None:
        idx = state["slice"]
        sample, atlas = align_longitudinal_pair(
            data,
            idx,
            atlas_plane=_current_atlas_plane(),
        )
        sample_layer.data = sample
        atlas_layer.data = atlas
        _layout_panels()
        _update_status()
        if state["_view_shape"] != sample.shape:
            viewer.reset_view()
            state["_view_shape"] = sample.shape

    def _sync_navigation_widget() -> None:
        state["_nav_syncing"] = True
        try:
            navigation.slice_index.max = _n_slices()
            navigation.slice_index.value = int(state["slice"])
            pmin, pmax = _atlas_plane_limits()
            navigation.atlas_plane.min = pmin
            navigation.atlas_plane.max = pmax
            navigation.atlas_plane.value = _current_atlas_plane()
        finally:
            state["_nav_syncing"] = False

    def _refocus_canvas() -> None:
        try:
            qt_viewer = viewer.window._qt_viewer
        except AttributeError:
            try:
                qt_viewer = viewer.window.qt_viewer
            except AttributeError:
                return
        canvas_native = getattr(getattr(qt_viewer, "canvas", None), "native", None)
        for target in (canvas_native, qt_viewer):
            if target is not None and hasattr(target, "setFocus"):
                target.setFocus()
                return

    def _navigate_to(slice_index: int | None = None, *, refocus_canvas: bool = False) -> None:
        if slice_index is not None:
            state["slice"] = max(1, min(_n_slices(), int(slice_index)))
            state["_atlas_plane"] = _default_plane(state["slice"])
        if refocus_canvas:
            _refresh()

            def _after_keyboard_nav() -> None:
                _sync_navigation_widget()
                _refocus_canvas()

            QTimer.singleShot(0, _after_keyboard_nav)
        else:
            _refresh()
            _sync_navigation_widget()

    def _confirm_and_advance() -> None:
        idx = state["slice"]
        anchor = _current_anchor()
        anchor.confirmed = True
        anchor.sample_index = int(data.chooselist[idx - 1, 0])
        anchor.atlas_plane = _current_atlas_plane()
        show_info(
            f"Confirmed anchor {idx}: sample z {anchor.sample_index} ↔ atlas z {anchor.atlas_plane}"
        )
        if idx < _n_slices():
            _navigate_to(idx + 1, refocus_canvas=True)
        else:
            _refresh()

    @magicgui(
        slice_index={
            "min": 1,
            "max": 20,
            "step": 1,
            "label": "Anchor #",
        },
        atlas_plane={
            "min": 1,
            "max": 100,
            "step": 1,
            "label": "Atlas z",
        },
        call_button="Show slice",
    )
    def navigation(slice_index: int = 1, atlas_plane: int = 1) -> None:
        _navigate_to(slice_index)
        _set_atlas_plane(atlas_plane)

    @navigation.atlas_plane.changed.connect
    def _atlas_plane_widget_changed() -> None:
        if state["_nav_syncing"]:
            return
        _set_atlas_plane(navigation.atlas_plane.value)

    @navigation.slice_index.changed.connect
    def _slice_index_widget_changed() -> None:
        if state["_nav_syncing"]:
            return
        _navigate_to(navigation.slice_index.value)

    @magicgui(call_button="◀  Previous anchor")
    def previous_slice() -> None:
        _navigate_to(state["slice"] - 1)

    @magicgui(call_button="Next anchor  ▶")
    def next_slice() -> None:
        _navigate_to(state["slice"] + 1)

    @magicgui(call_button="Confirm anchor (Enter)")
    def confirm_slice() -> None:
        _confirm_and_advance()

    @magicgui(call_button="Save && Close")
    def save_controls() -> None:
        data.correspondence.source = "manual"
        data.correspondence.save(data.correspondence_path)
        show_info(f"Saved {data.correspondence_path}")
        QTimer.singleShot(0, lambda: close_stage_or_viewer(viewer))

    def _previous_slice_key(_viewer) -> None:
        _navigate_to(state["slice"] - 1, refocus_canvas=True)

    def _next_slice_key(_viewer) -> None:
        _navigate_to(state["slice"] + 1, refocus_canvas=True)

    def _confirm_key(_viewer) -> None:
        _confirm_and_advance()

    def _atlas_plane_step(delta: int) -> None:
        _set_atlas_plane(_current_atlas_plane() + int(delta))
        _sync_navigation_widget()
        _refocus_canvas()

    def _atlas_plane_up(_viewer) -> None:
        _atlas_plane_step(1)

    def _atlas_plane_down(_viewer) -> None:
        _atlas_plane_step(-1)

    viewer.bind_key("Left", _previous_slice_key, overwrite=True)
    viewer.bind_key("Right", _next_slice_key, overwrite=True)
    viewer.bind_key("Enter", _confirm_key, overwrite=True)
    viewer.bind_key("PageUp", _atlas_plane_up, overwrite=True)
    viewer.bind_key("PageDown", _atlas_plane_down, overwrite=True)

    @viewer.window._qt_viewer.canvas.events.mouse_wheel.connect
    def _scroll_atlas_plane(event) -> None:
        if getattr(event, "delta", None) is None:
            return
        dy = event.delta[1] if len(event.delta) > 1 else event.delta[0]
        if dy == 0:
            return
        pos = getattr(viewer.window._qt_viewer.cursor, "position", None)
        if pos is None:
            return
        w = float(sample_layer.data.shape[1])
        col = float(pos[-1]) if len(pos) >= 2 else float(pos[0])
        if col < w + PANEL_GAP_X * 0.5:
            return
        step = 1 if dy < 0 else -1
        _atlas_plane_step(step)
        _sync_navigation_widget()

    def _initial_refresh() -> None:
        _refresh()
        _sync_navigation_widget()

    return DockStageController(
        dock_widgets=[
            (navigation, "Navigation"),
            (previous_slice, "Previous anchor"),
            (next_slice, "Next anchor"),
            (confirm_slice, "Confirm"),
            (save_controls, "Save"),
        ],
        _refresh_fn=_initial_refresh,
        result=data.correspondence_path,
    )


def run_spinal_align_longitudinal(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
) -> Path:
    """Launch Napari longitudinal alignment tool; returns correspondence JSON path."""
    if headless:
        return prepare_cord_align_longitudinal_session(config)

    title = f"LightSuite align-longitudinal — {config.sample.name}"

    def _attach(viewer: Any) -> DockStageController:
        return attach_spinal_align_longitudinal(viewer, config)

    run_attached_stage(
        title,
        _attach,
        before_run=lambda: console.print(
            "[bold]Napari align-longitudinal GUI[/bold] — straightened sample (left), atlas (right). "
            "Confirm ~20 rostrocaudal anchors: scroll atlas z with PgUp/PgDn or wheel over the atlas "
            "panel, then [bold]Enter[/bold] or Confirm. Save before running init-registration."
        ),
    )
    data = load_cord_align_longitudinal_data(config)
    return data.correspondence_path
