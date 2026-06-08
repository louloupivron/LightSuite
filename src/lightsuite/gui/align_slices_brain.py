"""Napari GUI for sample-to-atlas slice correspondence before control-point matching."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.brain_data import (
    BrainAlignSlicesData,
    _normalize_display,
    atlas_cut_axis_size,
    chooserow_with_atlas_plane,
    estimate_atlas_plane_index,
    load_brain_align_slices_data,
    prepare_brain_align_slices_session,
)
from lightsuite.gui.match_points_brain import PANEL_GAP_X, _chooselist_slice_label
from lightsuite.gui.slices import volume_index_to_image

console = Console()

_AXIS_NAMES = {1: "Y", 2: "X", 3: "Z"}


def _align_slices_pair(
    data: BrainAlignSlicesData,
    slice_idx: int,
    *,
    atlas_plane: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample and atlas images for one align-slices chooselist entry."""
    row = np.asarray(data.chooselist[slice_idx - 1], dtype=int)
    sample = _normalize_display(volume_index_to_image(data.sample_volume, row))
    atlas_row = chooserow_with_atlas_plane(row, atlas_plane)
    atlas = _normalize_display(volume_index_to_image(data.atlas_template, atlas_row))
    return sample, atlas


def _current_anchor(data: BrainAlignSlicesData, slice_idx: int):
    return data.correspondence.anchors[slice_idx - 1]


def _sync_anchor_plane(data: BrainAlignSlicesData, slice_idx: int, plane: int) -> None:
    anchor = _current_anchor(data, slice_idx)
    anchor.atlas_plane = int(plane)


def run_brain_align_slices(config: BrainPipelineConfig, *, headless: bool = False) -> Path:
    """Launch napari slice-alignment tool; returns path to slice_correspondence.json."""
    if headless:
        return prepare_brain_align_slices_session(config)

    try:
        import napari
        from magicgui import magicgui
        from napari.utils.notifications import show_info
        from qtpy.QtCore import QTimer
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    data = load_brain_align_slices_data(config)
    n_slices = int(data.chooselist.shape[0])
    state = {
        "slice": 1,
        "_nav_syncing": False,
        "_view_shape": None,
        "_atlas_plane": None,
    }

    viewer = napari.Viewer(title=f"LightSuite align-slices — {config.sample.name}")
    viewer.dims.ndisplay = 2
    sample_layer = viewer.add_image(np.zeros((10, 10)), name="sample", colormap="gray")
    atlas_layer = viewer.add_image(np.zeros((10, 10)), name="atlas", colormap="gray")

    def _layout_panels() -> None:
        _h, w = sample_layer.data.shape
        sample_layer.translate = (0.0, 0.0)
        atlas_layer.translate = (0.0, float(w + PANEL_GAP_X))

    def _atlas_plane_limits() -> tuple[int, int]:
        row = np.asarray(data.chooselist[state["slice"] - 1], dtype=int)
        return 1, atlas_cut_axis_size(data.atlas_template.shape, row)

    def _default_plane(slice_idx: int) -> int:
        anchor = _current_anchor(data, slice_idx)
        if anchor.atlas_plane > 0:
            return int(anchor.atlas_plane)
        row = np.asarray(data.chooselist[slice_idx - 1], dtype=int)
        return estimate_atlas_plane_index(
            data.sample_volume,
            row,
            data.original_trans,
            data.atlas_template.shape,
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
        _sync_anchor_plane(data, state["slice"], plane)
        _refresh()

    def _update_status() -> None:
        idx = state["slice"]
        plane = _current_atlas_plane()
        _pmin, pmax = _atlas_plane_limits()
        anchor = _current_anchor(data, idx)
        confirmed = sum(anchor.confirmed for anchor in data.correspondence.anchors)
        caption = _chooselist_slice_label(data.chooselist, idx)
        flag = "confirmed" if anchor.confirmed else "pending"
        axis = _AXIS_NAMES.get(data.cut_axis, str(data.cut_axis))
        viewer.status = (
            f"Slice {idx}/{n_slices} ({caption}) | {flag} "
            f"| atlas plane {plane}/{pmax} | confirmed {confirmed}/{n_slices} "
            f"| AP axis {axis}"
        )

    def _refresh() -> None:
        idx = state["slice"]
        sample, atlas = _align_slices_pair(data, idx, atlas_plane=_current_atlas_plane())
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
            state["slice"] = max(1, min(n_slices, int(slice_index)))
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
        anchor = _current_anchor(data, idx)
        anchor.confirmed = True
        anchor.sample_index = int(data.chooselist[idx - 1, 0])
        anchor.atlas_plane = _current_atlas_plane()
        show_info(f"Confirmed slice {idx}: sample {anchor.sample_index} ↔ atlas plane {anchor.atlas_plane}")
        if idx < n_slices:
            _navigate_to(idx + 1, refocus_canvas=True)
        else:
            _refresh()

    @magicgui(
        slice_index={
            "min": 1,
            "max": n_slices,
            "step": 1,
            "label": "Slice # (AP anchor index)",
        },
        atlas_plane={
            "min": 1,
            "max": 100,
            "step": 1,
            "label": "Atlas plane along cut axis",
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

    @magicgui(call_button="◀  Previous slice")
    def previous_slice() -> None:
        _navigate_to(state["slice"] - 1)

    @magicgui(call_button="Next slice  ▶")
    def next_slice() -> None:
        _navigate_to(state["slice"] + 1)

    @magicgui(call_button="Confirm slice (Enter)")
    def confirm_slice() -> None:
        _confirm_and_advance()

    @magicgui(call_button="Save && Close")
    def save_controls() -> None:
        data.correspondence.source = "manual"
        data.correspondence.save(data.correspondence_path)
        show_info(f"Saved {data.correspondence_path}")
        QTimer.singleShot(0, viewer.close)

    def _previous_slice_key(_viewer) -> None:
        _navigate_to(state["slice"] - 1, refocus_canvas=True)

    def _next_slice_key(_viewer) -> None:
        _navigate_to(state["slice"] + 1, refocus_canvas=True)

    def _confirm_key(_viewer) -> None:
        _confirm_and_advance()

    def _atlas_plane_step(delta: int) -> None:
        plane = _current_atlas_plane() + int(delta)
        _set_atlas_plane(plane)
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

    viewer.window.add_dock_widget(navigation, area="right", name="Navigation")
    viewer.window.add_dock_widget(previous_slice, area="right", name="Previous slice")
    viewer.window.add_dock_widget(next_slice, area="right", name="Next slice")
    viewer.window.add_dock_widget(confirm_slice, area="right", name="Confirm")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")
    _refresh()
    _sync_navigation_widget()

    console.print(
        "[bold]Napari align-slices GUI[/bold] — sample (left), atlas (right). "
        "Scroll the atlas plane (PgUp/PgDn or wheel over atlas) until anatomy matches, "
        "then press [bold]Enter[/bold] or Confirm. "
        "Shortcuts: [bold]←[/bold]/[bold]→[/bold] change AP anchor slice."
    )
    napari.run()
    return data.correspondence_path
