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
from lightsuite.gui.slice_correspondence import VOLUME_AXES
from lightsuite.gui.slices import volume_index_to_image

console = Console()

_AXIS_NAMES = {1: "Y", 2: "X", 3: "Z"}


def _align_slices_pair(
    data: BrainAlignSlicesData,
    cut_axis: int,
    slice_idx: int,
    *,
    atlas_plane: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample and atlas images for one align-slices chooselist entry."""
    row = np.asarray(data.chooselist_for_axis(cut_axis)[slice_idx - 1], dtype=int)
    sample = _normalize_display(volume_index_to_image(data.sample_volume, row))
    atlas_row = chooserow_with_atlas_plane(row, atlas_plane)
    atlas = _normalize_display(volume_index_to_image(data.atlas_template, atlas_row))
    return sample, atlas


def _current_anchor(data: BrainAlignSlicesData, cut_axis: int, slice_idx: int):
    return data.correspondence.anchors_for_axis(cut_axis)[slice_idx - 1]


def _sync_anchor_plane(data: BrainAlignSlicesData, cut_axis: int, slice_idx: int, plane: int) -> None:
    anchor = _current_anchor(data, cut_axis, slice_idx)
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
    state = {
        "axis": int(data.axis_order[0]),
        "slice": 1,
        "_nav_syncing": False,
        "_view_shape": None,
        "_atlas_plane": None,
    }

    def _n_slices() -> int:
        return int(data.chooselist_for_axis(state["axis"]).shape[0])

    viewer = napari.Viewer(title=f"LightSuite align-slices — {config.sample.name}")
    viewer.dims.ndisplay = 2
    sample_layer = viewer.add_image(np.zeros((10, 10)), name="sample", colormap="gray")
    atlas_layer = viewer.add_image(np.zeros((10, 10)), name="atlas", colormap="gray")

    def _layout_panels() -> None:
        _h, w = sample_layer.data.shape
        sample_layer.translate = (0.0, 0.0)
        atlas_layer.translate = (0.0, float(w + PANEL_GAP_X))

    def _atlas_plane_limits() -> tuple[int, int]:
        row = np.asarray(data.chooselist_for_axis(state["axis"])[state["slice"] - 1], dtype=int)
        return 1, atlas_cut_axis_size(data.atlas_template.shape, row)

    def _default_plane(slice_idx: int) -> int:
        axis = state["axis"]
        anchor = _current_anchor(data, axis, slice_idx)
        if anchor.atlas_plane > 0:
            return int(anchor.atlas_plane)
        row = np.asarray(data.chooselist_for_axis(axis)[slice_idx - 1], dtype=int)
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
        _sync_anchor_plane(data, state["axis"], state["slice"], plane)
        _refresh()

    def _update_status() -> None:
        idx = state["slice"]
        axis = state["axis"]
        plane = _current_atlas_plane()
        _pmin, pmax = _atlas_plane_limits()
        anchors = data.correspondence.anchors_for_axis(axis)
        confirmed = sum(anchor.confirmed for anchor in anchors)
        caption = _chooselist_slice_label(data.chooselist_for_axis(axis), idx)
        flag = "confirmed" if anchors[idx - 1].confirmed else "pending"
        axis_name = _AXIS_NAMES.get(axis, str(axis))
        total_confirmed = data.correspondence.confirmed_axis_count()
        viewer.status = (
            f"Axis {axis_name} ({axis}/3) | slice {idx}/{_n_slices()} ({caption}) | {flag} "
            f"| atlas plane {plane}/{pmax} | axis confirmed {confirmed}/{_n_slices()} "
            f"| axes done {total_confirmed}/3"
        )

    def _refresh() -> None:
        idx = state["slice"]
        axis = state["axis"]
        sample, atlas = _align_slices_pair(
            data,
            axis,
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
            navigation.cut_axis.value = int(state["axis"])
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

    def _navigate_to(
        slice_index: int | None = None,
        *,
        cut_axis: int | None = None,
        refocus_canvas: bool = False,
    ) -> None:
        if cut_axis is not None:
            state["axis"] = int(cut_axis)
            state["_atlas_plane"] = None
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
        axis = state["axis"]
        anchor = _current_anchor(data, axis, idx)
        anchor.confirmed = True
        anchor.sample_index = int(data.chooselist_for_axis(axis)[idx - 1, 0])
        anchor.atlas_plane = _current_atlas_plane()
        axis_name = _AXIS_NAMES.get(axis, str(axis))
        show_info(
            f"Confirmed {axis_name} slice {idx}: "
            f"sample {anchor.sample_index} ↔ atlas plane {anchor.atlas_plane}"
        )
        if idx < _n_slices():
            _navigate_to(idx + 1, refocus_canvas=True)
        else:
            next_axis = _next_incomplete_axis(axis)
            if next_axis is not None:
                show_info(f"Axis {axis_name} complete — switch to axis {_AXIS_NAMES.get(next_axis, next_axis)}")
                _navigate_to(1, cut_axis=next_axis, refocus_canvas=True)
            else:
                _refresh()

    def _next_incomplete_axis(after_axis: int) -> int | None:
        order = list(data.axis_order)
        start = order.index(after_axis) + 1 if after_axis in order else 0
        for axis in order[start:] + order[:start]:
            anchors = data.correspondence.anchors_for_axis(axis)
            if anchors and not all(anchor.confirmed for anchor in anchors):
                return axis
        return None

    @magicgui(
        cut_axis={"choices": list(VOLUME_AXES), "label": "Volume axis (1=Y, 2=X, 3=Z)"},
        slice_index={
            "min": 1,
            "max": 20,
            "step": 1,
            "label": "Anchor slice #",
        },
        atlas_plane={
            "min": 1,
            "max": 100,
            "step": 1,
            "label": "Atlas plane along cut axis",
        },
        call_button="Show slice",
    )
    def navigation(cut_axis: int = 1, slice_index: int = 1, atlas_plane: int = 1) -> None:
        _navigate_to(slice_index, cut_axis=cut_axis)
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

    @navigation.cut_axis.changed.connect
    def _cut_axis_widget_changed() -> None:
        if state["_nav_syncing"]:
            return
        _navigate_to(1, cut_axis=navigation.cut_axis.value)

    @magicgui(call_button="◀  Previous slice")
    def previous_slice() -> None:
        _navigate_to(state["slice"] - 1)

    @magicgui(call_button="Next slice  ▶")
    def next_slice() -> None:
        _navigate_to(state["slice"] + 1)

    @magicgui(call_button="Confirm slice (Enter)")
    def confirm_slice() -> None:
        _confirm_and_advance()

    @magicgui(call_button="Next axis ▶")
    def next_axis() -> None:
        order = list(data.axis_order)
        idx = order.index(state["axis"])
        next_ax = order[(idx + 1) % len(order)]
        _navigate_to(1, cut_axis=next_ax, refocus_canvas=True)

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
    viewer.window.add_dock_widget(next_axis, area="right", name="Next axis")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")
    _refresh()
    _sync_navigation_widget()

    console.print(
        "[bold]Napari align-slices GUI[/bold] — sample (left), atlas (right). "
        "Align all three volume axes (Y, X, Z): ~20 anchors per axis. "
        "Use the [bold]Volume axis[/bold] control to switch axes. "
        "Scroll the atlas plane (PgUp/PgDn or wheel over atlas), then "
        "[bold]Enter[/bold] or Confirm. After the last slice on an axis, "
        "you are prompted to continue on the next axis."
    )
    napari.run()
    return data.correspondence_path
