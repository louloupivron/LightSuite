"""Spinal cord control-point matching GUI (matchControlPointsSpine.m port)."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.cord_data import (
    CORD_ATLAS_PROVIDER,
    default_cord_session_path,
    load_cord_match_points_data,
    prepare_cord_match_points_session,
    resolve_atlas_plane_index,
    set_atlas_plane_index,
    slice_pair,
)
from lightsuite.gui.control_points import ControlPointSession, mark_session_saved_from_napari
from lightsuite.gui.slices import (
    layer_xy_from_slice_pixels,
    slice_pixels_from_layer_xy,
    volume_index_to_image,
)
from lightsuite.registration.warp import warp_volume_affine

console = Console()

CORRESPONDING_POINTS_JSON = "corresponding_points.json"
PANEL_GAP_X = 24
MIN_AFFINE_PAIRS = 16
TEXT_LABEL_COLOR = "white"
TEXT_LABEL_OFFSET = (0.0, 8.0)


def default_session_path(save_path: Path) -> Path:
    """Backward-compatible alias for ``default_cord_session_path``."""
    return default_cord_session_path(save_path)


def _plot_axes_for_row(chooserow: np.ndarray) -> list[int]:
    cut_axis = int(chooserow[1]) - 1
    return [d for d in range(3) if d != cut_axis]


def _raw_slice_shape(volume: np.ndarray, chooserow: np.ndarray) -> tuple[int, int]:
    return tuple(int(v) for v in volume_index_to_image(volume, chooserow).shape)


def _volume_points_to_layer_xy(
    points: list[list[float]],
    chooserow: np.ndarray,
    *,
    slice_shape: tuple[int, int],
) -> np.ndarray:
    if not points:
        return np.zeros((0, 2))
    pts = np.asarray(points, dtype=float)
    plot_axes = _plot_axes_for_row(chooserow)
    cut_axis = int(chooserow[1])
    return layer_xy_from_slice_pixels(
        pts[:, plot_axes[0]],
        pts[:, plot_axes[1]],
        slice_shape,
        cut_axis,
        CORD_ATLAS_PROVIDER,
    )


def _layer_xy_to_volume_point(
    xy: tuple[float, float],
    chooserow: np.ndarray,
    *,
    slice_shape: tuple[int, int],
    timestamp: float,
    plane_along_cut_axis: int | None = None,
) -> list[float]:
    plot_axes = _plot_axes_for_row(chooserow)
    cut_axis = int(chooserow[1]) - 1
    row, col = slice_pixels_from_layer_xy(
        np.asarray([xy], dtype=float),
        slice_shape,
        int(chooserow[1]),
        CORD_ATLAS_PROVIDER,
    )
    point = np.zeros(4, dtype=float)
    point[plot_axes[0]] = row[0]
    point[plot_axes[1]] = col[0]
    point[cut_axis] = float(
        plane_along_cut_axis if plane_along_cut_axis is not None else int(chooserow[0])
    )
    point[3] = timestamp
    return point.tolist()


def _configure_point_text(layer) -> None:
    if not hasattr(layer, "text"):
        return
    try:
        layer.text.color = TEXT_LABEL_COLOR
        layer.text.translation = np.array(TEXT_LABEL_OFFSET, dtype=float)
        layer.text.anchor = "center"
        layer.text.size = 12
    except (TypeError, ValueError, AttributeError):
        pass


def _slice_point_store(session: ControlPointSession, panel: str) -> list[list[list[float]]]:
    return (
        session.histology_control_points
        if panel == "sample"
        else session.atlas_control_points
    )


def _sync_store_from_layer(
    session: ControlPointSession,
    slice_idx: int,
    panel: str,
    layer_xy: np.ndarray,
    *,
    slice_shape: tuple[int, int],
    atlas_plane: int | None = None,
) -> None:
    chooserow = np.asarray(session.chooselist[slice_idx - 1], dtype=int)
    store = _slice_point_store(session, panel)
    existing = store[slice_idx - 1]
    plane = atlas_plane if panel == "atlas" else None
    updated: list[list[float]] = []
    for i, pt in enumerate(layer_xy):
        ts = float(existing[i][3]) if i < len(existing) else time.time()
        updated.append(
            _layer_xy_to_volume_point(
                (float(pt[0]), float(pt[1])),
                chooserow,
                slice_shape=slice_shape,
                timestamp=ts,
                plane_along_cut_axis=plane,
            )
        )
    store[slice_idx - 1] = updated


def _pair_labels(n_points: int) -> list[str]:
    return [str(i + 1) for i in range(n_points)]


def _pair_status(n_sample: int, n_atlas: int) -> str:
    n_pairs = min(n_sample, n_atlas)
    if n_sample == n_atlas:
        return f"pairs={n_pairs} (matched)"
    if n_sample > n_atlas:
        return f"pairs={n_pairs} | place atlas point #{n_atlas + 1}"
    return f"pairs={n_pairs} | place sample point #{n_sample + 1}"


def _boundary_overlay(
    warped_annotation: np.ndarray,
    chooserow: np.ndarray,
) -> np.ndarray:
    from scipy.ndimage import convolve

    from lightsuite.gui.slices import prepare_display_slice

    ann_slice = volume_index_to_image(warped_annotation, chooserow)
    edges = ann_slice.astype(float)
    kernel = np.ones((3, 3)) / 9.0
    blurred = convolve(edges, kernel, mode="constant")
    overlay = (np.round(blurred) != edges).astype(float)
    return prepare_display_slice(overlay, int(chooserow[1]), CORD_ATLAS_PROVIDER)


def _chooselist_slice_label(chooselist: np.ndarray, slice_idx: int) -> str:
    row = np.asarray(chooselist[slice_idx - 1], dtype=int)
    return f"sample plane {row[0]} / {row[1]}-axis cut"


def run_spinal_match_points(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
) -> Path:
    """Launch Napari control-point matcher; returns saved session path."""
    if headless:
        return prepare_cord_match_points_session(config)

    try:
        import napari
        from magicgui import magicgui
        from magicgui.widgets import Label
        from napari.utils.notifications import show_info
        from qtpy.QtCore import QTimer
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    data = load_cord_match_points_data(config)
    n_slices = int(data.chooselist.shape[0])
    state = {
        "slice": 1,
        "show_overlay": True,
        "_nav_syncing": False,
        "_view_shape": None,
        "_warped_annotation_key": None,
        "_warped_annotation": None,
    }

    viewer = napari.Viewer(title=f"LightSuite spinal — {config.sample.name}")
    viewer.dims.ndisplay = 2
    sample_layer = viewer.add_image(np.zeros((10, 10)), name="sample", colormap="gray")
    atlas_layer = viewer.add_image(np.zeros((10, 10)), name="atlas", colormap="gray")
    overlay_layer = viewer.add_image(
        np.zeros((10, 10)),
        name="overlay",
        colormap="red",
        opacity=0.35,
        visible=True,
        blending="additive",
    )
    sample_pts = viewer.add_points(
        np.zeros((0, 2)),
        name="sample_points",
        face_color="yellow",
        size=3,
        ndim=2,
    )
    atlas_pts = viewer.add_points(
        np.zeros((0, 2)),
        name="atlas_points",
        face_color="cyan",
        size=3,
        ndim=2,
    )
    _configure_point_text(sample_pts)
    _configure_point_text(atlas_pts)

    points_summary = Label(
        label="Session totals",
        value="Total points: 0 matched pairs (0 on sample, 0 on atlas)",
    )

    def _apply_layer_points(layer, xy: np.ndarray) -> None:
        labels = _pair_labels(int(xy.shape[0]))
        with layer.events.data.blocker():
            layer.data = xy
            if hasattr(layer, "text"):
                try:
                    layer.text = labels
                except (TypeError, ValueError, AttributeError):
                    pass
            _configure_point_text(layer)

    def _warped_annotation_volume(matrix: np.ndarray) -> np.ndarray:
        key = np.asarray(matrix, dtype=float).tobytes()
        if state["_warped_annotation_key"] == key and state["_warped_annotation"] is not None:
            return state["_warped_annotation"]
        warped = warp_volume_affine(
            data.atlas_annotation.astype(np.float32),
            matrix,
            data.sample_volume.shape,
            order=0,
            point_coords="xyz",
        )
        state["_warped_annotation_key"] = key
        state["_warped_annotation"] = warped
        return warped

    def _layout_panels() -> None:
        _h, w = sample_layer.data.shape
        sample_layer.translate = (0.0, 0.0)
        overlay_layer.translate = (0.0, 0.0)
        sample_pts.translate = (0.0, 0.0)
        atlas_offset = (0.0, float(w + PANEL_GAP_X))
        atlas_layer.translate = atlas_offset
        atlas_pts.translate = atlas_offset

    def _atlas_plane_limits() -> tuple[int, int]:
        row = np.asarray(data.chooselist[state["slice"] - 1], dtype=int)
        from lightsuite.gui.cord_data import atlas_cut_axis_size

        nmax = atlas_cut_axis_size(data.atlas_template.shape, row)
        return 1, nmax

    def _resolved_atlas_plane(plane: int | None = None) -> int:
        pmin, pmax = _atlas_plane_limits()
        if plane is None:
            plane = state.get("_atlas_plane")
        if plane is None:
            plane = resolve_atlas_plane_index(data, state["slice"])
        return int(np.clip(int(plane), pmin, pmax))

    def _set_atlas_plane(plane: int, *, persist: bool = True) -> None:
        plane = _resolved_atlas_plane(plane)
        if persist:
            set_atlas_plane_index(data.session, state["slice"], plane)
        _refresh(atlas_plane=plane)

    def _update_point_count_display() -> None:
        matched, total_sample, total_atlas = data.session.point_counts()
        points_summary.value = (
            f"Total points: {matched} matched pairs "
            f"({total_sample} on sample, {total_atlas} on atlas)"
        )

    def _update_status() -> None:
        idx = state["slice"]
        plane = _resolved_atlas_plane()
        _pmin, pmax = _atlas_plane_limits()
        n_s = len(data.session.histology_control_points[idx - 1])
        n_a = len(data.session.atlas_control_points[idx - 1])
        matched, total_sample, total_atlas = data.session.point_counts()
        caption = _chooselist_slice_label(data.chooselist, idx)
        overlay_flag = "on" if state["show_overlay"] else "off"
        pair_info = _pair_status(n_s, n_a)
        _update_point_count_display()
        fit_state = (
            "affine fit active"
            if matched >= MIN_AFFINE_PAIRS
            else f"coarse align (need {MIN_AFFINE_PAIRS} pairs for affine, have {matched})"
        )
        viewer.status = (
            f"Slice {idx}/{n_slices} ({caption}) | {pair_info} "
            f"| all slices: {matched} pairs ({total_sample} sample / {total_atlas} atlas) "
            f"| {fit_state} | atlas plane {plane}/{pmax} | overlay {overlay_flag} (O)"
        )

    def _refresh(*, atlas_plane: int | None = None) -> None:
        idx = state["slice"]
        state["_atlas_plane"] = _resolved_atlas_plane(atlas_plane)
        sample, atlas = slice_pair(data, idx, atlas_plane=state["_atlas_plane"])
        sample_layer.data = sample
        atlas_layer.data = atlas
        _layout_panels()
        chooserow = np.asarray(data.chooselist[idx - 1], dtype=int)
        slice_shape = _raw_slice_shape(data.sample_volume, chooserow)
        with sample_pts.events.data.blocker(), atlas_pts.events.data.blocker():
            _apply_layer_points(
                sample_pts,
                _volume_points_to_layer_xy(
                    data.session.histology_control_points[idx - 1],
                    chooserow,
                    slice_shape=slice_shape,
                ),
            )
            _apply_layer_points(
                atlas_pts,
                _volume_points_to_layer_xy(
                    data.session.atlas_control_points[idx - 1],
                    chooserow,
                    slice_shape=slice_shape,
                ),
            )
        matrix = np.asarray(data.session.atlas2histology_tform, dtype=float)
        if state["show_overlay"]:
            overlay_layer.data = _boundary_overlay(
                _warped_annotation_volume(matrix),
                chooserow,
            )
            overlay_layer.visible = True
        else:
            overlay_layer.visible = False
        _update_status()
        if state["_view_shape"] != sample.shape:
            viewer.reset_view()
            state["_view_shape"] = sample.shape

    def _try_align() -> None:
        mse = data.session.update_manual_alignment(
            min_pairs=MIN_AFFINE_PAIRS,
            fallback_tform=np.eye(4),
        )
        if mse is not None:
            show_info(f"Updated alignment fit (MSE={mse:.2f})")
        else:
            matched, _total_s, _total_a = data.session.point_counts()
            remaining = MIN_AFFINE_PAIRS - matched
            if remaining > 0:
                show_info(
                    f"{matched}/{MIN_AFFINE_PAIRS} matched pairs — "
                    f"add {remaining} more (across slices) before the affine fit runs"
                )
        _refresh()

    def _on_panel_points_changed(panel: str) -> None:
        idx = state["slice"]
        layer = sample_pts if panel == "sample" else atlas_pts
        plane = _resolved_atlas_plane() if panel == "atlas" else None
        chooserow = np.asarray(data.chooselist[idx - 1], dtype=int)
        _sync_store_from_layer(
            data.session,
            idx,
            panel,
            np.asarray(layer.data, dtype=float),
            slice_shape=_raw_slice_shape(data.sample_volume, chooserow),
            atlas_plane=plane,
        )
        n_s = len(data.session.histology_control_points[idx - 1])
        n_a = len(data.session.atlas_control_points[idx - 1])
        if n_s == n_a:
            _try_align()

    @sample_pts.events.data.connect
    def _sample_changed(_event=None) -> None:
        _on_panel_points_changed("sample")

    @atlas_pts.events.data.connect
    def _atlas_changed(_event=None) -> None:
        _on_panel_points_changed("atlas")

    def _sync_navigation_widget() -> None:
        state["_nav_syncing"] = True
        try:
            navigation.slice_index.value = int(state["slice"])
            navigation.show_overlay.value = bool(state["show_overlay"])
            pmin, pmax = _atlas_plane_limits()
            navigation.atlas_plane.min = pmin
            navigation.atlas_plane.max = pmax
            navigation.atlas_plane.value = _resolved_atlas_plane()
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
        show_overlay: bool | None = None,
        refocus_canvas: bool = False,
    ) -> None:
        if show_overlay is not None:
            state["show_overlay"] = bool(show_overlay)
        if slice_index is not None:
            new_slice = max(1, min(n_slices, int(slice_index)))
            if new_slice != state["slice"]:
                state["_atlas_plane"] = None
            state["slice"] = new_slice
        if refocus_canvas:
            _refresh()

            def _after_keyboard_nav() -> None:
                _sync_navigation_widget()
                _refocus_canvas()

            QTimer.singleShot(0, _after_keyboard_nav)
        else:
            _refresh()
            _sync_navigation_widget()

    @magicgui(
        slice_index={
            "min": 1,
            "max": n_slices,
            "step": 1,
            "label": "Slice # (chooselist index)",
        },
        atlas_plane={
            "min": 1,
            "max": 100,
            "step": 1,
            "label": "Atlas plane along cut axis",
        },
        show_overlay={"label": "Atlas boundary overlay on sample (shortcut: O)"},
        call_button="Show slice",
    )
    def navigation(slice_index: int = 1, atlas_plane: int = 1, show_overlay: bool = True) -> None:
        _navigate_to(slice_index, show_overlay=show_overlay)
        _set_atlas_plane(atlas_plane)

    @navigation.atlas_plane.changed.connect
    def _atlas_plane_widget_changed() -> None:
        if state["_nav_syncing"]:
            return
        _set_atlas_plane(navigation.atlas_plane.value)

    @magicgui(call_button="◀  Previous slice")
    def previous_slice() -> None:
        _navigate_to(state["slice"] - 1)

    @magicgui(call_button="Next slice  ▶")
    def next_slice() -> None:
        _navigate_to(state["slice"] + 1)

    @navigation.slice_index.changed.connect
    def _slice_index_widget_changed() -> None:
        if state["_nav_syncing"]:
            return
        _navigate_to(navigation.slice_index.value)

    @navigation.show_overlay.changed.connect
    def _overlay_widget_changed() -> None:
        if state["_nav_syncing"]:
            return
        _navigate_to(show_overlay=navigation.show_overlay.value)

    @magicgui(call_button="Save && Close")
    def save_controls() -> None:
        mark_session_saved_from_napari(data.session)
        data.session.save(data.session_path)
        show_info(f"Saved {data.session_path}")
        QTimer.singleShot(0, viewer.close)

    @magicgui(call_button="Clear slice points")
    def clear_slice() -> None:
        idx = state["slice"]
        data.session.histology_control_points[idx - 1] = []
        data.session.atlas_control_points[idx - 1] = []
        _refresh()

    def _toggle_overlay(_viewer) -> None:
        _navigate_to(show_overlay=not state["show_overlay"], refocus_canvas=True)

    def _previous_slice_key(_viewer) -> None:
        _navigate_to(state["slice"] - 1, refocus_canvas=True)

    def _next_slice_key(_viewer) -> None:
        _navigate_to(state["slice"] + 1, refocus_canvas=True)

    def _point_timestamp(point: list[float]) -> float:
        ts = point[3]
        if ts is None or (isinstance(ts, float) and np.isnan(ts)):
            return float("-inf")
        return float(ts)

    def _atlas_plane_step(delta: int) -> None:
        plane = _resolved_atlas_plane() + int(delta)
        _set_atlas_plane(plane)
        _sync_navigation_widget()
        _refocus_canvas()

    def _atlas_plane_up(_viewer) -> None:
        _atlas_plane_step(1)

    def _atlas_plane_down(_viewer) -> None:
        _atlas_plane_step(-1)

    def _delete_last_point_key(_viewer) -> None:
        idx = state["slice"]
        hist = data.session.histology_control_points[idx - 1]
        atlas = data.session.atlas_control_points[idx - 1]
        t_hist = _point_timestamp(hist[-1]) if hist else float("-inf")
        t_atlas = _point_timestamp(atlas[-1]) if atlas else float("-inf")
        if t_hist >= t_atlas and hist:
            hist.pop()
        elif atlas:
            atlas.pop()
        _navigate_to(refocus_canvas=True)

    viewer.bind_key("O", _toggle_overlay)
    viewer.bind_key("Left", _previous_slice_key, overwrite=True)
    viewer.bind_key("Right", _next_slice_key, overwrite=True)
    viewer.bind_key("Backspace", _delete_last_point_key, overwrite=True)
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

    viewer.window.add_dock_widget(points_summary, area="right", name="Point counts")
    viewer.window.add_dock_widget(navigation, area="right", name="Navigation")
    viewer.window.add_dock_widget(previous_slice, area="right", name="Previous slice")
    viewer.window.add_dock_widget(next_slice, area="right", name="Next slice")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")
    viewer.window.add_dock_widget(clear_slice, area="right", name="Edit")
    _refresh()
    _sync_navigation_widget()

    console.print(
        "[bold]Napari spinal control-point GUI[/bold] — sample (left), atlas (right). "
        "Numbered markers are matched pairs (1↔1, 2↔2, …). "
        "Scroll the wheel over the atlas (or PgUp/PgDn) to adjust the atlas plane. "
        "Shortcuts: [bold]←[/bold]/[bold]→[/bold] chooselist slices, [bold]O[/bold] overlay, "
        "[bold]Backspace[/bold] undo last point."
    )
    napari.run()
    if not data.session_path.is_file():
        mark_session_saved_from_napari(data.session)
        data.session.save(data.session_path)
    return data.session_path
