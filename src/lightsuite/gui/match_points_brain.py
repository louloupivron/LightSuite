"""Napari GUI for brain control-point matching (matchControlPoints_unified.m)."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.affine import fit_affine_transform, transform_points
from lightsuite.gui.brain_data import (
    BrainMatchPointsData,
    load_brain_match_points_data,
    prepare_brain_match_points_session,
    slice_pair,
)
from lightsuite.gui.control_points import ControlPointSession
from lightsuite.gui.slices import volume_index_to_image

console = Console()

# Horizontal gap between sample+overlay (left) and atlas (right), in pixels (dim 1 / X).
PANEL_GAP_X = 24
# Pair number labels: white text, nudged right of each marker (row, col) in layer data coords.
TEXT_LABEL_COLOR = "white"
TEXT_LABEL_OFFSET = (0.0, 8.0)


def _plot_axes_for_row(chooserow: np.ndarray) -> list[int]:
    """0-based volume axes shown in the 2D viewer (MATLAB toplot)."""
    cut_axis = int(chooserow[1]) - 1
    return [d for d in range(3) if d != cut_axis]


def _volume_points_to_layer_xy(points: list[list[float]], chooserow: np.ndarray) -> np.ndarray:
    """Map stored 4-column points to napari layer (row, col) coordinates."""
    if not points:
        return np.zeros((0, 2))
    pts = np.asarray(points, dtype=float)
    plot_axes = _plot_axes_for_row(chooserow)
    return pts[:, [plot_axes[1], plot_axes[0]]]


def _layer_xy_to_volume_point(
    xy: tuple[float, float],
    chooserow: np.ndarray,
    *,
    timestamp: float,
) -> list[float]:
    """Map a 2D layer click/drag position back to MATLAB-style [x, y, z, t] storage."""
    plot_axes = _plot_axes_for_row(chooserow)
    cut_axis = int(chooserow[1]) - 1
    point = np.zeros(4, dtype=float)
    point[plot_axes[1]] = xy[0]
    point[plot_axes[0]] = xy[1]
    point[cut_axis] = int(chooserow[0])
    point[3] = timestamp
    return point.tolist()


def _pair_labels(n_points: int) -> list[str]:
    return [str(i + 1) for i in range(n_points)]


def _configure_point_text(layer) -> None:
    """Style numbered pair labels (white, offset right of the marker)."""
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
) -> None:
    """Persist napari layer coordinates; preserve timestamps when dragging existing points."""
    chooserow = np.asarray(session.chooselist[slice_idx - 1], dtype=int)
    store = _slice_point_store(session, panel)
    existing = store[slice_idx - 1]
    updated: list[list[float]] = []
    for i, pt in enumerate(layer_xy):
        ts = float(existing[i][3]) if i < len(existing) else time.time()
        updated.append(
            _layer_xy_to_volume_point((float(pt[0]), float(pt[1])), chooserow, timestamp=ts)
        )
    store[slice_idx - 1] = updated


def _pair_status(n_sample: int, n_atlas: int) -> str:
    n_pairs = min(n_sample, n_atlas)
    if n_sample == n_atlas:
        return f"pairs={n_pairs} (matched)"
    if n_sample > n_atlas:
        return f"pairs={n_pairs} | place atlas point #{n_atlas + 1}"
    return f"pairs={n_pairs} | place sample point #{n_sample + 1}"


def _boundary_overlay(
    annotation: np.ndarray,
    chooserow: np.ndarray,
    matrix: np.ndarray,
    sample_shape: tuple[int, int],
) -> np.ndarray:
    """Warp annotation slice boundaries onto sample slice (approximate overlay)."""
    from scipy.ndimage import convolve

    ann_slice = volume_index_to_image(annotation, chooserow)
    warped_vol = transform_points(
        np.column_stack(
            [
                np.repeat(np.arange(ann_slice.shape[0]), ann_slice.shape[1]),
                np.tile(np.arange(ann_slice.shape[1]), ann_slice.shape[0]),
                np.zeros(ann_slice.size),
            ]
        ),
        matrix,
    )
    # Fallback: show atlas edge map on sample grid
    edges = ann_slice.astype(float)
    kernel = np.ones((3, 3)) / 9.0
    blurred = convolve(edges, kernel, mode="constant")
    boundary = (np.round(blurred) != edges).astype(float)
    if boundary.shape != sample_shape:
        from skimage.transform import resize

        boundary = resize(boundary, sample_shape, order=0, preserve_range=True, anti_aliasing=False)
    return boundary


def _chooselist_slice_label(chooselist: np.ndarray, slice_idx: int) -> str:
    """Human-readable description of one chooselist row."""
    row = np.asarray(chooselist[slice_idx - 1], dtype=int)
    axis_names = {1: "Y", 2: "X", 3: "Z"}
    axis = axis_names.get(int(row[1]), str(row[1]))
    return f"volume index {row[0]}, cut along axis {axis}"


def run_brain_match_points(config: BrainPipelineConfig, *, headless: bool = False) -> Path:
    """Launch napari control-point matcher; returns saved session path."""
    if headless:
        return prepare_brain_match_points_session(config)

    try:
        import napari
        from magicgui import magicgui
        from napari.utils.notifications import show_info
        from qtpy.QtCore import QTimer
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    data = load_brain_match_points_data(config)
    n_slices = int(data.chooselist.shape[0])
    state = {"slice": 1, "show_overlay": True, "_nav_syncing": False, "_view_shape": None}

    viewer = napari.Viewer(title=f"LightSuite — {config.sample.name}")
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
        size=10,
        ndim=2,
    )
    atlas_pts = viewer.add_points(
        np.zeros((0, 2)),
        name="atlas_points",
        face_color="cyan",
        size=10,
        ndim=2,
    )
    _configure_point_text(sample_pts)
    _configure_point_text(atlas_pts)

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

    def _layout_panels() -> None:
        """Place sample + overlay on the left, atlas + atlas points on the right."""
        _h, w = sample_layer.data.shape
        sample_layer.translate = (0.0, 0.0)
        overlay_layer.translate = (0.0, 0.0)
        sample_pts.translate = (0.0, 0.0)
        atlas_offset = (0.0, float(w + PANEL_GAP_X))
        atlas_layer.translate = atlas_offset
        atlas_pts.translate = atlas_offset

    def _refresh() -> None:
        idx = state["slice"]
        sample, atlas = slice_pair(data, idx)
        sample_layer.data = sample
        atlas_layer.data = atlas
        _layout_panels()
        # Block data events: programmatic updates must not re-enter the point
        # change handlers (_sample_changed / _atlas_changed), which call _try_align
        # and would otherwise recurse until vispy transform updates overflow.
        chooserow = np.asarray(data.chooselist[idx - 1], dtype=int)
        with sample_pts.events.data.blocker(), atlas_pts.events.data.blocker():
            _apply_layer_points(
                sample_pts,
                _volume_points_to_layer_xy(data.session.histology_control_points[idx - 1], chooserow),
            )
            _apply_layer_points(
                atlas_pts,
                _volume_points_to_layer_xy(data.session.atlas_control_points[idx - 1], chooserow),
            )
        matrix = np.asarray(data.session.atlas2histology_tform, dtype=float)
        if state["show_overlay"]:
            overlay_layer.data = _boundary_overlay(
                data.atlas_annotation,
                data.chooselist[idx - 1],
                matrix,
                sample.shape,
            )
            overlay_layer.visible = True
        else:
            overlay_layer.visible = False
        n_s = len(data.session.histology_control_points[idx - 1])
        n_a = len(data.session.atlas_control_points[idx - 1])
        caption = _chooselist_slice_label(data.chooselist, idx)
        overlay_flag = "on" if state["show_overlay"] else "off"
        pair_info = _pair_status(n_s, n_a)
        viewer.status = (
            f"Slice {idx}/{n_slices} ({caption}) | {pair_info} "
            f"| overlay {overlay_flag} (O) | ←/→ slice | Backspace undo"
        )
        if state["_view_shape"] != sample.shape:
            viewer.reset_view()
            state["_view_shape"] = sample.shape

    def _try_align() -> None:
        mse = data.session.update_manual_alignment(min_pairs=4)
        if mse is not None:
            show_info(f"Updated alignment fit (MSE={mse:.2f})")
        _refresh()

    def _on_panel_points_changed(panel: str) -> None:
        idx = state["slice"]
        layer = sample_pts if panel == "sample" else atlas_pts
        _sync_store_from_layer(data.session, idx, panel, np.asarray(layer.data, dtype=float))
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
        """Keep the spinbox in sync with state without re-firing navigation callbacks."""
        state["_nav_syncing"] = True
        try:
            navigation.slice_index.value = int(state["slice"])
            navigation.show_overlay.value = bool(state["show_overlay"])
        finally:
            state["_nav_syncing"] = False

    def _refocus_canvas() -> None:
        """Keep arrow/O shortcuts working after dock widgets take focus."""
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
            state["slice"] = max(1, min(n_slices, int(slice_index)))
        if refocus_canvas:
            # Refresh first; sync the spinbox then return focus to the canvas so
            # repeated ←/→ work without clicking back into the image.
            _refresh()

            def _after_keyboard_nav() -> None:
                _sync_navigation_widget()
                _refocus_canvas()

            QTimer.singleShot(0, _after_keyboard_nav)
        else:
            _sync_navigation_widget()
            _refresh()

    @magicgui(
        slice_index={
            "min": 1,
            "max": n_slices,
            "step": 1,
            "label": "Slice # (chooselist index)",
        },
        show_overlay={"label": "Atlas boundary overlay on sample (shortcut: O)"},
        call_button="Show slice",
    )
    def navigation(slice_index: int = 1, show_overlay: bool = True) -> None:
        _navigate_to(slice_index, show_overlay=show_overlay)

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

    def _delete_last_point_key(_viewer) -> None:
        """Remove the most recently placed point (MATLAB Backspace)."""
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

    # Napari normalizes key names; O and o are the same binding.
    viewer.bind_key("O", _toggle_overlay)
    viewer.bind_key("Left", _previous_slice_key, overwrite=True)
    viewer.bind_key("Right", _next_slice_key, overwrite=True)
    viewer.bind_key("Backspace", _delete_last_point_key, overwrite=True)

    viewer.window.add_dock_widget(navigation, area="right", name="Navigation")
    viewer.window.add_dock_widget(previous_slice, area="right", name="Previous slice")
    viewer.window.add_dock_widget(next_slice, area="right", name="Next slice")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")
    viewer.window.add_dock_widget(clear_slice, area="right", name="Edit")
    _sync_navigation_widget()
    _refresh()

    console.print(
        "[bold]Napari control-point GUI[/bold] — sample (left), atlas (right). "
        "Numbered markers are matched pairs (1↔1, 2↔2, …). "
        "Add points in order on each side; shortcuts: "
        "[bold]←[/bold]/[bold]→[/bold] slices, [bold]O[/bold] overlay, [bold]Backspace[/bold] undo last point."
    )
    napari.run()
    return data.session_path
