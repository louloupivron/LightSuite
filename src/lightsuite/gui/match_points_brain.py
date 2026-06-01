"""Napari GUI for brain control-point matching (matchControlPoints_unified.m)."""

from __future__ import annotations

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


def _points_for_slice(session: ControlPointSession, slice_idx: int, panel: str) -> np.ndarray:
    store = (
        session.histology_control_points
        if panel == "sample"
        else session.atlas_control_points
    )
    pts = np.asarray(store[slice_idx - 1], dtype=float)
    if pts.size == 0:
        return np.zeros((0, 2))
    row = session.chooselist[slice_idx - 1] if session.chooselist else None
    if row is None:
        return pts[:, :2]
    axis = int(row[1])
    dims = [0, 1, 2]
    plot_axes = [d for d in dims if d != axis - 1]
    return pts[:, [plot_axes[1], plot_axes[0]]]


def _append_point(session: ControlPointSession, slice_idx: int, panel: str, xy: tuple[float, float]) -> None:
    row = np.asarray(session.chooselist[slice_idx - 1], dtype=int)
    axis = int(row[1])
    slice_index = int(row[0])
    point = np.zeros(4, dtype=float)
    dims = [0, 1, 2]
    plot_axes = [d for d in dims if d != axis - 1]
    point[plot_axes[1]] = xy[0]
    point[plot_axes[0]] = xy[1]
    point[axis - 1] = slice_index
    point[3] = np.nan
    target = (
        session.histology_control_points if panel == "sample" else session.atlas_control_points
    )
    target[slice_idx - 1] = [*target[slice_idx - 1], point.tolist()]


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
    state = {"slice": 1, "show_overlay": True, "_nav_syncing": False}

    viewer = napari.Viewer(title=f"LightSuite — {config.sample.name}")
    sample_layer = viewer.add_image(np.zeros((10, 10)), name="sample", colormap="gray")
    atlas_layer = viewer.add_image(np.zeros((10, 10)), name="atlas", colormap="gray")
    overlay_layer = viewer.add_image(np.zeros((10, 10)), name="overlay", colormap="red", opacity=0.35, visible=True)
    sample_pts = viewer.add_points(
        np.zeros((0, 2)),
        name="sample_points",
        face_color="yellow",
        size=8,
        ndim=2,
    )
    atlas_pts = viewer.add_points(
        np.zeros((0, 2)),
        name="atlas_points",
        face_color="cyan",
        size=8,
        ndim=2,
    )

    def _layout_atlas() -> None:
        h, w = sample_layer.data.shape
        atlas_layer.translate = (w + 20, 0)
        atlas_pts.translate = (w + 20, 0)

    def _refresh() -> None:
        idx = state["slice"]
        sample, atlas = slice_pair(data, idx)
        sample_layer.data = sample
        atlas_layer.data = atlas
        _layout_atlas()
        # Block data events: programmatic updates must not re-enter the point
        # change handlers (_sample_changed / _atlas_changed), which call _try_align
        # and would otherwise recurse until vispy transform updates overflow.
        with sample_pts.events.data.blocker(), atlas_pts.events.data.blocker():
            sample_pts.data = _points_for_slice(data.session, idx, "sample")
            atlas_pts.data = _points_for_slice(data.session, idx, "atlas")
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
        viewer.status = (
            f"Slice {idx}/{n_slices} ({caption}) | sample pts={n_s} atlas pts={n_a}"
        )
        viewer.reset_view()

    def _try_align() -> None:
        mse = data.session.update_manual_alignment(min_pairs=4)
        if mse is not None:
            show_info(f"Updated alignment fit (MSE={mse:.2f})")
        _refresh()

    @sample_pts.events.data.connect
    def _sample_changed(_event=None) -> None:
        idx = state["slice"]
        data.session.histology_control_points[idx - 1] = []
        for pt in sample_pts.data:
            _append_point(data.session, idx, "sample", (float(pt[0]), float(pt[1])))
        if len(data.session.histology_control_points[idx - 1]) == len(
            data.session.atlas_control_points[idx - 1]
        ):
            _try_align()

    @atlas_pts.events.data.connect
    def _atlas_changed(_event=None) -> None:
        idx = state["slice"]
        data.session.atlas_control_points[idx - 1] = []
        offset = atlas_layer.translate[0]
        for pt in atlas_pts.data:
            _append_point(data.session, idx, "atlas", (float(pt[0] - offset), float(pt[1])))
        if len(data.session.histology_control_points[idx - 1]) == len(
            data.session.atlas_control_points[idx - 1]
        ):
            _try_align()

    def _sync_navigation_widget() -> None:
        """Keep the spinbox in sync with state without re-firing navigation callbacks."""
        state["_nav_syncing"] = True
        try:
            navigation.slice_index.value = int(state["slice"])
            navigation.show_overlay.value = bool(state["show_overlay"])
        finally:
            state["_nav_syncing"] = False

    def _navigate_to(slice_index: int | None = None, *, show_overlay: bool | None = None) -> None:
        if show_overlay is not None:
            state["show_overlay"] = bool(show_overlay)
        if slice_index is not None:
            state["slice"] = max(1, min(n_slices, int(slice_index)))
        _sync_navigation_widget()
        _refresh()

    @magicgui(
        slice_index={
            "min": 1,
            "max": n_slices,
            "step": 1,
            "label": "Slice # (chooselist index)",
        },
        show_overlay={"label": "Show atlas boundary overlay on sample"},
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

    viewer.window.add_dock_widget(navigation, area="right", name="Navigation")
    viewer.window.add_dock_widget(previous_slice, area="right", name="Previous slice")
    viewer.window.add_dock_widget(next_slice, area="right", name="Next slice")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")
    viewer.window.add_dock_widget(clear_slice, area="right", name="Edit")
    _sync_navigation_widget()
    _refresh()

    console.print(
        "[bold]Napari control-point GUI[/bold] — sample (left), atlas (right). "
        "Use the Navigation panel (Previous / Next / Show slice); "
        "the Napari dimension slider does not change slices here."
    )
    napari.run()
    return data.session_path
