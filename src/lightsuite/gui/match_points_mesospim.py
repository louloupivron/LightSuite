"""Napari GUI for mesoSPIM overview / ROI landmark matching."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.gui.mesospim_data import (
    load_mesospim_match_points_data,
    prepare_mesospim_match_points_session,
)
from lightsuite.mesospim.config_models import MesospimPipelineConfig
from lightsuite.mesospim.landmark_geometry import (
    fit_landmark_transform,
    update_landmark_session_fit,
)

console = Console()

PANEL_GAP_X = 24
TEXT_LABEL_COLOR = "white"
TEXT_LABEL_OFFSET = (0.0, 8.0)


def _pair_labels(n_points: int) -> list[str]:
    return [str(i + 1) for i in range(n_points)]


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


def _layer_xy_from_zyx(points_zyx: list[list[float]], z_index: int) -> np.ndarray:
    """Map stored ZYX points on the current Z slice to napari (row, col) = (Y, X)."""
    if not points_zyx:
        return np.zeros((0, 2))
    rows: list[list[float]] = []
    for pt in points_zyx:
        z, y, x = (float(v) for v in pt[:3])
        if int(round(z)) == int(z_index):
            rows.append([y, x])
    if not rows:
        return np.zeros((0, 2))
    return np.asarray(rows, dtype=float)


def _zyx_from_layer_xy(
    layer_xy: np.ndarray,
    z_index: int,
    existing_points: list[list[float]],
) -> list[list[float]]:
    """Replace points on the current Z slice; keep points on other slices."""
    kept = [pt for pt in existing_points if int(round(float(pt[0]))) != int(z_index)]
    new_points = [
        [float(z_index), float(row), float(col)] for row, col in np.asarray(layer_xy, dtype=float)
    ]
    return kept + new_points


def _pair_status(n_overview: int, n_roi: int) -> str:
    n_pairs = min(n_overview, n_roi)
    if n_overview == n_roi:
        return f"pairs={n_pairs} (matched)"
    if n_overview > n_roi:
        return f"pairs={n_pairs} | place ROI point #{n_roi + 1}"
    return f"pairs={n_pairs} | place overview point #{n_overview + 1}"


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


def run_mesospim_match_points(cfg: MesospimPipelineConfig, *, headless: bool = False) -> Path:
    """Launch Napari landmark matcher; returns saved session path."""
    if headless:
        return prepare_mesospim_match_points_session(cfg)

    try:
        import napari
        from magicgui import magicgui
        from magicgui.widgets import Label
        from napari.utils.notifications import show_info
        from qtpy.QtCore import QTimer
    except ImportError as exc:
        msg = "Napari GUI requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    data = load_mesospim_match_points_data(cfg)
    state = {
        "overview_z": data.overview_shape_zyx[0] // 2,
        "roi_z": data.roi_shape_zyx[0] // 2,
        "link_z": False,
        "_nav_syncing": False,
        "_view_shape": None,
        "_pending_overview_z": None,
        "_pending_roi_z": None,
    }

    z_nav_timer = QTimer()
    z_nav_timer.setSingleShot(True)
    z_nav_timer.setInterval(32)

    viewer = napari.Viewer(title=f"LightSuite mesospim — {cfg.sample.name}")
    viewer.dims.ndisplay = 2

    overview_layer = viewer.add_image(np.zeros((10, 10)), name="overview", colormap="gray")
    roi_layer = viewer.add_image(np.zeros((10, 10)), name="roi", colormap="gray")
    overview_pts = viewer.add_points(
        np.zeros((0, 2)),
        name="overview_points",
        face_color="yellow",
        size=10,
        ndim=2,
    )
    roi_pts = viewer.add_points(
        np.zeros((0, 2)),
        name="roi_points",
        face_color="cyan",
        size=10,
        ndim=2,
    )
    _configure_point_text(overview_pts)
    _configure_point_text(roi_pts)

    points_summary = Label(
        label="Landmark pairs",
        value="Total points: 0 matched pairs (0 overview, 0 ROI)",
    )
    fit_summary = Label(label="Fit preview", value="Add landmark pairs to preview the transform.")

    def _layout_panels() -> None:
        _h, w = overview_layer.data.shape
        overview_layer.translate = (0.0, 0.0)
        overview_pts.translate = (0.0, 0.0)
        offset = (0.0, float(w + PANEL_GAP_X))
        roi_layer.translate = offset
        roi_pts.translate = offset

    def _update_summaries(fit_message: str | None = None) -> None:
        matched, n_overview, n_roi = data.session.point_counts()
        points_summary.value = (
            f"Total points: {matched} matched pairs ({n_overview} overview, {n_roi} ROI)"
        )
        if fit_message is not None:
            fit_summary.value = fit_message

    def _try_fit_preview() -> None:
        matched, _n_o, _n_r = data.session.point_counts()
        if matched < data.min_pairs:
            _update_summaries(
                f"Need at least {data.min_pairs} matched pairs for fit preview "
                f"(have {matched})."
            )
            return
        try:
            import SimpleITK as sitk

            from lightsuite.mesospim.prepare import apply_volume_geometry

            overview_img = sitk.GetImageFromArray(
                np.zeros(data.overview_shape_zyx, dtype=np.float32)
            )
            roi_img = sitk.GetImageFromArray(np.zeros(data.roi_shape_zyx, dtype=np.float32))
            apply_volume_geometry(overview_img, roi_img, cfg.mesospim)
            fit = fit_landmark_transform(
                overview=overview_img,
                roi=roi_img,
                session=data.session,
                fit_mode=data.session.fit_mode,
                min_pairs=data.min_pairs,
            )
            stats = fit.fit_stats
            _update_summaries(
                f"Preview fit ({data.session.fit_mode}): RMS={fit.rms_error_um:.2f} µm, "
                f"median={stats['median']:.2f} µm, p95={stats['p95']:.2f} µm"
            )
        except (ImportError, ValueError) as exc:
            _update_summaries(f"Fit preview unavailable: {exc}")

    def _update_status() -> None:
        matched, n_overview, n_roi = data.session.point_counts()
        _update_summaries()
        viewer.status = (
            f"Overview Z={state['overview_z']}/{data.overview_shape_zyx[0] - 1} | "
            f"ROI Z={state['roi_z']}/{data.roi_shape_zyx[0] - 1} | "
            f"{_pair_status(n_overview, n_roi)} | "
            f"all: {matched} pairs ({n_overview} overview / {n_roi} ROI)"
        )

    def _refresh(*, fit_preview: bool = False) -> None:
        overview_layer.data = data.overview.read_display_slice(state["overview_z"])
        roi_layer.data = data.roi.read_display_slice(state["roi_z"])
        _layout_panels()
        with overview_pts.events.data.blocker(), roi_pts.events.data.blocker():
            _apply_layer_points(
                overview_pts,
                _layer_xy_from_zyx(data.session.overview_points_zyx, state["overview_z"]),
            )
            _apply_layer_points(
                roi_pts,
                _layer_xy_from_zyx(data.session.roi_points_zyx, state["roi_z"]),
            )
        if fit_preview:
            _try_fit_preview()
        else:
            _update_summaries()
        _update_status()
        if state["_view_shape"] != overview_layer.data.shape:
            viewer.reset_view()
            state["_view_shape"] = overview_layer.data.shape

    def _sync_store_from_layer(panel: str) -> None:
        layer = overview_pts if panel == "overview" else roi_pts
        z_index = state["overview_z"] if panel == "overview" else state["roi_z"]
        store = (
            data.session.overview_points_zyx
            if panel == "overview"
            else data.session.roi_points_zyx
        )
        updated = _zyx_from_layer_xy(np.asarray(layer.data, dtype=float), z_index, store)
        if panel == "overview":
            data.session.overview_points_zyx = updated
        else:
            data.session.roi_points_zyx = updated

    def _on_panel_points_changed(panel: str) -> None:
        _sync_store_from_layer(panel)
        matched, n_overview, n_roi = data.session.point_counts()
        fit_preview = n_overview == n_roi and matched >= data.min_pairs
        _refresh(fit_preview=fit_preview)

    @overview_pts.events.data.connect
    def _overview_changed(_event=None) -> None:
        _on_panel_points_changed("overview")

    @roi_pts.events.data.connect
    def _roi_changed(_event=None) -> None:
        _on_panel_points_changed("roi")

    def _sync_navigation_widget() -> None:
        state["_nav_syncing"] = True
        try:
            navigation.overview_z.value = int(state["overview_z"])
            navigation.roi_z.value = int(state["roi_z"])
            navigation.link_z.value = bool(state["link_z"])
        finally:
            state["_nav_syncing"] = False

    def _set_z(
        *,
        overview_z: int | None = None,
        roi_z: int | None = None,
        refocus_canvas: bool = False,
    ) -> None:
        if overview_z is not None:
            state["overview_z"] = int(
                np.clip(overview_z, 0, data.overview_shape_zyx[0] - 1)
            )
        if roi_z is not None:
            state["roi_z"] = int(np.clip(roi_z, 0, data.roi_shape_zyx[0] - 1))
        if state["link_z"] and overview_z is not None:
            state["roi_z"] = state["overview_z"]
        elif state["link_z"] and roi_z is not None:
            state["overview_z"] = state["roi_z"]
        _refresh(fit_preview=False)
        _sync_navigation_widget()
        if refocus_canvas:
            QTimer.singleShot(0, _refocus_canvas)

    def _schedule_z_from_widget(
        *,
        overview_z: int | None = None,
        roi_z: int | None = None,
    ) -> None:
        if overview_z is not None:
            state["_pending_overview_z"] = int(overview_z)
        if roi_z is not None:
            state["_pending_roi_z"] = int(roi_z)
        z_nav_timer.start()

    @z_nav_timer.timeout.connect
    def _apply_pending_z_from_widget() -> None:
        overview_z = state.pop("_pending_overview_z", None)
        roi_z = state.pop("_pending_roi_z", None)
        _set_z(
            overview_z=overview_z,
            roi_z=roi_z,
        )

    def _refocus_canvas() -> None:
        try:
            qt_viewer = viewer.window._qt_viewer
        except AttributeError:
            return
        canvas_native = getattr(getattr(qt_viewer, "canvas", None), "native", None)
        for target in (canvas_native, qt_viewer):
            if target is not None and hasattr(target, "setFocus"):
                target.setFocus()
                return

    @magicgui(
        overview_z={
            "min": 0,
            "max": max(data.overview_shape_zyx[0] - 1, 0),
            "step": 1,
            "label": "Overview Z index",
        },
        roi_z={
            "min": 0,
            "max": max(data.roi_shape_zyx[0] - 1, 0),
            "step": 1,
            "label": "ROI Z index",
        },
        link_z={"label": "Link overview / ROI Z"},
        call_button="Show slices",
    )
    def navigation(overview_z: int = 0, roi_z: int = 0, link_z: bool = False) -> None:
        state["link_z"] = bool(link_z)
        _set_z(overview_z=overview_z, roi_z=roi_z)

    @navigation.overview_z.changed.connect
    def _overview_z_changed() -> None:
        if state["_nav_syncing"]:
            return
        _schedule_z_from_widget(overview_z=navigation.overview_z.value)

    @navigation.roi_z.changed.connect
    def _roi_z_changed() -> None:
        if state["_nav_syncing"]:
            return
        _schedule_z_from_widget(roi_z=navigation.roi_z.value)

    @navigation.link_z.changed.connect
    def _link_z_changed() -> None:
        if state["_nav_syncing"]:
            return
        state["link_z"] = bool(navigation.link_z.value)
        if state["link_z"]:
            _set_z(overview_z=state["overview_z"])

    @magicgui(call_button="Save && Close")
    def save_controls() -> None:
        matched, n_overview, n_roi = data.session.point_counts()
        if matched >= data.min_pairs and n_overview == n_roi:
            try:
                import SimpleITK as sitk

                from lightsuite.mesospim.prepare import apply_volume_geometry

                overview_img = sitk.GetImageFromArray(
                    np.zeros(data.overview_shape_zyx, dtype=np.float32)
                )
                roi_img = sitk.GetImageFromArray(np.zeros(data.roi_shape_zyx, dtype=np.float32))
                apply_volume_geometry(overview_img, roi_img, cfg.mesospim)
                fit = fit_landmark_transform(
                    overview=overview_img,
                    roi=roi_img,
                    session=data.session,
                    fit_mode=data.session.fit_mode,
                    min_pairs=data.min_pairs,
                )
                update_landmark_session_fit(data.session, fit)
            except ValueError as exc:
                show_info(f"Saved landmarks without fit: {exc}")
        elif matched < data.min_pairs:
            show_info(
                f"Saved {matched} pairs (minimum recommended: {data.min_pairs}). "
                "Add more landmarks before check-geometry / register."
            )
        else:
            show_info(
                "Saved landmarks with unmatched overview / ROI counts. "
                "Add the missing counterpart points before register."
            )
        data.session.save(data.session_path)
        show_info(f"Saved {data.session_path}")
        QTimer.singleShot(0, viewer.close)

    @magicgui(call_button="Clear current-slice points")
    def clear_current_slice() -> None:
        oz = int(state["overview_z"])
        rz = int(state["roi_z"])
        data.session.overview_points_zyx = [
            pt
            for pt in data.session.overview_points_zyx
            if int(round(pt[0])) != oz
        ]
        data.session.roi_points_zyx = [
            pt for pt in data.session.roi_points_zyx if int(round(pt[0])) != rz
        ]
        _refresh(fit_preview=False)

    @magicgui(call_button="Delete last point")
    def delete_last_point() -> None:
        n_o = len(data.session.overview_points_zyx)
        n_r = len(data.session.roi_points_zyx)
        if n_o == 0 and n_r == 0:
            return
        if n_o >= n_r and n_o > 0:
            data.session.overview_points_zyx.pop()
        if n_r > 0:
            data.session.roi_points_zyx.pop()
        matched, n_overview, n_roi = data.session.point_counts()
        fit_preview = n_overview == n_roi and matched >= data.min_pairs
        _refresh(fit_preview=fit_preview)

    def _overview_z_step(delta: int) -> None:
        _set_z(overview_z=state["overview_z"] + delta, refocus_canvas=True)

    def _roi_z_step(delta: int) -> None:
        _set_z(roi_z=state["roi_z"] + delta, refocus_canvas=True)

    viewer.bind_key("Left", lambda _v: _overview_z_step(-1), overwrite=True)
    viewer.bind_key("Right", lambda _v: _overview_z_step(1), overwrite=True)
    viewer.bind_key("PageUp", lambda _v: _roi_z_step(1), overwrite=True)
    viewer.bind_key("PageDown", lambda _v: _roi_z_step(-1), overwrite=True)
    viewer.bind_key("Backspace", lambda _v: delete_last_point(), overwrite=True)

    viewer.window.add_dock_widget(points_summary, area="right", name="Point counts")
    viewer.window.add_dock_widget(fit_summary, area="right", name="Fit preview")
    viewer.window.add_dock_widget(navigation, area="right", name="Navigation")
    viewer.window.add_dock_widget(save_controls, area="right", name="Save")
    viewer.window.add_dock_widget(clear_current_slice, area="right", name="Edit")
    viewer.window.add_dock_widget(delete_last_point, area="right", name="Undo")

    matched, _n_o, _n_r = data.session.point_counts()
    _refresh(fit_preview=matched >= data.min_pairs)
    _sync_navigation_widget()

    console.print(
        "[bold]Napari mesoSPIM landmark GUI[/bold] — overview (left), ROI (right). "
        "Place numbered corresponding landmarks on salient structures. "
        "Points are stored as [Z, Y, X] voxel indices. "
        "Shortcuts: [bold]←[/bold]/[bold]→[/bold] overview Z, "
        "[bold]PgUp[/bold]/[bold]PgDn[/bold] ROI Z, [bold]Backspace[/bold] undo last point."
    )
    napari.run()
    return data.session_path
