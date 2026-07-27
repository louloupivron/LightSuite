"""Napari GUI for multiresolution overview / ROI landmark matching."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from rich.console import Console

from lightsuite.gui.match_points_shared import (
    PANEL_GAP_X,
    apply_layer_points,
    configure_point_text,
    pair_status,
)
from lightsuite.gui.multires_data import (
    load_multires_match_points_data,
    prepare_multires_match_points_session,
)
from lightsuite.multires.config_models import MultiresPipelineConfig

console = Console()

_SITK_HINT = "uv sync --extra gui --extra registration"


def run_multires_match_points(cfg: MultiresPipelineConfig, *, headless: bool = False) -> Path:
    """Launch Napari landmark matcher; returns saved session path."""
    try:
        import SimpleITK as sitk  # noqa: F401
    except ImportError as exc:
        msg = (
            "SimpleITK is required for multires match-points.\n"
            f"Install with: {_SITK_HINT}"
        )
        raise RuntimeError(msg) from exc

    if headless:
        return prepare_multires_match_points_session(cfg)

    try:
        import napari
        from magicgui import magicgui
        from magicgui.widgets import Label
        from napari.utils.notifications import show_info
        from qtpy.QtCore import QTimer
    except ImportError as exc:
        msg = f"Napari GUI requires: {_SITK_HINT}"
        raise RuntimeError(msg) from exc

    from lightsuite.multires.landmarks import fit_landmark_transform, update_landmark_session_fit
    from lightsuite.multires.spec_geometry import sitk_geometry_from_spec

    data = load_multires_match_points_data(cfg)
    state = {
        "overview_z": data.initial_overview_z,
        "roi_z": data.initial_roi_z,
        "link_z": True if data.crop_mode else False,
        "_nav_syncing": False,
        "_view_shape": None,
        "_pending_overview_z": None,
        "_pending_roi_z": None,
    }

    z_nav_timer = QTimer()
    z_nav_timer.setSingleShot(True)
    z_nav_timer.setInterval(32)

    mode_label = "hybrid crop" if data.crop_mode else "full volume"
    viewer = napari.Viewer(
        title=f"LightSuite multires — {cfg.sample.name} ({mode_label})"
    )
    viewer.dims.ndisplay = 2

    overview_layer = viewer.add_image(np.zeros((10, 10)), name="overview_crop", colormap="gray")
    roi_layer = viewer.add_image(np.zeros((10, 10)), name="roi_crop", colormap="gray")
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
    configure_point_text(overview_pts)
    configure_point_text(roi_pts)

    points_summary = Label(
        label="Landmark pairs",
        value="Total points: 0 matched pairs (0 overview, 0 ROI)",
    )
    fit_summary = Label(label="Fit preview", value="Add landmark pairs to preview the transform.")
    crop_summary = Label(
        label="Display",
        value=(
            f"Metadata overlap crop (±{data.margin_um:.0f} µm)"
            if data.crop_mode
            else "Full volumes (no metadata overlap — fallback)"
        ),
    )

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
            overview_geo = sitk_geometry_from_spec(data.overview_spec)
            roi_geo = sitk_geometry_from_spec(data.roi_spec)
            fit = fit_landmark_transform(
                overview=overview_geo,
                roi=roi_geo,
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
        oy0, ox0 = data.overview.xy_origin_yx
        ry0, rx0 = data.roi.xy_origin_yx
        viewer.status = (
            f"Overview Z={state['overview_z']} "
            f"[{data.overview.z_min}–{data.overview.z_max}] "
            f"cropYX@({oy0},{ox0}) | "
            f"ROI Z={state['roi_z']} "
            f"[{data.roi.z_min}–{data.roi.z_max}] "
            f"cropYX@({ry0},{rx0}) | "
            f"{pair_status(n_overview, n_roi)} | "
            f"all: {matched} pairs"
        )

    def _refresh(*, fit_preview: bool = False) -> None:
        overview_layer.data = data.overview.read_display_slice(state["overview_z"])
        roi_layer.data = data.roi.read_display_slice(state["roi_z"])
        _layout_panels()
        with overview_pts.events.data.blocker(), roi_pts.events.data.blocker():
            apply_layer_points(
                overview_pts,
                data.overview.display_xy_from_volume_zyx(
                    data.session.overview_points_zyx,
                    state["overview_z"],
                ),
            )
            apply_layer_points(
                roi_pts,
                data.roi.display_xy_from_volume_zyx(
                    data.session.roi_points_zyx,
                    state["roi_z"],
                ),
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
        if panel == "overview":
            data.session.overview_points_zyx = data.overview.volume_zyx_from_display_xy(
                np.asarray(overview_pts.data, dtype=float),
                state["overview_z"],
                data.session.overview_points_zyx,
            )
        else:
            data.session.roi_points_zyx = data.roi.volume_zyx_from_display_xy(
                np.asarray(roi_pts.data, dtype=float),
                state["roi_z"],
                data.session.roi_points_zyx,
            )

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
            state["overview_z"] = data.overview.clip_z(overview_z)
        if roi_z is not None:
            state["roi_z"] = data.roi.clip_z(roi_z)
        if state["link_z"] and overview_z is not None and roi_z is None:
            state["roi_z"] = data.roi.z_index_from_physical_z(
                data.overview.physical_z_um(state["overview_z"])
            )
        elif state["link_z"] and roi_z is not None and overview_z is None:
            state["overview_z"] = data.overview.z_index_from_physical_z(
                data.roi.physical_z_um(state["roi_z"])
            )
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

    _init_oz = int(data.initial_overview_z)
    _init_rz = int(data.initial_roi_z)
    _init_link = bool(state["link_z"])

    @magicgui(
        overview_z={
            "min": int(data.overview.z_min),
            "max": int(data.overview.z_max),
            "step": 1,
            "label": "Overview Z index",
        },
        roi_z={
            "min": int(data.roi.z_min),
            "max": int(data.roi.z_max),
            "step": 1,
            "label": "ROI Z index",
        },
        link_z={"label": "Link Z (physical)"},
        call_button="Show slices",
    )
    def navigation(
        overview_z: int = _init_oz,
        roi_z: int = _init_rz,
        link_z: bool = _init_link,
    ) -> None:
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
                overview_geo = sitk_geometry_from_spec(data.overview_spec)
                roi_geo = sitk_geometry_from_spec(data.roi_spec)
                fit = fit_landmark_transform(
                    overview=overview_geo,
                    roi=roi_geo,
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

    viewer.window.add_dock_widget(crop_summary, area="right", name="Crop mode")
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
        "[bold]Napari multires landmark GUI[/bold] — "
        f"{'metadata overlap crop' if data.crop_mode else 'full volumes'} "
        "(overview left, ROI right). "
        "Landmarks are stored as full-volume [Z, Y, X] indices. "
        "After saving, set [bold]geometry_mode: hybrid[/bold] then "
        "run check-geometry / register. "
        "Shortcuts: [bold]←[/bold]/[bold]→[/bold] overview Z, "
        "[bold]PgUp[/bold]/[bold]PgDn[/bold] ROI Z, [bold]Backspace[/bold] undo last point."
    )
    napari.run()
    return data.session_path
