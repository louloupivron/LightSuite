"""Napari GUI for multiresolution overview / ROI landmark matching."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from rich.console import Console

from lightsuite.gui.match_points_shared import (
    PANEL_GAP_X,
    apply_layer_points,
    configure_point_text,
    configure_z_index_spinbox,
    pair_status,
    read_spinbox_int,
    sync_z_index_spinboxes,
)
from lightsuite.gui.stage_controller import (
    DockStageController,
    close_stage_or_viewer,
    require_magicgui,
    run_attached_stage,
)
from lightsuite.multires.config_models import MultiresPipelineConfig

if TYPE_CHECKING:
    from lightsuite.gui.multires_data import MultiresMatchPointsData

console = Console()

_SITK_HINT = "uv sync --extra gui --extra registration"


def _require_sitk() -> None:
    try:
        import SimpleITK as sitk  # noqa: F401
    except ImportError as exc:
        msg = f"SimpleITK is required for multires match-points.\nInstall with: {_SITK_HINT}"
        raise RuntimeError(msg) from exc


def attach_multires_match_points(
    viewer: Any,
    cfg: MultiresPipelineConfig,
    *,
    data: MultiresMatchPointsData | None = None,
) -> DockStageController:
    """Attach multires landmark-matching controls to an existing napari viewer."""
    _require_sitk()
    from lightsuite.gui.multires_data import load_multires_match_points_data

    from magicgui.widgets import Container, Label
    from napari.utils.notifications import show_info
    from qtpy.QtCore import QTimer

    from lightsuite.multires.landmarks import fit_landmark_transform, update_landmark_session_fit
    from lightsuite.multires.spec_geometry import sitk_geometry_from_spec

    magicgui = require_magicgui()

    if data is None:
        data = load_multires_match_points_data(cfg)

    state = {
        "overview_z": data.initial_overview_z,
        "roi_z": data.initial_roi_z,
        "link_z": True if data.crop_mode else False,
        "_nav_syncing": False,
        "_view_shape": None,
    }

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
                f"Need at least {data.min_pairs} matched pairs for fit preview (have {matched})."
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
        keep_focus = sync_z_index_spinboxes(
            navigation.overview_z,
            navigation.roi_z,
            overview_z=int(state["overview_z"]),
            roi_z=int(state["roi_z"]),
            link_z=navigation.link_z,
            link_value=bool(state["link_z"]),
        )
        if keep_focus is not None:
            keep_focus.setFocus()

    def _set_z(
        *,
        overview_z: int | None = None,
        roi_z: int | None = None,
        refocus_canvas: bool = False,
    ) -> None:
        next_oz = state["overview_z"] if overview_z is None else data.overview.clip_z(overview_z)
        next_rz = state["roi_z"] if roi_z is None else data.roi.clip_z(roi_z)
        if state["link_z"] and overview_z is not None and roi_z is None:
            next_rz = data.roi.z_index_from_physical_z(data.overview.physical_z_um(next_oz))
        elif state["link_z"] and roi_z is not None and overview_z is None:
            next_oz = data.overview.z_index_from_physical_z(data.roi.physical_z_um(next_rz))
        if next_oz == state["overview_z"] and next_rz == state["roi_z"]:
            if refocus_canvas:
                QTimer.singleShot(0, _refocus_canvas)
            return
        state["_nav_syncing"] = True
        try:
            state["overview_z"] = next_oz
            state["roi_z"] = next_rz
            _refresh(fit_preview=False)
            _sync_navigation_widget()
        finally:
            state["_nav_syncing"] = False
        if refocus_canvas:
            QTimer.singleShot(0, _refocus_canvas)

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
        call_button=False,
        auto_call=False,
    )
    def navigation(
        overview_z: int = _init_oz,
        roi_z: int = _init_rz,
        link_z: bool = _init_link,
    ) -> None:
        state["link_z"] = bool(link_z)
        oz = read_spinbox_int(navigation.overview_z, fallback=int(state["overview_z"]))
        rz = read_spinbox_int(navigation.roi_z, fallback=int(state["roi_z"]))
        if state["link_z"]:
            _set_z(overview_z=oz)
        else:
            _set_z(overview_z=oz, roi_z=rz)

    configure_z_index_spinbox(
        navigation.overview_z,
        lambda: _set_z(
            overview_z=read_spinbox_int(
                navigation.overview_z,
                fallback=int(state["overview_z"]),
            )
        ),
        blocked=lambda: state["_nav_syncing"],
    )
    configure_z_index_spinbox(
        navigation.roi_z,
        lambda: _set_z(
            roi_z=read_spinbox_int(
                navigation.roi_z,
                fallback=int(state["roi_z"]),
            )
        ),
        blocked=lambda: state["_nav_syncing"],
    )

    @navigation.overview_z.changed.connect
    def _overview_z_widget_changed() -> None:
        if state["_nav_syncing"]:
            return
        _set_z(
            overview_z=read_spinbox_int(
                navigation.overview_z,
                fallback=int(state["overview_z"]),
            )
        )

    @navigation.roi_z.changed.connect
    def _roi_z_widget_changed() -> None:
        if state["_nav_syncing"]:
            return
        _set_z(
            roi_z=read_spinbox_int(
                navigation.roi_z,
                fallback=int(state["roi_z"]),
            )
        )

    @navigation.link_z.changed.connect
    def _link_z_changed() -> None:
        if state["_nav_syncing"]:
            return
        new_link = bool(navigation.link_z.value)
        if new_link == state["link_z"]:
            return
        state["link_z"] = new_link
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
        QTimer.singleShot(0, lambda: close_stage_or_viewer(viewer))

    @magicgui(call_button="Clear current-slice points")
    def clear_current_slice() -> None:
        oz = int(state["overview_z"])
        rz = int(state["roi_z"])
        data.session.overview_points_zyx = [
            pt for pt in data.session.overview_points_zyx if int(round(pt[0])) != oz
        ]
        data.session.roi_points_zyx = [
            pt for pt in data.session.roi_points_zyx if int(round(pt[0])) != rz
        ]
        _refresh(fit_preview=False)

    @magicgui(call_button="Delete all points")
    def delete_all_points() -> None:
        data.session.overview_points_zyx = []
        data.session.roi_points_zyx = []
        data.session.roi_to_overview_tform = None
        data.session.rms_error_um = None
        data.session.fit_point_errors_um = None
        _refresh(fit_preview=False)
        show_info("Cleared all landmark points")

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

    edit_controls = Container(
        widgets=[clear_current_slice, delete_all_points, delete_last_point],
    )

    def _overview_z_step(delta: int) -> None:
        _set_z(overview_z=state["overview_z"] + delta, refocus_canvas=True)

    def _roi_z_step(delta: int) -> None:
        _set_z(roi_z=state["roi_z"] + delta, refocus_canvas=True)

    viewer.bind_key("Left", lambda _v: _overview_z_step(-1), overwrite=True)
    viewer.bind_key("Right", lambda _v: _overview_z_step(1), overwrite=True)
    viewer.bind_key("PageUp", lambda _v: _roi_z_step(1), overwrite=True)
    viewer.bind_key("PageDown", lambda _v: _roi_z_step(-1), overwrite=True)
    viewer.bind_key("Backspace", lambda _v: delete_last_point(), overwrite=True)

    matched, _n_o, _n_r = data.session.point_counts()
    fit_preview = matched >= data.min_pairs

    def _initial_refresh() -> None:
        _refresh(fit_preview=fit_preview)
        _sync_navigation_widget()

    return DockStageController(
        dock_widgets=[
            (points_summary, "Point counts"),
            (fit_summary, "Fit preview"),
            (navigation, "Navigation"),
            (save_controls, "Save"),
            (edit_controls, "Edit"),
        ],
        _refresh_fn=_initial_refresh,
        result=data.session_path,
    )


def run_multires_match_points(cfg: MultiresPipelineConfig, *, headless: bool = False) -> Path:
    """Launch Napari landmark matcher; returns saved session path."""
    _require_sitk()
    from lightsuite.gui.multires_data import (
        load_multires_match_points_data,
        prepare_multires_match_points_session,
    )

    if headless:
        return prepare_multires_match_points_session(cfg)

    data = load_multires_match_points_data(cfg)
    mode_label = "hybrid crop" if data.crop_mode else "full volume"
    title = f"LightSuite multires — {cfg.sample.name} ({mode_label})"

    def _attach(viewer: Any) -> DockStageController:
        return attach_multires_match_points(viewer, cfg, data=data)

    final = run_attached_stage(
        title,
        _attach,
        before_run=lambda: console.print(
            "[bold]Napari multires landmark GUI[/bold] — "
            f"{'metadata overlap crop' if data.crop_mode else 'full volumes'} "
            "(overview left, ROI right). "
            "Landmarks are stored as full-volume [Z, Y, X] indices. "
            "After saving, set [bold]geometry_mode: hybrid[/bold] then "
            "run check-geometry / register. "
            "Shortcuts: [bold]←[/bold]/[bold]→[/bold] overview Z, "
            "[bold]PgUp[/bold]/[bold]PgDn[/bold] ROI Z, [bold]Backspace[/bold] undo last point."
        ),
    )
    return final if final is not None else data.session_path
