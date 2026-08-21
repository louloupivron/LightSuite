"""Simple Napari viewer to trial mesoSPIM lateral_flip before multires registration."""

from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import Any

from rich.console import Console

from lightsuite.gui.match_points_shared import (
    PANEL_GAP_X,
    configure_z_index_spinbox,
    read_spinbox_int,
    sync_z_index_spinboxes,
)
from lightsuite.gui.stage_controller import (
    DockStageController,
    require_magicgui,
    run_attached_stage,
)
from lightsuite.multires.channels import (
    default_apply_transform_to,
    multires_channel_names,
    resolved_reference_channel,
)
from lightsuite.multires.config_models import MultiresPipelineConfig
from lightsuite.multires.geometry_qc import (
    GeometryQcSlice,
    build_reference_specs_with_geometry,
    compute_geometry_qc_slice,
    lateral_flip_from_config,
    lateral_flip_tuple,
    mesospim_geometry_yaml_snippet,
)
from lightsuite.multires.spec_geometry import index_xyz_to_physical, physical_to_continuous_index_xyz

console = Console()

_SITK_HINT = "uv sync --extra gui --extra registration"


def evaluate_initial_slice(
    cfg: MultiresPipelineConfig,
    *,
    lateral_flip: tuple[int, int],
    channel: str | None = None,
) -> GeometryQcSlice:
    overview_spec, roi_spec, manifest_dir = build_reference_specs_with_geometry(
        cfg,
        lateral_flip,
        channel=channel,
    )
    return compute_geometry_qc_slice(
        overview_spec,
        roi_spec,
        manifest_dir=manifest_dir,
        lateral_flip=lateral_flip,
    )


def _roi_z_for_overview_z(
    overview_spec,
    roi_spec,
    overview_z: int,
) -> int:
    point = index_xyz_to_physical(overview_spec, (0.0, 0.0, float(overview_z)))
    roi_z = int(round(physical_to_continuous_index_xyz(roi_spec, tuple(point))[2]))
    nz_roi = int(roi_spec.shape_zyx[0])
    return int(np.clip(roi_z, 0, max(nz_roi - 1, 0)))


def _overview_z_for_roi_z(
    overview_spec,
    roi_spec,
    roi_z: int,
) -> int:
    point = index_xyz_to_physical(roi_spec, (0.0, 0.0, float(roi_z)))
    overview_z = int(round(physical_to_continuous_index_xyz(overview_spec, tuple(point))[2]))
    nz_overview = int(overview_spec.shape_zyx[0])
    return int(np.clip(overview_z, 0, max(nz_overview - 1, 0)))


def attach_multires_inspect_geometry(
    viewer: Any,
    cfg: MultiresPipelineConfig,
    *,
    config_path: str | Path | None = None,
    flip: tuple[int, int] | None = None,
    initial_result: GeometryQcSlice | None = None,
) -> DockStageController:
    """Attach multires geometry QC controls to an existing napari viewer."""
    from magicgui.widgets import Label
    from qtpy.QtCore import QTimer

    magicgui = require_magicgui()

    if flip is None:
        flip = lateral_flip_from_config(cfg)
    channel_names = multires_channel_names(cfg)
    ref_channel = resolved_reference_channel(cfg) or (channel_names[0] if channel_names else None)
    preview_channel = ref_channel
    result = initial_result or evaluate_initial_slice(
        cfg,
        lateral_flip=flip,
        channel=preview_channel,
    )

    overview_spec, roi_spec, manifest_dir = build_reference_specs_with_geometry(
        cfg,
        flip,
        channel=preview_channel,
    )
    state = {
        "flip_x": flip[0] == -1,
        "flip_y": flip[1] == -1,
        "overview_z": result.overview_z,
        "roi_z": result.roi_z,
        "link_z": True,
        "_refreshing": False,
        "_nav_syncing": False,
        "preview_channel": preview_channel,
        "reference_channel": ref_channel,
        "_view_shape": None,
        "overview_spec": overview_spec,
        "roi_spec": roi_spec,
        "manifest_dir": manifest_dir,
    }

    refresh_timer = QTimer()
    refresh_timer.setSingleShot(True)
    refresh_timer.setInterval(120)

    viewer.dims.ndisplay = 2

    overview_layer = viewer.add_image(np.zeros((10, 10)), name="overview", colormap="gray")
    roi_layer = viewer.add_image(np.zeros((10, 10)), name="roi resampled", colormap="gray")
    ncc_label = Label(label="Physical NCC", value="—")
    yaml_label = Label(label="YAML (paste into multires:)", value=mesospim_geometry_yaml_snippet(flip))

    def _layout_panels() -> None:
        _h, w = overview_layer.data.shape
        overview_layer.translate = (0.0, 0.0)
        roi_layer.translate = (0.0, float(w + PANEL_GAP_X))

    def _current_flip() -> tuple[int, int]:
        return lateral_flip_tuple(flip_x=state["flip_x"], flip_y=state["flip_y"])

    def _reload_specs() -> None:
        lateral_flip = _current_flip()
        overview_spec, roi_spec, manifest_dir = build_reference_specs_with_geometry(
            cfg,
            lateral_flip,
            channel=state["preview_channel"],
        )
        state["overview_spec"] = overview_spec
        state["roi_spec"] = roi_spec
        state["manifest_dir"] = manifest_dir
        navigation.overview_z.max = max(int(overview_spec.shape_zyx[0]) - 1, 0)
        navigation.roi_z.max = max(int(roi_spec.shape_zyx[0]) - 1, 0)

    def _evaluate() -> GeometryQcSlice:
        return compute_geometry_qc_slice(
            state["overview_spec"],
            state["roi_spec"],
            manifest_dir=state["manifest_dir"],
            overview_z=state["overview_z"],
            roi_z=state["roi_z"],
            lateral_flip=_current_flip(),
            link_z=state["link_z"],
        )

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
    ) -> None:
        state["_nav_syncing"] = True
        try:
            if overview_z is not None:
                state["overview_z"] = int(overview_z)
            if roi_z is not None:
                state["roi_z"] = int(roi_z)
            if state["link_z"] and overview_z is not None and roi_z is None:
                state["roi_z"] = _roi_z_for_overview_z(
                    state["overview_spec"],
                    state["roi_spec"],
                    state["overview_z"],
                )
            elif state["link_z"] and roi_z is not None and overview_z is None:
                state["overview_z"] = _overview_z_for_roi_z(
                    state["overview_spec"],
                    state["roi_spec"],
                    state["roi_z"],
                )
            _refresh()
        finally:
            state["_nav_syncing"] = False

    def _refresh() -> GeometryQcSlice:
        state["_refreshing"] = True
        try:
            qc = _evaluate()
            overview_layer.data = qc.overview_display
            roi_layer.data = qc.roi_resampled_display
            _layout_panels()
            if not qc.has_overlap:
                ncc_label.value = "No FOV overlap for this lateral_flip"
            elif not state["link_z"]:
                ncc_label.value = (
                    f"Unlinked Z — overview Z={qc.overview_z}, ROI Z={qc.roi_z} "
                    "(physical NCC disabled)"
                )
            else:
                ncc_label.value = (
                    f"Physical NCC = {qc.physical_ncc:.3f}  "
                    f"(overview Z={qc.overview_z}, ROI Z={qc.roi_z})"
                )
            yaml_label.value = mesospim_geometry_yaml_snippet(qc.lateral_flip)
            state["overview_z"] = qc.overview_z
            state["roi_z"] = qc.roi_z
            _sync_navigation_widget()
            if state["_view_shape"] != overview_layer.data.shape:
                viewer.reset_view()
                state["_view_shape"] = overview_layer.data.shape
            return qc
        finally:
            state["_refreshing"] = False

    def _schedule_refresh(*, reload_specs: bool = False) -> None:
        if reload_specs:
            _reload_specs()
        refresh_timer.start()

    @refresh_timer.timeout.connect
    def _on_refresh_timer() -> None:
        _refresh()

    flip_x_init, flip_y_init = (flip[0] == -1, flip[1] == -1)

    @magicgui(
        flip_x={"label": "Flip X (both volumes)"},
        flip_y={"label": "Flip Y (both volumes)"},
    )
    def geometry_controls(flip_x: bool = flip_x_init, flip_y: bool = flip_y_init) -> None:
        state["flip_x"] = bool(flip_x)
        state["flip_y"] = bool(flip_y)
        _schedule_refresh(reload_specs=True)

    @geometry_controls.flip_x.changed.connect
    def _flip_x_changed() -> None:
        if state["_refreshing"]:
            return
        state["flip_x"] = bool(geometry_controls.flip_x.value)
        _schedule_refresh(reload_specs=True)

    @geometry_controls.flip_y.changed.connect
    def _flip_y_changed() -> None:
        if state["_refreshing"]:
            return
        state["flip_y"] = bool(geometry_controls.flip_y.value)
        _schedule_refresh(reload_specs=True)

    channel_panel = None
    if len(channel_names) > 1 and ref_channel is not None:
        from qtpy.QtWidgets import QHBoxLayout, QPushButton, QWidget

        def _use_preview_for_coregistration() -> None:
            from lightsuite.config.loader import save_reference_channel_to_multires_config
            from napari.utils.notifications import show_info, show_warning

            if config_path is None:
                show_warning("No config path — pass -c /path/to/config.yaml on the CLI.")
                return
            channel = str(channel_panel.preview_channel.value)
            try:
                saved = save_reference_channel_to_multires_config(config_path, channel)
            except (OSError, ValueError) as exc:
                show_warning(f"Could not update config: {exc}")
                return
            state["reference_channel"] = channel
            state["preview_channel"] = channel
            channel_panel.preview_channel.value = channel
            cfg.multires.registration.reference_channel = channel
            cfg.multires.registration.apply_transform_to = default_apply_transform_to(cfg, channel)
            patch_status.value = (
                f"Updated {saved.name} → reference_channel {channel!r}"
            )
            show_info(f"Set co-registration reference to {channel}")

        @magicgui(
            preview_channel={"label": "Channel", "choices": channel_names},
            call_button=False,
            auto_call=False,
        )
        def channel_panel(preview_channel: str = ref_channel) -> None:
            state["preview_channel"] = str(preview_channel)

        @channel_panel.preview_channel.changed.connect
        def _preview_channel_changed(value: str) -> None:
            state["preview_channel"] = str(value)

        def _show_preview_channel() -> None:
            state["preview_channel"] = str(channel_panel.preview_channel.value)
            _schedule_refresh(reload_specs=True)

        button_row = QWidget()
        button_layout = QHBoxLayout(button_row)
        button_layout.setContentsMargins(0, 0, 0, 0)
        show_btn = QPushButton("Show channel")
        coreg_btn = QPushButton("Set reference")
        coreg_btn.setToolTip("Save the selected channel as multires.registration.reference_channel")
        show_btn.clicked.connect(_show_preview_channel)
        coreg_btn.clicked.connect(_use_preview_for_coregistration)
        button_layout.addWidget(show_btn)
        button_layout.addWidget(coreg_btn)
        channel_panel.native.layout().addWidget(button_row)
    @magicgui(
        overview_z={"min": 0, "max": max(int(overview_spec.shape_zyx[0]) - 1, 0), "step": 1},
        roi_z={"min": 0, "max": max(int(roi_spec.shape_zyx[0]) - 1, 0), "step": 1},
        link_z={"label": "Link Z (physical)"},
        call_button="Show slices",
    )
    def navigation(
        overview_z: int = result.overview_z,
        roi_z: int = result.roi_z,
        link_z: bool = True,
    ) -> None:
        state["link_z"] = bool(navigation.link_z.value)
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
        blocked=lambda: state["_nav_syncing"] or state["_refreshing"],
    )
    configure_z_index_spinbox(
        navigation.roi_z,
        lambda: _set_z(
            roi_z=read_spinbox_int(
                navigation.roi_z,
                fallback=int(state["roi_z"]),
            )
        ),
        blocked=lambda: state["_nav_syncing"] or state["_refreshing"],
    )

    @navigation.link_z.changed.connect
    def _link_z_changed() -> None:
        if state["_refreshing"] or state["_nav_syncing"]:
            return
        new_link = bool(navigation.link_z.value)
        if new_link == state["link_z"]:
            return
        state["link_z"] = new_link
        if state["link_z"]:
            _set_z(overview_z=state["overview_z"])

    def _set_overview_z(overview_z: int) -> None:
        state["link_z"] = bool(navigation.link_z.value)
        _set_z(overview_z=overview_z)

    viewer.bind_key("Left", lambda _v: _set_overview_z(state["overview_z"] - 1), overwrite=True)
    viewer.bind_key("Right", lambda _v: _set_overview_z(state["overview_z"] + 1), overwrite=True)

    patch_status = Label(
        label="Config file",
        value=str(Path(config_path).resolve()) if config_path is not None else "(no config path)",
    )

    @magicgui(call_button="Apply to YAML")
    def apply_to_config() -> None:
        from lightsuite.config.loader import save_mesospim_lateral_flip_to_multires_config
        from napari.utils.notifications import show_info, show_warning

        if config_path is None:
            show_warning("No config path — pass -c /path/to/config.yaml on the CLI.")
            return
        lateral_flip = _current_flip()
        try:
            saved = save_mesospim_lateral_flip_to_multires_config(config_path, lateral_flip)
        except (OSError, ValueError) as exc:
            show_warning(f"Could not update config: {exc}")
            return
        patch_status.value = f"Updated {saved.name} → lateral_flip {list(lateral_flip)}"
        show_info(f"Patched {saved}")

    dock_widgets: list[tuple[Any, str]] = []
    if channel_panel is not None:
        dock_widgets.append((channel_panel, "Channel"))
    dock_widgets.extend(
        [
            (geometry_controls, "Geometry"),
            (ncc_label, "Score"),
            (yaml_label, "Config snippet"),
            (patch_status, "Config file"),
            (apply_to_config, "Save"),
            (navigation, "Z navigation"),
        ]
    )

    controller = DockStageController(
        dock_widgets=dock_widgets,
        result=result,
    )

    def _refresh_with_result() -> None:
        controller.result = _refresh()

    controller._refresh_fn = _refresh_with_result
    return controller


def run_multires_inspect_geometry(
    cfg: MultiresPipelineConfig,
    *,
    config_path: str | Path | None = None,
    headless: bool = False,
    write_config: bool = False,
) -> GeometryQcSlice:
    """Open Napari to compare overview vs physically resampled ROI overlap."""
    if not cfg.multires.supports_inspect_geometry():
        suite = cfg.multires.vendor.suite.value
        msg = (
            "inspect-geometry is only for mesoSPIM configs (lateral_flip tuning). "
            f"This config uses vendor {suite!r}; run check-geometry instead."
        )
        raise RuntimeError(msg)
    try:
        import SimpleITK as sitk  # noqa: F401
    except ImportError as exc:
        msg = (
            "SimpleITK is required for multires inspect-geometry.\n"
            f"Install with: {_SITK_HINT}"
        )
        raise RuntimeError(msg) from exc

    flip = lateral_flip_from_config(cfg)
    ref_channel = resolved_reference_channel(cfg)
    result = evaluate_initial_slice(cfg, lateral_flip=flip, channel=ref_channel)
    if headless:
        if write_config:
            if config_path is None:
                msg = "config_path is required when write_config=True"
                raise ValueError(msg)
            from lightsuite.config.loader import save_mesospim_lateral_flip_to_multires_config

            saved = save_mesospim_lateral_flip_to_multires_config(config_path, result.lateral_flip)
            console.print(f"[green]Updated[/green] {saved}")
        return result

    title = f"LightSuite multires geometry — {cfg.sample.name}"

    def _attach(viewer: Any) -> DockStageController:
        return attach_multires_inspect_geometry(
            viewer,
            cfg,
            config_path=config_path,
            flip=flip,
            initial_result=result,
        )

    final = run_attached_stage(
        title,
        _attach,
        before_run=lambda: console.print(
            "[bold]Multires geometry inspect[/bold] — overview (left) vs ROI physically resampled "
            "(right). Toggle [bold]Flip X/Y[/bold], then [bold]Apply to YAML[/bold] to patch "
            "[bold]multires.mesospim_geometry[/bold]. Use [bold]←[/bold]/[bold]→[/bold] for overview Z."
        ),
    )
    return final if final is not None else result
