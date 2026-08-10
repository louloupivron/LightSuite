"""Simple Napari viewer to trial mesoSPIM lateral_flip before multires registration."""

from __future__ import annotations

import numpy as np
from pathlib import Path
from rich.console import Console

from lightsuite.gui.match_points_shared import PANEL_GAP_X
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
) -> GeometryQcSlice:
    overview_spec, roi_spec, manifest_dir = build_reference_specs_with_geometry(cfg, lateral_flip)
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


def run_multires_inspect_geometry(
    cfg: MultiresPipelineConfig,
    *,
    config_path: str | Path | None = None,
    headless: bool = False,
    write_config: bool = False,
) -> GeometryQcSlice:
    """Open Napari to compare overview vs physically resampled ROI overlap."""
    try:
        import SimpleITK as sitk  # noqa: F401
    except ImportError as exc:
        msg = (
            "SimpleITK is required for multires inspect-geometry.\n"
            f"Install with: {_SITK_HINT}"
        )
        raise RuntimeError(msg) from exc

    flip = lateral_flip_from_config(cfg)
    result = evaluate_initial_slice(cfg, lateral_flip=flip)
    if headless:
        if write_config:
            if config_path is None:
                msg = "config_path is required when write_config=True"
                raise ValueError(msg)
            from lightsuite.config.loader import save_mesospim_lateral_flip_to_multires_config

            saved = save_mesospim_lateral_flip_to_multires_config(config_path, result.lateral_flip)
            console.print(f"[green]Updated[/green] {saved}")
        return result

    try:
        import napari
        from magicgui import magicgui
        from magicgui.widgets import Label
        from qtpy.QtCore import QTimer
    except ImportError as exc:
        msg = f"Napari GUI requires: {_SITK_HINT}"
        raise RuntimeError(msg) from exc

    overview_spec, roi_spec, manifest_dir = build_reference_specs_with_geometry(cfg, flip)
    state = {
        "flip_x": flip[0] == -1,
        "flip_y": flip[1] == -1,
        "overview_z": result.overview_z,
        "roi_z": result.roi_z,
        "link_z": True,
        "_refreshing": False,
        "_view_shape": None,
        "overview_spec": overview_spec,
        "roi_spec": roi_spec,
        "manifest_dir": manifest_dir,
    }

    refresh_timer = QTimer()
    refresh_timer.setSingleShot(True)
    refresh_timer.setInterval(120)

    viewer = napari.Viewer(title=f"LightSuite multires geometry — {cfg.sample.name}")
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
        )

    def _refresh() -> GeometryQcSlice:
        state["_refreshing"] = True
        try:
            qc = _evaluate()
            overview_layer.data = qc.overview_display
            roi_layer.data = qc.roi_resampled_display
            _layout_panels()
            if not qc.has_overlap:
                ncc_label.value = "No FOV overlap for this lateral_flip"
            else:
                ncc_label.value = (
                    f"Physical NCC = {qc.physical_ncc:.3f}  "
                    f"(overview Z={qc.overview_z}, ROI Z={qc.roi_z})"
                )
            yaml_label.value = mesospim_geometry_yaml_snippet(qc.lateral_flip)
            state["overview_z"] = qc.overview_z
            state["roi_z"] = qc.roi_z
            navigation.overview_z.value = int(qc.overview_z)
            navigation.roi_z.value = int(qc.roi_z)
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
        state["link_z"] = bool(link_z)
        state["overview_z"] = int(overview_z)
        state["roi_z"] = int(roi_z)
        if state["link_z"]:
            state["roi_z"] = _roi_z_for_overview_z(
                state["overview_spec"],
                state["roi_spec"],
                state["overview_z"],
            )
        _refresh()

    def _set_overview_z(overview_z: int) -> None:
        state["overview_z"] = int(overview_z)
        if state["link_z"]:
            state["roi_z"] = _roi_z_for_overview_z(
                state["overview_spec"],
                state["roi_spec"],
                state["overview_z"],
            )
        _refresh()

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

    viewer.window.add_dock_widget(geometry_controls, area="right", name="Geometry")
    viewer.window.add_dock_widget(ncc_label, area="right", name="Score")
    viewer.window.add_dock_widget(yaml_label, area="right", name="Config snippet")
    viewer.window.add_dock_widget(patch_status, area="right", name="Config file")
    viewer.window.add_dock_widget(apply_to_config, area="right", name="Save")
    viewer.window.add_dock_widget(navigation, area="right", name="Z navigation")

    _refresh()

    console.print(
        "[bold]Multires geometry inspect[/bold] — overview (left) vs ROI physically resampled "
        "(right). Toggle [bold]Flip X/Y[/bold], then [bold]Apply to YAML[/bold] to patch "
        "[bold]multires.mesospim_geometry[/bold]. Use [bold]←[/bold]/[bold]→[/bold] for overview Z."
    )
    napari.run()
    return _evaluate()
