"""Napari GUI for spinal cord straightening (spinal_cord_aligner.m port)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
from rich.console import Console

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.stage_controller import (
    DockStageController,
    require_napari,
    run_attached_stage,
)
from lightsuite.preprocess.cord_checkpoint import CordRegOptsCheckpoint, SpinalAlignmentCheckpoint
from lightsuite.registration.straightening_optimizer import run_straightening_optimizer

console = Console()

POINT_MARKER_SIZE = 2


def fit_overlay_geometry(fit: dict[str, np.ndarray], slice_index: int) -> dict[str, np.ndarray] | None:
    """Compute napari overlay geometry for one slice (spinal_cord_aligner.m update_view)."""
    cx = float(fit["fit_x"][slice_index])
    cy = float(fit["fit_y"][slice_index])
    th = float(fit["fit_theta"][slice_index])
    if not np.isfinite(cx) or not np.isfinite(cy) or not np.isfinite(th):
        return None
    r = float(fit["fit_rad"][slice_index])
    if not np.isfinite(r):
        r = 10.0
    dx = r * np.cos(th)
    dy = r * np.sin(th)
    row_c, col_c = cy, cx
    row_ant, col_ant = cy + dy, cx + dx
    row_pos, col_pos = cy - dy, cx - dx
    return {
        "center": np.array([[row_c, col_c]], dtype=float),
        "anterior": np.array([[row_ant, col_ant]], dtype=float),
        "posterior": np.array([[row_pos, col_pos]], dtype=float),
        "axis": np.array([[row_pos, col_pos], [row_ant, col_ant]], dtype=float),
    }


def clear_all_alignment_points(
    user_cen: np.ndarray,
    user_ant: np.ndarray,
    user_pos: np.ndarray,
) -> None:
    user_cen[:] = np.nan
    user_ant[:] = np.nan
    user_pos[:] = np.nan


def pop_last_alignment_edit(
    history: list[tuple[int, str, tuple[float, float] | None]],
    user_cen: np.ndarray,
    user_ant: np.ndarray,
    user_pos: np.ndarray,
) -> int | None:
    """Undo the most recent point edit; return affected slice index."""
    if not history:
        return None
    slice_idx, kind, before = history.pop()
    arrays = {"cen": user_cen, "ant": user_ant, "pos": user_pos}
    arr = arrays[kind]
    if before is None:
        arr[slice_idx] = [np.nan, np.nan]
    else:
        arr[slice_idx] = list(before)
    return slice_idx


@dataclass
class StraightenCordData:
    regvol: np.ndarray
    display_vol: np.ndarray
    n_slices: int
    user_cen: np.ndarray
    user_ant: np.ndarray
    user_pos: np.ndarray
    lambda_pos: float
    lambda_ang: float
    fit: dict[str, np.ndarray]


def _normalize_display(regvol: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(1)
    idx = rng.choice(regvol.size, size=min(regvol.size, 20_000), replace=False)
    sample = regvol.ravel()[idx].astype(np.float32)
    vmax = float(np.quantile(sample, 0.999))
    vmin = float(np.quantile(sample, 0.001))
    scaled = (regvol.astype(np.float32) - vmin) / max(vmax - vmin, 1e-6)
    return np.clip(scaled, 0, 1)


def load_straighten_data(config: SpinalCordPipelineConfig) -> StraightenCordData:
    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        msg = f"Missing {regopts_path}. Run 'lightsuite spinal preprocess' first."
        raise FileNotFoundError(msg)
    checkpoint = CordRegOptsCheckpoint.load(regopts_path)
    regvol = tifffile.imread(checkpoint.regvol_path).astype(np.uint16)
    n_slices = regvol.shape[2]
    align_path = save_path / "spinal_alignment_opt.json"
    user_cen = np.full((n_slices, 2), np.nan)
    user_ant = np.full((n_slices, 2), np.nan)
    user_pos = np.full((n_slices, 2), np.nan)
    lambda_pos = config.registration.straightening_lambda_pos
    lambda_ang = config.registration.straightening_lambda_ang
    if align_path.is_file():
        prev = SpinalAlignmentCheckpoint.load(align_path)
        user_cen = _restore_points(prev.user_cen)
        user_ant = _restore_points(prev.user_ant)
        user_pos = _restore_points(prev.user_pos)
        lambda_pos = prev.lambda_pos
        lambda_ang = prev.lambda_ang
    fit = run_straightening_optimizer(
        user_cen,
        user_ant,
        user_pos,
        lambda_pos=lambda_pos,
        lambda_ang=lambda_ang,
    )
    return StraightenCordData(
        regvol=regvol,
        display_vol=_normalize_display(regvol),
        n_slices=n_slices,
        user_cen=user_cen,
        user_ant=user_ant,
        user_pos=user_pos,
        lambda_pos=lambda_pos,
        lambda_ang=lambda_ang,
        fit=fit,
    )


def _restore_points(raw: list[list[float | None]]) -> np.ndarray:
    arr = np.full((len(raw), 2), np.nan)
    for i, pair in enumerate(raw):
        if pair[0] is not None and pair[1] is not None:
            arr[i, 0] = pair[0]
            arr[i, 1] = pair[1]
    return arr


def _serialize_points(points: np.ndarray) -> list[list[float | None]]:
    out: list[list[float | None]] = []
    for row in points:
        if np.isnan(row[0]):
            out.append([None, None])
        else:
            out.append([float(row[0]), float(row[1])])
    return out


def save_alignment_checkpoint(
    data: StraightenCordData,
    save_path: Path,
) -> Path:
    fit = run_straightening_optimizer(
        data.user_cen,
        data.user_ant,
        data.user_pos,
        lambda_pos=data.lambda_pos,
        lambda_ang=data.lambda_ang,
    )
    checkpoint = SpinalAlignmentCheckpoint(
        user_cen=_serialize_points(data.user_cen),
        user_ant=_serialize_points(data.user_ant),
        user_pos=_serialize_points(data.user_pos),
        fit_x=fit["fit_x"].tolist(),
        fit_y=fit["fit_y"].tolist(),
        fit_theta=fit["fit_theta"].tolist(),
        lambda_pos=data.lambda_pos,
        lambda_ang=data.lambda_ang,
    )
    out = save_path / "spinal_alignment_opt.json"
    checkpoint.save(out)
    return out


def _headless_default_alignment(data: StraightenCordData) -> StraightenCordData:
    """Place synthetic center/ant/pos clicks for tests without Napari."""
    cy, cx = data.regvol.shape[0] // 2, data.regvol.shape[1] // 2
    for z in range(data.n_slices):
        data.user_cen[z] = [cx + 1.0, cy + 1.0]
        data.user_ant[z] = [cx + 1.0, cy - 8.0]
        data.user_pos[z] = [cx + 1.0, cy + 8.0]
    data.fit = run_straightening_optimizer(
        data.user_cen,
        data.user_ant,
        data.user_pos,
        lambda_pos=data.lambda_pos,
        lambda_ang=data.lambda_ang,
    )
    return data


def attach_spinal_straighten(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    data: StraightenCordData | None = None,
) -> DockStageController:
    """Attach spinal cord straightening controls to an existing napari viewer."""
    from napari.utils.notifications import show_info
    from qtpy.QtCore import Qt
    from qtpy.QtGui import QFont, QKeySequence, QShortcut
    from qtpy.QtWidgets import (
        QHBoxLayout,
        QLabel,
        QPushButton,
        QSlider,
        QSpinBox,
        QVBoxLayout,
        QWidget,
    )

    require_napari()
    save_path = config.sample.save_path.expanduser()
    if data is None:
        data = load_straighten_data(config)

    viewer.dims.ndisplay = 2
    state = {"slice": 0, "_syncing_layers": False, "show_pred": True, "edit_history": []}

    def _stored_point(arr: np.ndarray, slice_idx: int) -> tuple[float, float] | None:
        if np.isnan(arr[slice_idx, 0]):
            return None
        return float(arr[slice_idx, 0]), float(arr[slice_idx, 1])

    def _register_shortcut(key: str, callback) -> None:
        shortcut = QShortcut(QKeySequence(key), viewer.window._qt_viewer)
        shortcut.setContext(Qt.ApplicationShortcut)
        shortcut.activated.connect(callback)

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

    image_layer = viewer.add_image(
        data.display_vol,
        name="cord",
        colormap="gray",
        contrast_limits=(0, 1),
    )

    fit_axis_layer = viewer.add_shapes(
        [],
        shape_type="line",
        edge_color="yellow",
        edge_width=1,
        name="fit axis",
    )
    fit_cen_layer = viewer.add_points(
        np.zeros((0, 2)),
        name="fit center",
        face_color="cyan",
        symbol="+",
        size=POINT_MARKER_SIZE,
        ndim=2,
    )
    fit_ant_layer = viewer.add_points(
        np.zeros((0, 2)),
        name="fit anterior",
        face_color="green",
        symbol="square",
        size=POINT_MARKER_SIZE,
        ndim=2,
    )
    fit_pos_layer = viewer.add_points(
        np.zeros((0, 2)),
        name="fit posterior",
        face_color="red",
        symbol="square",
        size=POINT_MARKER_SIZE,
        ndim=2,
    )
    fit_layers = (fit_axis_layer, fit_cen_layer, fit_ant_layer, fit_pos_layer)
    for layer in fit_layers:
        layer.editable = False

    cen_layer = viewer.add_points(
        np.zeros((0, 2)),
        name="center",
        face_color="blue",
        size=POINT_MARKER_SIZE,
        ndim=2,
    )
    ant_layer = viewer.add_points(
        np.zeros((0, 2)),
        name="anterior",
        face_color="green",
        size=POINT_MARKER_SIZE,
        ndim=2,
    )
    pos_layer = viewer.add_points(
        np.zeros((0, 2)),
        name="posterior",
        face_color="red",
        size=POINT_MARKER_SIZE,
        ndim=2,
    )
    point_layers = (cen_layer, ant_layer, pos_layer)
    layer_kind = {
        cen_layer: "cen",
        ant_layer: "ant",
        pos_layer: "pos",
    }

    def _point_to_layer(row: float, col: float) -> np.ndarray:
        """Map stored 1-based MATLAB-style (x=col, y=row) to napari (row, col)."""
        return np.array([[row, col]], dtype=float)

    def current_slice_points() -> None:
        state["_syncing_layers"] = True
        try:
            s = state["slice"]
            for layer, arr in (
                (cen_layer, data.user_cen),
                (ant_layer, data.user_ant),
                (pos_layer, data.user_pos),
            ):
                pt = arr[s]
                with layer.events.data.blocker():
                    if np.isnan(pt[0]):
                        layer.data = np.zeros((0, 2))
                    else:
                        layer.data = _point_to_layer(pt[1], pt[0])
        finally:
            state["_syncing_layers"] = False

    def _active_point_layer_name() -> str:
        active = viewer.layers.selection.active
        if active in point_layers:
            return str(active.name)
        return "—"

    def _store_from_layer(layer, arr: np.ndarray) -> None:
        if state["_syncing_layers"]:
            return
        s = state["slice"]
        before = _stored_point(arr, s)
        pts = np.asarray(layer.data, dtype=float)
        if pts.size == 0:
            arr[s] = [np.nan, np.nan]
        else:
            if pts.ndim == 1:
                row, col = float(pts[0]), float(pts[1])
            else:
                row, col = float(pts[-1, 0]), float(pts[-1, 1])
            arr[s] = [col + 1.0, row + 1.0]
            if pts.shape[0] > 1:
                with layer.events.data.blocker():
                    layer.data = np.array([[row, col]], dtype=float)
        after = _stored_point(arr, s)
        if before != after:
            state["edit_history"].append((s, layer_kind[layer], before))
        refresh_fit()
        update_status()

    def refresh_fit() -> None:
        data.fit = run_straightening_optimizer(
            data.user_cen,
            data.user_ant,
            data.user_pos,
            lambda_pos=data.lambda_pos,
            lambda_ang=data.lambda_ang,
        )
        update_fit_overlay()

    def update_fit_overlay() -> None:
        geom = fit_overlay_geometry(data.fit, state["slice"])
        show = state["show_pred"] and geom is not None
        for layer in fit_layers:
            layer.visible = show
        if not show or geom is None:
            fit_cen_layer.data = np.zeros((0, 2))
            fit_ant_layer.data = np.zeros((0, 2))
            fit_pos_layer.data = np.zeros((0, 2))
            fit_axis_layer.data = []
            return
        fit_cen_layer.data = geom["center"]
        fit_ant_layer.data = geom["anterior"]
        fit_pos_layer.data = geom["posterior"]
        fit_axis_layer.data = [geom["axis"]]

    status_label = QLabel()

    def _click_counts() -> tuple[int, int, int]:
        n_cen = int(np.sum(~np.isnan(data.user_cen[:, 0])))
        n_ant = int(np.sum(~np.isnan(data.user_ant[:, 0])))
        n_pos = int(np.sum(~np.isnan(data.user_pos[:, 0])))
        return n_cen, n_ant, n_pos

    def update_status() -> None:
        n_cen, n_ant, n_pos = _click_counts()
        hidden = "" if state["show_pred"] else " [fit hidden]"
        viewer.status = (
            f"Slice {state['slice'] + 1}/{data.n_slices} | "
            f"layer={_active_point_layer_name()} | clicks C={n_cen} A={n_ant} P={n_pos}{hidden} | "
            "select a point layer and click; Save/Toggle fit/Clear buttons or S/P keys"
        )
        status_label.setText(viewer.status)

    def clear_all_points() -> None:
        clear_all_alignment_points(data.user_cen, data.user_ant, data.user_pos)
        state["edit_history"].clear()
        refresh_fit()
        current_slice_points()
        update_fit_overlay()
        update_status()

    def clear_last_point() -> None:
        slice_idx = pop_last_alignment_edit(
            state["edit_history"],
            data.user_cen,
            data.user_ant,
            data.user_pos,
        )
        if slice_idx is None:
            show_info("No points to remove.")
            return
        refresh_fit()
        show_slice(slice_idx)
        _refocus_canvas()

    def toggle_fit_preview() -> None:
        state["show_pred"] = not state["show_pred"]
        update_fit_overlay()
        update_status()

    def show_slice(index: int) -> None:
        state["slice"] = int(np.clip(index, 0, data.n_slices - 1))
        image_layer.data = data.display_vol[:, :, state["slice"]]
        current_slice_points()
        update_fit_overlay()
        sync_nav_widgets()
        update_status()

    slice_spin = QSpinBox()
    slice_spin.setRange(1, data.n_slices)
    slice_spin.setAccelerated(True)
    slice_spin.setToolTip("Jump to slice (1-based)")

    slice_slider = QSlider(Qt.Horizontal)
    slice_slider.setRange(1, data.n_slices)
    slice_slider.setPageStep(max(1, data.n_slices // 20))
    slice_slider.setToolTip("Drag to scroll through slices")

    def sync_nav_widgets() -> None:
        one_based = state["slice"] + 1
        for widget in (slice_spin, slice_slider):
            widget.blockSignals(True)
        try:
            slice_spin.setValue(one_based)
            slice_slider.setValue(one_based)
        finally:
            for widget in (slice_spin, slice_slider):
                widget.blockSignals(False)

    def on_nav_widget_changed(value: int) -> None:
        show_slice(int(value) - 1)

    slice_spin.valueChanged.connect(on_nav_widget_changed)
    slice_slider.valueChanged.connect(on_nav_widget_changed)

    def _compact_font() -> QFont:
        font = QFont()
        font.setPointSize(9)
        return font

    controls = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(6, 6, 6, 6)
    layout.setSpacing(4)
    action_row = QHBoxLayout()
    save_btn = QPushButton("Save (s)")
    fit_btn = QPushButton("Toggle fit (p)")
    clear_last_btn = QPushButton("Clear last")
    clear_all_btn = QPushButton("Clear all")
    slice_label = QLabel("Slice")

    compact = _compact_font()
    for widget in (
        save_btn,
        fit_btn,
        clear_last_btn,
        clear_all_btn,
        slice_label,
        slice_spin,
        status_label,
    ):
        widget.setFont(compact)

    def do_save() -> None:
        path = save_alignment_checkpoint(data, save_path)
        console.print(f"[green]Saved[/green] {path}")
        show_info(f"Saved {path}")

    save_btn.clicked.connect(do_save)
    fit_btn.clicked.connect(toggle_fit_preview)
    clear_last_btn.clicked.connect(clear_last_point)
    clear_all_btn.clicked.connect(clear_all_points)
    action_row.addWidget(save_btn)
    action_row.addWidget(fit_btn)
    action_row.addWidget(clear_last_btn)
    action_row.addWidget(clear_all_btn)

    slice_row = QHBoxLayout()
    slice_row.addWidget(slice_label)
    slice_row.addWidget(slice_spin)
    slice_row.addWidget(slice_slider, stretch=1)

    layout.addLayout(action_row)
    layout.addLayout(slice_row)
    status_label.setWordWrap(True)
    layout.addWidget(status_label)
    controls.setLayout(layout)

    cen_layer.events.data.connect(lambda _event=None: _store_from_layer(cen_layer, data.user_cen))
    ant_layer.events.data.connect(lambda _event=None: _store_from_layer(ant_layer, data.user_ant))
    pos_layer.events.data.connect(lambda _event=None: _store_from_layer(pos_layer, data.user_pos))

    @viewer.layers.selection.events.active.connect
    def _on_active_layer_changed(_event=None) -> None:
        update_status()

    def _previous_slice_key(_viewer) -> None:
        show_slice(state["slice"] - 1)

    def _next_slice_key(_viewer) -> None:
        show_slice(state["slice"] + 1)

    _register_shortcut("S", do_save)
    _register_shortcut("P", toggle_fit_preview)
    viewer.bind_key("s", lambda _v: do_save(), overwrite=True)
    viewer.bind_key("p", lambda _v: toggle_fit_preview(), overwrite=True)

    viewer.bind_key("Left", _previous_slice_key, overwrite=True)
    viewer.bind_key("Right", _next_slice_key, overwrite=True)
    viewer.bind_key("PageUp", _next_slice_key, overwrite=True)
    viewer.bind_key("PageDown", _previous_slice_key, overwrite=True)

    @viewer.window._qt_viewer.canvas.events.mouse_wheel.connect
    def _scroll_slice(event) -> None:
        if getattr(event, "delta", None) is None:
            return
        dy = event.delta[1] if len(event.delta) > 1 else event.delta[0]
        if dy == 0:
            return
        step = 1 if dy < 0 else -1
        show_slice(state["slice"] + step)

    out = save_path / "spinal_alignment_opt.json"

    return DockStageController(
        dock_widgets=[(controls, "Controls")],
        _refresh_fn=lambda: show_slice(0),
        result=out,
    )


def run_spinal_straighten(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
) -> Path:
    """Interactive or headless straightening GUI."""
    save_path = config.sample.save_path.expanduser()
    data = load_straighten_data(config)
    if headless:
        data = _headless_default_alignment(data)
        return save_alignment_checkpoint(data, save_path)

    def _attach(viewer: Any) -> DockStageController:
        return attach_spinal_straighten(viewer, config, data=data)

    run_attached_stage("Spinal cord straightening", _attach)
    out = save_path / "spinal_alignment_opt.json"
    if not out.is_file():
        msg = f"Straightening closed without saving; expected {out}"
        console.print(f"[yellow]{msg}[/yellow]")
    return out
