"""Shared Napari helpers for overview / ROI landmark match-points GUIs."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

PANEL_GAP_X = 24
TEXT_LABEL_COLOR = "white"
TEXT_LABEL_OFFSET = (0.0, 8.0)


def pair_labels(n_points: int) -> list[str]:
    return [str(i + 1) for i in range(n_points)]


def configure_point_text(layer) -> None:
    if not hasattr(layer, "text"):
        return
    try:
        layer.text.color = TEXT_LABEL_COLOR
        layer.text.translation = np.array(TEXT_LABEL_OFFSET, dtype=float)
        layer.text.anchor = "center"
        layer.text.size = 12
    except (TypeError, ValueError, AttributeError):
        pass


def layer_xy_from_zyx(points_zyx: list[list[float]], z_index: int) -> np.ndarray:
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


def zyx_from_layer_xy(
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


def pair_status(n_overview: int, n_roi: int) -> str:
    n_pairs = min(n_overview, n_roi)
    if n_overview == n_roi:
        return f"pairs={n_pairs} (matched)"
    if n_overview > n_roi:
        return f"pairs={n_pairs} | place ROI point #{n_roi + 1}"
    return f"pairs={n_pairs} | place overview point #{n_overview + 1}"


def apply_layer_points(layer, xy: np.ndarray) -> None:
    labels = pair_labels(int(xy.shape[0]))
    with layer.events.data.blocker():
        layer.data = xy
        if hasattr(layer, "text"):
            try:
                layer.text = labels
            except (TypeError, ValueError, AttributeError):
                pass
        configure_point_text(layer)


def spinbox_native(spinbox: object) -> object | None:
    return getattr(spinbox, "native", None)


def focused_spinbox_native(*spinboxes: object) -> object | None:
    """Return the native Qt spinbox that currently has keyboard focus, if any."""
    try:
        from qtpy.QtWidgets import QApplication
    except ImportError:
        return None
    app = QApplication.instance()
    if app is None:
        return None
    focused = app.focusWidget()
    natives = [spinbox_native(widget) for widget in spinboxes]
    return focused if focused in natives else None


def read_spinbox_int(spinbox: object, *, fallback: int) -> int:
    """Read the text the user typed, falling back to the committed spinbox value."""
    native = spinbox_native(spinbox)
    if native is not None and hasattr(native, "lineEdit"):
        line_edit = native.lineEdit()
        if line_edit is not None:
            text = line_edit.text().strip()
            if text not in ("", "-", "+"):
                try:
                    return int(text)
                except ValueError:
                    pass
        if hasattr(native, "value"):
            try:
                return int(native.value())
            except (TypeError, ValueError):
                pass
    try:
        return int(spinbox.value)  # type: ignore[attr-defined]
    except (TypeError, ValueError):
        return fallback


def set_spinbox_int_value(spinbox: object, value: int) -> None:
    """Set a spinbox value without emitting intermediate Qt change signals."""
    native = spinbox_native(spinbox)
    if native is not None and hasattr(native, "blockSignals"):
        native.blockSignals(True)
        try:
            if hasattr(native, "setValue"):
                native.setValue(int(value))
            else:
                spinbox.value = int(value)  # type: ignore[attr-defined]
        finally:
            native.blockSignals(False)
        return
    spinbox.value = int(value)  # type: ignore[attr-defined]


def configure_z_index_spinbox(
    spinbox: object,
    on_commit: Callable[[], None],
    *,
    blocked: Callable[[], bool] | None = None,
) -> None:
    """Apply Z navigation when the committed plane changes.

    ``keyboardTracking=False`` avoids a refresh on every typed digit. Arrow
    buttons, Enter, and focus-out emit ``valueChanged`` and refresh the slices.
    """
    native = spinbox_native(spinbox)
    if native is None:
        return
    if hasattr(native, "setKeyboardTracking"):
        native.setKeyboardTracking(False)

    def _commit(*_args: object) -> None:
        if blocked is not None and blocked():
            return
        on_commit()

    if hasattr(native, "valueChanged"):
        try:
            native.valueChanged[int].connect(_commit)
        except (TypeError, KeyError, AttributeError):
            native.valueChanged.connect(_commit)
    if hasattr(native, "editingFinished"):
        native.editingFinished.connect(_commit)
    line_edit = native.lineEdit() if hasattr(native, "lineEdit") else None
    if line_edit is not None and hasattr(line_edit, "returnPressed"):
        line_edit.returnPressed.connect(_commit)


def sync_z_index_spinboxes(
    overview_spinbox: object,
    roi_spinbox: object,
    *,
    overview_z: int,
    roi_z: int,
    link_z: object | None = None,
    link_value: bool | None = None,
) -> object | None:
    """Update navigation spinboxes without interrupting in-progress text entry."""
    keep_focus = focused_spinbox_native(overview_spinbox, roi_spinbox)
    overview_native = spinbox_native(overview_spinbox)
    roi_native = spinbox_native(roi_spinbox)

    if keep_focus is not overview_native:
        set_spinbox_int_value(overview_spinbox, int(overview_z))
    if keep_focus is not roi_native:
        set_spinbox_int_value(roi_spinbox, int(roi_z))
    if link_z is not None and link_value is not None:
        if bool(link_z.value) != bool(link_value):  # type: ignore[attr-defined]
            link_z.value = bool(link_value)  # type: ignore[attr-defined]

    return keep_focus


# Backwards-compatible alias used by multires GUI modules.
connect_spinbox_apply_on_finish = configure_z_index_spinbox
