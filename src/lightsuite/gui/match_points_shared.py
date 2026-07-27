"""Shared Napari helpers for overview / ROI landmark match-points GUIs."""

from __future__ import annotations

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
