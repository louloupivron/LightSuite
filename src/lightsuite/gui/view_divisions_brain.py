"""Napari viewer: toggle atlas divisions on registered channel volumes."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from lightsuite.analysis.division_map import ensure_division_map
from lightsuite.atlas.registry import resolve_brain_atlas_from_config, resolve_brain_atlas_with_config
from lightsuite.config.models import BrainPipelineConfig
from lightsuite.export.brain_export import _load_transform_params
from lightsuite.gui.inspect_brain_imports import (
    _contrast_limits,
    discover_brain_import_inspect_paths,
    volume_yxz_to_napari_zyx,
)
from lightsuite.registration.volume import load_registration_volume


@dataclass(frozen=True)
class DivisionViewerPaths:
    volume_registered_dir: Path
    division_labels: Path
    division_legend: Path
    channel_paths: dict[str, Path] = field(default_factory=dict)


def discover_division_viewer_paths(config: BrainPipelineConfig) -> DivisionViewerPaths:
    """Resolve registered channels and division label assets for the viewer."""
    inspect_paths = discover_brain_import_inspect_paths(config)
    transform_params = _load_transform_params(config.sample.save_path.expanduser())
    atlas = resolve_brain_atlas_with_config(transform_params.brain_atlas, config.atlas)
    division = ensure_division_map(atlas)

    channel_paths: dict[str, Path] = {}
    for ichan, path in sorted(inspect_paths.registered_channels.items()):
        channel_paths[f"channel {ichan}"] = path

    if not channel_paths:
        msg = (
            "No chan_*_registered_atlas.tif found. "
            "Run 'lightsuite brain export --save-volume' first."
        )
        raise FileNotFoundError(msg)

    return DivisionViewerPaths(
        volume_registered_dir=inspect_paths.volume_registered_dir,
        division_labels=division.paths.labels_tiff,
        division_legend=division.paths.legend_csv,
        channel_paths=channel_paths,
    )


def _load_channel_volume(path: Path) -> np.ndarray:
    return load_registration_volume(path).astype(np.float32, copy=False)


def load_division_viewer_volumes(
    config: BrainPipelineConfig,
    *,
    paths: DivisionViewerPaths | None = None,
    stride: int = 1,
) -> tuple[dict[str, np.ndarray], np.ndarray, dict[str, tuple[float, float]]]:
    """Load division labels and channel volumes, optionally stride-downsampled."""
    paths = paths or discover_division_viewer_paths(config)
    transform_params = _load_transform_params(config.sample.save_path.expanduser())
    expected = tuple(int(v) for v in transform_params.atlassize)

    step = max(1, stride)
    s = slice(None, None, step)
    labels = np.asarray(tifffile.imread(paths.division_labels), dtype=np.int32)
    if step > 1:
        labels = labels[s, s, s]
    elif tuple(labels.shape) != expected:
        msg = f"Division labels shape {labels.shape} != atlas shape {expected}"
        raise ValueError(msg)

    channels: dict[str, np.ndarray] = {}
    contrast: dict[str, tuple[float, float]] = {}
    for name, path in paths.channel_paths.items():
        vol = _load_channel_volume(path)
        if step > 1:
            vol = vol[s, s, s]
        if tuple(vol.shape) != labels.shape:
            msg = f"{name!r} shape {vol.shape} != division labels {labels.shape}"
            raise ValueError(msg)
        channels[name] = vol
        contrast[name] = _contrast_limits(vol)

    return channels, labels, contrast


def _build_division_panel(
    viewer,
    channel_volumes: dict[str, np.ndarray],
    division_labels: np.ndarray,
    legend_table: pd.DataFrame,
    contrast_limits: dict[str, tuple[float, float]],
):
    """Build the Qt division checkbox dock widget."""
    from qtpy.QtCore import Qt
    from qtpy.QtWidgets import (
        QCheckBox,
        QFrame,
        QHBoxLayout,
        QLabel,
        QPushButton,
        QScrollArea,
        QSizePolicy,
        QVBoxLayout,
        QWidget,
    )

    class DivisionCheckboxPanel(QWidget):
        def __init__(self) -> None:
            super().__init__()
            self._viewer = viewer
            self._channels = {k: np.asarray(v, dtype=np.float32) for k, v in channel_volumes.items()}
            self._labels = np.asarray(division_labels, dtype=np.int32)
            self._contrast = contrast_limits
            self._layers: dict[str, object] = {}
            self._ids: list[int] = []
            self._boxes: list[QCheckBox] = []

            outer = QVBoxLayout(self)
            title = QLabel("Divisions (checked = visible on all channel layers)")
            title.setWordWrap(True)
            outer.addWidget(title)

            btn_row = QHBoxLayout()
            self._btn_all = QPushButton("Select all")
            self._btn_none = QPushButton("Clear")
            self._btn_inv = QPushButton("Invert")
            for btn in (self._btn_all, self._btn_none, self._btn_inv):
                btn_row.addWidget(btn)
            outer.addLayout(btn_row)

            self._btn_all.clicked.connect(self._select_all)
            self._btn_none.clicked.connect(self._select_none)
            self._btn_inv.clicked.connect(self._invert)

            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
            inner = QWidget()
            inner_layout = QVBoxLayout(inner)
            inner_layout.setContentsMargins(4, 4, 4, 4)

            rows = legend_table.sort_values("division_id")
            for _, row in rows.iterrows():
                did = int(row["division_id"])
                acr = str(row["division_acronym"])
                name = str(row["division_name"])
                cb = QCheckBox(f"{acr}  (id {did})")
                cb.setToolTip(name)
                cb.setChecked(True)
                cb.stateChanged.connect(lambda _=None: self._apply_mask())
                inner_layout.addWidget(cb)
                self._ids.append(did)
                self._boxes.append(cb)

            inner_layout.addStretch()
            scroll.setWidget(inner)
            scroll.setMinimumHeight(400)
            scroll.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            outer.addWidget(scroll)

            line = QFrame()
            line.setFrameShape(QFrame.HLine)
            outer.addWidget(line)
            self._status = QLabel("")
            self._status.setWordWrap(True)
            outer.addWidget(self._status)

            ref_shape = next(iter(self._channels.values())).shape
            for layer_name, vol in self._channels.items():
                lims = self._contrast[layer_name]
                masked = self._masked_volume(vol)
                self._layers[layer_name] = viewer.add_image(
                    volume_yxz_to_napari_zyx(masked),
                    name=layer_name,
                    colormap="gray",
                    contrast_limits=lims,
                    blending="translucent",
                    depiction="volume",
                )
            n = len(self._selected_ids())
            self._status.setText(
                f"Visible divisions: {n}  |  {len(self._layers)} layer(s)  |  shape {ref_shape}"
            )

        def _selected_ids(self) -> list[int]:
            return [i for i, cb in zip(self._ids, self._boxes) if cb.isChecked()]

        def _masked_volume(self, channel: np.ndarray) -> np.ndarray:
            ids = self._selected_ids()
            if not ids:
                return np.full_like(channel, np.nan, dtype=np.float32)
            mask = np.isin(self._labels, np.asarray(ids, dtype=np.int32))
            return np.where(mask, channel, np.nan).astype(np.float32)

        def _apply_mask(self) -> None:
            for name, layer in self._layers.items():
                layer.data = volume_yxz_to_napari_zyx(self._masked_volume(self._channels[name]))
            n = len(self._selected_ids())
            ref_shape = next(iter(self._channels.values())).shape
            self._status.setText(
                f"Visible divisions: {n}  |  {len(self._layers)} layer(s)  |  shape {ref_shape}"
            )

        def _bulk_set(self, fn) -> None:
            for cb in self._boxes:
                cb.blockSignals(True)
            try:
                for cb in self._boxes:
                    fn(cb)
            finally:
                for cb in self._boxes:
                    cb.blockSignals(False)
            self._apply_mask()

        def _select_all(self) -> None:
            self._bulk_set(lambda cb: cb.setChecked(True))

        def _select_none(self) -> None:
            self._bulk_set(lambda cb: cb.setChecked(False))

        def _invert(self) -> None:
            self._bulk_set(lambda cb: cb.setChecked(not cb.isChecked()))

    return DivisionCheckboxPanel()


def run_brain_division_viewer(
    config: BrainPipelineConfig,
    *,
    headless: bool = False,
    stride: int = 1,
    force_division_rebuild: bool = False,
) -> DivisionViewerPaths:
    """Open Napari with division checkboxes masking all registered channel layers."""
    transform_params = _load_transform_params(config.sample.save_path.expanduser())
    atlas = resolve_brain_atlas_with_config(transform_params.brain_atlas, config.atlas)
    if force_division_rebuild:
        ensure_division_map(atlas, force=True)

    paths = discover_division_viewer_paths(config)
    channels, labels, contrast = load_division_viewer_volumes(config, paths=paths, stride=stride)

    if headless:
        return paths

    try:
        import napari
        from napari.utils.notifications import show_info
    except ImportError as exc:
        msg = "Napari division viewer requires: uv sync --extra gui"
        raise RuntimeError(msg) from exc

    legend = pd.read_csv(paths.division_legend)
    viewer = napari.Viewer(title=f"LightSuite divisions — {config.sample.name}")
    panel = _build_division_panel(viewer, channels, labels, legend, contrast)
    viewer.window.add_dock_widget(panel, name="Divisions", area="right")
    viewer.add_labels(
        volume_yxz_to_napari_zyx(labels),
        name="division labels",
        opacity=0.2,
        visible=False,
    )
    show_info(
        f"Loaded {len(channels)} channel(s). Toggle divisions in the dock panel. "
        f"Labels: {paths.division_labels.name}"
    )
    napari.run()
    return paths
