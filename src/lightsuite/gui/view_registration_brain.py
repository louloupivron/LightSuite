"""Napari GUI for post-export brain registration review."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from lightsuite.config.models import BrainPipelineConfig
from lightsuite.gui.brain_view_data import (
    BrainViewPaths,
    BrainViewVolumes,
    ViewSpace,
    brain_points_to_napari_zyx,
    brain_volume_to_napari_zyx,
    contrast_limits,
    discover_brain_view_paths,
    load_brain_view_volumes,
)
from lightsuite.gui.stage_controller import DockStageController, run_attached_stage
from lightsuite.atlas.registry import atlas_display_provider_from_config
from lightsuite.export.brain_export import _load_transform_params

DisplayMode = Literal["2d", "3d"]


@dataclass
class _ViewState:
    space: ViewSpace = "sample"
    display_mode: DisplayMode = "2d"
    channel_layers: dict[int, Any] = field(default_factory=dict)
    channel_volumes: dict[int, np.ndarray] = field(default_factory=dict)
    division_labels: np.ndarray | None = None
    atlas_provider: str = "allen"
    permute_sample_to_atlas: list[int] | None = None
    reference_shape: tuple[int, int, int] = (1, 1, 1)


def _masked_channel(
    channel: np.ndarray,
    labels: np.ndarray | None,
    selected_ids: list[int],
) -> np.ndarray:
    if labels is None or not selected_ids:
        return channel.astype(np.float32, copy=False)
    mask = np.isin(labels, np.asarray(selected_ids, dtype=np.int32))
    return np.where(mask, channel, np.nan).astype(np.float32)


def _build_division_panel(
    viewer: Any,
    state: _ViewState,
    legend_table: pd.DataFrame | None,
    *,
    on_change,
):
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

    class DivisionPanel(QWidget):
        def __init__(self) -> None:
            super().__init__()
            self._ids: list[int] = []
            self._boxes: list[QCheckBox] = []
            outer = QVBoxLayout(self)
            outer.addWidget(QLabel("Divisions (checked = visible in channels)"))
            btn_row = QHBoxLayout()
            for label, slot in (
                ("Select all", self._select_all),
                ("Clear", self._select_none),
                ("Invert", self._invert),
            ):
                btn = QPushButton(label)
                btn.clicked.connect(slot)
                btn_row.addWidget(btn)
            outer.addLayout(btn_row)

            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
            inner = QWidget()
            inner_layout = QVBoxLayout(inner)
            if legend_table is not None and not legend_table.empty:
                for _, row in legend_table.sort_values("division_id").iterrows():
                    did = int(row["division_id"])
                    acr = str(row.get("division_acronym", did))
                    name = str(row.get("division_name", acr))
                    cb = QCheckBox(f"{acr}  (id {did})")
                    cb.setToolTip(name)
                    cb.setChecked(True)
                    cb.stateChanged.connect(lambda _=None: self._apply())
                    inner_layout.addWidget(cb)
                    self._ids.append(did)
                    self._boxes.append(cb)
            else:
                inner_layout.addWidget(QLabel("No division legend available."))
            inner_layout.addStretch()
            scroll.setWidget(inner)
            scroll.setMinimumHeight(280)
            scroll.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            outer.addWidget(scroll)
            line = QFrame()
            line.setFrameShape(QFrame.HLine)
            outer.addWidget(line)
            self._status = QLabel("")
            self._status.setWordWrap(True)
            outer.addWidget(self._status)
            self._apply()

        def _selected_ids(self) -> list[int]:
            return [i for i, cb in zip(self._ids, self._boxes) if cb.isChecked()]

        def _apply(self) -> None:
            on_change(self._selected_ids())
            n = len(self._selected_ids())
            self._status.setText(f"Visible divisions: {n} / {len(self._ids)}")

        def _bulk_set(self, fn) -> None:
            for cb in self._boxes:
                cb.blockSignals(True)
            try:
                for cb in self._boxes:
                    fn(cb)
            finally:
                for cb in self._boxes:
                    cb.blockSignals(False)
            self._apply()

        def _select_all(self) -> None:
            self._bulk_set(lambda cb: cb.setChecked(True))

        def _select_none(self) -> None:
            self._bulk_set(lambda cb: cb.setChecked(False))

        def _invert(self) -> None:
            self._bulk_set(lambda cb: cb.setChecked(not cb.isChecked()))

    return DivisionPanel()


def attach_brain_view_registration(
    viewer: Any,
    config: BrainPipelineConfig,
    *,
    space: ViewSpace = "sample",
    paths: BrainViewPaths | None = None,
    volumes: BrainViewVolumes | None = None,
) -> DockStageController:
    """Attach registration review layers to an existing napari viewer."""
    from napari.utils.notifications import show_info

    paths = paths or discover_brain_view_paths(config, space=space)
    volumes = volumes or load_brain_view_volumes(
        config,
        paths=paths,
        space=space,
        load_multires_roi=False,
    )
    atlas_provider = atlas_display_provider_from_config(config.atlas)
    permute_sample_to_atlas: list[int] | None = None
    if space == "sample":
        transform_params = _load_transform_params(config.sample.save_path.expanduser())
        permute_sample_to_atlas = transform_params.permute_sample_to_atlas or [1, 2, 3]

    reference_shape = (
        next(iter(volumes.registered_channels.values())).shape
        if volumes.registered_channels
        else (1, 1, 1)
    )

    state = _ViewState(
        space=space,
        division_labels=volumes.division_labels,
        atlas_provider=atlas_provider,
        permute_sample_to_atlas=permute_sample_to_atlas,
        reference_shape=reference_shape,
        channel_volumes=dict(volumes.registered_channels),
    )

    viewer.dims.ndisplay = 2

    def _to_napari(volume: np.ndarray) -> np.ndarray:
        return brain_volume_to_napari_zyx(
            volume,
            space=state.space,
            atlas_provider=state.atlas_provider,
            permute_sample_to_atlas=state.permute_sample_to_atlas,
        )

    def _apply_division_mask(selected_ids: list[int]) -> None:
        depiction = "volume" if state.display_mode == "3d" else "plane"
        for ichan, layer in state.channel_layers.items():
            base = state.channel_volumes[ichan]
            masked = _masked_channel(base, state.division_labels, selected_ids)
            layer.data = _to_napari(masked)
            layer.depiction = depiction

    channel_cmaps = ["gray", "magenta", "cyan", "yellow", "red"]
    for idx, (ichan, vol) in enumerate(sorted(volumes.registered_channels.items())):
        cmap = channel_cmaps[idx % len(channel_cmaps)]
        layer = viewer.add_image(
            _to_napari(vol),
            name=f"channel {ichan}",
            colormap=cmap,
            blending="additive" if len(volumes.registered_channels) > 1 else "opaque",
            opacity=0.85 if len(volumes.registered_channels) > 1 else 1.0,
            contrast_limits=contrast_limits(vol),
        )
        state.channel_layers[ichan] = layer
        state.channel_volumes[ichan] = vol

    roi_cmaps = ["orange", "lime", "violet", "pink"]
    background_workers: list[Any] = []
    load_state = {"cancelled": False}
    _teardown_multires_load = None

    for idx, (channel, vol) in enumerate(sorted(volumes.multires_roi_channels.items())):
        cmap = roi_cmaps[idx % len(roi_cmaps)]
        viewer.add_image(
            _to_napari(vol),
            name=f"ROI registered — ch {channel}",
            colormap=cmap,
            blending="additive",
            opacity=0.55,
            contrast_limits=contrast_limits(vol),
        )

    if config.multires_link is not None and not volumes.multires_roi_channels:

        def _teardown_multires_load() -> None:
            load_state["cancelled"] = True
            from lightsuite.gui.qt_workers import wait_background_workers

            wait_background_workers(background_workers)

        _schedule_multires_roi_layers(
            viewer,
            config,
            state=state,
            expected_shape=reference_shape,
            to_napari=_to_napari,
            roi_cmaps=roi_cmaps,
            background_workers=background_workers,
            load_state=load_state,
        )

    if volumes.annotation is not None:
        viewer.add_labels(
            _to_napari(volumes.annotation).astype(np.int64, copy=False),
            name="atlas annotation",
            opacity=0.45,
        )

    if volumes.boundary is not None:
        viewer.add_image(
            _to_napari(volumes.boundary),
            name="atlas boundary",
            colormap="red",
            blending="additive",
            opacity=0.25,
            contrast_limits=(0.0, 1.0),
            visible=False,
        )

    if volumes.division_labels is not None:
        viewer.add_labels(
            _to_napari(volumes.division_labels).astype(np.int64, copy=False),
            name="division labels",
            opacity=0.15,
            visible=False,
        )

    for label, mask in volumes.mask_layers.items():
        viewer.add_image(
            _to_napari(mask),
            name=f"mask: {label}",
            colormap="red",
            blending="additive",
            opacity=0.35,
            contrast_limits=(0.0, 1.0),
        )

    for label, mask in volumes.resampled_roi_masks.items():
        viewer.add_image(
            _to_napari(mask),
            name=f"{label} (resampled)",
            colormap="yellow",
            blending="additive",
            opacity=0.4,
            contrast_limits=(0.0, 1.0),
        )

    for label, coords in volumes.point_layers.items():
        viewer.add_points(
            brain_points_to_napari_zyx(
                coords,
                space=state.space,
                atlas_provider=state.atlas_provider,
                volume_shape_yxz=reference_shape,
                permute_sample_to_atlas=state.permute_sample_to_atlas,
            ),
            name=f"points: {label}",
            size=4,
            face_color="red",
            border_color="white",
        )

    for label, coords in volumes.resampled_roi_points.items():
        viewer.add_points(
            brain_points_to_napari_zyx(
                coords,
                space=state.space,
                atlas_provider=state.atlas_provider,
                volume_shape_yxz=reference_shape,
                permute_sample_to_atlas=state.permute_sample_to_atlas,
            ),
            name=f"{label} (resampled)",
            size=5,
            face_color="yellow",
            border_color="black",
        )

    division_panel = _build_division_panel(
        viewer,
        state,
        volumes.division_legend,
        on_change=_apply_division_mask,
    )

    viewer.window.add_dock_widget(division_panel, area="right", name="Divisions")

    n_roi = len(volumes.resampled_roi_masks) + len(volumes.resampled_roi_points)
    n_multires = len(volumes.multires_roi_channels)
    extras: list[str] = []
    if n_roi:
        extras.append(f"{n_roi} resampled ROI preview layer(s)")
    if n_multires:
        extras.append(f"{n_multires} multires ROI channel(s) on 20 µm grid")
    elif config.multires_link is not None:
        extras.append("multires ROI channel(s) loading in background")

    open_log_message: str | None = None
    if config.multires_link is not None and not volumes.multires_roi_channels:
        open_log_message = "multires ROI channel(s) loading in background"

    show_info(
        f"Loaded {len(volumes.registered_channels)} channel(s), "
        f"{len(volumes.point_layers)} import point layer(s), "
        f"{len(volumes.mask_layers)} import mask layer(s)"
        + (", " + ", ".join(extras) if extras else "")
        + "."
    )

    return DockStageController(
        result=paths,
        _teardown_fn=_teardown_multires_load,
        open_log_message=open_log_message,
    )


def _schedule_multires_roi_layers(
    viewer: Any,
    config: BrainPipelineConfig,
    *,
    state: _ViewState,
    expected_shape: tuple[int, int, int],
    to_napari,
    roi_cmaps: list[str],
    background_workers: list[Any],
    load_state: dict[str, bool],
) -> None:
    """Load multires ROI overlays off the Qt main thread, then add Napari layers."""
    from napari.utils.notifications import show_info, show_warning
    from lightsuite.gui.brain_multires_link import load_multires_roi_channels_on_registration_grid
    from lightsuite.gui.qt_workers import start_background_task
    from lightsuite.preprocess.checkpoint import RegOptsCheckpoint

    save_path = config.sample.save_path.expanduser()
    regopts_path = save_path / "regopts.json"
    if not regopts_path.is_file():
        return

    def _work() -> dict[str, np.ndarray]:
        regopts = RegOptsCheckpoint.load(regopts_path)
        transform_params = _load_transform_params(save_path)
        return load_multires_roi_channels_on_registration_grid(
            config,
            checkpoint=regopts,
            transform_params=transform_params,
            expected_shape=expected_shape,
        )

    def _on_success(channels: dict[str, np.ndarray]) -> None:
        if load_state.get("cancelled"):
            return
        if not channels:
            show_warning("No multires ROI channels could be loaded for view-registration.")
            return
        for idx, (channel, vol) in enumerate(sorted(channels.items())):
            cmap = roi_cmaps[idx % len(roi_cmaps)]
            viewer.add_image(
                to_napari(vol),
                name=f"ROI registered — ch {channel}",
                colormap=cmap,
                blending="additive",
                opacity=0.55,
                contrast_limits=contrast_limits(vol),
            )
        show_info(f"Added {len(channels)} multires ROI channel(s) on the 20 µm grid.")

    def _on_failure(exc: BaseException) -> None:
        if load_state.get("cancelled"):
            return
        show_warning(f"Multires ROI load failed: {exc}")

    background_workers.append(
        start_background_task(_work, on_success=_on_success, on_failure=_on_failure)
    )


def run_brain_view_registration(
    config: BrainPipelineConfig,
    *,
    space: ViewSpace = "sample",
    headless: bool = False,
) -> BrainViewPaths:
    """Open Napari to review registration quality after export."""
    paths = discover_brain_view_paths(config, space=space)
    if headless:
        load_brain_view_volumes(config, paths=paths, space=space)
        return paths

    title = f"LightSuite — {config.sample.name} (registration review)"

    def _attach(viewer: Any) -> DockStageController:
        return attach_brain_view_registration(
            viewer,
            config,
            space=space,
            paths=paths,
        )

    run_attached_stage(title, _attach)
    return paths
