"""Napari GUI for post-export brain registration review."""

from __future__ import annotations

from collections.abc import Callable
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
    brain_view_load_summary,
    brain_view_spaces_available,
    brain_volume_to_napari_zyx,
    contrast_limits,
    discover_brain_view_paths,
    load_brain_view_volumes,
    resolve_brain_view_space,
)
from lightsuite.gui.stage_controller import (
    DockStageController,
    clear_viewer_layers_safely,
    remove_dock_widget,
    run_attached_stage,
)
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


def _build_brain_space_switch_panel(
    *,
    available: dict[ViewSpace, bool],
    current: ViewSpace,
    on_change: Callable[[ViewSpace], None],
) -> Any:
    from qtpy.QtWidgets import (
        QButtonGroup,
        QLabel,
        QRadioButton,
        QVBoxLayout,
        QWidget,
    )

    class SpaceSwitchPanel(QWidget):
        def __init__(self) -> None:
            super().__init__()
            layout = QVBoxLayout(self)
            self._title = QLabel("Coordinate space")
            layout.addWidget(self._title)
            self._group = QButtonGroup(self)
            self._buttons: dict[ViewSpace, QRadioButton] = {}
            for space, label in (
                ("sample", "Sample (registration grid)"),
                ("atlas", "Atlas (template grid)"),
            ):
                button = QRadioButton(label)
                button.setEnabled(available.get(space, False))
                if not available.get(space, False):
                    button.setToolTip(
                        "Export not found for this space. "
                        "Enable it under Export spaces and run export."
                    )
                button.toggled.connect(
                    lambda checked, selected=space: checked and on_change(selected)
                )
                self._group.addButton(button)
                layout.addWidget(button)
                self._buttons[space] = button
            self._buttons[current].setChecked(True)

        def set_current(self, space: ViewSpace) -> None:
            button = self._buttons[space]
            button.blockSignals(True)
            button.setChecked(True)
            button.blockSignals(False)

        def set_loading(self, loading: bool) -> None:
            """Disable controls and update title while a space switch is loading."""
            self._title.setText(
                "Coordinate space  [loading…]" if loading else "Coordinate space"
            )
            for button in self._buttons.values():
                if loading:
                    button.setEnabled(False)
                else:
                    space = next(k for k, v in self._buttons.items() if v is button)
                    button.setEnabled(available.get(space, False))

    return SpaceSwitchPanel()


def add_brain_view_layers(
    viewer: Any,
    config: BrainPipelineConfig,
    *,
    paths: BrainViewPaths,
    volumes: BrainViewVolumes,
    space: ViewSpace,
    background_workers: list[Any] | None = None,
    load_state: dict[str, bool] | None = None,
) -> tuple[_ViewState, Callable[[list[int]], None]]:
    """Add Napari layers for one brain view-registration coordinate space."""
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

    def apply_division_mask(selected_ids: list[int]) -> None:
        depiction = "volume" if state.display_mode == "3d" else "plane"
        for ichan, layer in state.channel_layers.items():
            base = state.channel_volumes[ichan]
            masked = _masked_channel(base, state.division_labels, selected_ids)
            layer.data = _to_napari(masked)
            layer.depiction = depiction

    if volumes.template is not None:
        template_name = (
            "atlas template" if space == "atlas" else "template in sample"
        )
        viewer.add_image(
            _to_napari(volumes.template),
            name=template_name,
            colormap="gray",
            blending="additive",
            opacity=0.35 if space == "atlas" else 0.25,
            contrast_limits=contrast_limits(volumes.template),
        )

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

    roi_cmaps = ["orange", "lime", "violet", "pink"]
    workers = background_workers if background_workers is not None else []
    roi_load_state = load_state if load_state is not None else {"cancelled": False}

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

    if (
        space == "sample"
        and config.multires_link is not None
        and not volumes.multires_roi_channels
    ):
        _schedule_multires_roi_layers(
            viewer,
            config,
            state=state,
            expected_shape=reference_shape,
            to_napari=_to_napari,
            roi_cmaps=roi_cmaps,
            background_workers=workers,
            load_state=roi_load_state,
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

    return state, apply_division_mask


def attach_brain_view_registration(
    viewer: Any,
    config: BrainPipelineConfig,
    *,
    space: ViewSpace = "sample",
    paths: BrainViewPaths | None = None,
    volumes: BrainViewVolumes | None = None,
    auto_fallback: bool = True,
) -> DockStageController:
    """Attach registration review layers to an existing napari viewer."""
    from napari.utils.notifications import show_info, show_warning

    available = brain_view_spaces_available(config)
    if not any(available.values()):
        msg = (
            "No view-registration exports found. "
            "Run 'lightsuite brain export' (atlas and/or sample space) first."
        )
        raise FileNotFoundError(msg)

    if auto_fallback:
        resolved_space = resolve_brain_view_space(config, preferred=space)
    elif not available.get(space, False):
        msg = (
            f"{space.capitalize()}-space export is not available. "
            "Run export for that space first."
        )
        raise FileNotFoundError(msg)
    else:
        resolved_space = space

    fallback_note: str | None = None
    if auto_fallback and resolved_space != space and not available.get(space, False):
        fallback_note = (
            f"{space.capitalize()}-space export not found; showing {resolved_space} space. "
            "Use the Coordinate space panel to switch after exporting."
        )

    holder: dict[str, Any] = {
        "space": resolved_space,
        "paths": paths,
        "volumes": volumes,
        "view_state": None,
        "switching": False,
        "division_panel": None,
        "division_handle": None,
    }
    background_workers: list[Any] = []
    load_state = {"cancelled": False}

    def _teardown_multires_load() -> None:
        load_state["cancelled"] = True
        from lightsuite.gui.qt_workers import wait_background_workers

        wait_background_workers(background_workers)

    def _load_space(target_space: ViewSpace) -> tuple[BrainViewPaths, BrainViewVolumes]:
        loaded_paths = discover_brain_view_paths(config, space=target_space)
        loaded_volumes = load_brain_view_volumes(
            config,
            paths=loaded_paths,
            space=target_space,
            load_multires_roi=target_space == "sample",
        )
        return loaded_paths, loaded_volumes

    def _mount_division_panel(
        apply_division_mask: Callable[[list[int]], None],
        loaded_volumes: BrainViewVolumes,
    ) -> None:
        if holder["division_handle"] is not None:
            remove_dock_widget(
                viewer,
                holder["division_panel"],
                dock_handle=holder["division_handle"],
            )
            holder["division_panel"] = None
            holder["division_handle"] = None

        panel = _build_division_panel(
            loaded_volumes.division_legend,
            on_change=apply_division_mask,
        )
        holder["division_panel"] = panel
        holder["division_handle"] = viewer.window.add_dock_widget(
            panel,
            area="right",
            name="Divisions",
        )

    if holder["paths"] is None or holder["volumes"] is None:
        holder["paths"], holder["volumes"] = _load_space(resolved_space)
    elif holder["paths"].space != resolved_space:
        holder["paths"], holder["volumes"] = _load_space(resolved_space)

    holder["view_state"], apply_division_mask = add_brain_view_layers(
        viewer,
        config,
        paths=holder["paths"],
        volumes=holder["volumes"],
        space=resolved_space,
        background_workers=background_workers,
        load_state=load_state,
    )
    _mount_division_panel(apply_division_mask, holder["volumes"])

    def _full_teardown() -> None:
        if holder["division_handle"] is not None:
            remove_dock_widget(
                viewer,
                holder["division_panel"],
                dock_handle=holder["division_handle"],
            )
            holder["division_panel"] = None
            holder["division_handle"] = None
        _teardown_multires_load()

    controller = DockStageController(
        result=holder["paths"],
        _teardown_fn=_full_teardown,
        open_log_message=fallback_note,
    )

    def _switch_space(target_space: ViewSpace) -> None:
        if holder["switching"] or target_space == holder["space"]:
            return
        if not available.get(target_space, False):
            show_warning(
                f"{target_space.capitalize()}-space export is not available. "
                "Run export for that space first."
            )
            space_panel.set_current(holder["space"])
            return

        holder["switching"] = True
        space_panel.set_loading(True)

        # Signal any in-flight background jobs to exit early, then clear the list.
        # cancel_background_workers is non-blocking (short timeout only); the
        # workers check load_state["cancelled"] and bail before doing heavy I/O.
        load_state["cancelled"] = True
        from lightsuite.gui.qt_workers import cancel_background_workers, start_background_task

        cancel_background_workers(background_workers)
        background_workers.clear()
        load_state["cancelled"] = False

        def _work() -> tuple[BrainViewPaths, BrainViewVolumes]:
            return _load_space(target_space)

        def _on_success(result: tuple[BrainViewPaths, BrainViewVolumes]) -> None:
            loaded_paths, loaded_volumes = result
            clear_viewer_layers_safely(viewer)
            vs, adm = add_brain_view_layers(
                viewer,
                config,
                paths=loaded_paths,
                volumes=loaded_volumes,
                space=target_space,
                background_workers=background_workers,
                load_state=load_state,
            )
            holder["view_state"] = vs
            holder["space"] = target_space
            holder["paths"] = loaded_paths
            holder["volumes"] = loaded_volumes
            controller.result = loaded_paths
            _mount_division_panel(adm, loaded_volumes)
            holder["switching"] = False
            space_panel.set_loading(False)
            show_info(
                brain_view_load_summary(
                    loaded_paths,
                    loaded_volumes,
                    space=target_space,
                    config=config,
                )
            )

        def _on_failure(exc: BaseException) -> None:
            show_warning(f"Could not switch to {target_space} space: {exc}")
            space_panel.set_current(holder["space"])
            holder["switching"] = False
            space_panel.set_loading(False)

        background_workers.append(
            start_background_task(_work, on_success=_on_success, on_failure=_on_failure)
        )

    space_panel = _build_brain_space_switch_panel(
        available=available,
        current=resolved_space,
        on_change=_switch_space,
    )
    controller.dock_widgets = [(space_panel, "Coordinate space")]

    def _notify() -> None:
        show_info(
            brain_view_load_summary(
                holder["paths"],
                holder["volumes"],
                space=holder["space"],
                config=config,
            )
        )

    controller._refresh_fn = _notify
    if fallback_note is not None:
        show_info(fallback_note)

    return controller


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
