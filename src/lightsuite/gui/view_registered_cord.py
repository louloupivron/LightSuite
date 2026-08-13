"""Napari viewer for spinal cord registration review and imported annotations."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, Literal

from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.inspect_cord_imports import (
    CordImportInspectPaths,
    CordImportInspectVolumes,
    add_cord_view_layers,
    cord_view_load_summary,
    cord_view_spaces_available,
    discover_cord_import_inspect_paths,
    load_cord_import_inspect_volumes,
    resolve_cord_view_space,
)
from lightsuite.gui.stage_controller import (
    DockStageController,
    clear_viewer_layers_safely,
    run_attached_stage,
)

ViewSpace = Literal["atlas", "sample"]

CordViewPaths = CordImportInspectPaths
CordViewVolumes = CordImportInspectVolumes
discover_cord_view_paths = discover_cord_import_inspect_paths
load_cord_view_volumes = load_cord_import_inspect_volumes


def _build_cord_space_switch_panel(
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
            layout.addWidget(QLabel("Coordinate space"))
            self._group = QButtonGroup(self)
            self._buttons: dict[ViewSpace, QRadioButton] = {}
            for space, label in (
                ("sample", "Sample (straightened)"),
                ("atlas", "Atlas (registered)"),
            ):
                button = QRadioButton(label)
                button.setEnabled(available.get(space, False))
                if not available.get(space, False):
                    button.setToolTip(
                        "Export not found for this space. "
                        "Run 'lightsuite spinal export --space both'."
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

    return SpaceSwitchPanel()


def attach_spinal_view_registration(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    space: ViewSpace = "sample",
    paths: CordViewPaths | None = None,
    volumes: CordViewVolumes | None = None,
    recompute_annotation: bool = False,
    auto_fallback: bool = True,
) -> DockStageController:
    """Attach registration channels, atlas labels, and imported point layers."""
    from napari.utils.notifications import show_info, show_warning

    available = cord_view_spaces_available(config)
    if not any(available.values()):
        msg = (
            "No view-registration exports found. "
            "Run 'lightsuite spinal export --space both' first."
        )
        raise FileNotFoundError(msg)

    if auto_fallback:
        resolved_space = resolve_cord_view_space(config, preferred=space)
    elif not available.get(space, False):
        msg = (
            f"{space.capitalize()}-space export is not available. "
            "Run 'lightsuite spinal export --space both' first."
        )
        raise FileNotFoundError(msg)
    else:
        resolved_space = space

    fallback_note: str | None = None
    if auto_fallback and resolved_space != space and available.get(space, False) is False:
        fallback_note = (
            "Sample-space export not found; showing atlas-space registration. "
            "Use the Coordinate space panel to switch after exporting sample space."
        )

    state: dict[str, Any] = {
        "space": resolved_space,
        "paths": paths,
        "volumes": volumes,
        "switching": False,
    }

    def _load_space(target_space: ViewSpace) -> tuple[CordViewPaths, CordViewVolumes]:
        loaded_paths = discover_cord_view_paths(config, space=target_space)
        loaded_volumes = load_cord_view_volumes(
            config,
            paths=loaded_paths,
            space=target_space,
            recompute_annotation=recompute_annotation,
        )
        return loaded_paths, loaded_volumes

    if state["paths"] is None or state["volumes"] is None:
        state["paths"], state["volumes"] = _load_space(resolved_space)
    elif state["paths"].space != resolved_space:
        state["paths"], state["volumes"] = _load_space(resolved_space)

    add_cord_view_layers(
        viewer,
        config,
        paths=state["paths"],
        volumes=state["volumes"],
        space=resolved_space,
    )

    controller = DockStageController(result=state["paths"])

    def _switch_space(target_space: ViewSpace) -> None:
        if state["switching"] or target_space == state["space"]:
            return
        if not available.get(target_space, False):
            show_warning(
                f"{target_space.capitalize()}-space export is not available. "
                "Run 'lightsuite spinal export --space both' first."
            )
            space_panel.set_current(state["space"])
            return
        state["switching"] = True
        try:
            loaded_paths, loaded_volumes = _load_space(target_space)
            clear_viewer_layers_safely(viewer)
            add_cord_view_layers(
                viewer,
                config,
                paths=loaded_paths,
                volumes=loaded_volumes,
                space=target_space,
            )
            state["space"] = target_space
            state["paths"] = loaded_paths
            state["volumes"] = loaded_volumes
            controller.result = loaded_paths
            show_info(
                cord_view_load_summary(
                    loaded_paths,
                    loaded_volumes,
                    space=target_space,
                    config=config,
                )
            )
        except (FileNotFoundError, OSError, ValueError) as exc:
            show_warning(f"Could not switch to {target_space} space: {exc}")
            space_panel.set_current(state["space"])
        finally:
            state["switching"] = False

    space_panel = _build_cord_space_switch_panel(
        available=available,
        current=resolved_space,
        on_change=_switch_space,
    )
    controller.dock_widgets = [(space_panel, "Coordinate space")]

    def _notify() -> None:
        show_info(
            cord_view_load_summary(
                state["paths"],
                state["volumes"],
                space=state["space"],
                config=config,
            )
        )

    controller._refresh_fn = _notify
    if fallback_note is not None:
        show_info(fallback_note)
        controller.open_log_message = fallback_note

    return controller


def run_spinal_view_registration(
    config: SpinalCordPipelineConfig,
    *,
    space: ViewSpace = "sample",
    headless: bool = False,
    recompute_annotation: bool = False,
    auto_fallback: bool = True,
) -> CordViewPaths:
    """Open Napari to review registration quality and imported annotations."""
    resolved_space = (
        resolve_cord_view_space(config, preferred=space) if auto_fallback else space
    )
    paths = discover_cord_view_paths(config, space=resolved_space)
    if headless:
        load_cord_view_volumes(
            config,
            paths=paths,
            space=resolved_space,
            recompute_annotation=recompute_annotation,
        )
        return paths

    title = f"LightSuite — {config.sample.name} (registration review)"

    def _attach(viewer: Any) -> DockStageController:
        return attach_spinal_view_registration(
            viewer,
            config,
            space=resolved_space,
            paths=paths,
            recompute_annotation=recompute_annotation,
            auto_fallback=False,
        )

    run_attached_stage(title, _attach)
    return paths


def run_spinal_registered_view(
    config: SpinalCordPipelineConfig,
    *,
    space: ViewSpace = "atlas",
    headless: bool = False,
    recompute_annotation: bool = False,
) -> CordViewPaths:
    """Backward-compatible alias for ``lightsuite spinal view``."""
    return run_spinal_view_registration(
        config,
        space=space,
        headless=headless,
        recompute_annotation=recompute_annotation,
    )
