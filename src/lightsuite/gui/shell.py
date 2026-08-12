"""Unified Napari shell for LightSuite pipeline workflows."""

from __future__ import annotations

import re
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lightsuite.cli.stage_registry import (
    StageContext,
    StageKind,
    get_workflow,
    run_stage,
    stage_kind,
)
from lightsuite.cli.stages import StageState, StageStatus
from lightsuite.config.workflow import load_project
from lightsuite.gui.qt_workers import start_background_task
from lightsuite.gui.stage_attach import get_stage_attach
from lightsuite.gui.stage_controller import StageController, require_napari
from lightsuite.reporter import CallbackReporter, capture_pipeline_output

_STATE_ICONS = {
    StageState.DONE: "✓",
    StageState.PENDING: "○",
    StageState.OPTIONAL: "◌",
    StageState.SKIPPED: "—",
}


def _strip_rich_markup(text: str) -> str:
    return re.sub(r"\[/?[^\]]+\]", "", text)


@dataclass
class PipelineProject:
    """Loaded config and workflow metadata for the GUI shell."""

    workflow: str
    config_path: Path
    config: Any


class LightsuiteShell:
    """Single-window Napari host for pipeline stages."""

    def __init__(self, *, config_path: Path | None = None) -> None:
        from qtpy.QtWidgets import (
            QComboBox,
            QFileDialog,
            QHBoxLayout,
            QLabel,
            QListWidget,
            QListWidgetItem,
            QPushButton,
            QTextEdit,
            QVBoxLayout,
            QWidget,
        )

        from qtpy.QtCore import Qt, QObject, Signal

        self._Qt = Qt

        class _LogEmitter(QObject):
            message = Signal(str)

        self._log_emitter = _LogEmitter()
        self._log_emitter.message.connect(self.log)

        napari = require_napari()
        self.viewer = napari.Viewer(title="LightSuite")
        self._log_emitter.setParent(self.viewer.window._qt_window)
        self._QWidget = QWidget
        self._QListWidget = QListWidget
        self._QListWidgetItem = QListWidgetItem
        self._QVBoxLayout = QVBoxLayout
        self._QHBoxLayout = QHBoxLayout
        self._QPushButton = QPushButton
        self._QLabel = QLabel
        self._QTextEdit = QTextEdit
        self._QComboBox = QComboBox
        self._QFileDialog = QFileDialog
        self.project: PipelineProject | None = None
        self._statuses: list[StageStatus] = []
        self._active_controller: StageController | None = None
        self._active_stage_id: str | None = None
        self._worker = None

        self._panel = self._build_panel()
        self.viewer.window.add_dock_widget(self._panel, area="right", name="LightSuite")

        if config_path is not None:
            self.load_project(config_path)

    def _build_panel(self) -> Any:
        QWidget = self._QWidget
        QVBoxLayout = self._QVBoxLayout
        QHBoxLayout = self._QHBoxLayout
        QPushButton = self._QPushButton
        QLabel = self._QLabel
        QTextEdit = self._QTextEdit
        QComboBox = self._QComboBox
        QListWidget = self._QListWidget

        panel = QWidget()
        layout = QVBoxLayout(panel)

        self._title_label = QLabel("No project loaded")
        self._title_label.setWordWrap(True)
        layout.addWidget(self._title_label)

        open_row = QHBoxLayout()
        self._open_button = QPushButton("Open config…")
        self._open_button.clicked.connect(self._on_open_config)
        self._refresh_button = QPushButton("Refresh")
        self._refresh_button.clicked.connect(self.refresh_statuses)
        open_row.addWidget(self._open_button)
        open_row.addWidget(self._refresh_button)
        layout.addLayout(open_row)

        self._channel_row = QWidget()
        channel_layout = QHBoxLayout(self._channel_row)
        channel_layout.setContentsMargins(0, 0, 0, 0)
        self._channel_label = QLabel("Reference channel:")
        self._channel_combo = QComboBox()
        self._channel_apply_button = QPushButton("Apply")
        self._channel_apply_button.clicked.connect(self._on_apply_reference_channel)
        channel_layout.addWidget(self._channel_label)
        channel_layout.addWidget(self._channel_combo, stretch=1)
        channel_layout.addWidget(self._channel_apply_button)
        self._channel_row.setVisible(False)
        layout.addWidget(self._channel_row)

        self._stage_list = QListWidget()
        self._stage_list.itemDoubleClicked.connect(self._on_stage_activated)
        layout.addWidget(self._stage_list, stretch=3)

        self._action_button = QPushButton("Run")
        self._action_button.clicked.connect(self._on_action)
        layout.addWidget(self._action_button)
        self._stage_list.currentItemChanged.connect(self._on_stage_selection_changed)

        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.setPlaceholderText("Pipeline log…")
        layout.addWidget(self._log, stretch=2)

        self._set_actions_enabled(False)
        return panel

    def _set_actions_enabled(self, enabled: bool) -> None:
        self._open_button.setEnabled(True)
        self._refresh_button.setEnabled(enabled)
        self._stage_list.setEnabled(enabled)
        self._action_button.setEnabled(enabled)
        if self._channel_row.isVisible():
            self._channel_combo.setEnabled(enabled)
            self._channel_apply_button.setEnabled(enabled)
        self._update_action_button()

    def _on_stage_selection_changed(self, _current: Any = None, _previous: Any = None) -> None:
        self._update_action_button()

    def _update_action_button(self) -> None:
        if self.project is None:
            self._action_button.setText("Run")
            self._action_button.setToolTip("Load a config, then select a stage")
            return
        selected = self._selected_stage()
        if selected is not None:
            stage_id, _status = selected
            spec = next(item.stage for item in self._statuses if item.stage.id == stage_id)
            if stage_kind(spec) == StageKind.INTERACTIVE:
                self._action_button.setText("Open stage")
                self._action_button.setToolTip(f"Open {spec.title} in the viewer")
            else:
                self._action_button.setText("Run stage")
                self._action_button.setToolTip(f"Run {spec.title} in the background")
            return
        pending = self._pending_stages()
        if pending:
            next_item = pending[0]
            if stage_kind(next_item.stage) == StageKind.INTERACTIVE:
                self._action_button.setText("Open next stage")
            else:
                self._action_button.setText("Resume pipeline")
            self._action_button.setToolTip(f"Continue with {next_item.stage.title}")
            return
        self._action_button.setText("Run")
        self._action_button.setToolTip("Select a stage from the list")

    def _pending_stages(self) -> list[StageStatus]:
        return [item for item in self._statuses if item.state == StageState.PENDING]

    def _update_channel_selector(self) -> None:
        if self.project is None or self.project.workflow != "multires":
            self._channel_row.setVisible(False)
            return
        from lightsuite.multires.channels import (
            multires_channel_names,
            resolved_apply_transform_to,
            resolved_reference_channel,
        )

        names = multires_channel_names(self.project.config)
        if len(names) <= 1:
            self._channel_row.setVisible(False)
            return

        ref = resolved_reference_channel(self.project.config) or names[0]
        targets = resolved_apply_transform_to(self.project.config)
        target_note = f" (apply to {', '.join(targets)})" if targets else ""
        self._channel_label.setText(f"Reference channel{target_note}:")
        self._channel_combo.blockSignals(True)
        self._channel_combo.clear()
        self._channel_combo.addItems(names)
        ref_index = names.index(ref) if ref in names else 0
        self._channel_combo.setCurrentIndex(ref_index)
        self._channel_combo.blockSignals(False)
        self._channel_row.setVisible(True)

    def _on_apply_reference_channel(self) -> None:
        if self.project is None:
            return
        from lightsuite.config.loader import save_reference_channel_to_multires_config
        from lightsuite.multires.channels import default_apply_transform_to

        channel = self._channel_combo.currentText()
        try:
            save_reference_channel_to_multires_config(self.project.config_path, channel)
        except (OSError, ValueError) as exc:
            self.log(f"Could not update reference channel: {exc}")
            return
        self._teardown_active_stage()
        workflow, config = load_project(self.project.config_path)
        self.project = PipelineProject(
            workflow=workflow,
            config_path=self.project.config_path,
            config=config,
        )
        self._update_channel_selector()
        targets = default_apply_transform_to(config, channel)
        target_note = f", apply transform to {', '.join(targets)}" if targets else ""
        self.log(f"Co-registration reference channel set to {channel}{target_note}")
        self.refresh_statuses()

    def log(self, text: str) -> None:
        self._log.append(_strip_rich_markup(text))

    def load_project(self, config_path: str | Path) -> None:
        """Load a YAML config and populate the stage checklist."""
        try:
            path = Path(config_path).expanduser().resolve()
            workflow, config = load_project(path)
        except Exception as exc:
            self.log(f"Failed to load config: {exc}")
            return
        self.project = PipelineProject(workflow=workflow, config_path=path, config=config)
        self._title_label.setText(
            f"{workflow.title()} — {config.sample.name}\n{path}"
        )
        self.log(f"Loaded {workflow} config for sample '{config.sample.name}'")
        self._set_actions_enabled(True)
        self._update_channel_selector()
        self.refresh_statuses()

    def refresh_statuses(self) -> None:
        """Re-read checkpoint artifacts and refresh the stage list."""
        if self.project is None:
            return
        workflow, config = load_project(self.project.config_path)
        self.project = PipelineProject(
            workflow=workflow,
            config_path=self.project.config_path,
            config=config,
        )
        spec = get_workflow(self.project.workflow)
        self._statuses = spec.stage_statuses(self.project.config)
        self._stage_list.clear()
        for item in self._statuses:
            stage = item.stage
            icon = _STATE_ICONS.get(item.state, "?")
            label = stage.title
            if stage.manual:
                label += " (GUI)"
            if stage.optional:
                label += " (optional)"
            text = f"{icon}  {label}"
            list_item = self._QListWidgetItem(text)
            list_item.setData(self._Qt.UserRole, stage.id)
            list_item.setToolTip(item.detail or stage.checkpoint_hint)
            self._stage_list.addItem(list_item)

        self._update_action_button()

    def _selected_stage(self) -> tuple[str, StageStatus] | None:
        row = self._stage_list.currentRow()
        if row < 0 or row >= len(self._statuses):
            return None
        status = self._statuses[row]
        return status.stage.id, status

    def _stage_context(self) -> StageContext:
        if self.project is None:
            msg = "No project loaded"
            raise RuntimeError(msg)
        return StageContext(
            config_path=self.project.config_path,
            headless=True,
            force_preprocess=False,
        )

    def _teardown_active_stage(self) -> None:
        if self._active_controller is not None:
            self._active_controller.teardown(self.viewer)
            self._active_controller = None
        self._active_stage_id = None
        self.viewer.layers.clear()

    def finish_interactive_stage(self, *, refresh: bool = True) -> None:
        """Tear down the active GUI stage without closing the shell window."""
        stage_id = self._active_stage_id
        self._teardown_active_stage()
        if refresh:
            self.refresh_statuses()
        if stage_id is not None:
            self.log(f"Saved and closed stage: {stage_id}")

    def _attach_interactive_stage(self, stage_id: str) -> None:
        if self.project is None:
            return
        factory = get_stage_attach(self.project.workflow, stage_id)
        if factory is None:
            self.log(f"Stage {stage_id!r} is not interactive.")
            return
        self._teardown_active_stage()
        ctx = StageContext(
            config_path=self.project.config_path,
            headless=False,
            force_preprocess=False,
        )
        try:
            controller = factory(self.viewer, self.project.config, ctx)
        except (ImportError, RuntimeError) as exc:
            self.log(f"Failed to open {stage_id}: {exc}")
            try:
                from napari.utils.notifications import show_warning

                show_warning(str(exc))
            except ImportError:
                pass
            return
        controller.mount(self.viewer)
        self._active_controller = controller
        self._active_stage_id = stage_id
        self.log(f"Opened interactive stage: {stage_id}")

    def _emit_log(self, text: str) -> None:
        self._log_emitter.message.emit(text)

    def _run_auto_stage(self, stage_id: str) -> None:
        if self.project is None or self._worker is not None:
            return

        ctx = self._stage_context()
        workflow = self.project.workflow
        config = self.project.config
        reporter = CallbackReporter(on_message=self._emit_log)

        def _work() -> Any:
            with capture_pipeline_output(reporter):
                return run_stage(workflow, stage_id, config, ctx)

        self._set_running(True)
        self.log(f"Running {stage_id}…")
        try:
            self._worker = start_background_task(
                _work,
                on_success=lambda result: self._on_stage_finished(stage_id, result),
                on_failure=self._on_stage_error,
            )
        except Exception as exc:
            self._worker = None
            self._set_running(False)
            self.log(f"Could not start {stage_id}: {exc}")
            self.log(traceback.format_exc())

    def _set_running(self, running: bool) -> None:
        enabled = not running and self.project is not None
        self._action_button.setEnabled(enabled)
        self._open_button.setEnabled(not running)
        self._refresh_button.setEnabled(enabled)
        if not running:
            self._update_action_button()

    def _on_stage_finished(self, stage_id: str, _result: Any) -> None:
        self._worker = None
        self._set_running(False)
        self.log(f"Completed {stage_id}")
        self.refresh_statuses()

    def _on_stage_error(self, error_info: Any) -> None:
        self._worker = None
        self._set_running(False)
        if isinstance(error_info, BaseException):
            exc = error_info
        elif isinstance(error_info, tuple) and len(error_info) >= 2:
            exc = error_info[1]
        else:
            exc = error_info
        self.log(f"Stage failed: {exc}")
        self.log(traceback.format_exc())

    def _on_open_config(self) -> None:
        path, _ = self._QFileDialog.getOpenFileName(
            self._panel,
            "Open LightSuite config",
            str(Path.home()),
            "YAML files (*.yaml *.yml)",
        )
        if path:
            self.load_project(path)

    def _on_stage_activated(self, _item: Any) -> None:
        selected = self._selected_stage()
        if selected is None:
            return
        stage_id, status = selected
        spec = next(item.stage for item in self._statuses if item.stage.id == stage_id)
        if stage_kind(spec) == StageKind.INTERACTIVE:
            self._attach_interactive_stage(stage_id)
        elif status.state != StageState.DONE:
            self._run_auto_stage(stage_id)

    def _on_action(self) -> None:
        if self.project is None or self._worker is not None:
            return
        selected = self._selected_stage()
        if selected is not None:
            stage_id, _status = selected
            spec = next(item.stage for item in self._statuses if item.stage.id == stage_id)
            if stage_kind(spec) == StageKind.INTERACTIVE:
                self._attach_interactive_stage(stage_id)
            else:
                self._run_auto_stage(stage_id)
            return
        self._run_next_pending_stage()

    def _run_next_pending_stage(self) -> None:
        if self.project is None or self._worker is not None:
            return
        pending = self._pending_stages()
        if not pending:
            self.log("No pending stages. Select a stage from the list.")
            return
        next_item = pending[0]
        if stage_kind(next_item.stage) == StageKind.INTERACTIVE:
            self.log(f"Next stage is interactive: {next_item.stage.title}.")
            self._attach_interactive_stage(next_item.stage.id)
            return
        self._run_auto_stage(next_item.stage.id)


def launch_lightsuite_gui(config_path: str | Path | None = None) -> None:
    """Open the unified LightSuite Napari shell."""
    napari = require_napari()
    path = Path(config_path).expanduser().resolve() if config_path is not None else None
    shell = LightsuiteShell(config_path=path)
    # Keep a strong reference for the Qt event loop; otherwise the shell is GC'd
    # while the dock widget remains visible and button slots stop firing.
    shell.viewer.window._lightsuite_shell = shell
    napari.run()
