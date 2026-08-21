"""Unified Napari shell for LightSuite pipeline workflows."""

from __future__ import annotations

import re
import threading
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lightsuite.cli.spaces import (
    default_export_space_checks,
    export_spaces_from_checks,
    format_export_spaces,
)
from lightsuite.cli.stage_registry import (
    StageContext,
    StageKind,
    get_workflow,
    run_stage,
    stage_kind,
)
from lightsuite.cli.stages import StageState, StageStatus
from lightsuite.config.workflow import load_project
from lightsuite.exceptions import StageCancelledError
from lightsuite.gui.config_editor import ConfigEditorDock
from lightsuite.gui.qt_workers import start_background_task
from lightsuite.gui.stage_attach import get_stage_attach
from lightsuite.gui.orientation_cord import cord_orientation_missing
from lightsuite.gui.stage_controller import (
    StageController,
    clear_viewer_layers_safely,
    defer_clear_viewer_layers,
    require_napari,
    validate_gui_dependencies,
)
from lightsuite.reporter import CallbackReporter, capture_pipeline_output, stage_cancellation

_STATE_ICONS = {
    StageState.DONE: "✓",
    StageState.PENDING: "○",
    StageState.OPTIONAL: "◌",
    StageState.SKIPPED: "—",
}

_INTERACTIVE_LOADING_HINTS: dict[str, str] = {
    "check-orientation": "reading volumes",
    "align-slices": "reading volumes and slice correspondence",
    "match-points": "reading volumes and control-point session",
    "straighten": "reading volumes",
    "align-longitudinal": "reading volumes and correspondence",
    "inspect-geometry": "reading multiresolution manifest",
    "inspect-registration": "reading registered ROI volumes",
    "view-registration": "reading registered volumes and annotations",
    "plot-stats": "loading region stats tables",
}


def interactive_loading_message(stage_title: str, stage_id: str) -> str:
    """Log line shown while an interactive stage is opening."""
    hint = _INTERACTIVE_LOADING_HINTS.get(stage_id)
    if hint:
        return f"Loading {stage_title}… ({hint})"
    return f"Loading {stage_title}…"


def _strip_rich_markup(text: str) -> str:
    return re.sub(r"\[/?[^\]]+\]", "", text)


def filter_stage_statuses(
    statuses: list[StageStatus],
    *,
    show_optional: bool,
) -> list[StageStatus]:
    """Return stage statuses visible in the GUI checklist."""
    if show_optional:
        return list(statuses)
    return [item for item in statuses if not item.stage.optional]


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
            QCheckBox,
            QComboBox,
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
        self._QCheckBox = QCheckBox
        self.project: PipelineProject | None = None
        self._workflow_preview: str | None = None
        self._statuses: list[StageStatus] = []
        self._show_optional_stages = True
        self._active_controller: StageController | None = None
        self._active_stage_id: str | None = None
        self._opening_stage = False
        self._worker = None
        self._stage_cancel_event: threading.Event | None = None
        self._pending_auto_stage: str | None = None

        self._panel = self._build_panel()
        self._lightsuite_dock = self.viewer.window.add_dock_widget(
            self._panel,
            area="left",
            name="LightSuite",
            tabify=True,
            add_vertical_stretch=False,
        )
        self._lightsuite_dock.show()
        self._lightsuite_dock.raise_()

        self._config_editor = ConfigEditorDock(
            on_saved=self._on_config_saved,
            on_log=self.log,
            on_template_loaded=self.preview_workflow,
        )
        self._config_dock = self.viewer.window.add_dock_widget(
            self._config_editor.widget,
            area="right",
            name="Config",
        )
        self._config_dock.show()

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
        QCheckBox = self._QCheckBox
        QListWidget = self._QListWidget

        panel = QWidget()
        layout = QVBoxLayout(panel)

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

        self._optional_toggle = QPushButton("Hide optional")
        self._optional_toggle.clicked.connect(self._toggle_optional_stages)
        self._optional_toggle.setEnabled(False)
        layout.addWidget(self._optional_toggle)

        self._stage_list = QListWidget()
        self._stage_list.itemDoubleClicked.connect(self._on_stage_activated)
        from qtpy.QtWidgets import QSizePolicy

        self._stage_list.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Maximum,
        )
        layout.addWidget(self._stage_list, stretch=0)

        self._export_row = QWidget()
        export_layout = QHBoxLayout(self._export_row)
        export_layout.setContentsMargins(0, 0, 0, 0)
        export_layout.addWidget(QLabel("Export spaces:"))
        self._export_atlas_check = QCheckBox("Atlas")
        self._export_atlas_check.setToolTip(
            "Warp channels to atlas / Fiederling native grid (volume_registered/*.tiff)"
        )
        self._export_sample_check = QCheckBox("Sample")
        self._export_sample_check.setToolTip(
            "Warp atlas labels onto the straightened registration grid "
            "(volume_registered/sample_space/)"
        )
        export_layout.addWidget(self._export_atlas_check)
        export_layout.addWidget(self._export_sample_check)
        export_layout.addStretch(1)
        self._export_atlas_check.stateChanged.connect(self._on_export_space_changed)
        self._export_sample_check.stateChanged.connect(self._on_export_space_changed)
        self._export_row.setVisible(False)
        layout.addWidget(self._export_row)

        self._action_button = QPushButton("Run")
        self._action_button.clicked.connect(self._on_action)
        layout.addWidget(self._action_button)
        self._cancel_button = QPushButton("Cancel stage")
        self._cancel_button.clicked.connect(self._on_cancel_stage)
        self._cancel_button.setVisible(False)
        layout.addWidget(self._cancel_button)
        self._stage_list.currentItemChanged.connect(self._on_stage_selection_changed)

        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.setPlaceholderText("Pipeline log…")
        layout.addWidget(self._log, stretch=1)

        self._set_actions_enabled(False)
        return panel

    def _set_actions_enabled(self, enabled: bool) -> None:
        self._optional_toggle.setEnabled(enabled)
        self._stage_list.setEnabled(enabled)
        self._action_button.setEnabled(enabled)
        if self._channel_row.isVisible():
            self._channel_combo.setEnabled(enabled)
            self._channel_apply_button.setEnabled(enabled)
        if self._export_row.isVisible():
            self._export_atlas_check.setEnabled(enabled)
            self._export_sample_check.setEnabled(enabled)
        self._update_action_button()

    def _on_stage_selection_changed(self, _current: Any = None, _previous: Any = None) -> None:
        self._update_export_row_visibility()
        self._update_action_button()

    def _update_action_button(self) -> None:
        if self.project is None:
            if self._workflow_preview is not None:
                self._action_button.setText("Save config first")
                self._action_button.setToolTip(
                    "Fill sample paths in the Config panel, then Save as… to enable running stages"
                )
            else:
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
            elif stage_id == "export" and self._export_row.isVisible():
                try:
                    spaces_label = format_export_spaces(self._selected_export_spaces())
                except ValueError:
                    spaces_label = "none selected"
                self._action_button.setText("Run stage")
                self._action_button.setToolTip(f"Run export ({spaces_label})")
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
        visible_ids = {item.stage.id for item in self._visible_statuses()}
        return [
            item
            for item in self._statuses
            if item.state == StageState.PENDING and item.stage.id in visible_ids
        ]

    def _visible_statuses(self) -> list[StageStatus]:
        return filter_stage_statuses(
            self._statuses,
            show_optional=self._show_optional_stages,
        )

    def _status_by_id(self, stage_id: str) -> StageStatus | None:
        for item in self._statuses:
            if item.stage.id == stage_id:
                return item
        return None

    def _toggle_optional_stages(self) -> None:
        self._show_optional_stages = not self._show_optional_stages
        self._update_optional_toggle_label()
        self._populate_stage_list()
        self._update_export_row_visibility()
        self._update_action_button()

    def _update_optional_toggle_label(self) -> None:
        if self._show_optional_stages:
            self._optional_toggle.setText("Hide optional")
            self._optional_toggle.setToolTip("Hide optional pipeline stages in the list")
        else:
            self._optional_toggle.setText("Show optional")
            self._optional_toggle.setToolTip("Show optional pipeline stages in the list")

    def _format_stage_list_label(self, item: StageStatus) -> str:
        stage = item.stage
        icon = _STATE_ICONS.get(item.state, "?")
        label = stage.title
        if stage.manual:
            label += " (GUI)"
        if stage.optional:
            label += " (optional)"
        return f"{icon}  {label}"

    def _populate_stage_list(self) -> None:
        current_id: str | None = None
        current_item = self._stage_list.currentItem()
        if current_item is not None:
            current_id = str(current_item.data(self._Qt.UserRole))

        self._stage_list.clear()
        for item in self._visible_statuses():
            stage = item.stage
            list_item = self._QListWidgetItem(self._format_stage_list_label(item))
            list_item.setData(self._Qt.UserRole, stage.id)
            list_item.setToolTip(item.detail or stage.checkpoint_hint)
            self._stage_list.addItem(list_item)
            if current_id is not None and stage.id == current_id:
                self._stage_list.setCurrentItem(list_item)
        self._fit_stage_list_height()

    def _fit_stage_list_height(self) -> None:
        """Size the stage list to its rows so the log can use remaining space."""
        count = self._stage_list.count()
        if count == 0:
            self._stage_list.setFixedHeight(0)
            return
        row_height = self._stage_list.sizeHintForRow(0)
        if row_height <= 0:
            row_height = self._stage_list.fontMetrics().height() + 6
        frame = self._stage_list.frameWidth() * 2
        visible_rows = min(count, 10)
        height = row_height * visible_rows + frame + 2
        self._stage_list.setFixedHeight(height)

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

    def _load_export_space_defaults(self) -> None:
        """Initialize export-space checkboxes from YAML (row stays hidden until export)."""
        if self.project is None or self.project.workflow not in {"brain", "spinal"}:
            return
        atlas_default, sample_default = default_export_space_checks(
            getattr(self.project.config.export, "spaces", None),
        )
        self._export_atlas_check.blockSignals(True)
        self._export_sample_check.blockSignals(True)
        self._export_atlas_check.setChecked(atlas_default)
        self._export_sample_check.setChecked(sample_default)
        self._export_atlas_check.blockSignals(False)
        self._export_sample_check.blockSignals(False)

    def _update_export_row_visibility(self) -> None:
        if self.project is None or self.project.workflow not in {"brain", "spinal"}:
            self._export_row.setVisible(False)
            return
        selected = self._selected_stage()
        show = selected is not None and selected[0] == "export"
        self._export_row.setVisible(show)

    def _export_spaces_for_stage(self, stage_id: str | None) -> list[str] | None:
        if stage_id != "export" or self.project is None:
            return None
        if self.project.workflow not in {"brain", "spinal"}:
            return None
        if self._export_row.isVisible():
            return self._selected_export_spaces()
        atlas_default, sample_default = default_export_space_checks(
            getattr(self.project.config.export, "spaces", None),
        )
        return export_spaces_from_checks(atlas=atlas_default, sample=sample_default)

    def _selected_export_spaces(self) -> list[str]:
        return export_spaces_from_checks(
            atlas=self._export_atlas_check.isChecked(),
            sample=self._export_sample_check.isChecked(),
        )

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

    def _flush_ui(self) -> None:
        """Paint pending log lines before a blocking stage open."""
        try:
            from qtpy.QtWidgets import QApplication
        except ImportError:
            return
        app = QApplication.instance()
        if app is not None:
            app.processEvents()

    def _is_busy(self) -> bool:
        return self._worker is not None or self._opening_stage

    def _stage_title(self, stage_id: str) -> str:
        for item in self._statuses:
            if item.stage.id == stage_id:
                return item.stage.title
        return stage_id

    def _loading_message(self, stage_id: str) -> str:
        return interactive_loading_message(self._stage_title(stage_id), stage_id)

    def preview_workflow(self, workflow: str) -> None:
        """Show the stage checklist for a workflow before a config file is saved."""
        from lightsuite.cli.stages import preview_stage_statuses

        self._teardown_active_stage()
        self.project = None
        self._workflow_preview = workflow
        self._statuses = preview_stage_statuses(workflow)
        self._populate_stage_list()
        self._optional_toggle.setEnabled(True)
        self._stage_list.setEnabled(True)
        self._action_button.setEnabled(False)
        self._channel_row.setVisible(False)
        self._update_export_row_visibility()
        self._update_action_button()
        self.log(
            f"Loaded {workflow} template — stages shown as preview. "
            "Fill paths in Config, Save as…, then stages become runnable."
        )

    def load_project(self, config_path: str | Path) -> None:
        """Load a YAML config and populate the stage checklist."""
        try:
            path = Path(config_path).expanduser().resolve()
            workflow, config = load_project(path)
        except Exception as exc:
            self.log(f"Failed to load config: {exc}")
            return
        self._workflow_preview = None
        self.project = PipelineProject(workflow=workflow, config_path=path, config=config)
        self.log(f"Loaded {workflow} config for sample '{config.sample.name}'")
        self._config_editor.load_from_path(path)
        self._set_actions_enabled(True)
        self._update_channel_selector()
        self._load_export_space_defaults()
        self._update_export_row_visibility()
        self.refresh_statuses()

    def _on_config_saved(self, config_path: Path) -> None:
        """Reload the shell after the YAML editor saves a valid config."""
        self.load_project(config_path)

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
        self._populate_stage_list()
        self._update_export_row_visibility()
        self._update_action_button()
        self._config_editor.maybe_reload_from_disk()

    def _selected_stage(self) -> tuple[str, StageStatus] | None:
        item = self._stage_list.currentItem()
        if item is None:
            return None
        stage_id = item.data(self._Qt.UserRole)
        if not stage_id:
            return None
        status = self._status_by_id(str(stage_id))
        if status is None:
            return None
        return str(stage_id), status

    def _stage_context(self, stage_id: str | None = None) -> StageContext:
        if self.project is None:
            msg = "No project loaded"
            raise RuntimeError(msg)
        export_spaces: list[str] | None = None
        if stage_id == "export":
            export_spaces = self._export_spaces_for_stage(stage_id)
        return StageContext(
            config_path=self.project.config_path,
            headless=True,
            force_preprocess=False,
            export_spaces=export_spaces,
        )

    def _teardown_active_stage(self, *, defer_layer_clear: bool = False) -> None:
        if self._active_controller is not None:
            self._active_controller.teardown(self.viewer)
            self._active_controller = None
        self._active_stage_id = None
        if defer_layer_clear:
            defer_clear_viewer_layers(self.viewer)
        else:
            clear_viewer_layers_safely(self.viewer)

    def finish_interactive_stage(self, *, refresh: bool = True) -> None:
        """Tear down the active GUI stage without closing the shell window."""
        stage_id = self._active_stage_id
        self._teardown_active_stage(defer_layer_clear=True)
        if refresh:
            self.refresh_statuses()
        if stage_id is not None:
            self.log(f"Saved and closed stage: {stage_id}")

        pending = self._pending_auto_stage
        self._pending_auto_stage = None
        if pending is None or self.project is None:
            return
        config = self.project.config
        if pending == "preprocess" and cord_orientation_missing(config):
            self.log("Preprocess paused: save cord orientation first.")
            return
        from qtpy.QtCore import QTimer

        QTimer.singleShot(0, lambda: self._run_auto_stage(pending))

    def _attach_interactive_stage(self, stage_id: str) -> None:
        if self.project is None or self._opening_stage:
            return
        factory = get_stage_attach(self.project.workflow, stage_id)
        if factory is None:
            self.log(f"Stage {stage_id!r} is not interactive.")
            return
        self._teardown_active_stage()
        self._opening_stage = True
        self._set_running(True)
        self.log(self._loading_message(stage_id))
        self._flush_ui()

        ctx = StageContext(
            config_path=self.project.config_path,
            headless=False,
            force_preprocess=False,
        )
        config = self.project.config

        def _open_stage() -> None:
            try:
                try:
                    validate_gui_dependencies()
                    controller = factory(self.viewer, config, ctx)
                except (ImportError, RuntimeError, FileNotFoundError) as exc:
                    message = str(exc)
                    if "default_en.txt" in message:
                        message = (
                            "The pint package in this environment is incomplete "
                            "(missing default_en.txt). Repair with: "
                            "uv sync --extra gui --reinstall-package pint, then restart the GUI."
                        )
                    self.log(f"Failed to open {stage_id}: {message}")
                    try:
                        from napari.utils.notifications import show_warning

                        show_warning(message)
                    except ImportError:
                        pass
                    return
                controller.mount(self.viewer)
                self._active_controller = controller
                self._active_stage_id = stage_id
                open_log = getattr(controller, "open_log_message", None)
                if open_log:
                    self.log(open_log)
                elif stage_id != "view-registration":
                    self.log(f"Opened interactive stage: {stage_id}")
            finally:
                self._opening_stage = False
                self._set_running(False)

        from qtpy.QtCore import QTimer

        QTimer.singleShot(0, _open_stage)

    def _emit_log(self, text: str) -> None:
        self._log_emitter.message.emit(text)

    def _run_auto_stage(self, stage_id: str) -> None:
        if self.project is None or self._is_busy():
            return

        if self._config_editor.is_dirty:
            if self._config_editor.config_path is None:
                self.log("Save the config before running stages.")
                return
            self._config_editor.save_if_dirty()

        self.refresh_statuses()
        if self.project is None:
            return

        workflow = self.project.workflow
        config = self.project.config
        if (
            workflow == "spinal"
            and stage_id == "preprocess"
            and cord_orientation_missing(config)
        ):
            self._pending_auto_stage = "preprocess"
            self.log("Cord orientation required — opening check-orientation before preprocess.")
            self._attach_interactive_stage("check-orientation")
            return

        try:
            ctx = self._stage_context(stage_id)
        except ValueError as exc:
            self.log(str(exc))
            return

        reporter = CallbackReporter(on_message=self._emit_log)

        def _work() -> Any:
            with stage_cancellation(self._stage_cancel_event):
                with capture_pipeline_output(reporter):
                    return run_stage(workflow, stage_id, config, ctx)

        self._stage_cancel_event = threading.Event()
        self._set_running(True)
        if stage_id == "export" and self._export_row.isVisible():
            try:
                spaces_label = format_export_spaces(ctx.export_spaces or [])
            except (TypeError, ValueError):
                spaces_label = "atlas + sample"
            self.log(f"Running {stage_id} ({spaces_label})…")
        else:
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
        self._action_button.setVisible(not running)
        self._cancel_button.setVisible(running)
        self._optional_toggle.setEnabled(enabled)
        self._stage_list.setEnabled(enabled)
        if not running:
            self._stage_cancel_event = None
        if not running:
            self._update_action_button()

    def _on_cancel_stage(self) -> None:
        if self._stage_cancel_event is not None:
            self._stage_cancel_event.set()
        self.log("Cancelling stage…")

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
        if isinstance(exc, StageCancelledError):
            self.log("Stage cancelled.")
            return
        self.log(f"Stage failed: {exc}")
        self.log(traceback.format_exc())

    def _on_stage_activated(self, _item: Any) -> None:
        if self.project is None:
            self.log("Save a valid config file before running stages.")
            return
        selected = self._selected_stage()
        if selected is None:
            return
        stage_id, status = selected
        spec = next(item.stage for item in self._statuses if item.stage.id == stage_id)
        if stage_kind(spec) == StageKind.INTERACTIVE:
            self._attach_interactive_stage(stage_id)
        elif status.state != StageState.DONE:
            self._run_auto_stage(stage_id)

    def _on_export_space_changed(self, _state: int = 0) -> None:
        self._update_action_button()

    def _on_action(self) -> None:
        if self.project is None or self._is_busy():
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
        if self.project is None or self._is_busy():
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
