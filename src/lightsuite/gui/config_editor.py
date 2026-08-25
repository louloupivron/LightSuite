"""Form-based config editor dock for the LightSuite Napari shell."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from lightsuite.config.loader import write_config_yaml
from lightsuite.exceptions import LightsuiteConfigError
from lightsuite.atlas.brainglobe_backend import brainglobe_name_key
from lightsuite.gui.config_form_tooltips import tooltips_for_workflow
from lightsuite.analysis.intensity_metrics import DEFAULT_INTENSITY_METRICS
from lightsuite.gui.config_form_data import (
    AnnotationImportRow,
    BrainFormState,
    ChannelPaths,
    MultiresFormState,
    SpinalFormState,
    brain_atlas_help_text,
    brain_atlas_provider_options,
    brain_brainglobe_catalog,
    default_local_atlas_resolution_um,
    brain_form_from_raw,
    brain_form_to_raw,
    dump_config_dict,
    load_raw_config,
    load_template_raw,
    merge_yaml_only_multires_fields,
    multires_form_from_raw,
    multires_form_to_raw,
    parse_orientation_text,
    resolve_brain_brainglobe_form_fields,
    spinal_form_from_raw,
    spinal_form_to_raw,
    try_validate_config_dict,
)


_VOLUME_FILE_FILTER = "TIFF images (*.tif *.tiff);;All files (*)"


def _browse_directory(parent: Any, title: str, start: str) -> str:
    from qtpy.QtWidgets import QFileDialog

    path = QFileDialog.getExistingDirectory(parent, title, start)
    return path or ""


def _browse_file(parent: Any, title: str, start: str, name_filter: str = "") -> str:
    from qtpy.QtWidgets import QFileDialog

    path, _ = QFileDialog.getOpenFileName(parent, title, start, name_filter)
    return path or ""


def _browse_file_or_directory(
    parent: Any,
    title: str,
    start: str,
    *,
    button: Any | None = None,
    name_filter: str = "",
) -> str:
    """Native file or folder picker via a short menu on the Browse button."""
    from qtpy.QtGui import QCursor
    from qtpy.QtWidgets import QMenu

    menu = QMenu(parent)
    file_action = menu.addAction("Select file…")
    folder_action = menu.addAction("Select folder…")
    if button is not None:
        pos = button.mapToGlobal(button.rect().bottomLeft())
    else:
        pos = QCursor.pos()
    chosen = menu.exec(pos)
    if chosen is file_action:
        return _browse_file(parent, title, start, name_filter)
    if chosen is folder_action:
        return _browse_directory(parent, title, start)
    return ""


_INTENSITY_METRIC_OPTIONS: list[tuple[str, str]] = [
    ("median_intensity", "Median"),
    ("mean_intensity", "Mean"),
    ("std", "Std"),
    ("variance", "Variance"),
    ("volume_mm3", "Volume mm³"),
]


@dataclass
class _ChannelRow:
    row_widget: Any
    field: "_PathField"


@dataclass
class _AnnotationRow:
    row_widget: Any
    format_combo: Any
    path_field: _PathField
    label_edit: Any


class _PathField:
    """Line edit with a browse button for directories, files, or either."""

    def __init__(
        self,
        *,
        parent: Any,
        browse_label: str,
        browse_mode: str = "directory",
        on_change: Callable[[], None] | None = None,
        name_filter: str = "",
    ) -> None:
        from qtpy.QtWidgets import QHBoxLayout, QLineEdit, QPushButton, QWidget

        self.widget = QWidget(parent)
        layout = QHBoxLayout(self.widget)
        layout.setContentsMargins(0, 0, 0, 0)
        self.line = QLineEdit(parent)
        if on_change is not None:
            self.line.textChanged.connect(on_change)
        self._browse_mode = browse_mode
        self._browse_label = browse_label
        self._name_filter = name_filter
        self._browse = QPushButton("Browse…", parent)
        self._browse.clicked.connect(self._on_browse)
        layout.addWidget(self.line, stretch=1)
        layout.addWidget(self._browse)

    def _on_browse(self) -> None:
        start = self.line.text().strip() or str(Path.home())
        if self._browse_mode == "file":
            path = _browse_file(self.widget, self._browse_label, start, self._name_filter)
        elif self._browse_mode == "file_or_directory":
            path = _browse_file_or_directory(
                self.widget,
                self._browse_label,
                start,
                button=self._browse,
                name_filter=self._name_filter,
            )
        else:
            path = _browse_directory(self.widget, self._browse_label, start)
        if path:
            self.line.setText(path)

    def text(self) -> str:
        return self.line.text().strip()

    def set_text(self, value: str) -> None:
        self.line.setText(value)

    def set_tooltip(self, text: str) -> None:
        self.widget.setToolTip(text)
        self.line.setToolTip(text)
        self._browse.setToolTip(text)


class ConfigEditorDock:
    """Dock widget for editing pipeline configs via form fields."""

    def __init__(
        self,
        *,
        on_saved: Callable[[Path], None],
        on_log: Callable[[str], None],
        on_template_loaded: Callable[[str], None] | None = None,
    ) -> None:
        from qtpy.QtCore import Qt
        from qtpy.QtWidgets import (
            QCheckBox,
            QComboBox,
            QDoubleSpinBox,
            QFormLayout,
            QHBoxLayout,
            QLabel,
            QPushButton,
            QScrollArea,
            QSpinBox,
            QVBoxLayout,
            QWidget,
        )

        self._on_saved = on_saved
        self._on_log = on_log
        self._on_template_loaded = on_template_loaded
        self._config_path: Path | None = None
        self._workflow: str | None = None
        self._raw: dict[str, Any] = {}
        self._dirty = False
        self._loading = False

        self.widget = QWidget()
        root = QVBoxLayout(self.widget)
        root.setContentsMargins(4, 4, 4, 4)

        file_row = QHBoxLayout()
        file_row.addWidget(QLabel("File:"))
        self._path_label = QLabel("No config loaded")
        self._path_label.setWordWrap(True)
        self._path_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        file_row.addWidget(self._path_label, stretch=1)
        self._open_button = QPushButton("Open…")
        self._open_button.clicked.connect(self._on_open_config)
        file_row.addWidget(self._open_button)
        root.addLayout(file_row)

        template_row = QHBoxLayout()
        template_row.addWidget(QLabel("New:"))
        self._template_combo = QComboBox()
        self._template_combo.addItem("brain", "brain")
        self._template_combo.addItem("spinal", "spinal")
        self._template_combo.addItem("multires", "multires")
        self._template_button = QPushButton("From template")
        self._template_button.clicked.connect(self._on_load_template)
        template_row.addWidget(self._template_combo, stretch=1)
        template_row.addWidget(self._template_button)
        root.addLayout(template_row)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        self._form_host = QWidget()
        self._form_layout = QVBoxLayout(self._form_host)
        scroll.setWidget(self._form_host)
        root.addWidget(scroll, stretch=1)

        self._placeholder = QLabel("Open a config or start from a template.")
        self._placeholder.setWordWrap(True)
        self._form_layout.addWidget(self._placeholder)

        self._config_block = QWidget()
        self._config_form = QFormLayout(self._config_block)
        self._config_form.setContentsMargins(0, 0, 0, 0)
        self._name_edit = self._line_edit(self._mark_dirty)
        self._config_form.addRow("Name", self._name_edit)
        self._tiff_type_combo = QComboBox()
        self._tiff_type_combo.addItem("channelperfile", "channelperfile")
        self._tiff_type_combo.addItem("planeperfile", "planeperfile")
        self._tiff_type_combo.currentIndexChanged.connect(self._on_tiff_type_changed)
        self._config_form.addRow("TIFF layout", self._tiff_type_combo)
        self._source_path = _PathField(
            parent=self._config_block,
            browse_label="Sample source folder",
            on_change=self._mark_dirty,
        )
        self._config_form.addRow("Source path", self._source_path.widget)
        self._channel_list_host = QWidget()
        self._channel_list_layout = QVBoxLayout(self._channel_list_host)
        self._channel_list_layout.setContentsMargins(0, 0, 0, 0)
        channel_buttons = QHBoxLayout()
        self._add_channel_button = QPushButton("Add channel folder")
        self._add_channel_button.clicked.connect(lambda _checked=False: self._add_channel_row())
        channel_buttons.addWidget(self._add_channel_button)
        channel_buttons.addStretch(1)
        channel_column = QVBoxLayout()
        channel_column.addWidget(self._channel_list_host)
        channel_column.addLayout(channel_buttons)
        self._channel_list_wrapper = QWidget()
        self._channel_list_wrapper.setLayout(channel_column)
        self._config_form.addRow("Channel folders", self._channel_list_wrapper)
        self._scratch_path = _PathField(
            parent=self._config_block,
            browse_label="Scratch directory",
            on_change=self._mark_dirty,
        )
        self._config_form.addRow("Scratch", self._scratch_path.widget)
        self._save_path = _PathField(
            parent=self._config_block,
            browse_label="Results directory",
            on_change=self._mark_dirty,
        )
        self._config_form.addRow("Save path", self._save_path.widget)
        self._voxel_x = self._float_spin(0.1, 100.0, self._mark_dirty)
        self._voxel_y = self._float_spin(0.1, 100.0, self._mark_dirty)
        self._voxel_z = self._float_spin(0.1, 100.0, self._mark_dirty)
        voxel_row = QHBoxLayout()
        voxel_row.addWidget(self._voxel_x)
        voxel_row.addWidget(self._voxel_y)
        voxel_row.addWidget(self._voxel_z)
        voxel_widget = QWidget()
        voxel_widget.setLayout(voxel_row)
        self._voxel_widget = voxel_widget
        self._config_form.addRow("Voxel µm (x, y, z)", voxel_widget)
        self._atlas_source_combo = QComboBox()
        self._atlas_source_combo.addItem("Local NIfTI files", "files")
        self._atlas_source_combo.addItem("BrainGlobe (auto-download)", "brainglobe")
        self._atlas_source_combo.currentIndexChanged.connect(self._on_atlas_fields_changed)
        self._config_form.addRow("Atlas source", self._atlas_source_combo)
        self._atlas_provider_combo = QComboBox()
        self._atlas_provider_combo.currentIndexChanged.connect(self._on_atlas_fields_changed)
        self._config_form.addRow("Atlas", self._atlas_provider_combo)
        self._brainglobe_atlas_combo = QComboBox()
        self._brainglobe_atlas_combo.currentIndexChanged.connect(self._on_atlas_fields_changed)
        self._config_form.addRow("BrainGlobe atlas", self._brainglobe_atlas_combo)
        self._atlas_help_label = QLabel("")
        self._atlas_help_label.setWordWrap(True)
        self._atlas_help_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self._atlas_help_label.setStyleSheet("color: palette(mid); font-size: 11px;")
        self._config_form.addRow(self._atlas_help_label)
        self._atlas_resolution = self._float_spin(1.0, 100.0, self._mark_dirty)
        self._config_form.addRow("Atlas resolution µm", self._atlas_resolution)
        self._atlas_dir = _PathField(
            parent=self._config_block,
            browse_label="Atlas directory",
            on_change=self._mark_dirty,
        )
        self._config_form.addRow("Atlas directory", self._atlas_dir.widget)
        self._channel_primary = self._int_spin(1, 16, self._mark_dirty)
        self._config_form.addRow("Primary channel", self._channel_primary)
        self._channel_secondary = self._int_spin(0, 16, self._mark_dirty)
        self._channel_secondary.setSpecialValueText("none")
        self._config_form.addRow("Secondary channel", self._channel_secondary)
        self._registration_resolution = self._float_spin(1.0, 100.0, self._mark_dirty)
        self._config_form.addRow("Registration resolution µm", self._registration_resolution)
        self._bspline_spatial_scale = self._float_spin(0.01, 10.0, self._mark_dirty)
        self._config_form.addRow("B-spline grid spacing mm", self._bspline_spatial_scale)
        self._control_point_weight = self._float_spin(0.0, 1.0, self._mark_dirty)
        self._control_point_weight.setSingleStep(0.05)
        self._config_form.addRow("Control point weight", self._control_point_weight)
        self._augment_points_check = QCheckBox("Add auto-landmarks to control points")
        self._augment_points_check.stateChanged.connect(self._mark_dirty)
        self._config_form.addRow(self._augment_points_check)
        self._dual_channel_mi_primary = self._float_spin(0.0, 1.0, self._mark_dirty)
        self._dual_channel_mi_primary.setSingleStep(0.05)
        self._config_form.addRow("Dual MI weight (primary)", self._dual_channel_mi_primary)
        self._dual_channel_mi_secondary = self._float_spin(0.0, 1.0, self._mark_dirty)
        self._dual_channel_mi_secondary.setSingleStep(0.05)
        self._config_form.addRow("Dual MI weight (secondary)", self._dual_channel_mi_secondary)
        self._orientation_edit = self._line_edit(self._mark_dirty)
        self._orientation_edit.setPlaceholderText("1, 2, 3 (empty = from check-orientation)")
        self._config_form.addRow("Orientation", self._orientation_edit)
        self._canvas_mode_combo = QComboBox()
        self._canvas_mode_combo.addItem("None", "off")
        self._canvas_mode_combo.addItem("Pad to atlas", "pad")
        self._canvas_mode_combo.addItem("Crop to sample", "crop")
        self._canvas_mode_combo.addItem("Union bbox", "union")
        self._canvas_mode_combo.currentIndexChanged.connect(self._mark_dirty)
        self._config_form.addRow("Registration canvas", self._canvas_mode_combo)
        self._geometry_mode_combo = QComboBox()
        self._geometry_mode_combo.addItem("Metadata (manifest)", "metadata")
        self._geometry_mode_combo.addItem("Hybrid (landmarks)", "hybrid")
        self._geometry_mode_combo.currentIndexChanged.connect(self._mark_dirty)
        self._config_form.addRow("Geometry mode", self._geometry_mode_combo)
        self._landmark_fit_mode_combo = QComboBox()
        self._landmark_fit_mode_combo.addItem("Similarity", "similarity")
        self._landmark_fit_mode_combo.addItem("Affine", "affine")
        self._landmark_fit_mode_combo.addItem("Rigid", "rigid")
        self._landmark_fit_mode_combo.currentIndexChanged.connect(self._mark_dirty)
        self._config_form.addRow("Landmark fit mode", self._landmark_fit_mode_combo)
        self._overlap_margin_um = self._float_spin(-500.0, 500.0, self._mark_dirty)
        self._config_form.addRow("Overlap margin µm", self._overlap_margin_um)
        self._write_full_overview_canvas = QCheckBox("Write full overview canvas")
        self._write_full_overview_canvas.stateChanged.connect(self._mark_dirty)
        self._config_form.addRow(self._write_full_overview_canvas)
        self._import_annotations_host = QWidget()
        self._import_annotations_layout = QVBoxLayout(self._import_annotations_host)
        self._import_annotations_layout.setContentsMargins(0, 0, 0, 0)
        import_annotation_buttons = QHBoxLayout()
        self._add_import_annotation_button = QPushButton("Add annotation")
        self._add_import_annotation_button.clicked.connect(
            lambda _checked=False: self._add_import_annotation_row()
        )
        import_annotation_buttons.addWidget(self._add_import_annotation_button)
        import_annotation_buttons.addStretch(1)
        import_annotation_column = QVBoxLayout()
        import_annotation_column.addWidget(self._import_annotations_host)
        import_annotation_column.addLayout(import_annotation_buttons)
        import_annotation_wrapper = QWidget()
        import_annotation_wrapper.setLayout(import_annotation_column)
        self._import_annotations_wrapper = import_annotation_wrapper
        self._config_form.addRow("Import annotations", import_annotation_wrapper)
        self._intensity_metric_checks: dict[str, Any] = {}
        metrics_host = QWidget()
        metrics_layout = QVBoxLayout(metrics_host)
        metrics_layout.setContentsMargins(0, 0, 0, 0)
        for metric_id, label in _INTENSITY_METRIC_OPTIONS:
            check = QCheckBox(label)
            check.stateChanged.connect(self._mark_dirty)
            self._intensity_metric_checks[metric_id] = check
            metrics_layout.addWidget(check)
        self._analysis_metrics_wrapper = metrics_host
        self._config_form.addRow("Intensity metrics", metrics_host)
        self._parcellate_intensities_check = QCheckBox("Parcellate channel intensities")
        self._parcellate_intensities_check.stateChanged.connect(self._mark_dirty)
        self._config_form.addRow(self._parcellate_intensities_check)
        self._workers = self._int_spin(1, 64, self._mark_dirty)
        self._config_form.addRow("Workers", self._workers)
        self._detection_check = QCheckBox("Enable cell detection")
        self._detection_check.stateChanged.connect(self._mark_dirty)
        self._config_form.addRow(self._detection_check)
        self._pair_label_edit = self._line_edit(self._mark_dirty)
        self._config_form.addRow("Pair label", self._pair_label_edit)
        self._pair_manifest = _PathField(
            parent=self._config_block,
            browse_label="Pair manifest JSON",
            browse_mode="file",
            on_change=self._mark_dirty,
        )
        self._config_form.addRow("Pair manifest", self._pair_manifest.widget)
        self._reference_channel_edit = self._line_edit(self._mark_dirty)
        self._config_form.addRow("Reference channel", self._reference_channel_edit)
        self._multires_channels_host = QWidget()
        self._multires_channels_layout = QVBoxLayout(self._multires_channels_host)
        self._multires_channels_layout.setContentsMargins(0, 0, 0, 0)
        multires_channel_buttons = QHBoxLayout()
        self._add_multires_channel_button = QPushButton("Add channel")
        self._add_multires_channel_button.clicked.connect(
            lambda _checked=False: self._add_multires_channel_row()
        )
        multires_channel_buttons.addWidget(self._add_multires_channel_button)
        multires_channel_buttons.addStretch(1)
        multires_channel_column = QVBoxLayout()
        multires_channel_column.addWidget(self._multires_channels_host)
        multires_channel_column.addLayout(multires_channel_buttons)
        multires_channel_wrapper = QWidget()
        multires_channel_wrapper.setLayout(multires_channel_column)
        self._multires_channels_wrapper = multires_channel_wrapper
        self._config_form.addRow("Multires channels", multires_channel_wrapper)

        self._form_layout.addWidget(self._config_block)
        self._config_block.hide()
        self._apply_field_tooltips()

        button_row = QHBoxLayout()
        self._reload_button = QPushButton("Reload")
        self._reload_button.clicked.connect(self.reload_from_disk)
        self._save_button = QPushButton("Save")
        self._save_button.clicked.connect(self._on_save)
        self._save_as_button = QPushButton("Save as…")
        self._save_as_button.clicked.connect(self._on_save_as)
        button_row.addWidget(self._reload_button)
        button_row.addWidget(self._save_button)
        button_row.addWidget(self._save_as_button)
        root.addLayout(button_row)

        self._status_label = QLabel("Edit fields, then Save.")
        self._status_label.setWordWrap(True)
        root.addWidget(self._status_label)

        self._channel_rows: list[_ChannelRow] = []
        self._multires_channel_rows: list[tuple[Any, _PathField, _PathField]] = []
        self._import_annotation_rows: list[_AnnotationRow] = []
        self._loading_metrics = False
        self._brainglobe_catalog: list[Any] = []
        self._update_buttons()

    def _line_edit(self, on_change: Callable[[], None]) -> Any:
        from qtpy.QtWidgets import QLineEdit

        edit = QLineEdit()
        edit.textChanged.connect(on_change)
        return edit

    def _float_spin(
        self,
        minimum: float,
        maximum: float,
        on_change: Callable[[], None],
    ) -> Any:
        from qtpy.QtWidgets import QDoubleSpinBox

        spin = QDoubleSpinBox()
        spin.setRange(minimum, maximum)
        spin.setDecimals(3)
        spin.setSingleStep(0.1)
        spin.valueChanged.connect(on_change)
        return spin

    def _int_spin(
        self,
        minimum: int,
        maximum: int,
        on_change: Callable[[], None],
    ) -> Any:
        from qtpy.QtWidgets import QSpinBox

        spin = QSpinBox()
        spin.setRange(minimum, maximum)
        spin.valueChanged.connect(on_change)
        return spin

    @property
    def config_path(self) -> Path | None:
        return self._config_path

    @property
    def is_dirty(self) -> bool:
        return self._dirty

    def set_config_path(self, path: Path | None) -> None:
        self._config_path = path.expanduser().resolve() if path is not None else None
        self._update_path_label()

    def reload_from_disk(self) -> None:
        if self._config_path is None or not self._config_path.is_file():
            self._on_log("No config file to reload.")
            return
        try:
            workflow, raw = load_raw_config(self._config_path)
        except LightsuiteConfigError as exc:
            self._on_log(str(exc))
            return
        self._populate_form(workflow, raw)
        self._dirty = False
        self._status_label.setText("Reloaded from disk.")
        self._update_buttons()

    def maybe_reload_from_disk(self) -> None:
        if not self._dirty:
            self.reload_from_disk()

    def load_from_path(self, path: Path) -> None:
        resolved = path.expanduser().resolve()
        if not resolved.is_file():
            self._on_log(f"Config file not found: {resolved}")
            return
        try:
            workflow, raw = load_raw_config(resolved)
        except LightsuiteConfigError as exc:
            self._on_log(str(exc))
            return
        self._config_path = resolved
        self._populate_form(workflow, raw)
        self._dirty = False
        self._update_path_label()
        self._status_label.setText("Loaded from disk.")
        self._update_buttons()

    def clear(self) -> None:
        self._config_path = None
        self._workflow = None
        self._raw = {}
        self._dirty = False
        self._clear_channel_rows()
        self._clear_multires_channel_rows()
        self._clear_import_annotation_rows()
        self._placeholder.show()
        self._config_block.hide()
        self._update_path_label()
        self._status_label.setText("Open a config or start from a template.")
        self._update_buttons()

    def _set_form_row_visible(self, form: Any, field: Any, visible: bool) -> None:
        field.setVisible(visible)
        label = form.labelForField(field)
        if label is not None:
            label.setVisible(visible)

    def _update_workflow_layout(self) -> None:
        workflow = self._workflow
        is_brain = workflow == "brain"
        is_spinal = workflow == "spinal"
        is_multires = workflow == "multires"
        form = self._config_form

        self._set_form_row_visible(form, self._tiff_type_combo, not is_multires)
        self._set_form_row_visible(form, self._source_path.widget, not is_multires)
        self._set_form_row_visible(form, self._channel_list_wrapper, not is_multires)
        voxel_host = self._voxel_x.parentWidget()
        if voxel_host is not None:
            self._set_form_row_visible(form, voxel_host, not is_multires)

        self._set_form_row_visible(form, self._atlas_source_combo, is_brain)
        self._set_form_row_visible(form, self._atlas_provider_combo, is_brain)
        self._set_form_row_visible(form, self._brainglobe_atlas_combo, is_brain)
        self._set_form_row_visible(form, self._atlas_help_label, is_brain)
        self._set_form_row_visible(form, self._atlas_resolution, is_brain)
        self._set_form_row_visible(form, self._atlas_dir.widget, is_brain or is_spinal)

        self._set_form_row_visible(form, self._channel_primary, not is_multires)
        self._set_form_row_visible(form, self._channel_secondary, is_brain)
        self._set_form_row_visible(form, self._registration_resolution, not is_multires)

        self._set_form_row_visible(form, self._bspline_spatial_scale, is_brain)
        self._set_form_row_visible(form, self._control_point_weight, is_brain or is_spinal)
        self._set_form_row_visible(form, self._augment_points_check, is_brain)
        self._set_form_row_visible(form, self._dual_channel_mi_primary, is_brain)
        self._set_form_row_visible(form, self._dual_channel_mi_secondary, is_brain)
        self._set_form_row_visible(form, self._orientation_edit, is_brain)
        self._set_form_row_visible(form, self._canvas_mode_combo, is_brain)
        self._set_form_row_visible(form, self._geometry_mode_combo, is_multires)
        self._set_form_row_visible(form, self._landmark_fit_mode_combo, is_multires)
        self._set_form_row_visible(form, self._overlap_margin_um, is_multires)
        self._set_form_row_visible(form, self._write_full_overview_canvas, is_multires)
        self._set_form_row_visible(
            form,
            self._import_annotations_wrapper,
            is_brain or is_spinal or is_multires,
        )
        self._set_form_row_visible(form, self._analysis_metrics_wrapper, is_brain or is_spinal)
        self._set_form_row_visible(form, self._parcellate_intensities_check, is_spinal)

        self._set_form_row_visible(form, self._workers, True)
        self._set_form_row_visible(form, self._detection_check, is_brain)

        self._set_form_row_visible(form, self._pair_label_edit, is_multires)
        self._set_form_row_visible(form, self._pair_manifest.widget, is_multires)
        self._set_form_row_visible(form, self._reference_channel_edit, is_multires)
        multires_channels_wrapper = self._multires_channels_host.parentWidget()
        if multires_channels_wrapper is not None:
            self._set_form_row_visible(form, multires_channels_wrapper, is_multires)

        if not is_multires:
            self._on_tiff_type_changed()
        if is_brain:
            self._update_atlas_field_state(mark_dirty=False)
        self._apply_field_tooltips()

    def _set_field_tooltip(self, field: Any, text: str) -> None:
        if not text:
            return
        form_field = field
        if isinstance(field, _PathField):
            field.set_tooltip(text)
            form_field = field.widget
        else:
            field.setToolTip(text)
        label = self._config_form.labelForField(form_field)
        if label is not None:
            label.setToolTip(text)

    def _apply_field_tooltips(self) -> None:
        tips = tooltips_for_workflow(self._workflow)
        shared_fields: list[tuple[Any, str]] = [
            (self._name_edit, "name"),
            (self._scratch_path, "scratch"),
            (self._save_path, "save_path"),
            (self._workers, "workers"),
        ]
        for field, key in shared_fields:
            self._set_field_tooltip(field, tips.get(key, ""))

        if self._workflow == "multires":
            multires_fields: list[tuple[Any, str]] = [
                (self._pair_label_edit, "pair_label"),
                (self._pair_manifest, "pair_manifest"),
                (self._reference_channel_edit, "reference_channel"),
                (self._multires_channels_wrapper, "multires_channels"),
                (self._geometry_mode_combo, "geometry_mode"),
                (self._landmark_fit_mode_combo, "landmark_fit_mode"),
                (self._overlap_margin_um, "overlap_margin_um"),
                (self._write_full_overview_canvas, "write_full_overview_canvas"),
                (self._import_annotations_wrapper, "import_annotations"),
            ]
            for field, key in multires_fields:
                self._set_field_tooltip(field, tips.get(key, ""))
            self._add_multires_channel_button.setToolTip(tips.get("multires_channels", ""))
            self._add_import_annotation_button.setToolTip(tips.get("import_annotations", ""))
            return

        sample_fields: list[tuple[Any, str]] = [
            (self._tiff_type_combo, "tiff_type"),
            (self._source_path, "source_path"),
            (self._channel_list_wrapper, "channel_folders"),
            (self._voxel_widget, "voxel_um"),
            (self._channel_primary, "channel_primary"),
            (self._registration_resolution, "registration_resolution"),
        ]
        for field, key in sample_fields:
            self._set_field_tooltip(field, tips.get(key, ""))
        for spin in (self._voxel_x, self._voxel_y, self._voxel_z):
            self._set_field_tooltip(spin, tips.get("voxel_um", ""))
        channel_tip = tips.get("channel_folders", "")
        self._add_channel_button.setToolTip(channel_tip)
        self._add_import_annotation_button.setToolTip(tips.get("import_annotations", ""))
        analysis_fields: list[tuple[Any, str]] = [
            (self._analysis_metrics_wrapper, "intensity_metrics"),
        ]
        for field, key in analysis_fields:
            self._set_field_tooltip(field, tips.get(key, ""))

        if self._workflow == "brain":
            brain_fields: list[tuple[Any, str]] = [
                (self._atlas_source_combo, "atlas_source"),
                (self._atlas_provider_combo, "atlas"),
                (self._brainglobe_atlas_combo, "brainglobe_atlas"),
                (self._atlas_resolution, "atlas_resolution"),
                (self._atlas_dir, "atlas_dir"),
                (self._channel_secondary, "channel_secondary"),
                (self._bspline_spatial_scale, "bspline_spatial_scale_mm"),
                (self._control_point_weight, "control_point_weight"),
                (self._augment_points_check, "augment_points"),
                (self._dual_channel_mi_primary, "dual_channel_mi_weight_primary"),
                (self._dual_channel_mi_secondary, "dual_channel_mi_weight_secondary"),
                (self._orientation_edit, "orientation"),
                (self._canvas_mode_combo, "canvas_mode"),
                (self._import_annotations_wrapper, "import_annotations"),
                (self._detection_check, "detection"),
            ]
            for field, key in brain_fields:
                self._set_field_tooltip(field, tips.get(key, ""))
            for spin in (self._converter_voxel_x, self._converter_voxel_y, self._converter_voxel_z):
                self._set_field_tooltip(spin, tips.get("converter_voxel_um", ""))
        elif self._workflow == "spinal":
            spinal_fields: list[tuple[Any, str]] = [
                (self._atlas_dir, "atlas_dir"),
                (self._control_point_weight, "control_point_weight"),
                (self._import_annotations_wrapper, "import_annotations"),
                (self._parcellate_intensities_check, "parcellate_intensities"),
            ]
            for field, key in spinal_fields:
                self._set_field_tooltip(field, tips.get(key, ""))
            for spin in (self._converter_voxel_x, self._converter_voxel_y, self._converter_voxel_z):
                self._set_field_tooltip(spin, tips.get("converter_voxel_um", ""))

    def _refresh_atlas_provider_combo(self, source: str, *, preferred_provider: str) -> None:
        current = str(self._atlas_provider_combo.currentData() or "")
        self._atlas_provider_combo.blockSignals(True)
        self._atlas_provider_combo.clear()
        options = brain_atlas_provider_options(source)
        provider_ids = [provider for provider, _label in options]
        selected = preferred_provider if preferred_provider in provider_ids else current
        if selected not in provider_ids and provider_ids:
            selected = provider_ids[0]
        for provider, label in options:
            self._atlas_provider_combo.addItem(label, provider)
        if selected:
            self._set_combo_value(self._atlas_provider_combo, selected)
        self._atlas_provider_combo.blockSignals(False)

    def _refresh_brainglobe_atlas_combo(self, *, preferred_name: str) -> None:
        self._brainglobe_catalog = brain_brainglobe_catalog()
        current = str(self._brainglobe_atlas_combo.currentData() or "")
        self._brainglobe_atlas_combo.blockSignals(True)
        self._brainglobe_atlas_combo.clear()
        selected = preferred_name or current
        selected_key = brainglobe_name_key(selected) if selected else ""
        for entry in self._brainglobe_catalog:
            self._brainglobe_atlas_combo.addItem(entry.label, entry.name)
        if selected_key:
            for index in range(self._brainglobe_atlas_combo.count()):
                data = str(self._brainglobe_atlas_combo.itemData(index) or "")
                if brainglobe_name_key(data) == selected_key:
                    self._brainglobe_atlas_combo.setCurrentIndex(index)
                    break
        elif self._brainglobe_atlas_combo.count() > 0:
            self._brainglobe_atlas_combo.setCurrentIndex(0)
        self._brainglobe_atlas_combo.blockSignals(False)

    def _selected_brainglobe_fields(self) -> tuple[str, str, float]:
        name = str(self._brainglobe_atlas_combo.currentData() or "")
        return resolve_brain_brainglobe_form_fields(
            brainglobe_name=name,
            catalog=self._brainglobe_catalog,
        )

    def _update_atlas_field_state(self, *, mark_dirty: bool) -> None:
        if self._workflow != "brain":
            return
        source = str(self._atlas_source_combo.currentData() or "files")
        use_local_files = source != "brainglobe"
        self._set_form_row_visible(self._config_form, self._atlas_provider_combo, use_local_files)
        self._set_form_row_visible(self._config_form, self._brainglobe_atlas_combo, not use_local_files)
        self._set_form_row_visible(self._config_form, self._atlas_resolution, use_local_files)
        self._set_form_row_visible(self._config_form, self._atlas_dir.widget, use_local_files)
        if use_local_files:
            provider = str(self._atlas_provider_combo.currentData() or "allen")
            self._atlas_help_label.setText(brain_atlas_help_text(source, provider))
        else:
            _name, provider, resolution_um = self._selected_brainglobe_fields()
            help_text = brain_atlas_help_text(source, provider)
            self._atlas_help_label.setText(
                f"{help_text} Selected: {resolution_um:g} µm native atlas voxels."
            )
        if mark_dirty and not self._loading:
            self._mark_dirty()

    def _on_atlas_fields_changed(self) -> None:
        if self._loading:
            return
        source = str(self._atlas_source_combo.currentData() or "files")
        if source == "brainglobe":
            if not self._brainglobe_catalog:
                self._refresh_brainglobe_atlas_combo(preferred_name="")
            self._update_atlas_field_state(mark_dirty=True)
            return
        provider = str(self._atlas_provider_combo.currentData() or "allen")
        self._refresh_atlas_provider_combo(source, preferred_provider=provider)
        self._atlas_resolution.setValue(default_local_atlas_resolution_um(provider))
        self._update_atlas_field_state(mark_dirty=True)

    def _populate_form(self, workflow: str, raw: dict[str, Any]) -> None:
        self._loading = True
        self._workflow = workflow
        self._raw = raw
        self._placeholder.hide()
        self._config_block.show()

        if workflow == "brain":
            state = brain_form_from_raw(raw)
            self._fill_sample_fields(
                state.sample_name,
                state.tiff_type,
                state.source_path,
                state.channel_paths,
                state.use_channel_list,
                state.scratch,
                state.save_path,
                state.voxel_um,
            )
            self._set_combo_value(self._atlas_source_combo, state.atlas_source)
            if state.atlas_source == "brainglobe":
                self._refresh_brainglobe_atlas_combo(preferred_name=state.brainglobe_name)
            else:
                self._refresh_atlas_provider_combo(
                    state.atlas_source,
                    preferred_provider=state.atlas_provider,
                )
                self._atlas_resolution.setValue(state.atlas_resolution_um)
            self._atlas_dir.set_text(state.atlas_dir)
            self._channel_primary.setValue(state.channel_primary)
            self._channel_secondary.setValue(state.channel_secondary or 0)
            self._registration_resolution.setValue(state.registration_resolution_um)
            self._bspline_spatial_scale.setValue(state.bspline_spatial_scale_mm)
            self._control_point_weight.setValue(state.control_point_weight)
            self._augment_points_check.setChecked(state.augment_points)
            self._dual_channel_mi_primary.setValue(state.dual_channel_mi_weight_primary)
            self._dual_channel_mi_secondary.setValue(state.dual_channel_mi_weight_secondary)
            self._orientation_edit.setText(
                "" if state.orientation is None else f"{state.orientation[0]}, {state.orientation[1]}, {state.orientation[2]}"
            )
            self._set_combo_value(self._canvas_mode_combo, state.canvas_mode)
            self._fill_import_annotations(state.import_annotations)
            self._fill_intensity_metrics(state.intensity_metrics)
            self._workers.setValue(state.workers)
            self._detection_check.setChecked(state.detection_enabled)
        elif workflow == "spinal":
            state = spinal_form_from_raw(raw)
            self._fill_sample_fields(
                state.sample_name,
                state.tiff_type,
                state.source_path,
                state.channel_paths,
                state.use_channel_list,
                state.scratch,
                state.save_path,
                state.voxel_um,
            )
            self._atlas_dir.set_text(state.atlas_dir)
            self._channel_primary.setValue(state.channel_primary)
            self._registration_resolution.setValue(state.registration_resolution_um)
            self._control_point_weight.setValue(state.control_point_weight)
            self._fill_import_annotations(state.import_annotations)
            self._fill_intensity_metrics(state.intensity_metrics)
            self._parcellate_intensities_check.setChecked(state.parcellate_intensities)
            self._workers.setValue(state.workers)
        else:
            state = multires_form_from_raw(raw)
            self._name_edit.setText(state.sample_name)
            self._save_path.set_text(state.save_path)
            self._scratch_path.set_text(state.scratch)
            self._pair_label_edit.setText(state.pair_label)
            self._pair_manifest.set_text(state.pair_manifest)
            self._reference_channel_edit.setText(state.reference_channel)
            self._set_combo_value(self._geometry_mode_combo, state.geometry_mode)
            self._set_combo_value(self._landmark_fit_mode_combo, state.landmark_fit_mode)
            self._overlap_margin_um.setValue(state.overlap_margin_um)
            self._write_full_overview_canvas.setChecked(state.write_full_overview_canvas)
            self._fill_import_annotations(state.import_annotations)
            compute = raw.get("compute") or {}
            self._workers.setValue(int(compute.get("workers") or 4))
            self._clear_multires_channel_rows()
            for item in state.channels:
                self._add_multires_channel_row(item.name, item.overview, item.roi)
            if not state.channels:
                self._add_multires_channel_row()

        self._loading = False
        self._update_workflow_layout()

    def _fill_intensity_metrics(self, metrics: list[str]) -> None:
        selected = set(metrics or DEFAULT_INTENSITY_METRICS)
        self._loading_metrics = True
        for metric_id, check in self._intensity_metric_checks.items():
            check.setChecked(metric_id in selected)
        self._loading_metrics = False

    def _collect_intensity_metrics(self) -> list[str]:
        metrics = [
            metric_id
            for metric_id, check in self._intensity_metric_checks.items()
            if check.isChecked()
        ]
        return metrics or list(DEFAULT_INTENSITY_METRICS)

    def _fill_import_annotations(self, rows: list[AnnotationImportRow]) -> None:
        self._clear_import_annotation_rows()
        if rows:
            for row in rows:
                self._add_import_annotation_row(row.format, row.path, row.label)
        else:
            self._add_import_annotation_row()

    def _collect_import_annotations(self) -> list[AnnotationImportRow]:
        rows: list[AnnotationImportRow] = []
        for row in self._import_annotation_rows:
            rows.append(
                AnnotationImportRow(
                    format=str(row.format_combo.currentData() or "points_csv"),
                    path=row.path_field.text(),
                    label=row.label_edit.text().strip(),
                )
            )
        return rows

    def _clear_import_annotation_rows(self) -> None:
        self._clear_layout_widgets(self._import_annotations_layout)
        self._import_annotation_rows.clear()

    def _add_import_annotation_row(
        self,
        format_name: str = "points_csv",
        path: str = "",
        label: str = "",
    ) -> None:
        from qtpy.QtWidgets import QComboBox, QHBoxLayout, QLineEdit, QPushButton, QWidget

        row_widget = QWidget(self._import_annotations_host)
        layout = QHBoxLayout(row_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        format_combo = QComboBox()
        format_combo.addItem("points_csv", "points_csv")
        format_combo.addItem("mask_tiff", "mask_tiff")
        format_combo.currentIndexChanged.connect(self._mark_dirty)
        self._set_combo_value(format_combo, format_name)
        path_field = _PathField(
            parent=row_widget,
            browse_label="Annotation file",
            browse_mode="file",
            on_change=self._mark_dirty,
        )
        path_field.set_text(path)
        label_edit = QLineEdit()
        label_edit.setPlaceholderText("label (optional)")
        label_edit.setText(label)
        label_edit.textChanged.connect(self._mark_dirty)
        remove = QPushButton("Remove")
        remove.clicked.connect(
            lambda _checked=False, row_path=path_field: self._remove_import_annotation_row(row_path)
        )
        layout.addWidget(format_combo)
        layout.addWidget(path_field.widget, stretch=1)
        layout.addWidget(label_edit)
        layout.addWidget(remove)
        self._import_annotations_layout.addWidget(row_widget)
        self._import_annotation_rows.append(
            _AnnotationRow(
                row_widget=row_widget,
                format_combo=format_combo,
                path_field=path_field,
                label_edit=label_edit,
            )
        )

    def _remove_import_annotation_row(self, path_field: _PathField) -> None:
        for index, row in enumerate(self._import_annotation_rows):
            if row.path_field is not path_field:
                continue
            row.row_widget.setParent(None)
            row.row_widget.deleteLater()
            self._import_annotation_rows.pop(index)
            self._mark_dirty()
            return

    def _fill_sample_fields(
        self,
        name: str,
        tiff_type: str,
        source_path: str,
        channel_paths: list[str],
        use_channel_list: bool,
        scratch: str,
        save_path: str,
        voxel_um: tuple[float, float, float],
    ) -> None:
        self._name_edit.setText(name)
        self._set_combo_value(self._tiff_type_combo, tiff_type)
        self._source_path.set_text(source_path)
        self._scratch_path.set_text(scratch)
        self._save_path.set_text(save_path)
        self._voxel_x.setValue(voxel_um[0])
        self._voxel_y.setValue(voxel_um[1])
        self._voxel_z.setValue(voxel_um[2])
        self._clear_channel_rows()
        if use_channel_list and channel_paths:
            for path in channel_paths:
                self._add_channel_row(path)
        elif use_channel_list:
            self._add_channel_row()

    def _set_combo_value(self, combo: Any, value: str) -> None:
        index = combo.findData(value)
        if index >= 0:
            combo.setCurrentIndex(index)

    def _on_tiff_type_changed(self) -> None:
        if self._workflow == "multires":
            return
        tiff_type = str(self._tiff_type_combo.currentData())
        multi = tiff_type == "planeperfile"
        self._set_form_row_visible(self._config_form, self._source_path.widget, not multi)
        self._set_form_row_visible(self._config_form, self._channel_list_wrapper, multi)
        if not self._loading:
            self._mark_dirty()

    def _clear_layout_widgets(self, layout: Any) -> None:
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()

    def _clear_channel_rows(self) -> None:
        self._clear_layout_widgets(self._channel_list_layout)
        self._channel_rows.clear()

    def _add_channel_row(self, path: str = "") -> None:
        from qtpy.QtWidgets import QHBoxLayout, QPushButton, QWidget

        row_widget = QWidget(self._channel_list_host)
        layout = QHBoxLayout(row_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        field = _PathField(
            parent=row_widget,
            browse_label="Channel folder",
            on_change=self._mark_dirty,
        )
        field.set_text(path)
        remove = QPushButton("Remove")
        remove.clicked.connect(lambda _checked=False, row_field=field: self._remove_channel_row(row_field))
        layout.addWidget(field.widget, stretch=1)
        layout.addWidget(remove)
        self._channel_list_layout.addWidget(row_widget)
        self._channel_rows.append(_ChannelRow(row_widget=row_widget, field=field))

    def _remove_channel_row(self, field: _PathField) -> None:
        for index, row in enumerate(self._channel_rows):
            if row.field is not field:
                continue
            row.row_widget.setParent(None)
            row.row_widget.deleteLater()
            self._channel_rows.pop(index)
            self._mark_dirty()
            return

    def _clear_multires_channel_rows(self) -> None:
        self._clear_layout_widgets(self._multires_channels_layout)
        self._multires_channel_rows.clear()

    def _add_multires_channel_row(
        self,
        name: str = "",
        overview: str = "",
        roi: str = "",
    ) -> None:
        from qtpy.QtWidgets import QHBoxLayout, QLabel, QLineEdit, QPushButton, QWidget

        row_widget = QWidget(self._multires_channels_host)
        layout = QHBoxLayout(row_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        name_edit = QLineEdit()
        name_edit.setPlaceholderText("488")
        name_edit.setText(name)
        name_edit.textChanged.connect(self._mark_dirty)
        overview_field = _PathField(
            parent=row_widget,
            browse_label="Overview TIFF or folder",
            browse_mode="file_or_directory",
            on_change=self._mark_dirty,
            name_filter=_VOLUME_FILE_FILTER,
        )
        overview_field.set_text(overview)
        roi_field = _PathField(
            parent=row_widget,
            browse_label="ROI TIFF or folder",
            browse_mode="file_or_directory",
            on_change=self._mark_dirty,
            name_filter=_VOLUME_FILE_FILTER,
        )
        roi_field.set_text(roi)
        remove = QPushButton("Remove")
        remove.clicked.connect(
            lambda _checked=False, name=name_edit: self._remove_multires_channel_row(name)
        )
        layout.addWidget(QLabel("Name"))
        layout.addWidget(name_edit)
        layout.addWidget(QLabel("Overview"))
        layout.addWidget(overview_field.widget, stretch=1)
        layout.addWidget(QLabel("ROI"))
        layout.addWidget(roi_field.widget, stretch=1)
        layout.addWidget(remove)
        self._multires_channels_layout.addWidget(row_widget)
        self._multires_channel_rows.append((name_edit, overview_field, roi_field))

    def _remove_multires_channel_row(self, name_edit: Any) -> None:
        for index, row in enumerate(self._multires_channel_rows):
            if row[0] is not name_edit:
                continue
            row_widget = name_edit.parentWidget()
            if row_widget is not None:
                row_widget.setParent(None)
                row_widget.deleteLater()
            self._multires_channel_rows.pop(index)
            self._mark_dirty()
            return

    def _collect_raw(self) -> dict[str, Any]:
        if self._workflow == "brain":
            atlas_source = str(self._atlas_source_combo.currentData() or "files")
            if atlas_source == "brainglobe":
                brainglobe_name, atlas_provider, atlas_resolution_um = self._selected_brainglobe_fields()
            else:
                brainglobe_name = ""
                atlas_provider = str(self._atlas_provider_combo.currentData())
                atlas_resolution_um = self._atlas_resolution.value()
            state = BrainFormState(
                sample_name=self._name_edit.text().strip(),
                source_path=self._source_path.text(),
                channel_paths=[row.field.text() for row in self._channel_rows if row.field.text()],
                use_channel_list=str(self._tiff_type_combo.currentData()) == "planeperfile",
                tiff_type=str(self._tiff_type_combo.currentData()),
                scratch=self._scratch_path.text(),
                save_path=self._save_path.text(),
                voxel_um=(
                    self._voxel_x.value(),
                    self._voxel_y.value(),
                    self._voxel_z.value(),
                ),
                atlas_source=atlas_source,
                atlas_provider=atlas_provider,
                brainglobe_name=brainglobe_name,
                atlas_resolution_um=atlas_resolution_um,
                atlas_dir=self._atlas_dir.text(),
                channel_primary=self._channel_primary.value(),
                channel_secondary=(
                    None if self._channel_secondary.value() == 0 else self._channel_secondary.value()
                ),
                registration_resolution_um=self._registration_resolution.value(),
                bspline_spatial_scale_mm=self._bspline_spatial_scale.value(),
                control_point_weight=self._control_point_weight.value(),
                augment_points=self._augment_points_check.isChecked(),
                dual_channel_mi_weight_primary=self._dual_channel_mi_primary.value(),
                dual_channel_mi_weight_secondary=self._dual_channel_mi_secondary.value(),
                orientation=parse_orientation_text(self._orientation_edit.text()),
                canvas_mode=str(self._canvas_mode_combo.currentData() or "off"),
                import_annotations=self._collect_import_annotations(),
                intensity_metrics=self._collect_intensity_metrics(),
                workers=self._workers.value(),
                detection_enabled=self._detection_check.isChecked(),
            )
            return brain_form_to_raw(state, self._raw)
        if self._workflow == "spinal":
            state = SpinalFormState(
                sample_name=self._name_edit.text().strip(),
                source_path=self._source_path.text(),
                channel_paths=[row.field.text() for row in self._channel_rows if row.field.text()],
                use_channel_list=str(self._tiff_type_combo.currentData()) == "planeperfile",
                tiff_type=str(self._tiff_type_combo.currentData()),
                scratch=self._scratch_path.text(),
                save_path=self._save_path.text(),
                voxel_um=(
                    self._voxel_x.value(),
                    self._voxel_y.value(),
                    self._voxel_z.value(),
                ),
                atlas_dir=self._atlas_dir.text(),
                channel_primary=self._channel_primary.value(),
                registration_resolution_um=self._registration_resolution.value(),
                control_point_weight=self._control_point_weight.value(),
                import_annotations=self._collect_import_annotations(),
                intensity_metrics=self._collect_intensity_metrics(),
                parcellate_intensities=self._parcellate_intensities_check.isChecked(),
                workers=self._workers.value(),
            )
            return spinal_form_to_raw(state, self._raw)
        channels = [
            ChannelPaths(
                name=name_edit.text().strip(),
                overview=overview_field.text(),
                roi=roi_field.text(),
            )
            for name_edit, overview_field, roi_field in self._multires_channel_rows
        ]
        state = MultiresFormState(
            sample_name=self._name_edit.text().strip(),
            save_path=self._save_path.text(),
            scratch=self._scratch_path.text(),
            pair_label=self._pair_label_edit.text().strip(),
            pair_manifest=self._pair_manifest.text(),
            reference_channel=self._reference_channel_edit.text().strip(),
            channels=channels,
            geometry_mode=str(self._geometry_mode_combo.currentData() or "metadata"),
            landmark_fit_mode=str(self._landmark_fit_mode_combo.currentData() or "similarity"),
            overlap_margin_um=self._overlap_margin_um.value(),
            write_full_overview_canvas=self._write_full_overview_canvas.isChecked(),
            import_annotations=self._collect_import_annotations(),
        )
        raw = multires_form_to_raw(state, self._raw)
        compute = dict(raw.get("compute") or {})
        compute["workers"] = self._workers.value()
        raw["compute"] = compute
        return raw

    def _update_path_label(self) -> None:
        if self._config_path is None:
            self._path_label.setText("Unsaved — use Save as…")
            return
        suffix = " *" if self._dirty else ""
        self._path_label.setText(f"{self._config_path}{suffix}")

    def _update_buttons(self) -> None:
        has_form = self._workflow is not None
        has_path = self._config_path is not None and self._config_path.is_file()
        self._reload_button.setEnabled(has_path and not self._loading)
        self._save_button.setEnabled(has_form)
        self._save_as_button.setEnabled(has_form)

    def _mark_dirty(self) -> None:
        if self._loading:
            return
        if getattr(self, "_loading_metrics", False):
            return
        self._dirty = True
        self._update_path_label()
        self._status_label.setText("Unsaved changes.")
        self._update_buttons()

    def _on_open_config(self) -> None:
        from qtpy.QtWidgets import QFileDialog

        start_dir = str(self._config_path.parent) if self._config_path else str(Path.home())
        path, _ = QFileDialog.getOpenFileName(
            self.widget,
            "Open LightSuite config",
            start_dir,
            "YAML files (*.yaml *.yml)",
        )
        if path:
            self.load_from_path(Path(path))
            self._on_saved(Path(path))

    def _on_load_template(self) -> None:
        workflow = str(self._template_combo.currentData())
        try:
            workflow, raw = load_template_raw(workflow)
        except ValueError as exc:
            self._on_log(str(exc))
            return
        self._config_path = None
        self._populate_form(workflow, raw)
        self._dirty = True
        self._update_path_label()
        self._status_label.setText(f"Loaded {workflow} template — edit fields, then Save as…")
        self._on_log(f"Loaded {workflow} config template.")
        if self._on_template_loaded is not None:
            self._on_template_loaded(workflow)
        self._update_buttons()

    def _on_save(self) -> None:
        if self._config_path is None:
            self._on_save_as()
            return
        self._save_to_path(self._config_path)

    def _on_save_as(self) -> None:
        from qtpy.QtWidgets import QFileDialog

        start_dir = str(self._config_path.parent) if self._config_path else str(Path.home())
        path, _ = QFileDialog.getSaveFileName(
            self.widget,
            "Save LightSuite config",
            start_dir,
            "YAML files (*.yaml *.yml)",
        )
        if not path:
            return
        save_path = Path(path)
        if save_path.suffix not in {".yaml", ".yml"}:
            save_path = save_path.with_suffix(".yaml")
        self._save_to_path(save_path)

    def _save_to_path(self, path: Path) -> None:
        if self._workflow is None:
            self._on_log("No config loaded in the editor.")
            return
        self._merge_yaml_only_fields_from_disk(path)
        raw = self._collect_raw()
        text = dump_config_dict(raw)
        try:
            write_config_yaml(path, text)
        except (LightsuiteConfigError, OSError, ValueError) as exc:
            self._on_log(f"Could not save config: {exc}")
            self._status_label.setText(str(exc))
            return

        self._config_path = path.expanduser().resolve()
        self._raw = raw
        self._dirty = False
        self._update_path_label()
        self._on_log(f"Saved config to {self._config_path}")

        validation = try_validate_config_dict(raw)
        if isinstance(validation, LightsuiteConfigError):
            self._status_label.setText(f"Saved (validation pending): {validation}")
            self._on_log(f"Config saved but not loaded — fix validation errors: {validation}")
            return

        workflow, config = validation
        self._status_label.setText(f"Saved and validated ({workflow}, {config.sample.name}).")
        self._on_saved(self._config_path)

    def _merge_yaml_only_fields_from_disk(self, path: Path) -> None:
        """Pull inspect-geometry YAML (lateral_flip) into the editor before Save."""
        if self._workflow != "multires":
            return
        source = path if path.is_file() else self._config_path
        if source is None:
            return
        source = Path(source)
        if not source.is_file():
            return
        try:
            _workflow, disk = load_raw_config(source)
        except LightsuiteConfigError:
            return
        self._raw = merge_yaml_only_multires_fields(self._raw, disk)
