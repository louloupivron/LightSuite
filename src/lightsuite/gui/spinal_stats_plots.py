"""Interactive stats / plots panel for spinal cord region tables."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from lightsuite.analysis.cord_heatmap import (
    DIFFERENCE_CMAP,
    build_cord_heatmap_comparison_figure,
    build_cord_heatmap_figure,
    cord_hemisphere_matrices,
    cord_stats_to_matrix,
    default_heatmap_output_path,
    default_segment_range_bounds,
    difference_color_limits,
    difference_heatmap_matrix,
    discover_cord_plot_options,
    ensure_cord_rollups_in_dataframe,
    load_segment_order,
    metric_label,
    parse_segment_range,
    resolve_region_stats_path,
    structure_acronym_to_label,
    write_cord_heatmap_compare_csv,
    write_cord_heatmap_matrix_csv,
)
from lightsuite.config.models import SpinalCordPipelineConfig
from lightsuite.gui.stage_controller import DockStageController, run_attached_stage

_X_LABELS_BY_ROLLUP = {
    "structure": "Lamina / white matter",
    "division": "Division",
    "horn": "Horn",
    "region": "Atlas region",
}

_STATS_PANEL_TOOLTIPS: dict[str, str] = {
    "plot_type": (
        "Value shown in the heatmap. Intensity metrics use registered channel volumes; "
        "cell count and density come from imported point clouds."
    ),
    "channel": (
        "Intensity channel index or import label to plot. The list updates when you "
        "change plot type."
    ),
    "hemisphere": (
        "Which side to plot. Whole cord, left, or right are single heatmaps. "
        "Compare shows L | R | whole side by side with a shared color scale. "
        "Difference maps left minus right (red = higher on the left). "
        "Compare and difference appear only when region_stats.csv has left/right rows."
    ),
    "rollup": (
        "How fine the X axis is: structure = laminae I–X plus dorsal/lateral/ventral "
        "funiculi; division = gray vs white matter; horn = dorsal/ventral/central; "
        "region = finest atlas parcels."
    ),
    "segments": (
        "Rostrocaudal range to show, inclusive, in atlas order (same as C1:Co2 on the CLI). "
        "Segments with no values stay grey."
    ),
    "normalize": (
        "Rescale colors before plotting. None keeps raw values. Per segment divides each "
        "row by its max; per region divides each column by its max. Ignored for the "
        "left − right difference map."
    ),
    "update": "Redraw the heatmap with the current controls.",
    "save_png": "Save the heatmap currently on screen as a PNG (default: <save_path>/plots/).",
    "save_csv": (
        "Save the displayed matrix as CSV. Compare mode writes columns left, right, and whole; "
        "other modes write a segment × region table."
    ),
}


def _parse_channel_value(text: str) -> int | str:
    if text.isdigit():
        return int(text)
    return text


def _default_cmap(metric: str) -> str:
    if "intensity" in metric or metric == "std":
        return "hot"
    return "viridis"


def _add_form_row(form, title: str, widget, tooltip: str) -> None:
    from qtpy.QtWidgets import QLabel

    label = QLabel(title)
    label.setToolTip(tooltip)
    widget.setToolTip(tooltip)
    form.addRow(label, widget)


class CordStatsPlotsPanel:
    """Qt panel to explore segment × lamina/WM heatmaps from region_stats.csv."""

    def __init__(
        self,
        config: SpinalCordPipelineConfig,
        *,
        stats_df: pd.DataFrame,
        stats_path: Path,
    ) -> None:
        from qtpy.QtWidgets import (
            QComboBox,
            QFormLayout,
            QHBoxLayout,
            QLabel,
            QPushButton,
            QScrollArea,
            QVBoxLayout,
            QWidget,
        )

        self._config = config
        self._stats_path = stats_path
        self._atlas_dir = config.atlas.atlas_dir
        self._segment_order = load_segment_order(self._atlas_dir)
        self._stats_df, self._regions_df = ensure_cord_rollups_in_dataframe(
            stats_df,
            atlas_dir=self._atlas_dir,
        )
        self._options = discover_cord_plot_options(self._stats_df)
        self._has_split = {"left", "right"}.issubset(set(self._options["hemispheres"]))
        self._figure = None
        self._canvas: FigureCanvasQTAgg | None = None
        self._export_kind = "matrix"
        self._export_matrices: dict[str, pd.DataFrame] = {}

        self.widget = QWidget()
        root = QVBoxLayout(self.widget)

        source_label = QLabel(f"Source: {stats_path}")
        source_label.setWordWrap(True)
        root.addWidget(source_label)

        form = QFormLayout()
        self._metric_combo = QComboBox()
        for metric in self._options["metrics"]:
            self._metric_combo.addItem(metric_label(str(metric)), str(metric))
        self._metric_combo.currentIndexChanged.connect(self._on_metric_changed)
        _add_form_row(form, "Plot type:", self._metric_combo, _STATS_PANEL_TOOLTIPS["plot_type"])

        self._channel_combo = QComboBox()
        _add_form_row(form, "Channel / label:", self._channel_combo, _STATS_PANEL_TOOLTIPS["channel"])

        self._hemisphere_combo = QComboBox()
        hemisphere_labels = {"whole": "Whole cord", "left": "Left", "right": "Right"}
        for hemisphere in self._options["hemispheres"]:
            if hemisphere in hemisphere_labels:
                self._hemisphere_combo.addItem(hemisphere_labels[hemisphere], hemisphere)
        if self._has_split:
            self._hemisphere_combo.addItem("Compare L | R | whole", "compare")
            self._hemisphere_combo.addItem("Difference L − R", "difference")
        whole_index = self._hemisphere_combo.findData("whole")
        if whole_index >= 0:
            self._hemisphere_combo.setCurrentIndex(whole_index)
        _add_form_row(form, "Hemisphere:", self._hemisphere_combo, _STATS_PANEL_TOOLTIPS["hemisphere"])

        self._rollup_combo = QComboBox()
        rollup_levels = self._options["rollup_levels"]
        preferred = "structure" if "structure" in rollup_levels else rollup_levels[0]
        for level in rollup_levels:
            self._rollup_combo.addItem(str(level))
        rollup_index = self._rollup_combo.findText(preferred)
        if rollup_index >= 0:
            self._rollup_combo.setCurrentIndex(rollup_index)
        _add_form_row(form, "Rollup:", self._rollup_combo, _STATS_PANEL_TOOLTIPS["rollup"])

        self._segment_start = QComboBox()
        self._segment_end = QComboBox()
        for label in self._segment_order:
            self._segment_start.addItem(label)
            self._segment_end.addItem(label)
        default_start, default_end = default_segment_range_bounds(self._atlas_dir)
        start_index = self._segment_start.findText(default_start)
        end_index = self._segment_end.findText(default_end)
        if start_index >= 0:
            self._segment_start.setCurrentIndex(start_index)
        if end_index >= 0:
            self._segment_end.setCurrentIndex(end_index)
        self._segment_start.setToolTip(_STATS_PANEL_TOOLTIPS["segments"])
        self._segment_end.setToolTip(_STATS_PANEL_TOOLTIPS["segments"])
        segment_row = QWidget()
        segment_layout = QHBoxLayout(segment_row)
        segment_layout.setContentsMargins(0, 0, 0, 0)
        to_label = QLabel("to")
        to_label.setToolTip(_STATS_PANEL_TOOLTIPS["segments"])
        segment_layout.addWidget(self._segment_start)
        segment_layout.addWidget(to_label)
        segment_layout.addWidget(self._segment_end)
        _add_form_row(form, "Segments:", segment_row, _STATS_PANEL_TOOLTIPS["segments"])

        self._normalize_combo = QComboBox()
        for label, value in (
            ("None", "none"),
            ("Per segment (row)", "row"),
            ("Per region (column)", "column"),
        ):
            self._normalize_combo.addItem(label, value)
        _add_form_row(form, "Normalize:", self._normalize_combo, _STATS_PANEL_TOOLTIPS["normalize"])

        root.addLayout(form)

        button_row = QHBoxLayout()
        self._update_button = QPushButton("Update plot")
        self._update_button.clicked.connect(self._update_plot)
        self._save_png_button = QPushButton("Save PNG…")
        self._save_png_button.clicked.connect(self._save_png)
        self._save_csv_button = QPushButton("Save CSV…")
        self._save_csv_button.clicked.connect(self._save_csv)
        self._update_button.setToolTip(_STATS_PANEL_TOOLTIPS["update"])
        self._save_png_button.setToolTip(_STATS_PANEL_TOOLTIPS["save_png"])
        self._save_csv_button.setToolTip(_STATS_PANEL_TOOLTIPS["save_csv"])
        button_row.addWidget(self._update_button)
        button_row.addWidget(self._save_png_button)
        button_row.addWidget(self._save_csv_button)
        root.addLayout(button_row)

        self._status_label = QLabel("")
        self._status_label.setWordWrap(True)
        root.addWidget(self._status_label)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        self._canvas_host = QWidget()
        self._canvas_layout = QVBoxLayout(self._canvas_host)
        scroll.setWidget(self._canvas_host)
        root.addWidget(scroll, stretch=1)

        self._on_metric_changed()
        self._update_plot()

    def _current_metric(self) -> str:
        return str(self._metric_combo.currentData())

    def _current_channel(self) -> int | str:
        return _parse_channel_value(str(self._channel_combo.currentData()))

    def _current_hemisphere_mode(self) -> str:
        return str(self._hemisphere_combo.currentData() or "whole")

    def _current_view(self) -> str:
        mode = self._current_hemisphere_mode()
        if mode in {"compare", "difference"}:
            return mode
        return "single"

    def _current_segments(self) -> list[str]:
        start = str(self._segment_start.currentText())
        end = str(self._segment_end.currentText())
        return parse_segment_range(f"{start}:{end}", atlas_dir=self._atlas_dir)

    def _on_metric_changed(self) -> None:
        metric = self._current_metric()
        channels = self._options["channels_by_metric"].get(metric, [])
        self._channel_combo.blockSignals(True)
        self._channel_combo.clear()
        for channel in channels:
            label = f"Channel {channel}" if channel.isdigit() else channel
            self._channel_combo.addItem(label, channel)
        self._channel_combo.blockSignals(False)

    def _rollup_x_label(self, rollup_level: str) -> str:
        return _X_LABELS_BY_ROLLUP.get(rollup_level.lower(), "Region")

    def _column_labels(self, matrix: pd.DataFrame) -> list[str]:
        return [structure_acronym_to_label(col) for col in matrix.columns]

    def _build_single_matrix(self, *, hemisphere: str, segments: list[str]) -> pd.DataFrame:
        return cord_stats_to_matrix(
            self._stats_df,
            metric=self._current_metric(),  # type: ignore[arg-type]
            channel=self._current_channel(),
            rollup_level=str(self._rollup_combo.currentText()),
            hemisphere=hemisphere,
            segments=segments,
            regions_df=self._regions_df,
        )

    def _plot_title(self, *, hemisphere: str | None = None) -> str:
        metric = self._current_metric()
        channel = self._current_channel()
        title = f"{self._config.sample.name} — {metric_label(metric)} (ch {channel})"
        if hemisphere and hemisphere != "whole":
            title = f"{title}, {hemisphere}"
        return title

    def _update_plot(self) -> None:
        from qtpy.QtWidgets import QMessageBox

        try:
            segments = self._current_segments()
            rollup_level = str(self._rollup_combo.currentText())
            x_label = self._rollup_x_label(rollup_level)
            metric = self._current_metric()
            cmap = _default_cmap(metric)
            normalize = str(self._normalize_combo.currentData())
            view = self._current_view()
            if view == "compare":
                matrices = cord_hemisphere_matrices(
                    self._stats_df,
                    metric=metric,  # type: ignore[arg-type]
                    channel=self._current_channel(),
                    rollup_level=rollup_level,
                    segments=segments,
                    regions_df=self._regions_df,
                )
                first = next(iter(matrices.values()))
                fig = build_cord_heatmap_comparison_figure(
                    matrices,
                    column_labels=self._column_labels(first),
                    title=self._plot_title(),
                    cmap=cmap,
                    normalize=normalize,
                    x_label=x_label,
                )
                self._export_kind = "compare"
                self._export_matrices = matrices
                status = (
                    f"{first.shape[0]} segments × {first.shape[1]} regions "
                    f"(left | right | whole)"
                )
            elif view == "difference":
                matrices = cord_hemisphere_matrices(
                    self._stats_df,
                    metric=metric,  # type: ignore[arg-type]
                    channel=self._current_channel(),
                    rollup_level=rollup_level,
                    segments=segments,
                    regions_df=self._regions_df,
                    hemispheres=("left", "right"),
                )
                diff = difference_heatmap_matrix(matrices["left"], matrices["right"])
                vmin, vmax = difference_color_limits(diff)
                fig = build_cord_heatmap_figure(
                    diff,
                    column_labels=self._column_labels(diff),
                    title=f"{self._plot_title()} — left − right",
                    cmap=DIFFERENCE_CMAP,
                    vmin=vmin,
                    vmax=vmax,
                    x_label=x_label,
                )
                self._export_kind = "matrix"
                self._export_matrices = {"difference": diff}
                status = f"{diff.shape[0]} segments × {diff.shape[1]} regions (left − right)"
            else:
                hemisphere = self._current_hemisphere_mode()
                matrix = self._build_single_matrix(hemisphere=hemisphere, segments=segments)
                fig = build_cord_heatmap_figure(
                    matrix,
                    column_labels=self._column_labels(matrix),
                    title=self._plot_title(hemisphere=hemisphere),
                    cmap=cmap,
                    normalize=normalize,
                    x_label=x_label,
                )
                self._export_kind = "matrix"
                self._export_matrices = {"data": matrix}
                status = f"{matrix.shape[0]} segments × {matrix.shape[1]} regions"
        except (ValueError, KeyError) as exc:
            self._status_label.setText(str(exc))
            QMessageBox.warning(self.widget, "Could not plot", str(exc))
            return

        if self._figure is not None:
            plt.close(self._figure)
        if self._canvas is not None:
            self._canvas_layout.removeWidget(self._canvas)
            self._canvas.deleteLater()
            self._canvas = None

        self._figure = fig
        self._canvas = FigureCanvasQTAgg(fig)
        self._canvas_layout.addWidget(self._canvas)
        self._status_label.setText(status)

    def _default_export_path(self, suffix: str) -> Path:
        return default_heatmap_output_path(
            self._config,
            metric=self._current_metric(),
            channel=self._current_channel(),
            view=self._current_hemisphere_mode() if self._current_hemisphere_mode() != "whole" else "single",
            suffix=suffix,
        )

    def _save_png(self) -> None:
        from qtpy.QtWidgets import QFileDialog, QMessageBox

        if self._figure is None:
            self._update_plot()
        if self._figure is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self.widget,
            "Save heatmap",
            str(self._default_export_path(".png")),
            "PNG images (*.png)",
        )
        if not path:
            return
        try:
            output = Path(path)
            output.parent.mkdir(parents=True, exist_ok=True)
            self._figure.savefig(output, bbox_inches="tight", pad_inches=0.08, dpi=300)
        except OSError as exc:
            QMessageBox.warning(self.widget, "Could not save", str(exc))
            return
        self._status_label.setText(f"Saved {output}")

    def _save_csv(self) -> None:
        from qtpy.QtWidgets import QFileDialog, QMessageBox

        if not self._export_matrices:
            self._update_plot()
        if not self._export_matrices:
            return
        path, _ = QFileDialog.getSaveFileName(
            self.widget,
            "Save heatmap matrix",
            str(self._default_export_path(".csv")),
            "CSV files (*.csv)",
        )
        if not path:
            return
        try:
            output = Path(path)
            if self._export_kind == "compare":
                saved = write_cord_heatmap_compare_csv(self._export_matrices, output)
            else:
                matrix = next(iter(self._export_matrices.values()))
                saved = write_cord_heatmap_matrix_csv(matrix, output)
        except (ValueError, OSError) as exc:
            QMessageBox.warning(self.widget, "Could not save", str(exc))
            return
        self._status_label.setText(f"Saved {saved}")

    def teardown(self) -> None:
        if self._figure is not None:
            plt.close(self._figure)
            self._figure = None
        self._canvas = None


def attach_spinal_stats_plots(
    viewer: Any,
    config: SpinalCordPipelineConfig,
    *,
    stats_path: Path | None = None,
) -> DockStageController:
    """Attach the stats/plots dock widgets to a napari viewer."""
    resolved = resolve_region_stats_path(config, input_path=stats_path)
    stats_df = pd.read_csv(resolved)
    panel = CordStatsPlotsPanel(config, stats_df=stats_df, stats_path=resolved)

    def _teardown() -> None:
        panel.teardown()

    return DockStageController(
        dock_widgets=[(panel.widget, "Stats / plots")],
        _teardown_fn=_teardown,
        open_log_message=(
            f"Stats / plots — {config.sample.name} "
            f"({len(stats_df)} rows from {resolved.name})"
        ),
    )


def run_spinal_stats_plots(
    config: SpinalCordPipelineConfig,
    *,
    headless: bool = False,
    stats_path: Path | None = None,
) -> Any:
    """Open the stats/plots panel in napari or the unified shell."""
    if headless:
        msg = "Stats / plots requires an interactive GUI session."
        raise RuntimeError(msg)

    title = f"LightSuite spinal — {config.sample.name} stats"
    return run_attached_stage(
        title,
        lambda viewer: attach_spinal_stats_plots(viewer, config, stats_path=stats_path),
    )
